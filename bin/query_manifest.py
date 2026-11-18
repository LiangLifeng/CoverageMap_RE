#!/usr/bin/env python3
"""Auditable multi-sample TSV queries against an existing native reference index."""
import argparse
from collections import Counter, defaultdict
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor, wait, FIRST_COMPLETED
from contextlib import nullcontext
import csv
import gzip
import hashlib
import io
import json
import math
import os
from pathlib import Path
import resource
import sys
import tempfile
import time
import traceback
from types import SimpleNamespace

# Set before importing NumPy, including spawned workers. One explicit jobs budget.
for _name in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ[_name] = '1'
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import numpy as np
from bin.bin_count_eval_math import UndefinedCorrelationError, distances, exact_edges
from bin.coverage_core import BlockPool, canonical_reference, compatible_policies, runtime_pool, extract_pool, layout_metadata, merge_references, reference_from_layout, select_units, window_sums
from bin.gcf_windows import COORDINATE_COLUMNS, GCF_CONCAT, PER_CONTIG, coordinate_rows, edges_for
from bin.native_concordance import NativeCache
from bin.native_reference import covered_sums, native_edges, validate_native
from bin.query_detection import DOMAINS, coverage_profile, read_rules, reference_center, select_rule
from bin.query_manifest_stats import similarity, aggregate, low_null, reference_row_means
from bin.jaccard_breadth import thresholds
from bin.query_decision_methods import COLUMNS as METHOD_COLUMNS, evaluate_saved, parse_methods, validate_options
from bin.ref_depth_catalog import atomic_json, digest, file_stat
from bin.reference_index import ReferenceIndex
from bin.stage_timing import StageTiming

REQUIRED = ('Sample', 'GCFID', 'Depth', 'CramPath')
DRAW_COLUMNS = ('jaccard_coverage_threshold', 'query_covered_bins', 'reference_covered_bins', 'query_id', 'Sample', 'GCFID', 'draw_id', 'replicate_id', 'reference_id',
                'reference_sample', 'source_root', 'manifest_path', 'source_sha256', 'same_sample',
                'same_source', 'reference_full_depth', 'requested_depth', 'actual_depth', 'units',
                'selected_ids_sha256', 'rng_sha256', 'status', 'error', 'jaccard', 'pearson',
                'intersection', 'union', 'jaccard_status', 'pearson_status',
                'sampling_unit', 'coverage_alignment_policy', 'target_bases', 'target_depth', 'actual_bases',
                'overshoot_bases', 'selected_units', 'selected_reads', 'selected_fragments',
                'selected_segments', 'selected_blocks')
SUMMARY_COLUMNS = ('jaccard_coverage_threshold', 'query_covered_bins', 'jaccard_zero_union_count', 'jaccard_null_zero_union_count', 'jaccard_zero_union_fraction', 'jaccard_null_zero_union_fraction', 'query_id', 'Sample', 'GCFID', 'status', 'error', 'arm', 'manifest_depth',
                   'actual_query_depth', 'manifest_coverage', 'actual_query_coverage', 'available_count',
                   'requested_draws', 'valid_draws', 'unique_references', 'duplicate_draws', 'windows',
                   'layout_sha256', 'gcf_reference_length', 'nominal_window_size', 'requested_bin_count', 'window_policy',
                   'high_null_policy', 'high_null_candidate_pairs', 'high_null_sampled_pairs',
                   'jaccard_mean', 'jaccard_max', 'jaccard',
                   'jaccard_p', 'jaccard_status', 'jaccard_null_count', 'jaccard_null_mean', 'jaccard_null_sd',
                   'pearson_mean', 'pearson_max', 'pearson', 'pearson_p', 'pearson_status',
                   'pearson_null_count', 'pearson_null_mean', 'pearson_null_sd')
NULL_COLUMNS = ('metric', 'left_draw', 'right_draw', 'value', 'rng_sha256',
                'left_reference_sample', 'right_reference_sample', 'left_reference_id',
                'right_reference_id', 'intersection', 'union')
DETECTION_COLUMNS = ('query_id', 'Sample', 'GCFID', 'status', 'error', 'domain', 'gcf_reference_length',
                     'min_genome_bp', 'max_genome_bp_exclusive', 'bin_count', 'distance', 'distance_value',
                     'max_distance', 'detection_status', 'available_count', 'center_reference_count',
                     'excluded_self_count', 'manifest_depth', 'actual_query_depth')
REFERENCE_POOL_SIZE = 100
DISTRIBUTION_COLUMNS = tuple(column for column in SUMMARY_COLUMNS if not column.startswith('high_null_')) + (
    'excluded_self_count', 'eligible_count', 'alpha', 'jaccard_decision', 'pearson_decision', 'detection_status')


def rng(seed, *parts):
    token = hashlib.sha256(json.dumps([seed, *parts], ensure_ascii=False, separators=(',', ':')).encode()).digest()
    return np.random.default_rng(np.frombuffer(token, dtype='<u4')), token.hex()


def draw_count(specification, available):
    if specification.endswith('%'):
        percent = float(specification[:-1])
        if not math.isfinite(percent) or not 0 < percent <= 100:
            raise ValueError('reference-draws percent must be in (0,100]')
        return math.ceil(available*percent/100)
    if not specification.isdecimal() or int(specification) < 1:
        raise ValueError('reference-draws must be a positive count or (0,100]%')
    return int(specification)


def windows(layout, size=None, count=None):
    """Legacy manifest windows for existing per-contig indexes."""
    edges = edges_for(layout, PER_CONTIG, size=size, count=count)
    width = size if size is not None else (layout.length+count-1)//count
    rows = [dict(contig=r['contig'], start=r['local_start'], end=r['local_end'],
                 width=r['width'], nominal_width=width)
            for r in coordinate_rows(layout, edges, width)]
    return edges, rows


def table(path, columns, rows):
    """Deterministic gzip headers, TSV LF, and blank missing values."""
    if str(path).endswith('.gz'):
        raw = open(path, 'wb')
        handle = io.TextIOWrapper(gzip.GzipFile(filename='', fileobj=raw, mode='wb', mtime=0), encoding='utf-8', newline='')
    else:
        raw = None
        handle = open(path, 'w', encoding='utf-8', newline='')
    try:
        writer = csv.DictWriter(handle, columns, delimiter='\t', lineterminator='\n', extrasaction='raise')
        writer.writeheader()
        writer.writerows(rows)
    finally:
        handle.close()
        if raw is not None:
            raw.close()


def read_manifest(path, coverage_unit=None, query_id=None):
    with Path(path).open(encoding='utf-8-sig', newline='') as handle:
        reader = csv.DictReader(handle, delimiter='\t')
        if not reader.fieldnames or len(set(reader.fieldnames)) != len(reader.fieldnames) or not set(REQUIRED) <= set(reader.fieldnames):
            raise ValueError('TSV requires unique headers: '+', '.join(REQUIRED))
        rows = []
        for number, row in enumerate(reader, 1):
            row['query_id'] = f'q{number:08d}'
            if query_id is not None and row['query_id'] != query_id:
                continue
            row['_error'] = ''
            extra_fields = row.pop(None, None)
            if extra_fields is not None:
                row['_extra_fields'] = extra_fields
            try:
                if extra_fields is not None or any(row.get(k) is None for k in reader.fieldnames):
                    raise ValueError('Wrong number of TSV fields')
                if any(not row[k].strip() for k in REQUIRED):
                    raise ValueError('Empty required value')
                depth = float(row['Depth'])
                if not math.isfinite(depth) or depth < 0:
                    raise ValueError('Depth must be finite and nonnegative')
                for key in ('CramPath', 'CraiPath', 'BedGraphPath'):
                    if row.get(key):
                        candidate = Path(row[key])
                        if not candidate.is_absolute() or not candidate.is_file() or not candidate.stat().st_size:
                            raise ValueError(f'{key} must be an absolute nonempty file: {candidate}')
                if row.get('Coverage'):
                    value = float(row['Coverage'])
                    if coverage_unit not in ('fraction', 'percent'):
                        raise ValueError('Nonempty Coverage requires explicit --coverage-unit fraction|percent')
                    maximum = 1 if coverage_unit == 'fraction' else 100
                    if not math.isfinite(value) or not 0 <= value <= maximum:
                        raise ValueError(f'Coverage outside [0,{maximum}] for {coverage_unit}')
            except (ValueError, OSError) as error:
                row['_error'] = str(error)
            rows.append(row)
            if query_id is not None:
                break
    if not rows:
        if query_id is not None:
            raise ValueError(f'Query ID not found in manifest: {query_id}')
        raise ValueError('Manifest contains no query rows')
    return rows


def identity(path):
    before = file_stat(path)
    sha = digest(path)
    if before != file_stat(path):
        raise ValueError(f'Input changed while hashing: {path}')
    return dict(path=str(Path(path).resolve()), stat=before, sha256=sha)


def check_identity(record):
    if identity(record['path']) != record:
        raise ValueError('Input identity changed: '+record['path'])


def code_identity():
    root = Path(__file__).resolve().parents[1]
    files = [root/'Tools.py', *sorted((root/'bin').glob('*.py'))]
    return {str(p.relative_to(root)): digest(p) for p in files}


def source_draws(rows, count, seed, query):
    generator, token = rng(seed, query['GCFID'], query['Sample'], query['query_id'], 'reference-samples')
    return generator.integers(0, len(rows), size=count), token


def require_v2_index(index):
    if index.metadata.get('storage_mode') != 'native':
        raise ValueError('v2 requires a native index with verified reference geometry')


def align_reference_geometry(pool, source, target):
    """Rebase legacy contig offsets, preserving unit identities and every block length."""
    if source.sort_policy == target.sort_policy:
        merge_references(target, source)
        return pool, source
    aligned = canonical_reference(source)
    merge_references(target, aligned)
    if source.sort_policy is not None or pool.genome_len != target.length:
        raise ValueError('Unsupported reference geometry coordinate conversion')
    if source.names == target.names:
        return pool, aligned
    offsets = np.asarray(list(source.offsets.values()), dtype=np.int64)
    ends = offsets + np.asarray(source.lengths, dtype=np.int64)
    contigs = np.searchsorted(offsets, pool.starts, side='right') - 1
    if np.any(pool.ends > ends[contigs]):
        raise ValueError('Reference aligned block crosses a contig boundary')
    target_offsets = target.offsets
    shifts = np.asarray([target_offsets[name] - offset
                         for name, offset in zip(source.names, offsets)], dtype=np.int64)
    shifted = shifts[contigs]
    result = BlockPool(pool.starts.astype(np.int64) + shifted,
                       pool.ends.astype(np.int64) + shifted, pool.indptr, pool.genome_len)
    if not np.array_equal(result.unit_bases, pool.unit_bases):
        raise AssertionError('Coordinate conversion changed sampling-unit aligned bases')
    # Unit identities belong to whole Reads/Fragments, not their contig offsets.
    # Keep the Geometry reference for exact Mate/Segment audit after rebinning.
    for name in ('alignment_geometry', 'geometry_unit_ids', 'sampling_unit', 'coverage_alignment_policy'):
        if hasattr(pool, name):
            setattr(result, name, getattr(pool, name))
    return result, aligned


def reference_vector(cache, row, layout, edges, timing, policies=None):
    if policies is not None and ('coverage_alignment_policy' in policies or policies.get('sampling_unit') != cache.index.metadata['config']['policies'].get('sampling_unit', 'read')):
        stamps = cache.stamps(row)
        pool, other, record, _ = cache.load(row, timing)
        pool = runtime_pool(pool, record['policies'], policies)
        pool, other = align_reference_geometry(pool, other, layout)
        merge_references(layout, other)
        return window_sums(pool, np.arange(len(pool)), edges), dict(record, units=len(pool)), stamps
    stamps = cache.stamps(row)
    key = ('query-full-windows', cache.index.window_policy, row['id'], layout.names,
           layout.lengths, layout.md5, hashlib.sha256(edges.tobytes()).hexdigest())
    value = cache.get(key, stamps)
    if value is None:
        attempt = (cache.index.root_for(row)/row['manifest_path']).parent
        with timing.measure('reference_metadata'):
            manifest = attempt/'result.json'
            if digest(manifest) != row['manifest_sha'] or digest(attempt/'checksums.json') != row['checksums_sha']:
                raise ValueError('Reference metadata differs from SQLite')
            record = json.loads(manifest.read_text())
            cache.check_window_policy(record)
            other = reference_from_layout(record['layout'])
            merge_references(layout, canonical_reference(other) if layout.sort_policy else other)
        if other.names == layout.names and np.array_equal(edges, native_edges(other, record)):
            # Exact existing window sums need no geometry decompression. Verify every
            # USED payload against the registered checksum ledger, not the whole pool.
            with timing.measure('reference_native_windows'):
                checks = json.loads((attempt/'checksums.json').read_text())
                for name in ('native_coverage.npz', 'edges.npy'):
                    if digest(attempt/name) != checks[name]:
                        raise ValueError('Native window payload checksum mismatch')
                validate_native(attempt, record)
                with np.load(attempt/'native_coverage.npz', allow_pickle=False) as saved:
                    value = saved['aligned_bases'].copy()
        else:
            pool, other, record, _ = cache.load(row, timing)
            pool, other = align_reference_geometry(pool, other, layout)
            merge_references(layout, other)
            value = window_sums(pool, np.arange(len(pool)), edges)
        if (record['gcf'] != row['gcf'] or record['sample'] != row['sample']
                or record['full_mean_depth'] != row['full_mean_depth'] or stamps != cache.stamps(row)):
            raise ValueError('Reference identity changed or differs from SQLite')
        cache.put(key, (value, record), value.nbytes+len(json.dumps(record)), stamps)
    else:
        value, record = value
    return value, record, stamps


def reference_coverage(cache, row, layout, edges, timing, policies):
    stamps = cache.stamps(row)
    key = ('query-breadth-v1', row['id'], layout.names, layout.lengths, layout.md5,
           hashlib.sha256(edges.tobytes()).hexdigest(), json.dumps(policies, sort_keys=True))
    cached = cache.get(key, stamps)
    if cached is None:
        pool, other, record, _ = cache.load(row, timing)
        pool = runtime_pool(pool, record['policies'], policies)
        pool, other = align_reference_geometry(pool, other, layout)
        merge_references(layout, other)
        if not compatible_policies(record['policies'], policies):
            raise ValueError('Reference/query read policy mismatch')
        vector = window_sums(pool, np.arange(len(pool)), edges)
        covered = covered_sums(pool, edges)
        record = dict(record, units=len(pool))
        if stamps != cache.stamps(row):
            raise ValueError('Reference changed during coverage computation')
        cached = vector, covered, record
        cache.put(key, cached, vector.nbytes+covered.nbytes+len(json.dumps(record)), stamps)
    return *cached, stamps


def pair_population_size(draws):
    """Unordered distinct-Sample draw pairs, retaining repeated draw weights."""
    counts = Counter(d['reference_sample'] for d in draws if d['status'] == 'READY')
    size = sum(counts.values())
    return (size*(size-1)-sum(n*(n-1) for n in counts.values()))//2


def high_null(matrix, draws, config, query, *, covered=None, lengths=None, threshold=0.):
    """User-defined null: up to 1000 unordered distinct-Sample draw pairs."""
    if config['high_null'] != 'draw-pairs':
        return {}, []
    # Reservoir sample pair identities without O(draws squared) memory.
    generator, token = rng(config['seed'], query['query_id'], 'high-null-pairs')
    pairs, seen = [], 0
    for i in range(len(draws)):
        for j in range(i+1, len(draws)):
            if (draws[i]['status'] != 'READY' or draws[j]['status'] != 'READY'
                    or draws[i]['reference_sample'] == draws[j]['reference_sample']):
                continue
            seen += 1
            if len(pairs) < config['null_limit']:
                pairs.append((i, j))
            else:
                place = int(generator.integers(0, seen))
                if place < config['null_limit']:
                    pairs[place] = (i, j)
    pairs.sort()
    values = {'jaccard': [], 'pearson': []}
    records = []
    for i, j in pairs:
        result = similarity(matrix[i], matrix[j], covered_x=None if covered is None else covered[i],
                            covered_y=None if covered is None else covered[j], lengths=lengths, threshold=threshold)
        for metric in values:
            values[metric].append(result[metric])
            records.append(dict(metric=metric, left_draw=i, right_draw=j, value=result[metric], rng_sha256=token,
                                left_reference_sample=draws[i]['reference_sample'],
                                right_reference_sample=draws[j]['reference_sample'],
                                left_reference_id=draws[i].get('reference_id', ''),
                                right_reference_id=draws[j].get('reference_id', ''),
                                intersection=result['intersection'], union=result['union']))
    return values, records


def extract_query(query, layout, output, policies, timing):
    fai = output/'layout.fai'
    fai.write_text(''.join(f'{n}\t{v}\n' for n, v in zip(layout.names, layout.lengths)))
    inputs = {key: identity(query[key]) for key in ('CramPath', 'CraiPath', 'BedGraphPath') if query.get(key)}
    with timing.measure('query_extract'):
        order = {'sort_contigs': True} if layout.sort_policy else {}
        capture = {'capture_geometry': output/'query_alignment_geometry'} if 'coverage_alignment_policy' in policies else {}
        pool, actual_layout = extract_pool(query['CramPath'], reference_index=fai, **order, **policies, **capture)
    merge_references(layout, actual_layout)
    return pool, actual_layout, inputs


def sample_reference_windows(ref, generator, target_bases, edges, audit=None, *, covered_output=None):
    ids = select_units(ref, generator, target_bases=target_bases, depth_sampler='permutation')
    vector = window_sums(ref, ids, edges)
    if covered_output is not None:
        covered_output[:] = covered_sums(ref, edges, ids)
    if audit is not None:
        from bin.alignment_geometry import selection_audit
        audit.update(selection_audit(ref, ids, target_bases))
    if int(vector.sum()) != int(ref.unit_bases[ids].sum()):
        raise AssertionError('Sampled aligned-base mass differs')
    return vector, len(ids), hashlib.sha256(ids.astype('<i8').tobytes()).hexdigest()


def fill_distribution_draws(groups, matrix, test, layout, edges, cache, timing, *,
                            arm, target_bases, policies, threads, covered_matrix=None, query_covered=None, threshold=0.):
    """Load bounded batches of unique sources and process their independent draw slots."""
    # StageTiming accounts for process-wide CPU. Measure each joined stage on the
    # controller rather than sharing its mutable state across the worker threads.
    worker_timing = SimpleNamespace(measure=lambda _stage: nullcontext())

    def prepare(group):
        row, _, _ = group[0]
        try:
            if arm == 'high':
                if covered_matrix is None:
                    data, record, stamps = reference_vector(cache, row, layout, edges, worker_timing, policies)
                else:
                    vector, covered, record, stamps = reference_coverage(cache, row, layout, edges, worker_timing, policies)
                    data = (vector, covered)
                units = record['units']
            else:
                stamps = cache.stamps(row)
                data, other, record, _ = cache.load(row, worker_timing)
                data = runtime_pool(data, record['policies'], policies)
                units = len(data) if hasattr(data, 'alignment_geometry') else record['units']
                data, other = align_reference_geometry(data, other, layout)
                merge_references(layout, other)
            if not compatible_policies(record['policies'], policies):
                raise ValueError('Reference/query read policy mismatch')
        except (ValueError, OSError) as error:
            for _, entry, _ in group:
                entry.update(status='ERROR', error=str(error))
            return None
        return group, data, units, (row, stamps)

    def fill(task):
        entry, generator, data, units = task
        try:
            if arm == 'low':
                extra = {} if covered_matrix is None else dict(covered_output=covered_matrix[entry['draw_id']])
                data, units, selected_sha = sample_reference_windows(data, generator, target_bases, edges, audit=entry, **extra)
                entry['selected_ids_sha256'] = selected_sha
            if arm == 'high' and covered_matrix is not None:
                data, covered = data
                covered_matrix[entry['draw_id']] = covered
            matrix[entry['draw_id']] = data
            entry.update(status='READY', actual_depth=int(data.sum())/layout.length, units=units,
                         **similarity(test, data, covered_x=query_covered,
                                      covered_y=None if covered_matrix is None else covered_matrix[entry['draw_id']],
                                      lengths=np.diff(edges), threshold=threshold))
        except (ValueError, OSError) as error:
            entry.update(status='ERROR', error=str(error))

    sources = list(groups.values())
    workers = min(threads, len(matrix))
    verified = {}
    with ThreadPoolExecutor(max_workers=workers) as executor:
        for start in range(0, len(sources), workers):
            with timing.measure('reference_prepare'):
                prepared = [result for result in executor.map(prepare, sources[start:start+workers])
                            if result is not None]
            with timing.measure('reference_draws'):
                tasks = []
                for group, data, units, (row, stamps) in prepared:
                    verified[row['id']] = (row, stamps)
                    tasks.extend((entry, generator, data, units) for _, entry, generator in group)
                for _ in executor.map(fill, tasks):
                    pass
            # Drop active reference buffers before preparing the next batch.
            tasks.clear()
            prepared.clear()
            data = None
    return verified


def calculate_detection(query, config, index, cache, output):
    summary = dict(query_id=query['query_id'], Sample=query.get('Sample', ''), GCFID=query.get('GCFID', ''),
                   status='ERROR', error=query['_error'], domain=config['domain'],
                   manifest_depth=query.get('Depth', ''), detection_status='NOT_EVALUATED')
    if query['_error']:
        return summary, {}
    group = index.group(query['GCFID'], missing_ok=True)
    if group is None:
        summary.update(status='GCF_NOT_INDEXED', available_count=0)
        return summary, {}
    layout = reference_from_layout(group)
    summary['gcf_reference_length'] = layout.length
    rule = select_rule(config['detection_rules'], layout.length)
    if rule is None:
        summary.update(status='NO_DETECTION_RULE')
        return summary, {}
    summary.update(rule)
    bins = rule['bin_count']
    if bins > layout.length:
        summary.update(status='INVALID_BIN_COUNT', error='Bin count exceeds the GCF reference length')
        return summary, {}
    candidates = index.candidates(query['GCFID'])
    summary['available_count'] = len(candidates)
    if len(candidates) <= config['min_samples']:
        summary.update(status='INSUFFICIENT_REFERENCE_SAMPLES', error='Requires available_count > min-samples')
        return summary, {}

    timing = StageTiming()
    pool, actual_layout, inputs = extract_query(query, layout, output, config['policies'], timing)
    summary['actual_query_depth'] = pool.total_bases / layout.length
    references = [row for row in candidates if row['sample'] != query['Sample']
                  and row['source_sha'] != inputs['CramPath']['sha256']]
    summary.update(center_reference_count=len(references), excluded_self_count=len(candidates)-len(references))
    audit = dict(inputs=inputs, policies=config['policies'],
                 missing_header_contigs=list(actual_layout.missing_header_contigs),
                 unknown_m5_contigs=[n for n, a, b in zip(layout.names, layout.md5, actual_layout.md5) if not a or not b],
                 center='arithmetic mean of normalized full reference profiles; exclude same Sample or source SHA256')
    if pool.total_bases == 0:
        summary.update(status='ZERO_QUERY_DEPTH', error='No aligned query bases')
    elif len(references) <= config['min_samples']:
        summary.update(status='INSUFFICIENT_REFERENCE_SAMPLES', error='Requires center_reference_count > min-samples after self exclusion')
    else:
        edges = exact_edges(layout.length, bins)  # Frozen grid-distance thresholds use legacy bins.
        widths = np.diff(edges)
        reference_stamps = []

        def profiles():
            for row in references:
                vector, record, stamps = reference_vector(cache, row, actual_layout, edges, timing)
                if not compatible_policies(record['policies'], config['policies']) or int(vector.sum()) != record['aligned_bases']:
                    raise ValueError('Reference read policy or aligned-base mass mismatch')
                reference_stamps.append(stamps)
                yield coverage_profile(vector, widths)

        center = reference_center(profiles(), len(references), bins, config['work_bytes'])
        mass = window_sums(pool, np.arange(len(pool)), edges)
        if int(mass.sum()) != pool.total_bases:
            raise ValueError('Query rebinning mass conservation failed')
        profile = coverage_profile(mass, widths)
        np.save(output/'bin_edges.npy', edges, allow_pickle=False)
        np.save(output/'query_profile.npy', profile, allow_pickle=False)
        np.save(output/'reference_center.npy', center, allow_pickle=False)
        audit['center_references'] = [dict(reference_id=r['id'], sample=r['sample'], source_sha256=r['source_sha'])
                                      for r in references]
        for row, stamps in zip(references, reference_stamps):
            if stamps != cache.stamps(row):
                raise ValueError('Reference changed during detection query')
        audit['reference_input_stamps'] = reference_stamps
        try:
            value = float(distances(profile, center, rule['distance'])[0, 0])
        except UndefinedCorrelationError as error:
            summary.update(status='DISTANCE_UNDEFINED', error=str(error))
            audit['distance_diagnostics'] = error.issues
        else:
            threshold = rule['max_distance']
            result = 'THRESHOLD_UNSET' if threshold is None else 'DETECTED' if value <= threshold else 'NOT_DETECTED'
            summary.update(status='COMPLETE', distance_value=value, detection_status=result)
    for item in inputs.values():
        check_identity(item)
    audit['stage_timing'] = timing.snapshot()
    return summary, audit


def calculate(query, config, index, cache, output):
    if config.get('detection_rules') is not None:
        return calculate_detection(query, config, index, cache, output)
    threshold = config.get('jaccard_coverage_threshold', 0.)
    alpha = config.get('alpha')
    distribution = alpha is not None
    summary = dict(query_id=query['query_id'], Sample=query.get('Sample', ''), GCFID=query.get('GCFID', ''),
                   status='ERROR', error=query['_error'], manifest_depth=query.get('Depth', ''),
                   manifest_coverage=query.get('Coverage', ''))
    summary['jaccard_coverage_threshold'] = threshold
    if distribution:
        summary.update(alpha=alpha, detection_status='NOT_EVALUATED', requested_draws=config.get('reference_draws', REFERENCE_POOL_SIZE))
    else:
        summary['high_null_policy'] = config['high_null']
    if query['_error']:
        return summary, {}
    layout_dict = index.group(query['GCFID'], missing_ok=True)
    candidates = index.candidates(query['GCFID']) if layout_dict is not None else []
    summary['available_count'] = len(candidates)
    if len(candidates) <= config['min_samples']:
        summary.update(status='INSUFFICIENT_REFERENCE_SAMPLES', error='Requires available_count > min-samples')
        return summary, {}
    layout = reference_from_layout(layout_dict)
    policy = GCF_CONCAT if distribution else index.window_policy
    if distribution:
        layout = canonical_reference(layout)
    size, count = config['window_size'], config['bin_count']
    if distribution:
        require_v2_index(index)
    if policy == GCF_CONCAT:
        edges = edges_for(layout, policy, size=size, count=count)
        nominal = size
        table(output/'windows.tsv', COORDINATE_COLUMNS, coordinate_rows(layout, edges, nominal))
        np.save(output/'edges.npy', edges, allow_pickle=False)
    else:
        edges, coordinates = windows(layout, size, count)
        nominal = coordinates[0]['nominal_width']
        table(output/'windows.tsv', ('contig', 'start', 'end', 'width', 'nominal_width'), coordinates)
    summary.update(windows=len(edges)-1, layout_sha256=digest(output/'windows.tsv'),
                   gcf_reference_length=layout.length, nominal_window_size=nominal,
                   requested_bin_count=count, window_policy=policy)
    timing = StageTiming()
    pool, actual_layout, inputs = extract_query(query, layout, output, config['policies'], timing)
    test = window_sums(pool, np.arange(len(pool)), edges)
    depth = pool.total_bases/layout.length
    query_covered = covered_sums(pool, edges)
    coverage = int(query_covered.sum())/layout.length
    arm = 'high' if float(query['Depth']) > config['depth_threshold'] else 'low'
    summary.update(actual_query_depth=depth, actual_query_coverage=coverage, arm=arm)
    audit = dict(inputs=inputs, missing_header_contigs=list(actual_layout.missing_header_contigs),
                 unknown_m5_contigs=[n for n, a, b in zip(layout.names, layout.md5, actual_layout.md5) if not a or not b],
                 depth_branch_basis='manifest Depth (strict > threshold)',
                 depth_target_basis='actual CRAM aligned bases / full layout length',
                 depth_difference=depth-float(query['Depth']), policies=config['policies'],
                 window_value='sum', window_policy=policy, layout=layout_metadata(actual_layout),
                 reference_window_policy=index.window_policy,
                 reference_geometry_rebinned=index.window_policy != policy,
                 statistics_calibrated=False, gcf_reference_length=layout.length,
                 nominal_window_size=nominal, requested_bin_count=count)
    audit.update(query_numeric_pool_bytes=pool.nbytes, coverage_unit=config['coverage_unit'],
                 reference_verification='registered metadata, status/stamps and used payload SHA256; geometry only for low arm or rebin')
    layout = merge_references(layout, actual_layout)
    if distribution:
        references = [row for row in candidates if row['sample'] != query['Sample']
                      and row['source_sha'] != inputs['CramPath']['sha256']]
        summary['excluded_self_count'] = len(candidates) - len(references)
        candidates = references if arm == 'high' else [
            row for row in references if row['aligned_bases'] >= pool.total_bases]
        summary['eligible_count'] = len(candidates)
    np.save(output/'query_windows.npy', test, allow_pickle=False)
    np.save(output/'query_covered_bases.npy', query_covered, allow_pickle=False)
    np.save(output/'window_lengths.npy', np.diff(edges), allow_pickle=False)
    if depth <= 0:
        summary.update(status='ZERO_QUERY_DEPTH', error='No aligned query bases')
        for item in inputs.values():
            check_identity(item)
        return summary, audit
    if distribution and len(candidates) <= config['min_samples']:
        status = ('INSUFFICIENT_REFERENCE_SAMPLES' if len(references) <= config['min_samples']
                  else 'INSUFFICIENT_REFERENCE_DEPTH')
        summary.update(status=status, error='Requires eligible_count > min-samples after reference filtering')
        for item in inputs.values():
            check_identity(item)
        return summary, audit
    # Actual CRAM mean depth defines matching simulation mass; supplied Depth defines routing.
    if distribution:
        requested = config.get('reference_draws', REFERENCE_POOL_SIZE)
    else:
        requested = config['simulations'] if arm == 'low' else draw_count(config['reference_draws'], len(candidates))
    summary['requested_draws'] = requested
    matrix_bytes = requested*(len(edges)-1)*8
    if config['work_bytes'] is not None and matrix_bytes*5 > config['work_bytes']:
        raise ValueError(f'Window matrix and statistics scratch exceed work budget: {matrix_bytes*5} bytes')
    matrix = np.zeros((requested, len(edges)-1), dtype=np.int64)
    covered_matrix = np.zeros_like(matrix)
    choices, selection_token = source_draws(candidates, requested, config['seed'], query)
    draws, verified = [], {}
    groups = defaultdict(list)
    for number, choice in enumerate(choices):
        row = candidates[int(choice)]
        generator, token = rng(config['seed'], query['query_id'], query['Sample'], query['GCFID'], 'read-units', number)
        entry = dict(query_id=query['query_id'], Sample=query['Sample'], GCFID=query['GCFID'], draw_id=number,
                     replicate_id=number if arm == 'low' else '', reference_id=row['id'], reference_sample=row['sample'],
                     source_root=str(index.root_for(row)), manifest_path=row['manifest_path'], source_sha256=row['source_sha'],
                     same_sample=row['sample'] == query['Sample'], same_source=row['source_sha'] == inputs['CramPath']['sha256'],
                     reference_full_depth=row['full_mean_depth'], requested_depth=depth if arm == 'low' else row['full_mean_depth'],
                     rng_sha256=token, status='ERROR')
        if distribution:
            groups[row['id']].append((row, entry, generator))
            draws.append(entry)
            continue
        try:
            if arm == 'low' and row['full_mean_depth'] < depth:
                entry.update(status='INSUFFICIENT_REFERENCE_DEPTH', error='No upsampling, substitution or redraw')
                draws.append(entry)
                continue
            if arm == 'high':
                vector, covered, record, stamps = reference_coverage(cache, row, layout, edges, timing, config['policies'])
                covered_matrix[number] = covered
                matrix[number] = vector
                units = record['units']
            else:
                ref, other, record, _ = cache.load(row, timing)
                ref = runtime_pool(ref, record['policies'], config['policies'])
                merge_references(layout, other)
                stamps = cache.stamps(row)
                with timing.measure('reference_sample_and_windows'):
                    vector, units, selected_sha = sample_reference_windows(ref, generator, pool.total_bases, edges, audit=entry, covered_output=covered_matrix[number])
                    matrix[number] = vector
                entry['selected_ids_sha256'] = selected_sha
            if not compatible_policies(record['policies'], config['policies']):
                raise ValueError('Reference/query read policy mismatch')
            verified[row['id']] = (row, stamps)
            entry.update(status='READY', actual_depth=int(matrix[number].sum())/layout.length, units=units,
                         **similarity(test, matrix[number], covered_x=query_covered, covered_y=covered_matrix[number], lengths=np.diff(edges), threshold=threshold))
        except (ValueError, OSError) as error:
            entry.update(status='ERROR', error=str(error))
        draws.append(entry)
    if distribution:
        verified = fill_distribution_draws(
            groups, matrix, test, layout, edges, cache, timing, arm=arm, target_bases=pool.total_bases,
            policies=config['policies'], threads=config['threads'], covered_matrix=covered_matrix,
            query_covered=query_covered, threshold=threshold)
    with timing.measure('statistics'):
        if (arm == 'low' or distribution) and all(d['status'] == 'READY' for d in draws):
            nulls = reference_row_means(matrix, config['threads'], covered=covered_matrix, lengths=np.diff(edges), threshold=threshold) if distribution else low_null(matrix, covered=covered_matrix, lengths=np.diff(edges), threshold=threshold)
            records = [dict(metric=k, left_draw=i, right_draw='', value=v, rng_sha256='')
                       for k, values in nulls.items() for i, v in enumerate(values)]
            if distribution:
                for record in records:
                    draw = draws[record['left_draw']]
                    record.update(left_reference_sample=draw['reference_sample'], left_reference_id=draw['reference_id'])
        elif arm == 'high' and not distribution:
            nulls, records = high_null(matrix, draws, config, query, covered=covered_matrix, lengths=np.diff(edges), threshold=threshold)
        else:
            nulls, records = {}, []
        results = aggregate(draws, nulls, 'low' if distribution else arm, requested, stable_tail=distribution)
    unique = len({d['reference_id'] for d in draws})
    summary.update(valid_draws=sum(d['status'] == 'READY' for d in draws),
                   unique_references=unique, duplicate_draws=requested-unique)
    if arm == 'high' and not distribution:
        summary.update(high_null_candidate_pairs=pair_population_size(draws),
                       high_null_sampled_pairs=len(records)//2)
    for metric, result in results.items():
        if distribution:
            p = result['p']
            summary[metric+'_decision'] = 'NOT_EVALUATED' if p is None else 'PASS' if p >= alpha else 'FAIL'
            if result['valid_draws'] != requested:
                result.update(value=None, mean=None)
        for source, suffix in (('value', ''), ('mean', '_mean'), ('max', '_max'), ('p', '_p'),
                               ('status', '_status'), ('null_count', '_null_count'), ('null_mean', '_null_mean'), ('null_sd', '_null_sd')):
            summary[metric+suffix] = result[source]
    summary['status'] = 'COMPLETE' if summary['valid_draws'] == requested else 'COMPLETE_WITH_GAPS'
    if any(d['status'] == 'ERROR' for d in draws):
        summary['status'] = 'ERROR'
    if distribution:
        decisions = {summary['jaccard_decision'], summary['pearson_decision']}
        if 'NOT_EVALUATED' in decisions:
            if summary['status'] == 'COMPLETE':
                summary['status'] = 'COMPLETE_WITH_GAPS'
        else:
            summary['detection_status'] = 'DETECTED' if decisions == {'PASS'} else 'NOT_DETECTED'
    elif arm == 'high' and config['high_null'] == 'none':
        for metric in results:
            if results[metric]['status'] == 'insufficient_null':
                summary[metric+'_status'] = 'HIGH_NULL_NOT_DEFINED'
    from bin.jaccard_sensitivity import write_sensitivity
    sensitivity = write_sensitivity(output, summary, draws, test, matrix, query_covered,
                                    covered_matrix, edges, config, query, nulls)
    summary.update(sensitivity[0])
    np.save(output/'reference_covered_bases.npy', covered_matrix, allow_pickle=False)
    audit['jaccard_breadth'] = dict(threshold=threshold, thresholds=thresholds(threshold, config.get('jaccard_coverage_thresholds')),
        definition='unique covered bases / actual window length; tau0 is C>0',
        covered_sources='query and identical selected reference units',
        sampling_unit=config['policies'].get('sampling_unit','read'),
        coverage_alignment_policy=config['policies'].get('coverage_alignment_policy','primary_only'))
    table(output/'draws.tsv.gz', DRAW_COLUMNS, draws)
    table(output/'null.tsv.gz', NULL_COLUMNS, records)
    np.save(output/'reference_windows.npy', matrix, allow_pickle=False)
    audit.update(reference_selection_rng_sha256=selection_token,
                 simulation='independent complete read units without replacement per draw; no truncation',
                 results=results, window_matrix_bytes=matrix.nbytes, stage_timing=timing.snapshot(),
                 depth_overshoot='whole last unit retained')
    if distribution:
        audit.update(selection=f'{requested} uniform draws with replacement after same Sample or source SHA256 exclusion',
                     reference_distribution=f'per-draw mean against the other {requested-1} draw slots, including repeated sources',
                     query_statistic=f'mean similarity against the {requested} draw slots',
                     test=f'lower-tail normal approximation using population SD of the {requested} reference means',
                     alpha=alpha, detection_rule='both metrics must have p >= alpha',
                     reference_input_stamps=[stamps for _, stamps in verified.values()])
    else:
        audit.update(selection='uniform sample draws with replacement; no depth band or self exclusion',
                     high_null=config['high_null'])
    if arm == 'high' and not distribution:
        audit['high_null_distribution'] = dict(
            candidate_pairs=summary['high_null_candidate_pairs'], sampled_pairs=summary['high_null_sampled_pairs'],
            pair_limit=config['null_limit'], selection='uniform without replacement over unordered draw pairs',
            identity_exclusion='equal reference_sample; repeated draws otherwise retain multiplicity',
            statistic='maximum query-vs-reference similarity', z_variance='population', tail='lower normal CDF')
    for row, stamps in verified.values():
        if stamps != cache.stamps(row):
            raise ValueError('Reference changed during query')
    for item in inputs.values():
        check_identity(item)
    return summary, audit


def finish_query(query, config, index, cache):
    target = Path(config['outdir'])/'queries'/query['query_id']
    if target.exists():
        checks = json.loads((target/'checksums.json').read_text())
        if not checks or any(digest(target/name) != sha for name, sha in checks.items()):
            raise ValueError('Resume output checksum mismatch: '+str(target))
        audit = json.loads((target/'audit.json').read_text())
        for item in audit.get('inputs', {}).values():
            check_identity(item)
        for reference in audit.get('reference_input_stamps', []):
            for path, stamp in reference:
                if file_stat(path) != stamp:
                    raise ValueError('Reference changed since detection query: '+path)
        return json.loads((target/'summary.json').read_text())
    temporary = Path(tempfile.mkdtemp(prefix='.'+query['query_id']+'-', dir=target.parent))
    started = time.perf_counter()
    before = resource.getrusage(resource.RUSAGE_SELF)
    try:
        summary, audit = calculate(query, config, index, cache, temporary)
    except Exception as error:
        summary = dict(query_id=query['query_id'], Sample=query.get('Sample', ''), GCFID=query.get('GCFID', ''),
                       manifest_depth=query.get('Depth', ''), manifest_coverage=query.get('Coverage', ''),
                       status='ERROR', error=str(error))
        if config.get('detection_rules') is not None:
            summary.pop('manifest_coverage')
            summary.update(domain=config['domain'], detection_status='NOT_EVALUATED')
        elif 'alpha' in config:
            summary.update(alpha=config['alpha'], detection_status='NOT_EVALUATED',
                           requested_draws=config.get('reference_draws', REFERENCE_POOL_SIZE))
        audit = dict(traceback=traceback.format_exc())
    after = resource.getrusage(resource.RUSAGE_SELF)
    audit.update(query=query, wall_seconds=time.perf_counter()-started,
                 user_cpu_seconds=after.ru_utime-before.ru_utime, system_cpu_seconds=after.ru_stime-before.ru_stime,
                 process_peak_rss_bytes=int(after.ru_maxrss*(1 if sys.platform == 'darwin' else 1024)), cache=cache.stats())
    if not (temporary/'sensitivity.tsv').exists() and config.get('detection_rules') is None:
        columns = DISTRIBUTION_COLUMNS if 'alpha' in config else SUMMARY_COLUMNS
        table(temporary/'sensitivity.tsv', columns, (dict(summary, jaccard_coverage_threshold=tau)
            for tau in thresholds(config.get('jaccard_coverage_threshold',0.), config.get('jaccard_coverage_thresholds'))))
    if 'alpha' in config:
        # Scores come only from the frozen per-query draw/null files. No CRAM
        # decoding, extra reference selection or changes to the normal call.
        candidates = evaluate_saved(temporary, summary, methods=config['jaccard_methods'],
            alpha=config['alpha'], bootstrap_replicates=config['bootstrap_replicates'],
            bootstrap_confidence=config['bootstrap_confidence'], seed=config['seed'])
        table(temporary/'decision_methods.tsv', METHOD_COLUMNS, candidates)
    atomic_json(temporary/'audit.json', audit)
    atomic_json(temporary/'summary.json', summary)
    atomic_json(temporary/'checksums.json', {str(p.relative_to(temporary)): digest(p) for p in sorted(temporary.rglob('*')) if p.is_file()})
    temporary.rename(target)
    return summary


def run_group(task):
    queries, config = task
    index = ReferenceIndex(config['index'])
    index.connection.execute('PRAGMA query_only=ON')
    index.connection.execute('BEGIN')
    cache = NativeCache(index, config['cache_bytes'])
    try:
        return [finish_query(query, config, index, cache) for query in queries]
    finally:
        cache.close()
        index.close()


def bounded_groups(tasks, jobs):
    if jobs == 1:
        for task in tasks:
            yield run_group(task)
        return
    # Spawn so no open SQLite connection or BLAS pool is inherited.
    import multiprocessing
    with ProcessPoolExecutor(jobs, mp_context=multiprocessing.get_context('spawn')) as executor:
        iterator = iter(tasks)
        pending = {executor.submit(run_group, task) for task in list_next(iterator, jobs)}
        while pending:
            done, pending = wait(pending, return_when=FIRST_COMPLETED)
            for future in done:
                yield future.result()
                task = next(iterator, None)
                if task is not None:
                    pending.add(executor.submit(run_group, task))


def list_next(iterator, size):
    for _ in range(size):
        value = next(iterator, None)
        if value is None:
            break
        yield value


def parser(*, reference_distribution=False):
    description = ('Batch detection using configurable-reference Jaccard and Pearson similarity distributions.'
                   if reference_distribution else __doc__)
    result = argparse.ArgumentParser(description=description)
    if reference_distribution:
        result.add_argument('--version', action='version', version='2.3.0')
    result.add_argument('--manifest', type=Path, required=True, help='TSV: Sample GCFID Depth CramPath; optional Coverage CraiPath BedGraphPath')
    result.add_argument('--index', type=Path, required=True)
    result.add_argument('--coverage-unit', choices=('fraction', 'percent'), help='Required when optional Coverage values are nonempty; never inferred')
    result.add_argument('--sampling-unit', choices=('read', 'fragment'))
    result.add_argument('--coverage-alignment-policy', choices=('primary_only', 'primary_plus_supplementary'))
    result.add_argument('--jaccard-coverage-threshold', type=float, default=0., help='Covered-bases/window-length; 0 means any covered base')
    result.add_argument('--jaccard-coverage-thresholds', help='Sensitivity grid: default or comma-separated fractions; reuses identical draws')
    result.add_argument('--jaccard-methods', default='normal',
                        help='Comma-separated normal,empirical,minimum,bootstrap; normal remains the official detection rule')
    result.add_argument('--bootstrap-replicates', type=int, default=1000)
    result.add_argument('--bootstrap-confidence', type=float, default=.95)
    result.add_argument('--outdir', type=Path, required=True)
    group = result.add_mutually_exclusive_group(required=True)
    group.add_argument('--window-size', '--bin-size', type=int,
                       help='Fixed window length in bp (v2: whole GCF)')
    group.add_argument('--bin-count', '--windows-per-contig', dest='bin_count', type=int,
                       help='Bin count (v2: exact number of balanced windows per GCF); --windows-per-contig is deprecated')
    if reference_distribution:
        result.add_argument('--reference-draws', type=int, default=REFERENCE_POOL_SIZE,
                            help='Reference draws with replacement for both arms (minimum 2; default 100)')
        result.add_argument('--alpha', type=float, required=True,
                            help='Lower-tail significance level (0 < alpha < 1); both p-values must be >= alpha')
        result.add_argument('--threads', type=int, default=2,
                            help='Reference preparation, draw and scoring threads per query (minimum 2)')
        result.add_argument('--query-id', help='Process one original manifest row, e.g. q00000001; requires --jobs 1')
    else:
        group.add_argument('--detection-config', type=Path, help='JSON genome-length rules using Bin Count grid distances')
        result.add_argument('--domain', choices=DOMAINS, help='Required with --detection-config; selects a top-level domain')
        result.add_argument('--simulations', type=int, default=100)
        result.add_argument('--reference-draws', default='100', help='High arm: count or percent, e.g. 100 or 20%%; WITH replacement')
        result.add_argument('--high-null', choices=('none', 'draw-pairs'), default='draw-pairs',
                            help='Default draw-pairs: user-defined unordered distinct-Sample draw-pair null; none explicitly requests descriptive-only replay')
        result.add_argument('--null-limit', type=int, default=1000, help='Maximum retained high-null pairs, 1..1000')
    result.add_argument('--min-samples', type=int, default=10, help='Strict available count > this value')
    result.add_argument('--depth-threshold', type=float, default=5, help='Strict manifest Depth > this value routes high; equality routes low')
    result.add_argument('--seed', type=int, default=401)
    result.add_argument('--jobs', type=int, default=1)
    result.add_argument('--cache-mib', type=int, default=131072, help='Total accounted cache divided among workers; 0 = unlimited; not an RSS limit')
    result.add_argument('--work-mib', type=int, default=131072, help='Total numeric window work budget divided among workers; 0 = unlimited')
    result.add_argument('--resume', action='store_true', help='Validate existing run fingerprint, input identities and query checksums')
    return result


def run(args):
    thresholds(args.jaccard_coverage_threshold, args.jaccard_coverage_thresholds)
    method_names = parse_methods(args.jaccard_methods)
    validate_options(getattr(args, 'alpha', None) or .05,
                     args.bootstrap_replicates, args.bootstrap_confidence)
    if getattr(args, 'alpha', None) is None and method_names != ('normal',):
        raise ValueError('Alternative Jaccard decision methods require --alpha / v2')
    started = time.perf_counter()
    before = resource.getrusage(resource.RUSAGE_SELF)
    child_before = resource.getrusage(resource.RUSAGE_CHILDREN)
    if (args.jobs < 1 or min(args.cache_mib, args.work_mib, args.min_samples, args.seed) < 0
            or not math.isfinite(args.depth_threshold) or args.depth_threshold <= 0):
        raise ValueError('Invalid numeric configuration')
    alpha = getattr(args, 'alpha', None)
    if alpha is None:
        if args.simulations < 1 or not 1 <= args.null_limit <= 1000:
            raise ValueError('Invalid numeric configuration')
        draw_count(args.reference_draws, 1)
    elif not 0 < alpha < 1:
        raise ValueError('--alpha must satisfy 0 < alpha < 1')
    elif type(args.reference_draws) is not int or args.reference_draws < 2:
        raise ValueError('--reference-draws must be an integer at least 2')
    elif args.threads < 2:
        raise ValueError('--threads must be at least 2')
    query_id = getattr(args, 'query_id', None)
    if query_id is not None and args.jobs != 1:
        raise ValueError('--query-id requires --jobs 1')
    detection_path = getattr(args, 'detection_config', None)
    domain = getattr(args, 'domain', None)
    if (detection_path is None) != (domain is None):
        raise ValueError('--detection-config and --domain must be supplied together')
    if detection_path is None and (args.window_size or args.bin_count or 0) <= 0:
        raise ValueError('Window parameter must be positive')
    detection_input = detection_rules = None
    if detection_path is not None:
        detection_input = identity(detection_path)
        detection_rules = read_rules(detection_path, domain)
        check_identity(detection_input)
        try:
            import scipy.spatial.distance
        except ImportError as error:
            raise ValueError('Detection mode requires SciPy; use the CoverageMap Bin Count runtime') from error
    rows = read_manifest(args.manifest, args.coverage_unit, query_id)
    index = ReferenceIndex(args.index)
    try:
        if index.metadata.get('storage_mode') != 'native':
            raise ValueError('Only existing native reference indexes are supported')
        if alpha is not None:
            require_v2_index(index)
        if detection_rules is not None and index.window_policy != PER_CONTIG:
            raise ValueError('Frozen grid-distance rules require a legacy per-contig index')
        policies = dict(index.metadata['config']['policies'])
        if getattr(args, 'sampling_unit', None) is not None or index.metadata['config'].get('capture_alignment'):
            policies['sampling_unit'] = getattr(args, 'sampling_unit', None) or 'read'
        if getattr(args, 'coverage_alignment_policy', None) is not None or index.metadata['config'].get('capture_alignment'):
            policies['coverage_alignment_policy'] = getattr(args, 'coverage_alignment_policy', None) or 'primary_only'
        capabilities = index.metadata.get('geometry_capabilities', dict(
            sampling_units=[index.metadata['config']['policies'].get('sampling_unit','read')], supplementary_available=False))
        if policies.get('sampling_unit','read') not in capabilities['sampling_units']:
            raise ValueError('Legacy index lacks requested runtime sampling capability; rebuild complete geometry')
        if policies.get('coverage_alignment_policy','primary_only') != 'primary_only' and not capabilities['supplementary_available']:
            raise ValueError('Legacy index has no supplementary alignment data')
        metadata = index.metadata
        indexed_domain = metadata.get('taxonomy_partition', {}).get('domain')
        if domain is not None and indexed_domain is not None and indexed_domain != domain:
            raise ValueError(f'Index domain is {indexed_domain}, not {domain}')
    finally:
        index.close()
    config = {key: str(value.resolve()) if isinstance(value, Path) else value for key, value in vars(args).items() if key != 'resume'}
    # None is the internal unlimited sentinel; positive budgets still divide by jobs.
    # Keep a divided zero-byte positive budget bounded (distinct from CLI zero).
    config.update(policies=policies,
                  cache_bytes=None if args.cache_mib == 0 else args.cache_mib*1024**2//args.jobs,
                  work_bytes=None if args.work_mib == 0 else args.work_mib*1024**2//args.jobs)
    # Worker count/cache size do not change results but remain part of resume identity.
    provenance = dict(config=config, manifest=identity(args.manifest), database=identity(args.index), code=code_identity(),
                      schema='query-manifest-v1', database_metadata=metadata)
    if detection_input is not None:
        config['detection_rules'] = detection_rules
        provenance.update(schema='query-manifest-grid-detection-v1', detection_config=detection_input)
    elif alpha is not None:
        provenance.update(schema='query-manifest-reference-distribution-v2', reference_pool_size=args.reference_draws,
                          test='lower-tail normal approximation', statistics_calibrated=False)
    fingerprint = hashlib.sha256(json.dumps(provenance, sort_keys=True).encode()).hexdigest()
    outdir = args.outdir.resolve()
    if args.resume:
        if json.loads((outdir/'run.json').read_text())['fingerprint'] != fingerprint:
            raise ValueError('Resume configuration/source/input/database identity differs')
    else:
        outdir.mkdir(parents=True, exist_ok=False)
        (outdir/'queries').mkdir()
        atomic_json(outdir/'run.json', dict(provenance, fingerprint=fingerprint))
    lock = outdir/'.running.lock'
    descriptor = os.open(lock, os.O_CREAT | os.O_EXCL | os.O_WRONLY, 0o600)
    os.close(descriptor)
    try:
        groups = defaultdict(list)
        for row in rows:
            groups[row.get('GCFID') or ''].append(row)
        summaries = []
        for group in bounded_groups(((queries, config) for _, queries in sorted(groups.items())), args.jobs):
            summaries.extend(group)
            print(f'Completed {len(summaries)}/{len(rows)} queries', file=sys.stderr, flush=True)
        summaries.sort(key=lambda row: row['query_id'])
        check_identity(provenance['database'])
        check_identity(provenance['manifest'])
        if detection_input is not None:
            check_identity(detection_input)
        temporary = outdir/'summary.tsv.tmp'
        if alpha is not None:
            columns = DISTRIBUTION_COLUMNS
        else:
            columns = DETECTION_COLUMNS if detection_input is not None else SUMMARY_COLUMNS
        table(temporary, columns, summaries)
        os.replace(temporary, outdir/'summary.tsv')
        if detection_input is None:
            def sensitivity_rows():
                for row in summaries:
                    with (outdir/'queries'/row['query_id']/'sensitivity.tsv').open() as handle:
                        yield from csv.DictReader(handle, delimiter='\t')
            table(outdir/'sensitivity.tsv', columns, sensitivity_rows())
        if alpha is not None:
            def decision_rows():
                for row in summaries:
                    with (outdir/'queries'/row['query_id']/'decision_methods.tsv').open() as handle:
                        yield from csv.DictReader(handle, delimiter='\t')
            table(outdir/'decision_methods.tsv', METHOD_COLUMNS, decision_rows())
        after = resource.getrusage(resource.RUSAGE_SELF)
        children = resource.getrusage(resource.RUSAGE_CHILDREN)
        state = 'FAILED' if any(r['status'] == 'ERROR' for r in summaries) else 'COMPLETE_WITH_GAPS' if any(r['status'] != 'COMPLETE' for r in summaries) else 'COMPLETE'
        report = dict(state=state, queries=len(rows), groups=len(groups), jobs=args.jobs,
                      wall_seconds=time.perf_counter()-started,
                      user_cpu_seconds=after.ru_utime-before.ru_utime+children.ru_utime-child_before.ru_utime,
                      system_cpu_seconds=after.ru_stime-before.ru_stime+children.ru_stime-child_before.ru_stime,
                      controller_peak_rss_bytes=int(after.ru_maxrss*(1 if sys.platform == 'darwin' else 1024)),
                      output_bytes=sum(p.stat().st_size for p in outdir.rglob('*') if p.is_file()),
                      summary_sha256=digest(outdir/'summary.tsv'), statistics_calibrated=False,
                      python=sys.version, numpy=np.__version__)
        if alpha is not None:
            report.update(threads=args.threads, requested_cpu_slots=args.jobs * args.threads)
        atomic_json(outdir/'completion.json', report)
        return 1 if state == 'FAILED' else 0
    finally:
        lock.unlink()


def main(*, reference_distribution=False):
    cli = parser(reference_distribution=reference_distribution)
    args = cli.parse_args()
    if any(word == '--windows-per-contig' or word.startswith('--windows-per-contig=') for word in sys.argv[1:]):
        print('WARNING: --windows-per-contig is deprecated; use --bin-count.', file=sys.stderr)
    try:
        status = run(args)
    except (ValueError, OSError, KeyError) as error:
        cli.exit(2, f'ERROR: {error}\n')
    raise SystemExit(status)


if __name__ == '__main__':
    main()
