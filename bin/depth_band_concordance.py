"""Query-depth reference pools and multi-layer low-depth simulations."""
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import re
import resource
import sys
import time
import traceback

import numpy as np
from Tools import jaccard_similarity, pearson_corr
from bin.coverage_core import compatible_policies, runtime_pool, merge_references, select_units, window_sums
from bin.concordance_core import task_rng
from bin.gcf_windows import GCF_CONCAT, edges_for
from bin.ref_depth_catalog import atomic_json, digest, file_stat, sample_folder
from bin.stage_timing import StageTiming


def selection_size(value):
    """Parse before any output is created; ceil percentages preserve a nonempty draw."""
    if value == 'all':
        return ('all', None)
    if re.fullmatch(r'[1-9][0-9]*', value):
        return ('count', int(value))
    if value.endswith('%'):
        percent = float(value[:-1])
        if math.isfinite(percent) and 0 < percent <= 100:
            return ('percent', percent)
    raise ValueError('simulation-references must be all, a positive integer, or (0,100]%')


def configure(args, depths):
    policy = getattr(args, 'comparison_policy', 'depth-band')
    if policy == 'nearest':
        if args.simulation_depths or args.simulation_references != 'all' or args.simulation_source_depth_range or args.simulation_match_depth_band:
            raise ValueError('Simulation layer/source controls require --comparison-policy depth-band')
        return dict(comparison_policy=policy)
    if depths:
        raise ValueError('Native depth-band uses --simulation-depths for low-depth layers; --depth/--depths/--depth-grid require --comparison-policy nearest')
    size, count = args.window_size, getattr(args, 'bin_count', None)
    if count is None and size is None:
        size = 1000
    if (size if count is None else count) < 1 or not math.isfinite(args.depth_tolerance) or args.depth_tolerance < 0:
        raise ValueError('window-size must be positive; depth-tolerance must be finite and nonnegative')
    if not math.isfinite(args.low_depth_threshold) or args.low_depth_threshold <= 0:
        raise ValueError('low-depth-threshold must be finite and positive')
    if args.replicates < 2 or args.seed < 0 or args.depth_sampler not in ('permutation', 'lazy'):
        raise ValueError('Native queries require N>=2, seed>=0 and an explicit permutation/lazy sampler')
    layers = [float(v) for v in args.simulation_depths.split(',')] if args.simulation_depths else []
    if len(layers) != len(set(layers)) or any(not math.isfinite(v) or v <= 0 for v in layers):
        raise ValueError('simulation-depths must be finite, positive and unique')
    bounds = None
    if args.simulation_source_depth_range:
        bounds = list(map(float, args.simulation_source_depth_range.split(':')))
        if len(bounds) != 2 or any(not math.isfinite(v) or v < 0 for v in bounds) or bounds[0] > bounds[1]:
            raise ValueError('simulation-source-depth-range must be finite 0<=MIN<=MAX')
    selection_size(args.simulation_references)
    return dict(comparison_policy=policy, include_self=True if args.include_self is None else args.include_self,
                depth_metric=args.depth_metric, window_size=size, bin_count=count,
                depth_tolerance=args.depth_tolerance, low_depth_threshold=args.low_depth_threshold,
                simulation_depths=layers, simulation_references=args.simulation_references,
                simulation_source_depth_range=bounds, simulation_match_depth_band=args.simulation_match_depth_band)


def select_references(rows, specification, seed, gcf):
    rows = sorted(rows, key=lambda row: row['sample'])
    kind, value = selection_size(specification)
    count = len(rows) if kind == 'all' else min(len(rows), value) if kind == 'count' else math.ceil(len(rows)*value/100)
    if count == len(rows):
        return rows
    words = np.frombuffer(hashlib.sha256(gcf.encode()).digest()[:16], dtype='<u4').tolist()
    rng = np.random.default_rng(np.random.SeedSequence([seed, *words, 2]))
    return [rows[i] for i in sorted(rng.choice(len(rows), size=count, replace=False))]


def high_coverage(test, references, alpha):
    """Tools.py high-coverage maximum; reference distribution is unordered pairs.

    Keep only O(R) scratch rather than a full R-by-R similarity matrix.
    Normal moments are sufficient for the original normal z calculation.
    """
    results = {}
    for metric, function in (('jaccard', jaccard_similarity), ('pearson', pearson_corr)):
        observed = np.atleast_1d(function(test, references))
        best = float(observed.max()) if np.all(np.isfinite(observed)) else None
        count = 0
        mean = m2 = 0.
        finite = True
        for i in range(len(references)-1):
            values = np.atleast_1d(function(references[i], references[i+1:]))
            if np.any(~np.isfinite(values)):
                finite = False
                continue
            # Combine population moments; avoids cancellation in E[x^2]-E[x]^2.
            n = len(values)
            delta = float(values.mean())-mean
            total = count+n
            m2 += float(np.square(values-values.mean()).sum())+delta*delta*count*n/total
            mean += delta*n/total
            count = total
        sd = math.sqrt(max(0., m2/count)) if count else None
        entry = dict(similarity=best, reference_mean=mean if count else None, reference_sd=sd,
                     reference_pairs=count, p_value=None, accepted=None)
        if best is None or not finite or not count:
            entry.update(status='non_estimable', reason='insufficient_or_constant_profiles')
        elif sd == 0:
            entry.update(status='non_estimable', reason='zero_variance_reference_distribution')
        else:
            from statistics import NormalDist
            p = NormalDist().cdf((best-mean)/sd)
            entry.update(status='estimable', p_value=p, accepted=bool(p >= alpha))
        results[metric] = entry
    return dict(results, model='tools-high-coverage-maximum-normal-not-calibrated', alpha=alpha,
                reference_distribution='all-unordered-distinct-reference-pair-similarities',
                reference_count=len(references), simulations=0, windows=len(test))


def result_row(base, result, output):
    row = dict(base, status='READY', result_path=str(output/'concordance.json'))
    for name in ('jaccard', 'pearson'):
        row[name] = result[name]['similarity']
        row[name+'_status'] = result[name]['status']
        row[name+'_p'] = result[name]['p_value']
        row[name+'_accepted'] = result[name]['accepted']
    return row


def run_comparisons(query, config, index, cache, directory, timing):
    """Extract once, choose sources before layers, and never upsample missing mass."""
    from bin import native_concordance as native
    group = index.group(query['gcf'], missing_ok=True)
    if group is None:
        exists = index.connection.execute('SELECT 1 FROM gcfs WHERE gcf=?',(query['gcf'],)).fetchone()
        reason = 'NO_SELECTED_REFERENCE_GCF' if exists else 'NO_REFERENCE_GCF'
        return [dict(query_sample=query['sample'], gcf=query['gcf'], status=reason)], dict(extraction_calls=0)
    layout_file = directory/'query_layout.fai'
    layout_file.write_text(''.join(f'{n}\t{length}\n' for n, length in zip(group['names'], group['lengths'])))
    with timing.measure('query_sha256'):
        source_stat = file_stat(query['alignment'])
        source_hash = digest(query['alignment'])
    with timing.measure('query_extract'):
        order = {'sort_contigs': True} if index.window_policy == GCF_CONCAT else {}
        pool, layout = native.extract_pool(query['alignment'], reference_index=layout_file,
                                           **order, **config['policies'])
    # The full FAI, not only covered positions or unweighted window means, is the denominator.
    depth = pool.total_bases/pool.genome_len
    count = config.get('bin_count')
    if count is not None and index.window_policy != GCF_CONCAT:
        raise ValueError('Depth-band bin-count requires a gcf-concat-v1 index')
    edges = edges_for(layout, index.window_policy, size=config['window_size'], count=count)
    widths = np.diff(edges)
    with timing.measure('query_full_coverage'):
        full = window_sums(pool, np.arange(len(pool)), edges)
    np.save(directory/'edges.npy', edges, allow_pickle=False)
    np.save(directory/'query_mean_depth.npy', full/widths, allow_pickle=False)
    low = depth < config['low_depth_threshold']
    minimum, maximum = (0., None) if low else (max(0., depth-config['depth_tolerance']), depth+config['depth_tolerance'])
    if low and config['simulation_source_depth_range']:
        minimum, maximum = config['simulation_source_depth_range']
    if low and config['simulation_match_depth_band']:
        minimum = max(minimum, depth-config['depth_tolerance'])
        maximum = min(maximum, depth+config['depth_tolerance']) if maximum is not None else depth+config['depth_tolerance']
    with timing.measure('reference_sql'):
        candidates = index.candidates(query['gcf'], minimum, maximum,
            exclude_sample=None if config['include_self'] else query['sample'],
            exclude_source=None if config['include_self'] else source_hash, samples=config['reference_samples'])
        band_count = len(candidates)
        selection_mode = 'simulation_sources' if low else 'depth_band'
        if not low and not candidates:
            available = index.candidates(query['gcf'],
                exclude_sample=None if config['include_self'] else query['sample'],
                exclude_source=None if config['include_self'] else source_hash, samples=config['reference_samples'])
            if available:
                difference = min(abs(row['full_mean_depth']-depth) for row in available)
                candidates = [row for row in available if abs(row['full_mean_depth']-depth) == difference]
                selection_mode = 'nearest_depth_fallback'
    selected = select_references(candidates, config['simulation_references'], config['seed'], query['gcf']) if low else candidates
    audit = dict(query_full_depth=depth, depth_metric=config['depth_metric'], window_size=config['window_size'],
                 bin_count=count, window_policy=index.window_policy,
                 route='low' if low else 'high', candidate_count=len(candidates), selected_samples=[r['sample'] for r in selected],
                 band_candidate_count=None if low else band_count, selection_mode=selection_mode,
                 depth_band_matched=None if low else selection_mode == 'depth_band', include_self=config['include_self'],
                 selected_references=[dict(sample=r['sample'], full_mean_depth=r['full_mean_depth'], source_sha=r['source_sha'],
                     manifest=str(index.root_for(r)/r['manifest_path'])) for r in selected],
                 source_depth_bounds=[minimum, maximum], extraction_calls=1, query_source_sha256=source_hash,
                 comparison_scope='self-consistency-included' if any(r['sample'] == query['sample'] or r['source_sha'] == source_hash for r in selected) else 'separate-input-pools',
                 test_missing_header_contigs=list(layout.missing_header_contigs),
                 test_unknown_m5_contigs=[n for n, m5 in zip(layout.names, layout.md5) if not m5])
    atomic_json(directory/'selection.json', audit)
    base = dict(query_sample=query['sample'], gcf=query['gcf'], query_full_depth=depth,
                depth_metric=config['depth_metric'], window_size=config['window_size'], candidate_count=len(candidates),
                reference_count=len(selected), comparison_scope=audit['comparison_scope'])
    base.update(selection_mode=selection_mode, depth_band_matched=audit['depth_band_matched'])
    verified = {}
    if depth <= 0:
        rows = [dict(base, status='ZERO_QUERY_DEPTH')]
    elif not selected:
        rows = [dict(base, status='NO_ELIGIBLE_REFERENCE' if low else 'NO_USABLE_REFERENCE')]
    elif low:
        rows = low_comparisons(base, selected, pool, layout, edges, full, config, index, cache, directory, timing, verified)
    else:
        vectors = []
        for row in selected:
            ref, ref_layout, record, _ = cache.load(row, timing)
            ref = runtime_pool(ref, record['policies'], config['policies'])
            verified[row['id']] = (row, cache.stamps(row))
            merge_references(ref_layout, layout)
            if not compatible_policies(record['policies'], config['policies']):
                raise ValueError('Reference/query policy mismatch')
            with timing.measure('reference_full_coverage'):
                vectors.append(window_sums(ref, np.arange(len(ref)), edges)/widths)
        with timing.measure('statistics'):
            result = high_coverage(full/widths, np.asarray(vectors), config['alpha'])
        if selection_mode == 'nearest_depth_fallback':
            for metric in ('jaccard','pearson'):
                result[metric].update(status='descriptive_only' if result[metric]['similarity'] is not None else 'non_estimable',
                    reason='outside_requested_depth_band',p_value=None,accepted=None)
        output = directory/'high-depth'
        output.mkdir()
        entry = dict(base, depth=depth, requested_depth=depth, comparison_depth=depth,
                     test_actual_depth=depth, comparison_mode='high_depth_reference_pool_maximum')
        result.update(**entry, reference_samples=audit['selected_samples'], statistics_calibrated=False)
        atomic_json(output/'concordance.json', result)
        rows = [result_row(entry, result, output)]
    with timing.measure('query_sha256'):
        if file_stat(query['alignment']) != source_stat or digest(query['alignment']) != source_hash:
            raise ValueError('Query changed during calculation')
    with timing.measure('reference_metadata'):
        for row, stamps in verified.values():
            if stamps != cache.stamps(row):
                raise ValueError('Reference changed during calculation')
    return rows, audit


def low_comparisons(base, selected, pool, layout, edges, full, config, index, cache, directory, timing, verified):
    from bin.native_concordance import compare
    rows = []
    widths = np.diff(edges)
    layers = config['simulation_depths'] or [base['query_full_depth']]
    for number, target in enumerate(layers):
        entry = dict(base, depth=target, requested_depth=target, comparison_depth=target)
        if target > base['query_full_depth']:
            rows.append(dict(entry, status='UNSUPPORTED_QUERY_DEPTH'))
            continue
        with timing.measure('query_sample_count'):
            ids = np.arange(len(pool)) if target == base['query_full_depth'] else select_units(
                pool, task_rng(config['seed'], target, 1, 0), target_bases=target*pool.genome_len, depth_sampler=config['sampler'])
            test = full if target == base['query_full_depth'] else window_sums(pool, ids, edges)
            if int(test.sum()) != int(pool.unit_bases[ids].sum()):
                raise AssertionError('Query mass mismatch')
        layer_dir = directory/f'layer-{number:05d}'
        layer_dir.mkdir()
        np.save(layer_dir/'test.npy', test, allow_pickle=False)
        np.save(layer_dir/'test_ids.npy', ids, allow_pickle=False)
        for row in selected:
            item = dict(entry, reference_sample=row['sample'], reference_full_depth=row['full_mean_depth'],
                        sample_manifest=str(index.root_for(row)/row['manifest_path']))
            if row['full_mean_depth'] < target:
                rows.append(dict(item, status='UNSUPPORTED_REFERENCE_DEPTH'))
                continue
            ref, ref_layout, record, _ = cache.load(row, timing)
            ref = runtime_pool(ref, record['policies'], config['policies'])
            current = cache.stamps(row)
            if row['id'] in verified and verified[row['id']][1] != current:
                raise ValueError('Reference changed between simulation layers')
            verified[row['id']] = (row, current)
            merge_references(ref_layout, layout)
            if not compatible_policies(record['policies'], config['policies']):
                raise ValueError('Reference/query policy mismatch')
            matrix, audits = cache.simulate(row, ref, edges, target, config, timing)
            with timing.measure('statistics'):
                result = compare(test, matrix, widths, config)
            output = layer_dir/sample_folder(row['sample'])
            output.mkdir()
            item.update(comparison_mode='low_depth_multilayer_simulation', test_actual_depth=int(test.sum())/pool.genome_len,
                        reference_actual_depth=float(np.mean([a['actual_depth'] for a in audits])),
                        reference_depth_difference=abs(row['full_mean_depth']-target), reference_depth_ratio=row['full_mean_depth']/target)
            result.update(**item, reference_replicates=audits, seed=config['seed'], sampler=config['sampler'],
                          test_selected_ids_sha256=hashlib.sha256(ids.astype('<i8').tobytes()).hexdigest(), statistics_calibrated=False,
                          reference_missing_header_contigs=record.get('missing_header_contigs',[]),
                          test_missing_header_contigs=list(layout.missing_header_contigs),
                          identity_unverified_contigs=[n for n,a,b in zip(layout.names,layout.md5,ref_layout.md5) if not a or not b])
            atomic_json(output/'concordance.json', result)
            rows.append(result_row(item, result, output))
    return rows


def perform_depth_band(query, config, index, cache):
    from bin.concordance_batch import COLUMNS
    directory = Path(config['outdir'])/'queries'/query['gcf']/sample_folder(query['sample'])
    status = directory/'status'
    status.mkdir(parents=True)
    (status/'RUNNING').touch()
    timing = StageTiming()
    started = time.perf_counter()
    before = resource.getrusage(resource.RUSAGE_SELF)
    old = cache.stats()
    audit = {}
    try:
        rows, audit = run_comparisons(query, config, index, cache, directory, timing)
    except BaseException as error:
        if isinstance(error, (KeyboardInterrupt, SystemExit)):
            (status/'RUNNING').unlink(missing_ok=True)
            (status/'FAILED').touch()
            (status/'exit_code.txt').write_text('1\n')
            raise
        rows = [dict(query_sample=query['sample'], gcf=query['gcf'], status='ERROR', error=str(error))]
        (directory/'error.txt').write_text(traceback.format_exc())
    ready = sum(row['status'] == 'READY' for row in rows)
    failed = sum(row['status'] == 'ERROR' for row in rows)
    gaps = len(rows)-ready-failed
    missing = len(rows) == 1 and rows[0]['status'] == 'NO_REFERENCE_GCF'
    state = 'FAILED' if failed else 'SKIPPED' if missing else 'COMPLETE_WITH_GAPS' if gaps else 'COMPLETE'
    reason = '' if ready else rows[0].get('error', rows[0]['status'])
    with (directory/'summary.tsv').open('w', newline='') as handle:
        writer = csv.DictWriter(handle, COLUMNS, delimiter='\t')
        writer.writeheader()
        writer.writerows(rows)
    after = resource.getrusage(resource.RUSAGE_SELF)
    now = cache.stats()
    report = dict(query=query, state=state, reason=reason, **audit, ready_comparisons=ready, failed_comparisons=failed,
                  depth_gaps=gaps, worker_pid=os.getpid(), elapsed_seconds=time.perf_counter()-started,
                  cache={k:now[k]-old[k] for k in ('hits', 'misses', 'evictions')}, cache_after=now,
                  stage_timing=timing.snapshot(), user_cpu_seconds=after.ru_utime-before.ru_utime,
                  system_cpu_seconds=after.ru_stime-before.ru_stime,
                  process_peak_rss_bytes=int(after.ru_maxrss if sys.platform == 'darwin' else after.ru_maxrss*1024),
                  statistics_calibrated=False)
    atomic_json(directory/'provenance.json', report)
    atomic_json(directory/'progress.json', dict(state=state))
    (status/'RUNNING').unlink()
    (status/('FAILED' if failed else 'SKIPPED' if missing else 'SUCCESS')).touch()
    (status/'exit_code.txt').write_text('1\n' if failed else '0\n')
    return dict(query=query, state=state, reason=reason, report=str(directory/'provenance.json'),
                summary=str(directory/'summary.tsv'), exit_code=int(bool(failed)))
