"""Frozen inputs, deterministic stratification and exact geometry matrix construction."""
from collections import Counter
from contextlib import closing, contextmanager, ExitStack
import csv
from functools import lru_cache
import hashlib
import json
import math
from pathlib import Path
import sqlite3
import sys
from threading import Lock
from time import perf_counter

import numpy as np
from bin.bin_count_eval_math import exact_edges, normalized
from bin.coverage_core import layout_metadata, reference_from_layout
from bin.gcf_windows import GCF_CONCAT, concat_edges, require_layout_policy, window_policy

LEGACY_MATRIX_WINDOWS = 'concatenated-contigs-exact-B-remainder-last'


def matrix_policy(plan):
    policy = plan.get('window_policy', LEGACY_MATRIX_WINDOWS)
    if policy not in (LEGACY_MATRIX_WINDOWS, GCF_CONCAT):
        raise ValueError(f'Unsupported matrix window policy: {policy}')
    return policy


def matrix_edges(length, count, policy):
    return concat_edges(length, count=count) if policy == GCF_CONCAT else exact_edges(length, count)


def sha(path):
    value = hashlib.sha256()
    with Path(path).open('rb') as h:
        for block in iter(lambda: h.read(1024**2), b''):
            value.update(block)
    return value.hexdigest()


def write_json(path, value):
    Path(path).write_text(json.dumps(value, ensure_ascii=False, sort_keys=True, indent=2,
                                    allow_nan=False) + '\n', encoding='utf-8')


def write_edges(path, lengths, bins, policy=LEGACY_MATRIX_WINDOWS):
    """Stream the canonical JSON layout without materializing boxed edge lists."""
    if lengths and not 1 <= bins <= min(lengths.values()):
        raise ValueError('Require integer 1 <= bin count <= genome length')
    with path.open('w', encoding='utf-8') as stream:
        if not lengths:
            stream.write('{}\n')
            return
        stream.write('{\n')
        for i, (gcf, length) in enumerate(sorted(lengths.items())):
            if i:
                stream.write(',\n')
            stream.write('  ' + json.dumps(gcf, ensure_ascii=False) + ': [\n')
            edges = matrix_edges(length, bins, policy)
            for start in range(0, bins, 8192):
                stop = min(start + 8192, bins)
                stream.write(''.join(f'    {edge},\n' for edge in edges[start:stop]))
            stream.write(f'    {length}\n  ]')
        stream.write('\n}\n')


def read_tsv(path):
    with Path(path).open(encoding='utf-8', newline='') as h:
        reader = csv.DictReader(h, delimiter='\t')
        if not reader.fieldnames or len(reader.fieldnames) != len(set(reader.fieldnames)):
            raise ValueError('Missing/duplicate TSV headers')
        rows = list(reader)
    if any(None in r or any(v is None for v in r.values()) for r in rows):
        raise ValueError('Malformed TSV row')
    return rows


def write_tsv(path, columns, rows):
    with Path(path).open('w', encoding='utf-8', newline='') as h:
        w = csv.DictWriter(h, columns, delimiter='\t', lineterminator='\n')
        w.writeheader()
        w.writerows(rows)


def child(root, relative):
    relative = Path(relative)
    if relative.is_absolute() or '..' in relative.parts or not relative.parts:
        raise ValueError('Unsafe relative path')
    result = (Path(root) / relative).resolve()
    if not result.is_relative_to(Path(root).resolve()):
        raise ValueError('Path escapes source root')
    return result


def new_directory(path):
    path = Path(path).absolute()
    path.mkdir(parents=True, exist_ok=False)
    return path


def code_hashes():
    root = Path(__file__).parent
    names = ['bin_count_evaluate.py', 'bin_count_eval_math.py', 'bin_count_eval_io.py',
             'bin_count_eval_render.py', 'bin_count_multinode.py', 'coverage_core.py',
             'reference_store.py', 'gcf_windows.py']
    return {n: sha(root / n) for n in names}


def seal(root, data, *, written_sha256=None):
    data.update(code_sha256=code_hashes(), python=sys.version, numpy=np.__version__)
    written_sha256 = written_sha256 or {}
    artifacts = {str(p.relative_to(root)): p for p in sorted(root.rglob('*'))
                 if p.is_file() and p not in (root/'result.json', root/'SUCCESS')}
    if written_sha256.keys() - artifacts.keys():
        raise ValueError('Written checksum does not belong to a published artifact')
    data['artifact_sha256'] = {name: written_sha256[name] if name in written_sha256 else sha(path)
                               for name, path in artifacts.items()}
    write_json(root / 'result.json', data)
    (root / 'SUCCESS').write_text(sha(root / 'result.json') + '\n')
    return data


def read_manifest(root, kind):
    root = Path(root).resolve()
    if (root / 'SUCCESS').read_text().strip() != sha(root / 'result.json'):
        raise ValueError('Incomplete/changed output manifest')
    data = json.loads((root / 'result.json').read_text())
    if data['kind'] != kind:
        raise ValueError('Wrong input artifact kind')
    return root, data


def load(root, kind, *, folders=None):
    root, data = read_manifest(root, kind)
    for name, expected in data['artifact_sha256'].items():
        if folders is not None and str(Path(name).parent) not in folders:
            continue
        if sha(child(root, name)) != expected:
            raise ValueError('Changed artifact: ' + name)
    return root, data


def sizes(explicit=None, start=None, stop=None, step=None):
    if explicit is not None:
        if any(v is not None for v in (start, stop, step)):
            raise ValueError('Bin-size list and range are mutually exclusive')
        result = explicit
    elif any(v is not None for v in (start, stop, step)):
        if any(v is None for v in (start, stop, step)) or step <= 0 or start <= 0 or stop < start:
            raise ValueError('Need valid start/stop/step for bin-size range')
        if (stop-start)//step + 1 > 10000:
            raise ValueError('Bin-size range exceeds 10000 values')
        result = list(range(start, stop+1, step))
    else:
        result = [1, 50, 100]
    if not result or any(type(v) is not int or v <= 0 for v in result):
        raise ValueError('Bin sizes must be positive integers')
    return sorted(set(result))


def read_strata_config(path):
    config = json.loads(Path(path).read_text(encoding='utf-8'))
    if not isinstance(config, dict) or set(config) != {'strata'}:
        raise ValueError('Strata config requires a single "strata" field')
    strata = config['strata']
    if not isinstance(strata, list) or not strata:
        raise ValueError('Strata config requires a nonempty "strata" list')
    fields = {'min_genome_bp', 'max_genome_bp_exclusive', 'nominal_bin_widths_bp'}
    resolved = []
    for index, item in enumerate(strata, 1):
        if not isinstance(item, dict) or set(item) != fields:
            raise ValueError(f'Stratum {index} requires fields: {", ".join(sorted(fields))}')
        lower = item['min_genome_bp']
        upper = item['max_genome_bp_exclusive']
        widths = item['nominal_bin_widths_bp']
        if type(lower) is not int or lower < 0:
            raise ValueError(f'Stratum {index} lower bound must be a nonnegative integer')
        if resolved and resolved[-1]['max_genome_bp_exclusive'] != lower:
            raise ValueError(f'Stratum {index} must start at the previous exclusive upper bound')
        if index == len(strata):
            if upper is not None:
                raise ValueError('Last stratum must have a null upper bound')
        elif type(upper) is not int or upper <= lower:
            raise ValueError(f'Stratum {index} needs an integer upper bound above its lower bound')
        if not isinstance(widths, list):
            raise ValueError(f'Stratum {index} nominal_bin_widths_bp must be a list')
        resolved.append(dict(min_genome_bp=lower, max_genome_bp_exclusive=upper,
                             nominal_bin_widths_bp=sizes(widths)))
    return resolved


def strata_plan(rows, layouts, breaks, bin_sizes, strata_config=None):
    if (not breaks or any(type(x) is not int or x < 0 for x in breaks) or
            breaks != sorted(set(breaks))):
        raise ValueError('Stratum lower bounds must be distinct increasing nonnegative integers')
    strata = []
    for i, lower in enumerate(breaks):
        upper = breaks[i+1] if i+1 < len(breaks) else None
        gcfs = sorted(g for g, layout in layouts.items()
                      if sum(layout['lengths']) >= lower and (upper is None or sum(layout['lengths']) < upper))
        selected = [j for j, r in enumerate(rows) if r['gcf'] in set(gcfs)]
        minimum = min((sum(layouts[g]['lengths']) for g in gcfs), default=None)
        candidates = {}
        widths = bin_sizes if strata_config is None else strata_config[i]['nominal_bin_widths_bp']
        for size in widths if minimum else []:
            b = minimum // size
            if b < 1:
                raise ValueError(f'Stratum starting at {lower} bp: bin size {size} exceeds shortest genome {minimum}')
            candidates.setdefault(b, []).append(size)
        strata.append(dict(id=f'stratum_{i+1:02d}', lower_bp=lower, upper_bp=upper,
                           min_genome_bp=minimum, gcfs=gcfs, row_indices=selected,
                           status='READY' if selected else 'SKIPPED_EMPTY',
                           candidates=[dict(bin_count=b, bin_sizes_bp=aliases,
                                            raw_plus_normalized_bytes=16*len(selected)*b)
                                       for b, aliases in sorted(candidates.items())]))
    if sum(len(s['row_indices']) for s in strata) != len(rows):
        raise ValueError('Some selected genomes fall below the first stratum boundary')
    return strata


def frozen_preflight(directory):
    directory = Path(directory).resolve()
    if (directory / 'SUCCESS').read_text().strip() != 'COMPLETE':
        raise ValueError('Preflight not complete/nonempty')
    ledger = {}
    for line in (directory / 'SHA256SUMS').read_text().splitlines():
        digest, name = line.split('  ', 1)
        if name in ledger or sha(child(directory, name)) != digest:
            raise ValueError('Preflight checksum mismatch or duplicate: ' + name)
        ledger[name] = digest
    required = {'summary.json', 'filtered_records.tsv', 'filtered_gcf_lengths.tsv', 'filtered_samples.txt'}
    if not required <= ledger.keys():
        raise ValueError('Preflight checksum ledger missing required payloads')
    summary = json.loads((directory / 'summary.json').read_text())
    rows = read_tsv(directory / 'filtered_records.tsv')
    gcfs = read_tsv(directory / 'filtered_gcf_lengths.tsv')
    keys = [(r['sample'], r['gcf']) for r in rows]
    if not rows or len(set(keys)) != len(keys) or any(not all(k) for k in keys):
        raise ValueError('Missing/duplicate Sample+GCF row key')
    counts = Counter(r['gcf'] for r in rows)
    if len(gcfs) != len({r['gcf'] for r in gcfs}) or counts != {r['gcf']: int(r['retained_records']) for r in gcfs}:
        raise ValueError('Preflight GCF support counts disagree')
    sample_ids = sorted({r['sample'] for r in rows})
    if sample_ids != (directory / 'filtered_samples.txt').read_text().splitlines():
        raise ValueError('Preflight sample list disagrees')
    c = summary['counts']
    if (len(rows), len(counts), len(sample_ids)) != (c['retained_records'], c['retained_gcfs'], c['retained_unique_samples']):
        raise ValueError('Preflight summary counts disagree')
    filters = summary['filters']
    if min(counts.values()) < filters.get('min_samples_per_gcf', 1):
        raise ValueError('Preflight support filter violated')
    for row in rows:
        if (not float(row['full_mean_depth']) > filters['depth_gt'] or
                not filters['cover_fraction_gt'] < float(row['covered_fraction']) <= 1 or
                not int(row['genome_length_bp']) > filters.get('genome_length_bp_gt', 0)):
            raise ValueError('Preflight record violates frozen filters')
    return directory, summary, rows, ledger


def validate_frozen_index(index, summary):
    index = Path(index).resolve()
    if not index.is_file() or any(Path(str(index)+s).exists() for s in ('-wal', '-journal')):
        raise ValueError('Need finalized source SQLite without WAL/journal')
    fingerprint = sha(index)
    if fingerprint != summary['input_sha256']:
        raise ValueError('SQLite hash differs from frozen preflight source')
    return index, fingerprint


def create_plan(preflight_dir, index, outdir, breaks, bin_sizes, root_map=None, *, strata_config=None):
    preflight_dir, summary, records, ledger = frozen_preflight(preflight_dir)
    index, fingerprint = validate_frozen_index(index, summary)
    con = sqlite3.connect(index.as_uri() + '?mode=ro', uri=True)
    con.row_factory = sqlite3.Row
    root_map = root_map or {}
    rows, layouts = [], {}
    try:
        con.execute('PRAGMA query_only=ON'); con.execute('BEGIN')
        meta = {r['key']: json.loads(r['value']) for r in con.execute('SELECT * FROM metadata')}
        policy = window_policy(meta.get('config', {}))
        if meta['schema'] == 'coverage-reference-sqlite-3':
            roots = {str(r['id']): r['root'] for r in con.execute('SELECT id,root FROM sources')}
        elif meta['schema'] == 'coverage-reference-sqlite-2':
            roots = {'1': meta['database_root']}
        else:
            raise ValueError('Need native reference SQLite v2/v3')
        if not set(root_map) <= roots.keys() or any(not Path(v).is_absolute() for v in root_map.values()):
            raise ValueError('Source-root-map needs known source IDs and absolute paths')
        roots.update(root_map)
        for row_number, record in enumerate(sorted(records, key=lambda r: (r['gcf'], r['sample']))):
            row = con.execute('SELECT * FROM references_table WHERE id=?', (int(record['reference_id']),)).fetchone()
            if row is None:
                raise ValueError('Reference ID absent from SQLite')
            row = dict(row)
            for name, source in [('gcf', 'gcf'), ('sample', 'sample'), ('genome_length_bp', 'genome_len'),
                                 ('full_mean_depth', 'full_mean_depth'), ('covered_fraction', 'covered_fraction')]:
                a, b = record[name], row[source]
                if (a != b if name in ('gcf', 'sample') else float(a) != b):
                    raise ValueError('Frozen record/index identity mismatch: ' + name)
            sid = str(row.get('source_id', 1))
            if record.get('source_id') and record['source_id'] != sid:
                raise ValueError('Frozen source ID mismatch')
            for field in ('manifest_path', 'geometry_path'):
                if not row[field] or record.get(field) != row[field]:
                    raise ValueError('Frozen/index path mismatch')
            root = Path(roots[sid])
            if not root.is_absolute():
                raise ValueError('Source root must be absolute')
            child(root, row['geometry_path'])
            manifest = child(root, row['manifest_path'])
            if manifest.parent / 'geometry' != child(root, row['geometry_path']):
                raise ValueError('Unexpected geometry/result relationship')
            gcf = row['gcf']
            if gcf not in layouts:
                item = con.execute('SELECT layout_json FROM gcfs WHERE gcf=?', (gcf,)).fetchone()
                layout = json.loads(item[0]) if item else {}
                from bin.reference_store import validate_layout
                ref = reference_from_layout(layout)
                validate_layout(ref)
                require_layout_policy(ref, policy)
                layouts[gcf] = layout_metadata(ref)
            if row['genome_len'] != sum(layouts[gcf]['lengths']):
                raise ValueError('Genome length/layout mismatch')
            selected = {k: row[k] for k in ('id', 'gcf', 'sample', 'genome_len', 'full_mean_depth',
                         'covered_fraction', 'aligned_bases', 'manifest_path', 'manifest_sha',
                         'checksums_sha', 'geometry_path')}
            selected.update(row_id=f'r{row_number+1:08d}',
                            row_key=json.dumps([row['sample'], gcf], ensure_ascii=False, separators=(',', ':')),
                            source_id=sid, source_root=str(root))
            rows.append(selected)
    finally:
        con.close()
    if sha(index) != fingerprint:
        raise ValueError('Source index changed during planning')
    strata = strata_plan(rows, layouts, breaks, bin_sizes, strata_config)
    out = new_directory(outdir)
    write_tsv(out/'rows.tsv', ['row_id', 'row_key', 'sample', 'gcf', 'reference_id'],
              [dict(row_id=r['row_id'], row_key=r['row_key'], sample=r['sample'], gcf=r['gcf'], reference_id=r['id']) for r in rows])
    write_tsv(out/'candidates.tsv', ['stratum', 'status', 'records', 'gcfs', 'minimum_genome_bp', 'bin_count', 'bin_sizes_bp', 'matrix_bytes'],
              [dict(stratum=s['id'], status=s['status'], records=len(s['row_indices']), gcfs=len(s['gcfs']),
                    minimum_genome_bp=s['min_genome_bp'], bin_count=c.get('bin_count'),
                    bin_sizes_bp=','.join(map(str,c.get('bin_sizes_bp',[]))), matrix_bytes=c.get('raw_plus_normalized_bytes',0))
               for s in strata for c in (s['candidates'] or [{}])])
    manifest = dict(kind='plan', schema=1, preflight_dir=str(preflight_dir), preflight_sha256=ledger,
                    index=str(index), index_sha256=fingerprint, rows=rows, layouts=layouts, strata=strata,
                    bin_sizes_bp=bin_sizes, genome_size_breaks_bp=breaks,
                    window_policy=GCF_CONCAT if policy == GCF_CONCAT else LEGACY_MATRIX_WINDOWS,
                    planned_matrix_bytes=sum(c['raw_plus_normalized_bytes'] for s in strata for c in s['candidates']))
    if strata_config is not None:
        manifest['strata_config'] = strata_config
    return seal(out, manifest)


def load_pool(row, layout):
    from bin.coverage_core import merge_references
    from bin.reference_store import load_reference_bundle
    manifest = child(row['source_root'], row['manifest_path'])
    status = manifest.parent / 'status'
    if (not (status/'SUCCESS').is_file() or (status/'FAILED').exists() or (status/'RUNNING').exists()
            or (status/'exit_code.txt').read_text().strip() != '0'):
        raise ValueError('Source attempt no longer successful')
    checksums = manifest.parent / 'checksums.json'
    if sha(manifest) != row['manifest_sha'] or sha(checksums) != row['checksums_sha']:
        raise ValueError('Source manifest/checksum changed')
    record = json.loads(manifest.read_text())
    if any(record[k] != row[k] for k in ('gcf', 'sample', 'genome_len', 'aligned_bases', 'full_mean_depth')):
        raise ValueError('Source record differs from index')
    geometry = child(row['source_root'], row['geometry_path'])
    outer = json.loads(checksums.read_text())
    for file in ('manifest.json', 'storage.json', 'checksums.json'):
        if outer.get('geometry/' + file) != sha(geometry/file):
            raise ValueError('Geometry metadata is not tied to source checksums')
    pool, observed, metadata = load_reference_bundle(geometry, workers=1)
    expected = reference_from_layout(layout)
    merge_references(expected, observed)
    if expected.sort_policy and window_policy(record['native']) != GCF_CONCAT:
        raise ValueError('Saved geometry record uses a different window policy')
    if record['layout'] != metadata['layout'] or pool.total_bases != row['aligned_bases']:
        raise ValueError('Geometry layout/mass mismatch')
    if not math.isclose(pool.total_bases/pool.genome_len, row['full_mean_depth'], rel_tol=1e-12):
        raise ValueError('Geometry full-depth mismatch')
    return pool


@contextmanager
def mapped_matrix(path, mode='r', shape=None):
    """Close the OS mapping even while Python references survive (important on NFS)."""
    matrix = (np.lib.format.open_memmap(path, mode=mode, dtype='float64', shape=shape)
              if shape is not None else np.load(path, mmap_mode=mode, allow_pickle=False))
    try:
        yield matrix
    finally:
        try:
            if mode != 'r':
                matrix.flush()
        finally:
            matrix._mmap.close()


def matrix_payload_bytes(raw_plus_normalized_bytes, normalized_only):
    return raw_plus_normalized_bytes // 2 if normalized_only else raw_plus_normalized_bytes


def check_sharded_matrix_budget(matrix_bytes, max_matrix_gib):
    required = 2 * matrix_bytes
    if (not math.isfinite(max_matrix_gib) or max_matrix_gib <= 0
            or required > max_matrix_gib * 1024**3):
        raise ValueError(f'Multi-node build needs {required/1024**3:.2f} GiB for both shards and '
                         'merged matrices within --max-matrix-gib')


def build(plan_dir, outdir, max_matrix_gib=64, *, jobs=1, row_shard=None, normalized_only=False):
    from bin.coverage_core import window_sums
    from bin.reference_store import bounded_map
    if type(jobs) is not int or jobs < 1:
        raise ValueError('Build jobs must be a positive integer')
    plan_dir, plan = load(plan_dir, 'plan')
    policy = matrix_policy(plan)
    if row_shard is not None and not 0 <= row_shard[0] < row_shard[1]:
        raise ValueError('Invalid matrix shard rank')
    payload_bytes = matrix_payload_bytes(plan['planned_matrix_bytes'], normalized_only)
    if not math.isfinite(max_matrix_gib) or max_matrix_gib <= 0 or payload_bytes > max_matrix_gib*1024**3:
        raise ValueError('Planned matrix disk payload exceeds --max-matrix-gib (or invalid limit)')
    if sha(Path(plan['index'])) != plan['index_sha256']:
        raise ValueError('Source index changed since planning')
    for name, expected in plan['preflight_sha256'].items():
        if sha(child(plan['preflight_dir'], name)) != expected:
            raise ValueError('Frozen preflight changed since planning')
    out = new_directory(outdir)
    results = []
    for stratum in plan['strata']:
        if stratum['status'] != 'READY':
            continue
        with ExitStack() as mappings:
            indices = stratum['row_indices']
            if row_shard is not None:
                rank, size = row_shard
                indices = indices[rank::size]
            if not indices:
                continue
            rows = [plan['rows'][i] for i in indices]
            arrays, entries = [], []
            lengths = {g: sum(plan['layouts'][g]['lengths']) for g in stratum['gcfs']}
            for candidate in stratum['candidates']:
                b = candidate['bin_count']
                folder = out / stratum['id'] / f'B{b}'
                folder.mkdir(parents=True)
                shape = (len(rows), b)
                raw = None
                if not normalized_only:
                    raw = mappings.enter_context(mapped_matrix(folder/'raw_mean_depth.npy', mode='w+', shape=shape))
                norm = mappings.enter_context(mapped_matrix(folder/'normalized.npy', mode='w+', shape=shape))
                write_edges(folder/'edges.json', lengths, b, policy)
                order = [dict(row_id=r['row_id'], row_key=r['row_key'], sample=r['sample'], gcf=r['gcf'], reference_id=r['id']) for r in rows]
                write_tsv(folder/'rows.tsv', ['row_id','row_key','sample','gcf','reference_id'], order)
                entry = dict(stratum=stratum['id'], bin_count=b, bin_sizes_bp=candidate['bin_sizes_bp'],
                             folder=str(folder.relative_to(out)), shape=[len(rows), b])
                entries.append(entry)
                arrays.append((raw, norm, b))
            print(f"BUILD {stratum['id']}: {len(rows)} records, {len(entries)} candidate counts", file=sys.stderr, flush=True)

            edge_lock = Lock()

            @lru_cache(maxsize=4)
            def boundaries(length):
                return tuple(matrix_edges(length, bins, policy) for _, _, bins in arrays)

            def build_row(i):
                row = rows[i]
                started = perf_counter()
                pool = load_pool(row, plan['layouts'][row['gcf']])
                loaded = perf_counter()
                ids = np.arange(len(pool), dtype=np.int64)
                with edge_lock:
                    row_edges = boundaries(lengths[row['gcf']])
                for (raw, norm, _), edges in zip(arrays, row_edges):
                    mass = window_sums(pool, ids, edges)
                    if int(mass.sum()) != pool.total_bases:
                        raise ValueError('Rebinning mass conservation failed')
                    mean = mass / np.diff(edges)
                    if raw is not None:
                        raw[i] = mean
                    norm[i] = normalized(mean)
                return loaded - started, perf_counter() - loaded

            worker_count = min(jobs, len(rows))
            indices = range(len(rows))
            # Close the worker iterator before the mappings, including on failure.
            times = (mappings.enter_context(closing(bounded_map(build_row, indices, worker_count)))
                     if worker_count > 1 else map(build_row, indices))
            load_seconds = bin_seconds = 0.0
            started = perf_counter()
            for completed, (load_time, bin_time) in enumerate(times, 1):
                load_seconds += load_time
                bin_seconds += bin_time
                if completed % 50 == 0 or completed == len(rows):
                    print(f"BUILD {stratum['id']} {completed}/{len(rows)}", file=sys.stderr, flush=True)
            row_wall_seconds = perf_counter() - started
            started = perf_counter()
            boundaries.cache_clear()
            for raw, norm, _ in arrays:
                if raw is not None:
                    raw.flush()
                norm.flush()
            print(f"TIMING build {stratum['id']} jobs={worker_count} rows_wall={row_wall_seconds:.2f}s "
                  f"geometry_load_verify_sum={load_seconds:.2f}s "
                  f"bin_normalize_write_sum={bin_seconds:.2f}s flush={perf_counter()-started:.2f}s",
                  file=sys.stderr, flush=True)
            arrays.clear()
            del raw, norm
            results.extend(entries)
    kind = 'matrices' if row_shard is None else 'matrix-shard'
    metadata = {} if row_shard is None else {'row_shard': list(row_shard)}
    if normalized_only:
        metadata['normalized_only'] = True
    return seal(out, dict(kind=kind, schema=1, plan_dir=str(plan_dir), plan_manifest_sha256=sha(plan_dir/'result.json'),
                         candidates=results, normalization_target=1e6,
                         build_jobs=jobs,
                         window_policy=policy,
                         source_geometry_verified=True, **metadata))
