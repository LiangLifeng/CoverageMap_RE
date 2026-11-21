#!/usr/bin/env python3
"""Plan, compute and merge a chunked CRAM index with query-compatible metrics."""
import argparse
from concurrent.futures import ProcessPoolExecutor
import csv
from contextlib import closing
from functools import lru_cache
import json
import multiprocessing
from pathlib import Path
import sqlite3
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from bin.build_cram_index import REQUIRED, add_catalog, load_catalog, quoted

VERSION = '2.0.0'
METRICS = {
    'genome_length': 'INTEGER', 'aligned_bases': 'INTEGER', 'covered_bases': 'INTEGER',
    'read_count': 'INTEGER', 'fragment_count': 'INTEGER',
    'read1_count': 'INTEGER', 'read2_count': 'INTEGER', 'unlabelled_count': 'INTEGER',
    'read1_bases': 'INTEGER', 'read2_bases': 'INTEGER', 'unlabelled_bases': 'INTEGER',
    'read1_depth': 'REAL', 'read2_depth': 'REAL',
    'read1_covered_bases': 'INTEGER', 'read2_covered_bases': 'INTEGER',
    'read1_coverage': 'REAL', 'read2_coverage': 'REAL',
    'Depth': 'REAL', 'Coverage': 'REAL', 'included_reference': 'INTEGER',
    'layout_json': 'TEXT', 'm5_audit_json': 'TEXT',
}
AUDIT_COLUMNS = ('id', 'Sample', 'GCFID', 'CramFilePath', 'status', 'included',
                 'read1_seen', 'read2_seen', 'unlabelled_seen', 'invalid_flags_seen', 'scan_complete')


def stamp(path):
    stat = Path(path).stat()
    return (stat.st_dev, stat.st_ino, stat.st_size, stat.st_mtime_ns)


def plan(records, fna, fai, outdir, chunk_size):
    if chunk_size < 1:
        raise ValueError('chunk-size must be positive')
    fna_map, _ = load_catalog(fna, 'GCF_path')
    fai_map, _ = load_catalog(fai, 'FAI_path')
    before = {str(p): stamp(p) for p in (records, fna, fai)}
    outdir.mkdir(exist_ok=False)
    (outdir / 'chunks').mkdir()
    with records.open(encoding='utf-8-sig', newline='') as handle:
        reader = csv.DictReader(handle)
        columns = reader.fieldnames
        if (not columns or not REQUIRED <= set(columns)
                or len({c.casefold() for c in columns}) != len(columns)
                or any(not c.strip() for c in columns) or 'id' in {c.casefold() for c in columns}):
            raise ValueError('Require unique columns, no id, and ' + ','.join(sorted(REQUIRED)))
        batch, count, chunk_count = [], 0, 0
        for count, row in enumerate(reader, 1):
            if None in row or any(v is None for v in row.values()):
                raise ValueError(f'CSV row {count}: wrong field count')
            for name in ('Sample', 'GCFID', 'CramFilePath'):
                if not row[name] or row[name] != row[name].strip():
                    raise ValueError(f'CSV row {count}: empty/padded {name}')
            if not Path(row['CramFilePath']).is_absolute():
                raise ValueError(f'CSV row {count}: alignment path must be absolute')
            gcf = row['GCFID']
            if gcf not in fna_map or gcf not in fai_map:
                raise ValueError(f'CSV row {count}: missing exact FNA/FAI mapping: {gcf}')
            # Source Depth/Coverage are opaque provenance, never parsed or used.
            batch.append(dict(id=count, row=row, fna=fna_map[gcf], fai=fai_map[gcf]))
            if len(batch) == chunk_size:
                chunk_count += 1
                write_chunk(outdir, chunk_count, columns, batch)
                batch = []
        if batch:
            chunk_count += 1
            write_chunk(outdir, chunk_count, columns, batch)
    if not count:
        raise ValueError('CSV has no data records')
    if before != {str(p): stamp(p) for p in (records, fna, fai)}:
        raise ValueError('Input CSV changed during planning')
    manifest = dict(columns=columns, source_records=count, chunks=chunk_count,
                    source_csv=str(records.resolve()), fna_csv=str(fna.resolve()), fai_csv=str(fai.resolve()))
    task_rows = []
    for chunk in sorted((outdir / 'chunks').glob('*.jsonl')):
        with chunk.open() as handle:
            next(handle)
            specs = [json.loads(line) for line in handle]
        dependencies = sorted({path for spec in specs for path in
                               (spec['row']['CramFilePath'], spec['fna'], spec['fai'])})
        if any(any(c in path for c in '|\t\r\n') for path in dependencies):
            raise ValueError('Dependency paths cannot contain pipe/tab/newline')
        task_rows.append((chunk.stem, str(chunk.relative_to(outdir)), '|'.join(dependencies)))
    with (outdir / 'tasks.tsv').open('x', newline='') as handle:
        writer = csv.writer(handle, delimiter='\t', lineterminator='\n')
        writer.writerow(('id', 'chunk', 'dependencies'))
        writer.writerows(task_rows)
    (outdir / 'plan.json').write_text(json.dumps(manifest, ensure_ascii=False, indent=2) + '\n')
    print(f'Planned {count} rows / {chunk_count} chunks of at most {chunk_size}', flush=True)
    return manifest


def write_chunk(outdir, number, columns, batch):
    with (outdir / 'chunks' / f'chunk-{number:06d}.jsonl').open('x') as handle:
        handle.write(json.dumps(dict(columns=columns), ensure_ascii=False) + '\n')
        for row in batch:
            handle.write(json.dumps(row, ensure_ascii=False) + '\n')


@lru_cache(maxsize=16)
def reference(fna, fai, fna_stamp, fai_stamp):
    from bin.coverage_core import canonical_reference, read_reference, read_reference_index
    ref = canonical_reference(read_reference(fna))
    index = canonical_reference(read_reference_index(fai))
    if (ref.names, ref.lengths) != (index.names, index.lengths):
        raise ValueError('FNA and FAI layouts differ')
    return ref


def measure(item):
    import numpy as np
    from bin.coverage_core import covered_sums, layout_metadata
    from bin.raw_cram_pool import extract_once
    spec, min_mapq, m5_policy = item
    row, fna, fai = spec['row'], spec['fna'], spec['fai']
    paths = (row['CramFilePath'], fna, fai)
    before = [stamp(p) for p in paths]
    ref = reference(fna, fai, before[1], before[2])
    pool, merged, missing = extract_once(paths[0], ref, sampling_unit='read',
        min_mapq=min_mapq, coverage_alignment_policy='primary_only',
        strict_reference=True, m5_policy=m5_policy)
    edges = np.array([0, ref.length], dtype=np.int64)
    metric = dict(genome_length=ref.length, aligned_bases=int(pool.total_bases),
                  covered_bases=int(covered_sums(pool, edges).sum()),
                  read_count=len(pool), fragment_count=int(len(np.unique(pool.unit_fragment_ids))),
                  layout_json=json.dumps(layout_metadata(merged), separators=(',', ':')),
                  m5_audit_json=json.dumps(getattr(pool, 'm5_audit', []), separators=(',', ':')))
    roles = pool.unit_mate_sides
    for name, role in [('read1', 1), ('read2', 2), ('unlabelled', 0)]:
        ids = np.flatnonzero(roles == role)
        metric[name + '_count'] = int(len(ids))
        metric[name + '_bases'] = int(pool.unit_bases[ids].sum(dtype=np.int64))
        if role:
            metric[name + '_depth'] = metric[name + '_bases'] / ref.length
            metric[name + '_covered_bases'] = int(covered_sums(pool, edges, ids).sum())
            metric[name + '_coverage'] = 100 * metric[name + '_covered_bases'] / ref.length
    if sum(metric[k + '_bases'] for k in ('read1', 'read2', 'unlabelled')) != metric['aligned_bases']:
        raise ValueError('Mate-side capacity does not conserve aligned bases')
    metric.update(Depth=metric['aligned_bases'] / ref.length,
                  Coverage=100 * metric['covered_bases'] / ref.length,
                  included_reference=int(bool(metric['read1_count'] or metric['read2_count'])))
    if before != [stamp(p) for p in paths]:
        raise ValueError(f'Input changed while measuring row {spec["id"]}')
    return spec, metric


def tables(connection, columns):
    definition = ', '.join(quoted(c) + ' TEXT NOT NULL' for c in columns)
    for table in ('reference_records', 'source_records'):
        connection.execute(f'CREATE TABLE {table} (id INTEGER PRIMARY KEY, {definition})')
    definition = ','.join(quoted(k) + ' ' + t + ' NOT NULL' for k, t in METRICS.items())
    connection.execute('CREATE TABLE reference_metrics (id INTEGER PRIMARY KEY,' + definition + ')')
    connection.execute('''CREATE TABLE reference_mate_audit (
        id INTEGER PRIMARY KEY, Sample TEXT NOT NULL, GCFID TEXT NOT NULL,
        CramFilePath TEXT NOT NULL, status TEXT NOT NULL, included INTEGER NOT NULL,
        read1_seen INTEGER NOT NULL, read2_seen INTEGER NOT NULL,
        unlabelled_seen INTEGER NOT NULL, invalid_flags_seen INTEGER NOT NULL,
        scan_complete INTEGER NOT NULL)''')


def insert_result(connection, columns, spec, metric):
    row, ident = spec['row'], spec['id']
    marks = ','.join('?' for _ in range(len(columns) + 1))
    connection.execute(f'INSERT INTO source_records VALUES ({marks})', (ident, *(row[c] for c in columns)))
    if metric['included_reference']:
        computed = dict(row, Depth=format(metric['Depth'], '.17g'), Coverage=format(metric['Coverage'], '.17g'))
        connection.execute(f'INSERT INTO reference_records VALUES ({marks})',
                           (ident, *(computed[c] for c in columns)))
    connection.execute('INSERT INTO reference_metrics VALUES (' + ','.join('?' for _ in range(len(METRICS)+1)) + ')',
                       (ident, *(metric[k] for k in METRICS)))
    r1, r2, unlabelled = (int(metric[n + '_count'] > 0) for n in ('read1', 'read2', 'unlabelled'))
    state = ('BOTH_LABELS_PRESENT' if r1 and r2 else 'ONLY_READ1' if r1 else
             'ONLY_READ2' if r2 else 'NO_MATE_LABELS' if unlabelled else 'NO_ELIGIBLE_PRIMARY')
    connection.execute('INSERT INTO reference_mate_audit VALUES (?,?,?,?,?,?,?,?,?,?,?)',
        (ident, row['Sample'], row['GCFID'], row['CramFilePath'], state,
         metric['included_reference'], r1, r2, unlabelled, 0, 1))


def compute(chunk, output, jobs, min_mapq, m5_policy):
    if jobs < 1 or min_mapq < 0:
        raise ValueError('jobs must be positive and min-mapq nonnegative')
    with chunk.open() as handle:
        columns = json.loads(next(handle))['columns']
        specs = [json.loads(line) for line in handle]
    if not specs:
        raise ValueError('Empty chunk')
    with output.open('xb'):
        pass
    try:
        with closing(sqlite3.connect(output, uri=True)) as connection, connection:
            tables(connection, columns)
            inputs = [(s, min_mapq, m5_policy) for s in specs]
            if jobs == 1:
                for number, result in enumerate(map(measure, inputs), 1):
                    insert_result(connection, columns, *result)
                    print(f'{chunk.stem}: {number}/{len(specs)}', flush=True)
            else:
                with ProcessPoolExecutor(max_workers=min(jobs, len(specs)),
                                         mp_context=multiprocessing.get_context('spawn')) as executor:
                    # Each process decodes one CRAM at a time; no nested decoder threads.
                    for number, result in enumerate(executor.map(measure, inputs), 1):
                        insert_result(connection, columns, *result)
                        print(f'{chunk.stem}: {number}/{len(specs)}', flush=True)
            connection.execute('CREATE TABLE chunk_info (columns_json TEXT, rows INTEGER, min_mapq INTEGER, m5_policy TEXT)')
            connection.execute('INSERT INTO chunk_info VALUES (?,?,?,?)',
                               (json.dumps(columns), len(specs), min_mapq, m5_policy))
    except BaseException:
        output.unlink(missing_ok=True)
        raise


def merge(plan_file, chunks, fna, fai, output, min_mapq, m5_policy):
    manifest = json.loads(plan_file.read_text())
    columns = manifest['columns']
    if len(chunks) != manifest['chunks'] or len({p.resolve() for p in chunks}) != len(chunks):
        raise ValueError('Missing or duplicate chunk databases')
    with output.open('xb'):
        pass
    try:
        with closing(sqlite3.connect(output, uri=True)) as connection, connection:
            tables(connection, columns)
            for chunk in sorted(chunks):
                # Read-only attach; databases contain only compact scalar metrics/provenance.
                connection.execute('ATTACH DATABASE ? AS part', (chunk.resolve().as_uri() + '?mode=ro',))
                try:
                    info = connection.execute('SELECT * FROM part.chunk_info').fetchone()
                    if not info or json.loads(info[0]) != columns or tuple(info[2:]) != (min_mapq, m5_policy):
                        raise ValueError('Chunk schema/calculation policy differs')
                    for table in ('source_records', 'reference_records', 'reference_metrics', 'reference_mate_audit'):
                        connection.execute(f'INSERT INTO main.{table} SELECT * FROM part.{table} ORDER BY id')
                    connection.commit()
                except BaseException:
                    connection.rollback()
                    raise
                finally:
                    connection.execute('DETACH DATABASE part')
            count, lo, hi = connection.execute('SELECT COUNT(*),MIN(id),MAX(id) FROM source_records').fetchone()
            if (count, lo, hi) != (manifest['source_records'], 1, manifest['source_records']):
                raise ValueError('Chunk rows do not exhaust original input IDs')
            for table in ('reference_metrics', 'reference_mate_audit'):
                if connection.execute(f'SELECT COUNT(*) FROM {table}').fetchone()[0] != count:
                    raise ValueError('Missing computed metrics/audit rows')
            connection.execute('CREATE INDEX reference_gcf_sample ON reference_records(GCFID, Sample, id)')
            records = connection.execute('SELECT COUNT(*) FROM reference_records').fetchone()[0]
            gcfs = connection.execute('SELECT COUNT(DISTINCT GCFID) FROM reference_records').fetchone()[0]
            meta = dict(schema='coverage-raw-cram-index-1', version=VERSION, columns=columns,
                source_csv=manifest['source_csv'], source_records=count, records=records, gcfs=gcfs,
                excluded_mate_records=count-records, mate_presence_audited=True, mate_audit_file=None,
                mate_filter='At least one identifiable mate label; no observed invalid flags',
                row_order='1-based original CSV data-record order',
                coverage_depth_basis='recalculated primary M/= /X aligned bases and covered union / full FAI length',
                coverage_unit='percent', alignment_files_checked=True,
                metrics_schema='coverage-reference-metrics-1', min_mapq=min_mapq,
                coverage_alignment_policy='primary_only', m5_policy=m5_policy,
                reference_metrics_complete=True)
            meta.update(add_catalog(connection, fna, fai))
            meta['reference_files_checked'] = True
            meta['reference_files_checked_scope'] = 'GCFs used by original CRAM rows; unused catalog paths not opened'
            connection.execute('CREATE TABLE metadata (key TEXT PRIMARY KEY,value TEXT NOT NULL)')
            connection.executemany('INSERT INTO metadata VALUES (?,?)',
                ((k, json.dumps(v, ensure_ascii=False)) for k, v in meta.items()))
            if connection.execute('PRAGMA main.quick_check').fetchone()[0] != 'ok':
                raise ValueError('SQLite quick_check failed')
            connection.commit()
        print(f'Complete: {records} eligible / {count} computed rows -> {output}', flush=True)
        return meta
    except BaseException:
        output.unlink(missing_ok=True)
        raise


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--version', action='version', version=VERSION)
    sub = parser.add_subparsers(dest='command', required=True)
    p = sub.add_parser('plan')
    for name in ('records', 'fna', 'fai', 'outdir'):
        p.add_argument('--' + name, type=Path, required=True)
    p.add_argument('--chunk-size', type=int, default=100)
    p = sub.add_parser('compute')
    p.add_argument('--chunk', type=Path, required=True)
    p.add_argument('--jobs', type=int, default=1)
    p.add_argument('--min-mapq', type=int, default=0)
    p.add_argument('--m5-policy', choices=('strict', 'geometry_only'), default='geometry_only')
    p.add_argument('--output', type=Path, required=True)
    p = sub.add_parser('merge')
    for name in ('plan', 'fna', 'fai', 'output'):
        p.add_argument('--' + name, type=Path, required=True)
    shards = p.add_mutually_exclusive_group(required=True)
    shards.add_argument('--chunks', nargs='+', type=Path)
    shards.add_argument('--chunk-dir', type=Path)
    p.add_argument('--min-mapq', type=int, default=0)
    p.add_argument('--m5-policy', choices=('strict', 'geometry_only'), default='geometry_only')
    args = parser.parse_args()
    try:
        if args.command == 'plan':
            plan(args.records, args.fna, args.fai, args.outdir, args.chunk_size)
        elif args.command == 'compute':
            compute(args.chunk, args.output, args.jobs, args.min_mapq, args.m5_policy)
        else:
            merge(args.plan, args.chunks if args.chunks is not None else sorted(args.chunk_dir.glob('*.sqlite')), args.fna, args.fai, args.output, args.min_mapq, args.m5_policy)
    except (OSError, ValueError, sqlite3.Error) as error:
        parser.exit(1, f'ERROR: {error}\n')


if __name__ == '__main__':
    main()
