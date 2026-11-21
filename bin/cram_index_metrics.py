#!/usr/bin/env python3
"""Compute and merge a chunked CRAM index with query-compatible metrics."""
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

VERSION = '2.1.0'
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
    'cram_size': 'INTEGER', 'cram_mtime_ns': 'INTEGER',
    'fna_size': 'INTEGER', 'fna_mtime_ns': 'INTEGER',
    'fai_size': 'INTEGER', 'fai_mtime_ns': 'INTEGER',
}
AUDIT_COLUMNS = ('id', 'Sample', 'GCFID', 'CramFilePath', 'status', 'included',
                 'read1_seen', 'read2_seen', 'unlabelled_seen', 'invalid_flags_seen', 'scan_complete')


def stamp(path):
    stat = Path(path).stat()
    return (stat.st_dev, stat.st_ino, stat.st_size, stat.st_mtime_ns)


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
    for name, observed in zip(('cram', 'fna', 'fai'), before):
        metric[name + '_size'] = observed[2]
        metric[name + '_mtime_ns'] = observed[3]
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
    rows = json.loads(chunk.read_text())
    columns = ['Sample', 'CramFilePath', 'GCFID', 'Coverage', 'Depth']
    specs = []
    seen = set()
    for number, values in enumerate(rows, 1):
        if (len(values) != 5 or any(not isinstance(v, str) or not v or v != v.strip() for v in values)
                or not all(Path(v).is_absolute() for v in (values[1], values[3], values[4]))):
            raise ValueError('Invalid row-channel identity/reference paths')
        sample, path, gcf, fna, fai = values
        if Path(fna).name != gcf + '.fna' or Path(fai).name not in (gcf+'.fai',gcf+'.fna.fai'):
            raise ValueError('Reference basename differs from exact GCF')
        key = sample, path, gcf
        if key in seen:
            raise ValueError('Duplicate identity within chunk')
        seen.add(key)
        specs.append(dict(id=number, row=dict(Sample=sample, CramFilePath=path, GCFID=gcf,
                                              Coverage='', Depth=''), fna=fna, fai=fai))
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


def merge(records, chunks, fna, fai, output, min_mapq, m5_policy):
    if not chunks or len({p.resolve() for p in chunks}) != len(chunks):
        raise ValueError('Missing or duplicate chunk databases')
    before = [stamp(p) for p in (records, fna, fai)]
    with records.open(encoding='utf-8-sig', newline='') as handle:
        reader = csv.DictReader(handle)
        columns = reader.fieldnames
        if (not columns or not REQUIRED <= set(columns)
                or len({c.casefold() for c in columns}) != len(columns)
                or any(not c.strip() for c in columns) or 'id' in {c.casefold() for c in columns}):
            raise ValueError('Require unique original CSV fields with Sample/CramFilePath/GCFID/Coverage/Depth')
        with output.open('xb'):
            pass
        try:
            with closing(sqlite3.connect(output, uri=True)) as connection, connection:
                connection.row_factory = sqlite3.Row
                tables(connection, columns)
                metric_fields = ','.join(quoted(k) for k in METRICS)
                metric_defs = ','.join(quoted(k)+' '+v+' NOT NULL' for k,v in METRICS.items())
                connection.execute('CREATE TEMP TABLE measured_sources (Sample TEXT,GCFID TEXT,CramFilePath TEXT,'
                                   + metric_defs + ',PRIMARY KEY(Sample,GCFID,CramFilePath))')
                connection.execute('CREATE TEMP TABLE wanted (Sample TEXT,GCFID TEXT,CramFilePath TEXT,'
                                   'PRIMARY KEY(Sample,GCFID,CramFilePath))')
                for chunk in sorted(chunks):
                    connection.execute('ATTACH DATABASE ? AS part', (chunk.resolve().as_uri()+'?mode=ro',))
                    try:
                        info = connection.execute('SELECT * FROM part.chunk_info').fetchone()
                        if not info or tuple(info[2:]) != (min_mapq,m5_policy):
                            raise ValueError('Chunk calculation policy differs')
                        connection.execute('INSERT INTO measured_sources SELECT s.Sample,s.GCFID,s.CramFilePath,'
                            + ','.join('m.'+quoted(k) for k in METRICS)
                            + ' FROM part.source_records s JOIN part.reference_metrics m ON s.id=m.id')
                        connection.commit()
                    except BaseException:
                        connection.rollback()
                        raise
                    finally:
                        connection.execute('DETACH DATABASE part')
                count = 0
                for count,row in enumerate(reader,1):
                    if None in row or any(v is None for v in row.values()):
                        raise ValueError(f'CSV row {count}: wrong field count')
                    key = (row['Sample'],row['GCFID'],row['CramFilePath'])
                    measured = connection.execute('SELECT '+metric_fields+' FROM measured_sources '
                        'WHERE Sample=? AND GCFID=? AND CramFilePath=?',key).fetchone()
                    if measured is None:
                        raise ValueError(f'CSV row {count}: no computed metrics; missing reference mapping or failed chunk')
                    insert_result(connection,columns,dict(id=count,row=row),dict(measured))
                    connection.execute('INSERT OR IGNORE INTO wanted VALUES (?,?,?)',key)
                if not count:
                    raise ValueError('Empty source CSV')
                if connection.execute('SELECT COUNT(*) FROM measured_sources').fetchone()[0] != connection.execute(
                        'SELECT COUNT(*) FROM wanted').fetchone()[0]:
                    raise ValueError('Extra computed identities do not belong to current source CSV')
                connection.execute('CREATE INDEX reference_gcf_sample ON reference_records(GCFID,Sample,id)')
                included = connection.execute('SELECT COUNT(*) FROM reference_records').fetchone()[0]
                gcfs = connection.execute('SELECT COUNT(DISTINCT GCFID) FROM reference_records').fetchone()[0]
                meta = dict(schema='coverage-raw-cram-index-1',version=VERSION,columns=columns,
                    source_csv=str(records.resolve()),source_records=count,records=included,gcfs=gcfs,
                    excluded_mate_records=count-included,mate_presence_audited=True,mate_audit_file=None,
                    mate_filter='At least one identifiable mate label; no observed invalid flags',
                    row_order='1-based current CSV data-record order; metrics cached by stable identity',
                    coverage_depth_basis='recalculated aligned M/= /X and union / full FAI length',
                    coverage_unit='percent',alignment_files_checked=True,
                    metrics_schema='coverage-reference-metrics-1',min_mapq=min_mapq,
                    coverage_alignment_policy='primary_only',m5_policy=m5_policy,reference_metrics_complete=True)
                meta.update(add_catalog(connection,fna,fai))
                meta.update(reference_files_checked=True,
                            reference_files_checked_scope='GCFs in source CRAM rows; unused catalog paths not opened')
                connection.execute('CREATE TABLE metadata (key TEXT PRIMARY KEY,value TEXT NOT NULL)')
                connection.executemany('INSERT INTO metadata VALUES (?,?)',
                                       ((k,json.dumps(v,ensure_ascii=False)) for k,v in meta.items()))
                if before != [stamp(p) for p in (records,fna,fai)]:
                    raise ValueError('Input CSV changed during merge')
                if connection.execute('PRAGMA main.quick_check').fetchone()[0] != 'ok':
                    raise ValueError('SQLite quick_check failed')
                connection.commit()
            print(f'Complete: {included} eligible / {count} computed rows -> {output}',flush=True)
            return meta
        except BaseException:
            output.unlink(missing_ok=True)
            raise


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--version', action='version', version=VERSION)
    sub = parser.add_subparsers(dest='command', required=True)
    p = sub.add_parser('compute')
    p.add_argument('--chunk', type=Path, required=True)
    p.add_argument('--jobs', type=int, default=1)
    p.add_argument('--min-mapq', type=int, default=0)
    p.add_argument('--m5-policy', choices=('strict', 'geometry_only'), default='geometry_only')
    p.add_argument('--output', type=Path, required=True)
    p = sub.add_parser('merge')
    for name in ('records', 'fna', 'fai', 'output'):
        p.add_argument('--' + name, type=Path, required=True)
    shards = p.add_mutually_exclusive_group(required=True)
    shards.add_argument('--chunks', nargs='+', type=Path)
    shards.add_argument('--chunk-dir', type=Path)
    p.add_argument('--min-mapq', type=int, default=0)
    p.add_argument('--m5-policy', choices=('strict', 'geometry_only'), default='geometry_only')
    args = parser.parse_args()
    try:
        if args.command == 'compute':
            compute(args.chunk, args.output, args.jobs, args.min_mapq, args.m5_policy)
        else:
            merge(args.records, args.chunks if args.chunks is not None else sorted(args.chunk_dir.glob('*.sqlite')), args.fna, args.fai, args.output, args.min_mapq, args.m5_policy)
    except (OSError, ValueError, sqlite3.Error) as error:
        parser.exit(1, f'ERROR: {error}\n')


if __name__ == '__main__':
    main()
