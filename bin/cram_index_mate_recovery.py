#!/usr/bin/env python3
"""Compute and merge a chunked CRAM index with query-compatible metrics."""
import argparse
from concurrent.futures import ProcessPoolExecutor
import csv
from contextlib import closing
import json
import multiprocessing
from pathlib import Path
import sqlite3
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from bin.build_cram_index import REQUIRED, add_catalog, quoted
from bin.cram_index_metrics import METRICS, insert_result, measure, stamp, tables

VERSION = '2.1.0'
MATE_ERROR = 'Ambiguous read identity or multiple primary alignments'


def mate_reason(metric):
    if not metric['read_count']:
        return 'NO_ELIGIBLE_PRIMARY'
    if metric['unlabelled_count']:
        return 'UNLABELLED_READS'
    if not metric['read1_count']:
        return 'MISSING_R1'
    if not metric['read2_count']:
        return 'MISSING_R2'
    return None


def measure_with_mate_filter(item):
    try:
        spec, metric = measure(item)
    except ValueError as error:
        if str(error) != MATE_ERROR:
            raise
        return item[0], None, 'AMBIGUOUS_READ_IDENTITY'
    return spec, metric, mate_reason(metric)


def compute(chunk, output, jobs, min_mapq, m5_policy, strict_mates=False):
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
            operation = measure_with_mate_filter if strict_mates else measure
            if strict_mates:
                connection.execute('CREATE TABLE mate_exclusions ('
                    'Sample TEXT,GCFID TEXT,CramFilePath TEXT,reason TEXT NOT NULL,'
                    'PRIMARY KEY(Sample,GCFID,CramFilePath))')

            def save(results):
                for number, result in enumerate(results, 1):
                    if strict_mates and result[2]:
                        row = result[0]['row']
                        connection.execute('INSERT INTO mate_exclusions VALUES (?,?,?,?)',
                            (row['Sample'],row['GCFID'],row['CramFilePath'],result[2]))
                    else:
                        insert_result(connection, columns, *result[:2])
                    print(f'{chunk.stem}: {number}/{len(specs)}', flush=True)

            if jobs == 1:
                save(map(operation, inputs))
            else:
                with ProcessPoolExecutor(max_workers=min(jobs, len(specs)),
                                         mp_context=multiprocessing.get_context('spawn')) as executor:
                    save(executor.map(operation, inputs))
            connection.execute('CREATE TABLE chunk_info (columns_json TEXT, rows INTEGER, min_mapq INTEGER, m5_policy TEXT)')
            connection.execute('INSERT INTO chunk_info VALUES (?,?,?,?)',
                               (json.dumps(columns), len(specs), min_mapq, m5_policy))
    except BaseException:
        output.unlink(missing_ok=True)
        raise


def merge(records, chunks, fna, fai, output, min_mapq, m5_policy, excluded_tsv=None):
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
                if excluded_tsv is not None:
                    connection.execute('CREATE TEMP TABLE omissions (Sample TEXT,GCFID TEXT,CramFilePath TEXT,'
                                       'reason TEXT NOT NULL,PRIMARY KEY(Sample,GCFID,CramFilePath))')
                for chunk in sorted(chunks):
                    connection.execute('ATTACH DATABASE ? AS part', (chunk.resolve().as_uri()+'?mode=ro',))
                    try:
                        info = connection.execute('SELECT * FROM part.chunk_info').fetchone()
                        if not info or tuple(info[2:]) != (min_mapq,m5_policy):
                            raise ValueError('Chunk calculation policy differs')
                        connection.execute('INSERT INTO measured_sources SELECT s.Sample,s.GCFID,s.CramFilePath,'
                            + ','.join('m.'+quoted(k) for k in METRICS)
                            + ' FROM part.source_records s JOIN part.reference_metrics m ON s.id=m.id')
                        if excluded_tsv is not None and connection.execute(
                                "SELECT 1 FROM part.sqlite_master WHERE name='mate_exclusions'").fetchone():
                            connection.execute('INSERT INTO omissions SELECT * FROM part.mate_exclusions')
                        connection.commit()
                    except BaseException:
                        connection.rollback()
                        raise
                    finally:
                        connection.execute('DETACH DATABASE part')
                count = 0
                excluded_rows = []
                if excluded_tsv is not None and connection.execute(
                        'SELECT 1 FROM measured_sources JOIN omissions USING(Sample,GCFID,CramFilePath) LIMIT 1').fetchone():
                    raise ValueError('Conflicting metrics and exclusion for the same identity')
                for count,row in enumerate(reader,1):
                    if None in row or any(v is None for v in row.values()):
                        raise ValueError(f'CSV row {count}: wrong field count')
                    key = (row['Sample'],row['GCFID'],row['CramFilePath'])
                    measured = connection.execute('SELECT '+metric_fields+' FROM measured_sources '
                        'WHERE Sample=? AND GCFID=? AND CramFilePath=?',key).fetchone()
                    reason = None
                    if excluded_tsv is not None:
                        if measured is not None:
                            reason = mate_reason(measured)
                        else:
                            omission = connection.execute('SELECT reason FROM omissions '
                                'WHERE Sample=? AND GCFID=? AND CramFilePath=?',key).fetchone()
                            reason = omission['reason'] if omission else None
                    if reason:
                        excluded_rows.append((count, key[0], key[1], key[2], reason))
                    elif measured is not None:
                        insert_result(connection,columns,dict(id=count,row=row),dict(measured))
                    else:
                        raise ValueError(f'CSV row {count}: no metrics or explicit exclusion for {key}')
                    connection.execute('INSERT OR IGNORE INTO wanted VALUES (?,?,?)',key)
                if not count:
                    raise ValueError('Empty source CSV')
                processed = connection.execute('SELECT COUNT(*) FROM measured_sources').fetchone()[0]
                if excluded_tsv is not None:
                    processed += connection.execute('SELECT COUNT(*) FROM omissions').fetchone()[0]
                if processed != connection.execute('SELECT COUNT(*) FROM wanted').fetchone()[0]:
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
                if excluded_tsv is not None:
                    if included + len(excluded_rows) != count:
                        raise ValueError('Input rows do not match included plus excluded rows')
                    meta.update(mate_filter='Both R1 and R2 present; no unlabelled eligible reads or ambiguous primary identity',
                                mate_audit_file=excluded_tsv.name)
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
            if excluded_tsv is not None:
                with excluded_tsv.open('x', newline='') as handle:
                    writer = csv.writer(handle, delimiter='\t', lineterminator='\n')
                    writer.writerow(('OriginalCSVRow','Sample','GCFID','CramFilePath','Reason'))
                    writer.writerows(excluded_rows)
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
    p.add_argument('--strict-mates', action='store_true')
    p = sub.add_parser('merge')
    for name in ('records', 'fna', 'fai', 'output'):
        p.add_argument('--' + name, type=Path, required=True)
    shards = p.add_mutually_exclusive_group(required=True)
    shards.add_argument('--chunks', nargs='+', type=Path)
    shards.add_argument('--chunk-dir', type=Path)
    p.add_argument('--min-mapq', type=int, default=0)
    p.add_argument('--m5-policy', choices=('strict', 'geometry_only'), default='geometry_only')
    p.add_argument('--excluded-tsv', type=Path)
    args = parser.parse_args()
    try:
        if args.command == 'compute':
            compute(args.chunk, args.output, args.jobs, args.min_mapq, args.m5_policy, args.strict_mates)
        else:
            merge(args.records, args.chunks if args.chunks is not None else sorted(args.chunk_dir.glob('*.sqlite')), args.fna, args.fai, args.output, args.min_mapq, args.m5_policy, args.excluded_tsv)
    except (OSError, ValueError, sqlite3.Error) as error:
        parser.exit(1, f'ERROR: {error}\n')


if __name__ == '__main__':
    main()
