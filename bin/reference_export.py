#!/usr/bin/env python3
"""Export Sample–GCF metadata from a schema-3 reference SQLite index, read-only."""
import argparse
import csv
import json
import math
import os
from pathlib import Path
import sqlite3
import tempfile


COLUMNS = ('Sample', 'GCFID', 'Depth', 'CramFilePath', 'DatabaseRoot',
           'Selected', 'EligibleDepth', 'ExclusionReason')


def validate_filters(samples, gcfs, minimum, maximum):
    for name, values in (('Sample', samples), ('GCF', gcfs)):
        if any(not isinstance(value, str) or not value.strip() for value in values):
            raise ValueError(f'{name} IDs must be nonempty strings')
    for value in (minimum, maximum):
        if value is not None and (not math.isfinite(value) or value < 0):
            raise ValueError('Depth bounds must be finite and nonnegative')
    if minimum is not None and maximum is not None and minimum > maximum:
        raise ValueError('min-depth must be <= max-depth')


def query_sql(samples, gcfs, minimum, maximum, selected_only):
    clauses, parameters = [], []
    for column, values in (('sample', samples), ('gcf', gcfs)):
        if values:
            clauses.append(f'r.{column} IN ({",".join("?" for _ in values)})')
            parameters.extend(values)
    for operator, value in (('>=', minimum), ('<=', maximum)):
        if value is not None:
            clauses.append(f'r.full_mean_depth {operator} ?')
            parameters.append(value)
    if selected_only:
        clauses.append('r.selected=1')
    where = ' WHERE ' + ' AND '.join(clauses) if clauses else ''
    sql = ('SELECT r.sample,r.gcf,r.full_mean_depth,r.record_json,s.root,'
           'r.selected,r.eligible_depth,r.exclusion_reason '
           'FROM references_table r LEFT JOIN sources s ON s.id=r.source_id' + where +
           ' ORDER BY r.sample COLLATE BINARY,r.gcf COLLATE BINARY')
    return sql, parameters


def export_record(row):
    sample, gcf, depth, record_json, root, selected, eligible, reason = row
    if root is None or not isinstance(depth, (int, float)) or not math.isfinite(depth) or depth < 0:
        raise ValueError(f'Invalid source/depth for {sample}/{gcf}')
    if selected not in (0, 1) or eligible not in (0, 1):
        raise ValueError(f'Invalid selection flags for {sample}/{gcf}')
    # The builder records the alignment first, followed by reference/FAI files.
    # Never probe or open these paths; missing provenance stays visibly empty.
    record = json.loads(record_json)
    if not isinstance(record, dict):
        raise ValueError(f'Invalid record provenance for {sample}/{gcf}')
    sources = record.get('sources', [])
    if not isinstance(sources, list):
        raise ValueError(f'Invalid sources provenance for {sample}/{gcf}')
    path = sources[0].get('path', '') if sources and isinstance(sources[0], dict) else ''
    alignment = path if isinstance(path, str) and Path(path).suffix.lower() in ('.cram', '.bam') else ''
    return (sample, gcf, depth, alignment, root, selected, eligible, reason)


def export_index(index, output, samples=(), gcfs=(), minimum=None, maximum=None, selected_only=False):
    """Stream stable TSV rows; atomically publish without replacing any output."""
    samples, gcfs = sorted(set(samples)), sorted(set(gcfs))
    validate_filters(samples, gcfs, minimum, maximum)
    index = Path(index).resolve(strict=True)
    output = Path(output).absolute()
    if not index.is_file():
        raise ValueError('index must be a SQLite file')
    if os.path.lexists(output):
        raise FileExistsError(f'Output already exists: {output}')
    sql, parameters = query_sql(samples, gcfs, minimum, maximum, selected_only)
    temporary = None
    connection = sqlite3.connect(index.as_uri() + '?mode=ro', uri=True)
    try:
        connection.execute('PRAGMA query_only=ON')
        connection.execute('BEGIN')  # A single consistent metadata snapshot.
        metadata = {key: json.loads(value) for key, value in connection.execute('SELECT key,value FROM metadata')}
        if metadata.get('schema') != 'coverage-reference-sqlite-3':
            raise ValueError('reference_export requires coverage-reference-sqlite-3')
        cursor = connection.execute(sql, parameters)
        with tempfile.NamedTemporaryFile(mode='w', encoding='utf-8', newline='',
                                         prefix='.reference-export-', suffix='.tmp',
                                         dir=output.parent, delete=False) as handle:
            temporary = Path(handle.name)
            writer = csv.writer(handle, delimiter='\t', lineterminator='\n')
            writer.writerow(COLUMNS)
            count, missing = 0, 0
            sample_ids, gcf_ids = set(), set()
            for row in cursor:
                values = export_record(row)
                writer.writerow(values)
                count += 1
                missing += not bool(values[3])
                sample_ids.add(values[0]); gcf_ids.add(values[1])
            handle.flush()
            os.fsync(handle.fileno())
        os.link(temporary, output)
        return dict(index=str(index), out=str(output), records=count, samples=len(sample_ids),
                    gcfs=len(gcf_ids), missing_alignment_paths=missing,
                    filters=dict(samples=samples, gcfs=gcfs, min_depth=minimum,
                                 max_depth=maximum, selected_only=selected_only),
                    index_selection=metadata.get('selection'), schema=metadata['schema'],
                    plan_sha256=metadata.get('plan_sha256'))
    finally:
        connection.close()
        if temporary is not None:
            temporary.unlink(missing_ok=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--index', required=True, type=Path, help='Existing schema-3 SQLite index')
    parser.add_argument('--out', required=True, type=Path, help='New TSV path; parent directory must exist')
    parser.add_argument('--sample', action='append', default=[], help='Exact Sample ID; repeat for multiple IDs (OR)')
    parser.add_argument('--gcf', action='append', default=[], help='Exact GCF/reference ID; repeat for multiple IDs (OR)')
    parser.add_argument('--min-depth', type=float, help='Original full-genome mean depth lower bound, inclusive')
    parser.add_argument('--max-depth', type=float, help='Original full-genome mean depth upper bound, inclusive')
    parser.add_argument('--selected-only', action='store_true', help='Only rows already selected=1 in this index')
    args = parser.parse_args()
    try:
        result = export_index(args.index, args.out, args.sample, args.gcf,
                              args.min_depth, args.max_depth, args.selected_only)
    except (OSError, ValueError, sqlite3.Error) as error:
        parser.exit(2, f'ERROR: {error}\n')
    print(json.dumps(result, ensure_ascii=False, indent=2, sort_keys=True))


if __name__ == '__main__':
    main()
