#!/usr/bin/env python3
"""Audit mapped primary mate flags from an explicit raw-CRAM SQL index range."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import csv
import json
from pathlib import Path
import sqlite3

VERSION = '1.0.0'
COLUMNS = ('id', 'Sample', 'GCFID', 'CramFilePath', 'status', 'read1_seen',
           'read2_seen', 'unlabelled_seen', 'invalid_flags_seen',
           'eligible_primary_observed', 'scan_complete', 'error')


def inspect_record(row, min_mapq=0):
    import pysam
    result = dict(row, read1_seen=False, read2_seen=False, unlabelled_seen=False,
                  invalid_flags_seen=False, eligible_primary_observed=0,
                  scan_complete=False, error='')
    path = Path(row['CramFilePath'])
    try:
        with path.open('rb') as handle:
            cram = handle.read(4) == b'CRAM'
        options = dict(format_options=[b'required_fields=4159', b'decode_md=0']) if cram else {}
        with pysam.AlignmentFile(str(path), 'rc' if cram else 'rb', **options) as bam:
            for read in bam.fetch(until_eof=True):
                if (read.is_unmapped or read.is_secondary or read.is_supplementary
                        or read.mapping_quality < min_mapq or not read.get_blocks()):
                    continue
                result['eligible_primary_observed'] += 1
                if read.is_read1 and not read.is_read2:
                    result['read1_seen'] = True
                elif read.is_read2 and not read.is_read1:
                    result['read2_seen'] = True
                elif read.is_paired or read.is_read1 or read.is_read2:
                    result['invalid_flags_seen'] = True
                else:
                    result['unlabelled_seen'] = True
                # Presence audit only: stop once both labels are proven. Counts
                # are explicitly observed lower bounds, not read distributions.
                if result['read1_seen'] and result['read2_seen']:
                    break
            else:
                result['scan_complete'] = True
        if result['read1_seen'] and result['read2_seen']:
            result['status'] = 'BOTH_LABELS_PRESENT'
        elif not result['eligible_primary_observed']:
            result['status'] = 'NO_ELIGIBLE_PRIMARY'
        elif result['invalid_flags_seen']:
            result['status'] = 'INVALID_MATE_FLAGS'
        elif result['read1_seen']:
            result['status'] = 'ONLY_READ1'
        elif result['read2_seen']:
            result['status'] = 'ONLY_READ2'
        else:
            result['status'] = 'NO_MATE_LABELS'
    except (OSError, ValueError) as error:
        result.update(status='ERROR', error=f'{type(error).__name__}: {error}')
    return result


def audit(index, output, start_id, end_id, jobs, min_mapq):
    if start_id < 1 or end_id < start_id or jobs < 1 or min_mapq < 0:
        raise ValueError('Invalid ID range, jobs or MAPQ')
    with sqlite3.connect(Path(index).resolve().as_uri()+'?mode=ro', uri=True) as connection:
        connection.row_factory = sqlite3.Row
        rows = [dict(r) for r in connection.execute(
            'SELECT id,Sample,GCFID,CramFilePath FROM reference_records '
            'WHERE id BETWEEN ? AND ? ORDER BY id', (start_id, end_id))]
    if not rows:
        raise ValueError('No indexed records in requested range')
    failures = 0
    with Path(output).open('x', newline='') as handle, ThreadPoolExecutor(max_workers=jobs) as executor:
        writer = csv.DictWriter(handle, COLUMNS, delimiter='\t', lineterminator='\n')
        writer.writeheader()
        for number, result in enumerate(executor.map(lambda r: inspect_record(r, min_mapq), rows), 1):
            writer.writerow(result)
            failures += result['status'] == 'ERROR'
            if number % 100 == 0 or number == len(rows):
                handle.flush()
                print(f'Audited {number}/{len(rows)}; errors={failures}', flush=True)
    return dict(records=len(rows), errors=failures, min_mapq=min_mapq,
                status='COMPLETE' if not failures else 'COMPLETE_WITH_INPUT_ERRORS')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--version', action='version', version=VERSION)
    parser.add_argument('--index', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--start-id', type=int, required=True)
    parser.add_argument('--end-id', type=int, required=True)
    parser.add_argument('--jobs', type=int, default=4)
    parser.add_argument('--min-mapq', type=int, default=0)
    args = parser.parse_args()
    try:
        result = audit(args.index, args.output, args.start_id, args.end_id, args.jobs, args.min_mapq)
    except (OSError, ValueError, sqlite3.Error) as error:
        parser.exit(1, f'ERROR: {error}\n')
    print(json.dumps(result, sort_keys=True))
    return 1 if result['errors'] else 0


if __name__ == '__main__':
    raise SystemExit(main())
