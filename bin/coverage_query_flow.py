#!/usr/bin/env python3
"""Plan, run and collect individual CoverageMap query tasks."""

import argparse
from collections import Counter
import csv
import hashlib
import json
from pathlib import Path
import subprocess
import sys

sys.path.insert(0, str(Path(__file__).absolute().parents[1]))
from bin.query_manifest import DISTRIBUTION_COLUMNS, read_manifest, require_v2_index, table
from bin.ref_depth_catalog import atomic_json, file_stat
from bin.jaccard_breadth import thresholds
from bin.query_decision_methods import COLUMNS as METHOD_COLUMNS, parse_methods
from bin.reference_index import ReferenceIndex


VERSION = '1.4.1'
OUTPUT_COLUMNS = (*DISTRIBUTION_COLUMNS, 'query_exit_code', 'result_source')


def query_metadata(row):
    return dict(query_id=row['query_id'], Sample=row.get('Sample') or '',
                GCFID=row.get('GCFID') or '', manifest_depth=row.get('Depth') or '')


def input_stamp(row):
    """Track manifest-declared files, including files that become readable on retry."""
    stamps = {}
    for name in ('CramPath', 'CraiPath', 'BedGraphPath'):
        path = row.get(name)
        if path:
            try:
                stamps[name] = file_stat(path)
            except (OSError, ValueError) as error:
                stamps[name] = dict(errno=getattr(error, 'errno', None))
    return hashlib.sha256(json.dumps(stamps, sort_keys=True).encode()).hexdigest()


def plan(args):
    index = ReferenceIndex(args.index)
    try:
        require_v2_index(index)
    finally:
        index.close()
    rows = read_manifest(args.manifest, args.coverage_unit)
    table(args.tasks, ('query_id', 'input_stamp'),
          (dict(query_id=row['query_id'], input_stamp=input_stamp(row)) for row in rows))
    atomic_json(args.roster, [query_metadata(row) for row in rows])


def read_result(path, query_id):
    with path.open(encoding='utf-8', newline='') as handle:
        reader = csv.DictReader(handle, delimiter='\t')
        rows = list(reader)
    if tuple(reader.fieldnames or ()) != DISTRIBUTION_COLUMNS or len(rows) != 1 or rows[0]['query_id'] != query_id:
        raise ValueError(f'Expected exactly one v2 summary row for {query_id}')
    return rows[0]


def validate_reference_draws(result, requested):
    """Reject successful output from a backend that ignored the requested N."""
    if result.get('status') == 'ERROR':
        return
    actual = result.get('requested_draws', '')
    if str(actual) != str(requested):
        raise ValueError(f'Reference draw mismatch: requested {requested}, backend reported {actual!r}')
    if result.get('status') == 'COMPLETE' and str(result.get('valid_draws', '')) != str(requested):
        raise ValueError(f'COMPLETE reference pool mismatch: requested {requested}, '
                         f"valid_draws={result.get('valid_draws', '')!r}")


def score(args):
    row = read_manifest(args.manifest, args.coverage_unit, args.query_id)[0]
    result = query_metadata(row)
    exit_code = 0
    result.update(result_source='TASK_RECEIPT', detection_status='NOT_EVALUATED')
    if args.stub:
        result.update(status='STUB', error='Stub run')
    elif input_stamp(row) != args.input_stamp:
        exit_code = 2
        result.update(status='ERROR', error='Query input identity changed since workflow planning')
    else:
        command = [
            str(Path(__file__).absolute().with_name('query_manifest_v2.py')),
            '--manifest', str(args.manifest), '--index', str(args.index),
            '--query-id', args.query_id, '--outdir', str(args.query_outdir),
            '--threads', str(args.threads), '--jobs', '1', '--alpha', str(args.alpha),
            '--seed', str(args.seed), '--min-samples', str(args.min_samples),
            '--depth-threshold', str(args.depth_threshold), '--cache-mib', str(args.cache_mib),
            '--work-mib', str(args.work_mib),
            '--reference-draws', str(getattr(args, 'reference_draws', 100)),
        ]
        command.extend(['--jaccard-coverage-threshold', str(args.jaccard_coverage_threshold)])
        command.extend(['--jaccard-methods', args.jaccard_methods,
                        '--bootstrap-replicates', str(args.bootstrap_replicates),
                        '--bootstrap-confidence', str(args.bootstrap_confidence)])
        if args.jaccard_coverage_thresholds:
            command.extend(['--jaccard-coverage-thresholds', args.jaccard_coverage_thresholds])
        if args.window_size is not None:
            command.extend(['--window-size', str(args.window_size)])
        else:
            command.extend(['--bin-count', str(args.bin_count)])
        for key in ('sampling_unit', 'coverage_alignment_policy'):
            if getattr(args, key, None):
                command.extend(['--'+key.replace('_', '-'), getattr(args, key)])
        if args.coverage_unit:
            command.extend(['--coverage-unit', args.coverage_unit])
        log = Path('query.stderr.log')
        with Path('query.stdout.log').open('wb') as stdout, log.open('wb') as stderr:
            exit_code = subprocess.run(command, stdout=stdout, stderr=stderr, check=False).returncode
        try:
            result.update(read_result(args.query_outdir / 'summary.tsv', args.query_id))
            if not exit_code:
                validate_reference_draws(result, args.reference_draws)
        except (OSError, ValueError) as error:
            exit_code = exit_code or 2
            result.update(status='ERROR', error=str(error))
        if exit_code:
            with log.open('rb') as handle:
                handle.seek(max(0, log.stat().st_size - 4096))
                stderr_tail = handle.read().decode(errors='replace').strip()
            result.update(status='ERROR', detection_status='NOT_EVALUATED',
                          error=' | '.join(filter(None, (result.get('error'), stderr_tail))))
    sensitivity = args.query_outdir / 'sensitivity.tsv'
    if not args.stub and sensitivity.exists():
        with sensitivity.open() as handle:
            result['threshold_results'] = list(csv.DictReader(handle, delimiter='\t'))
    else:
        result['threshold_results'] = [dict(result, jaccard_coverage_threshold=tau)
            for tau in thresholds(args.jaccard_coverage_threshold, args.jaccard_coverage_thresholds)]
    method_file = args.query_outdir / 'decision_methods.tsv'
    if not args.stub and method_file.is_file():
        with method_file.open() as handle:
            reader = csv.DictReader(handle, delimiter='\t')
            if tuple(reader.fieldnames or ()) != METHOD_COLUMNS:
                raise ValueError('Unexpected query decision method schema')
            method_rows = list(reader)
            expected_methods = set(parse_methods(args.jaccard_methods))
            if (len(method_rows) != len(expected_methods)
                    or {r.get('method') for r in method_rows} != expected_methods
                    or any(r.get('query_id') != args.query_id for r in method_rows)):
                raise ValueError('Query decision method rows differ from requested methods/ID')
            result['method_results'] = method_rows
    else:
        result['method_results'] = []
    if result.get('status') == 'ERROR' and not exit_code:
        exit_code = 2
    result['query_exit_code'] = 128 - exit_code if exit_code < 0 else exit_code
    atomic_json(args.receipt, result)
    if result.get('status') == 'ERROR' or exit_code:
        # Preserve diagnostics, but never turn a failed scorer into a successful
        # Nextflow cache entry. A subprocess signal uses the shell's 128+signal.
        status = 128 - exit_code if exit_code < 0 else exit_code or 2
        print(f"ERROR: {args.query_id}: {result.get('error', 'Query task failed')}; "
              f"exit={status}; receipt={args.receipt}; inspect query.stderr.log", file=sys.stderr)
        return status
    return 0


def collect(args):
    roster = json.loads(args.roster.read_text(encoding='utf-8'))
    expected = {row['query_id']: row for row in roster}
    if not roster or len(expected) != len(roster):
        raise ValueError('The query roster must contain unique query IDs')
    results = {}
    threshold_results = {}
    method_results = {}
    for path in sorted(args.receipts.glob('*.json')):
        row = json.loads(path.read_text(encoding='utf-8'))
        query_id = row['query_id']
        if query_id not in expected or query_id in results or path.stem != query_id:
            raise ValueError(f'Unexpected or duplicate task receipt: {path}')
        if any(row[key] != value for key, value in expected[query_id].items()):
            raise ValueError(f'Task receipt identity differs from roster: {query_id}')
        threshold_results[query_id] = row.pop('threshold_results', [])
        method_results[query_id] = row.pop('method_results', [])
        candidates = method_results[query_id]
        if candidates and (len(candidates) != len(parse_methods(args.jaccard_methods))
                or {r.get('method') for r in candidates} != set(parse_methods(args.jaccard_methods))
                or any(r.get('query_id') != query_id or r.get('Sample') != row.get('Sample')
                       or r.get('GCFID') != row.get('GCFID') for r in candidates)):
            raise ValueError(f'Unexpected method receipt identity for {query_id}')
        results[query_id] = row
    buckets = dict(negative=[], errors=[], unscored=[])
    counts = Counter()
    rows = []
    for query in roster:
        query_id = query['query_id']
        row = results.get(query_id)
        if row is None:
            row = dict(query, status='ERROR', detection_status='NOT_EVALUATED', result_source='NO_RECEIPT',
                       error='No task receipt after execution; inspect Nextflow trace and task logs')
        if row['status'] == 'ERROR':
            buckets['errors'].append(row)
            counts['ERROR'] += 1
        elif row['detection_status'] == 'NOT_DETECTED':
            buckets['negative'].append(row)
            counts['NOT_DETECTED'] += 1
        elif row['detection_status'] == 'NOT_EVALUATED':
            buckets['unscored'].append(row)
            counts['NOT_EVALUATED'] += 1
        elif row['detection_status'] == 'DETECTED':
            counts['DETECTED'] += 1
        else:
            raise ValueError(f'Unknown detection status for {query_id}')
        rows.append(row)
    args.outdir.mkdir(parents=True, exist_ok=False)
    table(args.outdir / 'summary.tsv', OUTPUT_COLUMNS, rows)
    def threshold_rows():
        for row in rows:
            candidates = threshold_results.get(row['query_id']) or [dict(row, jaccard_coverage_threshold=tau)
                for tau in thresholds(args.jaccard_coverage_threshold, args.jaccard_coverage_thresholds)]
            for candidate in candidates:
                candidate.update(query_exit_code=row.get('query_exit_code',''), result_source=row['result_source'])
                if row.get('status') == 'ERROR':
                    candidate.update(status='ERROR', detection_status='NOT_EVALUATED', error=row.get('error',''))
                yield candidate
    table(args.outdir / 'sensitivity.tsv', OUTPUT_COLUMNS, threshold_rows())
    def method_rows():
        for row in rows:
            candidates = method_results.get(row['query_id']) or [
                dict(query_id=row['query_id'], Sample=row.get('Sample',''), GCFID=row.get('GCFID',''),
                     status=row.get('status',''), method=method, detection_status='NOT_EVALUATED',
                     jaccard_decision='NOT_EVALUATED', pearson_decision='NOT_EVALUATED',
                     method_status='missing_query_method_result')
                for method in parse_methods(args.jaccard_methods)]
            for entry in candidates:
                if row.get('status') == 'ERROR':
                    entry.update(status='ERROR', detection_status='NOT_EVALUATED',
                                 jaccard_decision='NOT_EVALUATED', pearson_decision='NOT_EVALUATED',
                                 method_status='query_task_error')
                yield {name: entry.get(name) for name in METHOD_COLUMNS}
    table(args.outdir / 'decision_methods.tsv', METHOD_COLUMNS, method_rows())
    for label, subset in buckets.items():
        table(args.outdir / f'{label}.tsv', OUTPUT_COLUMNS, subset)
    atomic_json(args.outdir / 'completion.json', dict(
        queries=len(roster), receipts=len(results), missing_receipts=len(roster)-len(results),
        completion_status='PARTIAL_FAILURE' if counts['ERROR'] else 'COMPLETE',
        qc_status='EXECUTION_ERRORS' if counts['ERROR'] else
                  'HAS_UNRESOLVED' if counts['NOT_EVALUATED'] else 'COMPLETE',
        **{key: counts[key] for key in ('DETECTED', 'NOT_DETECTED', 'NOT_EVALUATED', 'ERROR')}))
    if counts['ERROR']:
        print(f"WARNING: {counts['ERROR']}/{len(roster)} queries failed; "
              f"partial results retained; inspect {args.outdir / 'errors.tsv'} "
              "and Nextflow trace/task logs before retrying", file=sys.stderr)


def parser():
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument('--version', action='version', version=VERSION)
    sub = cli.add_subparsers(dest='stage', required=True)
    planning = sub.add_parser('plan')
    scoring = sub.add_parser('score')
    collecting = sub.add_parser('collect')
    for stage in (planning, scoring):
        stage.add_argument('--manifest', type=Path, required=True)
        stage.add_argument('--coverage-unit', choices=('fraction', 'percent'))
    planning.add_argument('--tasks', type=Path, required=True)
    planning.add_argument('--index', type=Path, required=True)
    planning.add_argument('--roster', type=Path, required=True)
    for key in ('index', 'query-outdir', 'receipt'):
        scoring.add_argument('--'+key, type=Path, required=True)
    for key in ('query-id', 'input-stamp'):
        scoring.add_argument('--'+key, required=True)
    for key, kind in (('alpha', float), ('threads', int), ('seed', int), ('min-samples', int),
                      ('depth-threshold', float), ('cache-mib', int), ('work-mib', int)):
        scoring.add_argument('--'+key, type=kind, required=True)
    windows = scoring.add_mutually_exclusive_group(required=True)
    windows.add_argument('--window-size', type=int)
    windows.add_argument('--bin-count', type=int)
    scoring.add_argument('--sampling-unit', choices=('read', 'fragment'))
    scoring.add_argument('--coverage-alignment-policy', choices=('primary_only', 'primary_plus_supplementary'))
    scoring.add_argument('--reference-draws', type=int, default=100)
    scoring.add_argument('--bootstrap-replicates', type=int, default=1000)
    scoring.add_argument('--bootstrap-confidence', type=float, default=.95)
    for stage in (scoring, collecting):
        stage.add_argument('--jaccard-methods', default='normal')
        stage.add_argument('--jaccard-coverage-threshold',type=float,default=0.)
        stage.add_argument('--jaccard-coverage-thresholds')
    scoring.add_argument('--stub', action='store_true')
    for key in ('roster', 'receipts', 'outdir'):
        collecting.add_argument('--'+key, type=Path, required=True)
    planning.set_defaults(action=plan)
    scoring.set_defaults(action=score)
    collecting.set_defaults(action=collect)
    return cli


def main():
    cli = parser()
    args = cli.parse_args()
    try:
        status = args.action(args)
    except (ValueError, OSError, KeyError) as error:
        cli.exit(2, f'ERROR: {error}\n')
    raise SystemExit(status or 0)


if __name__ == '__main__':
    main()
