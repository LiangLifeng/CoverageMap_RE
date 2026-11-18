#!/usr/bin/env python3
"""Render CoverageMap v2 similarity diagnostics from saved scores and reference means."""

import argparse
from collections import Counter
from concurrent.futures import ProcessPoolExecutor
import csv
import gzip
import hashlib
import json
import math
import os
from pathlib import Path
import re
from statistics import fmean, pstdev
from tempfile import NamedTemporaryFile
from textwrap import shorten

VERSION = '1.3.0'
METRICS = {'jaccard': ('Jaccard similarity', 0.0), 'pearson': ('Pearson r', -1.0)}
IDENTITY = ('query_id', 'Sample', 'GCFID')
METRIC_FIELDS = ('mean', 'p', 'null_mean', 'null_sd', 'null_count', 'status', 'decision')
PLOT_FIELDS = (*IDENTITY, 'status', 'method', 'detection_status', 'plot_status', 'reason', 'plot_path',
               'plot_sha256', 'null_sha256', 'jaccard_reference_n', 'pearson_reference_n')
QUERY_ID = re.compile(r'q[0-9]{8,}')
TRACE_NAME = re.compile(r'(?:^|:)(?:RAW_)?QUERY_SCORE \((q[0-9]{8,})\)$')


def sha(path):
    with path.open('rb') as handle:
        return hashlib.file_digest(handle, 'sha256').hexdigest()


def read_table(path, fields, compressed=False):
    opener = gzip.open if compressed else open
    with opener(path, 'rt', encoding='utf-8', newline='') as handle:
        reader = csv.DictReader(handle, delimiter='\t', strict=True)
        if (not reader.fieldnames or len(set(reader.fieldnames)) != len(reader.fieldnames)
                or not set(fields) <= set(reader.fieldnames)):
            raise ValueError(f'Missing or duplicate TSV columns in {path}')
        rows = list(reader)
    if any(None in row or None in row.values() for row in rows):
        raise ValueError(f'Malformed TSV row in {path}')
    return rows


def trace_nulls(work, trace):
    """Resolve an explicit Nextflow trace and reject ambiguous successful tasks."""
    matches = {}
    for row in read_table(trace, ('name', 'hash', 'status')):
        name = TRACE_NAME.search(row['name'])
        if name is None or row['status'] not in ('COMPLETED', 'CACHED'):
            continue
        query_id = name.group(1)
        task_hash = row['hash']
        if not re.fullmatch(r'[a-f0-9]{2}/[a-f0-9]{6,30}', task_hash):
            raise ValueError(f'Invalid Nextflow task hash for {query_id}')
        prefix, stem = task_hash.split('/')
        paths = sorted((work/prefix).glob(stem+'*'))
        if len(paths) > 1 or query_id in matches:
            raise ValueError(f'Ambiguous completed Nextflow task for {query_id}')
        if paths:
            matches[query_id] = paths[0]/'query_result'
    return matches


def null_sources(args, rows):
    if args.work_dir is not None:
        if args.trace is None:
            raise ValueError('--work-dir requires --trace from the same run')
        matches = trace_nulls(args.work_dir, args.trace)
        return {qid: (path/'queries'/qid/'null.tsv.gz', path/'summary.tsv')
                for qid, path in matches.items()}
    if args.trace is not None:
        raise ValueError('--trace requires --work-dir')
    if args.null_root is not None:
        return {row['query_id']: (args.null_root/row['query_id']/'null.tsv.gz', None) for row in rows}
    if args.staged_ids is None:
        raise ValueError('--staged-nulls requires --staged-ids')
    ids = args.staged_ids.read_text().split()
    paths = sorted([*args.staged_nulls.glob('input*/null.tsv.gz'),
                    *args.staged_nulls.glob('input*/q*/null.tsv.gz')])
    if (len(paths) != len(ids) or len(set(ids)) != len(ids)
            or any(not QUERY_ID.fullmatch(qid) for qid in ids)):
        raise ValueError('Staged null files do not match the query ID list')
    return {qid: (path, None) for qid, path in zip(ids, paths, strict=True)}


def number(value, minimum, maximum):
    if value in ('', None):
        return None
    value = float(value)
    if not math.isfinite(value) or value < minimum - 1e-12 or value > maximum + 1e-12:
        raise ValueError('Numeric value outside its expected range')
    return value


def load_nulls(path, expected=100):
    result = {metric: {} for metric in METRICS}
    for row in read_table(path, ('metric', 'left_draw', 'right_draw', 'value'), compressed=True):
        metric = row['metric']
        if metric not in METRICS or row['right_draw']:
            raise ValueError('Expected per-draw reference means rather than pairwise similarities')
        slot = int(row['left_draw'])
        if not 0 <= slot < expected or slot in result[metric]:
            raise ValueError('Reference draw slot is invalid or duplicated')
        result[metric][slot] = number(row['value'], METRICS[metric][1], 1.0)
    return {metric: [slots[i] for i in sorted(slots) if slots[i] is not None]
            for metric, slots in result.items()}


def plot_data(row, source):
    path, task_summary = source
    if task_summary is not None and task_summary.is_file():
        saved, = read_table(task_summary, IDENTITY)
        keys = (*IDENTITY, 'actual_query_coverage', 'gcf_reference_length', 'windows', 'arm',
                *(f'{m}_{f}' for m in METRICS for f in METRIC_FIELDS))
        if any(saved.get(key, '') != row.get(key, '') for key in keys):
            raise ValueError('Collected and task summary values differ')
    have_null = path is not None and path.is_file()
    expected = int(row.get('requested_draws') or 100)
    references = load_nulls(path, expected) if have_null else {metric: [] for metric in METRICS}
    panels, reasons = {}, []
    for metric, (_, lower) in METRICS.items():
        values = references[metric]
        observed = number(row.get(metric+'_mean'), lower, 1.0)
        p = number(row.get(metric+'_p'), 0.0, 1.0)
        mean = number(row.get(metric+'_null_mean'), lower, 1.0)
        sd = number(row.get(metric+'_null_sd'), 0.0, 2.0)
        state = row.get(metric+'_status') or 'not_evaluated'
        count = row.get(metric+'_null_count')
        if have_null and count not in ('', None) and int(count) != len(values):
            reasons.append(f'{metric}: incomplete reference data ({len(values)}/{count} saved)')
        elif have_null and values:
            for actual, moment in ((fmean(values), mean), (pstdev(values), sd)):
                if moment is not None and not math.isclose(actual, moment, abs_tol=1e-12, rel_tol=1e-9):
                    raise ValueError(f'{metric} saved reference moments differ from summary')
        if len(values) != expected:
            reasons.append(f'{metric}: {len(values)}/{expected} reference means available')
        if state != 'estimable':
            reasons.append(f'{metric}: {state}')
        if observed is None:
            reasons.append(f'{metric}: query mean unavailable')
        if p is None:
            reasons.append(f'{metric}: p unavailable')
        panels[metric] = dict(values=values, observed=observed, p=p, mean=mean, sd=sd,
                              state=state, decision=row.get(metric+'_decision') or 'NOT_EVALUATED')
    if not have_null:
        reasons.append('null.tsv.gz missing')
    if row['status'] != 'COMPLETE':
        reasons.append(row.get('error') or row['status'])
    return panels, reasons, sha(path) if have_null else ''


def p_label(p):
    if p is None:
        return 'p unavailable'
    if p == 0:
        return 'p below float range'
    return f'p = {p:.3e}' if p < .001 else f'p = {p:.4g}'


def value_label(value):
    return 'NA' if value is None else f'{value:.5g}'


def depth_label(value):
    try:
        return value_label(number(value, 0, math.inf))
    except ValueError:
        return 'NA'


def coverage_data(row, source):
    """Read saved union-covered bases, never the thresholded occupied-bin count."""
    import numpy as np

    expected = int(row.get('requested_draws') or 100)
    query = number(row.get('actual_query_coverage'), 0, 1)
    result = dict(query=query, reference_mean=None, reference_min=None,
                  reference_max=None, reference_n=0, reference_expected=expected,
                  definition='unique covered bases / full GCF reference length')
    reasons = [] if query is not None else ['query coverage breadth unavailable']
    path = source[0]
    names = ('reference_covered_bases.npy', 'window_lengths.npy', 'draws.tsv.gz')
    if path is None or not all((path.parent/name).is_file() for name in names):
        return result, reasons + ['reference coverage breadth unavailable (saved arrays/draws missing)']
    length = number(row.get('gcf_reference_length'), 1, math.inf)
    if length is None or not length.is_integer():
        raise ValueError('Reference coverage breadth requires full GCF reference length')
    draws = read_table(path.parent/'draws.tsv.gz', (*IDENTITY, 'draw_id', 'status'), compressed=True)
    if any(any(d[key] != str(row[key]) for key in IDENTITY) for d in draws):
        raise ValueError('Coverage draw identity differs from summary')
    slots = [int(d['draw_id']) for d in draws]
    if len(set(slots)) != len(slots) or any(not 0 <= i < expected for i in slots):
        raise ValueError('Coverage draw slot is invalid or duplicated')
    arrays = []
    try:
        for name in names[:2]:
            arrays.append(np.load(path.parent/name, mmap_mode='r', allow_pickle=False))
        covered, lengths = arrays
        if (lengths.ndim != 1 or not lengths.size or lengths.dtype.kind not in 'iu'
                or np.any(lengths <= 0) or int(lengths.sum(dtype=np.uint64)) != int(length)
                or covered.ndim != 2 or covered.shape[1] != lengths.size
                or covered.dtype.kind not in 'iu' or covered.shape[0] != len(draws)
                or sorted(slots) != list(range(len(draws)))):
            raise ValueError('Coverage arrays differ from saved draw/window layout')
        values = []
        for draw in draws:
            if draw['status'] != 'READY':
                continue
            counts = covered[int(draw['draw_id'])]
            if np.any(counts < 0) or np.any(counts > lengths):
                raise ValueError('Reference covered bases outside window lengths')
            values.append(int(counts.sum(dtype=np.uint64))/length)
        if values:
            result.update(reference_mean=fmean(values), reference_min=min(values),
                          reference_max=max(values), reference_n=len(values))
        if len(values) != expected:
            reasons.append(f'reference coverage breadth: {len(values)}/{expected} READY draws')
    finally:
        for array in arrays:
            array._mmap.close()
    return result, reasons


def coverage_label(data):
    def percent(value):
        return 'NA' if value is None else f'{100*value:.2f}%'
    return (f"Coverage breadth (covered bp / full GCF bp): Query {percent(data['query'])}"
            f"  |  Ref mean {percent(data['reference_mean'])}"
            f" [min {percent(data['reference_min'])}, max {percent(data['reference_max'])}]"
            f"; n={data['reference_n']}/{data.get('reference_expected', 100)}")


def render(row, panels, reasons, target, coverage=None):
    from matplotlib.backends.backend_agg import FigureCanvasAgg
    from matplotlib.figure import Figure
    from matplotlib.lines import Line2D

    method = row.get('method', 'normal')
    blue, red = '#237cb5', '#c63e33'
    figure = Figure(figsize=(10, 4.8), dpi=140, facecolor='white')
    FigureCanvasAgg(figure)
    grid = figure.add_gridspec(2, 2, height_ratios=(4, 1), left=.09, right=.97,
                               top=.65, bottom=.23, wspace=.27, hspace=.07)
    figure.text(.5, .956, f"{row['query_id']}  |  {row['Sample']}  |  {row['GCFID']}",
                ha='center', va='top', fontsize=12, fontweight='bold', parse_math=False)
    window = (f"B={row['requested_bin_count']}" if row.get('requested_bin_count')
              else f"{row['nominal_window_size']} bp" if row.get('nominal_window_size') else 'windows unavailable')
    arm = {'low': 'matched-depth simulation', 'high': 'original references'}.get(row.get('arm'), 'route unavailable')
    actual_depth = depth_label(row.get('actual_query_depth'))
    manifest_depth = depth_label(row.get('manifest_depth'))
    figure.text(.5, .895, f"Depth {actual_depth}x"
                f" (manifest {manifest_depth}x)  |  {arm}  |  {window}"
                f"  |  {row.get('windows') or 'NA'} bins  |  alpha={row.get('alpha') or 'NA'}",
                ha='center', va='top', fontsize=9)
    figure.text(.5, .845, f"{method}  |  {row['detection_status']}  |  {row.get('window_policy') or 'saved legacy windows'}",
                ha='center', va='top', fontsize=9, color='#333333')
    coverage = coverage or dict(query=None, reference_mean=None, reference_min=None,
                                reference_max=None, reference_n=0)
    figure.text(.5, .787, coverage_label(coverage), ha='center', va='top', fontsize=8.5)
    for col, (metric, (label, lower)) in enumerate(METRICS.items()):
        data = panels[metric]
        values = data['values']
        observed = data['observed']
        ax = figure.add_subplot(grid[0, col])
        strip = figure.add_subplot(grid[1, col], sharex=ax)
        if values:
            if min(values) == max(values):
                ax.axvline(values[0], color=blue, linewidth=2)
                ax.set_yticks([])
            else:
                ax.hist(values, bins=min(18, max(5, round(math.sqrt(len(values))))),
                        density=True, color='#a8d4ec', edgecolor=blue, linewidth=.9)
            y = [.22 + ((i*37) % 101)/250 for i in range(len(values))]
            strip.scatter(values, y, s=11, color=blue, alpha=.54, linewidths=0, clip_on=False)
        else:
            ax.text(.5, .48, 'Reference means unavailable', transform=ax.transAxes,
                    ha='center', va='center', color='#555555', fontsize=9)
            ax.set_yticks([])
        if observed is not None:
            strip.scatter([observed], [.42], s=49, marker='D', color=red,
                          edgecolors='white', linewidths=.5, zorder=4, clip_on=False)
        plotted = values + ([observed] if observed is not None else [])
        if plotted:
            span = max(max(plotted)-min(plotted), .012)
            strip.set_xlim(max(lower, min(plotted)-.1*span), min(1, max(plotted)+.1*span))
        else:
            strip.set_xlim(lower, 1)
        strip.set_ylim(0, 1)
        strip.set_yticks([])
        strip.set_xlabel(label, fontsize=9)
        strip.tick_params(axis='x', labelsize=8)
        ax.tick_params(axis='x', labelbottom=False, bottom=False)
        ax.tick_params(axis='y', labelsize=8)
        if col == 0:
            ax.set_ylabel('Density', fontsize=9)
        criterion = p_label(data['p'])
        if metric == 'jaccard' and method != 'normal':
            criterion = f"score={value_label(observed)}, threshold={value_label(data['limit'])}"
            if method == 'bootstrap':
                criterion = f"score={value_label(observed)}, CI=[{value_label(data['ci_low'])}, {value_label(data['ci_high'])}]"
        ax.set_title(f"{label}  |  {criterion}  |  {data['decision']}\n"
                     f"Ref. n={len(values)}/{row.get('requested_draws') or 100}, mean={value_label(data['mean'])}, SD={value_label(data['sd'])}",
                     fontsize=9, loc='left', pad=6)
        if data['state'] != 'estimable':
            ax.text(.98, .94, data['state'], transform=ax.transAxes, ha='right', va='top',
                    fontsize=8, color='#555555', bbox=dict(facecolor='white', alpha=.85, edgecolor='none'))
    figure.legend(handles=[
        Line2D([], [], linestyle='none', marker='o', color=blue, label='Reference draw means'),
        Line2D([], [], linestyle='none', marker='D', color=red, label='Query mean')],
        loc='lower center', bbox_to_anchor=(.5, .079), ncol=2, frameon=False, fontsize=9)
    if reasons:
        figure.text(.5, .036, shorten(' | '.join(reasons), width=150, placeholder=' ...'),
                    ha='center', fontsize=8, color='#555555', parse_math=False)
    metadata = dict(Software=f'CoverageMap diagnostics {VERSION}',
                    Description=json.dumps(dict(summary=row, panels=panels, coverage=coverage), sort_keys=True, allow_nan=False))
    with NamedTemporaryFile(dir=target.parent, suffix='.png', delete=False) as temp:
        path = Path(temp.name)
    try:
        figure.savefig(path, format='png', metadata=metadata)
        os.replace(path, target)
    finally:
        path.unlink(missing_ok=True)
        figure.clear()


def plot_one(item):
    row, source, output, stub, method, method_row, subdir = item
    result = {key: row[key] for key in (*IDENTITY, 'status', 'detection_status')}
    result.update(method=method)
    if method_row is not None:
        result['detection_status'] = method_row['detection_status']
    if stub:
        return dict(result, plot_status='STUB', reason='Stub run')
    try:
        panels, reasons, null_sha = plot_data(row, source)
        row = dict(row, method=method, detection_status=result['detection_status'])
        if method_row is not None:
            jaccard = panels['jaccard']
            jaccard.update(observed=number(method_row['jaccard_score'], 0, 1),
                           decision=method_row['jaccard_decision'],
                           limit=number(method_row['reference_lower_limit'], -math.inf, math.inf),
                           ci_low=number(method_row['bootstrap_lower_ci'], 0, 1),
                           ci_high=number(method_row['bootstrap_upper_ci'], 0, 1))
            row.update(jaccard_decision=method_row['jaccard_decision'],
                       pearson_decision=method_row['pearson_decision'])
            if method != 'normal':
                row.update(jaccard_p=None, jaccard_mean=method_row['jaccard_score'])
                jaccard['p'] = None
                jaccard['state'] = method_row['method_status'] if jaccard['decision'] == 'NOT_EVALUATED' else 'estimable'
                reasons = [reason for reason in reasons if reason != 'jaccard: p unavailable'
                           and reason != f"jaccard: {row.get('jaccard_status')}"]
            panels['pearson'].update(p=number(method_row['pearson_p'], 0, 1),
                                     decision=method_row['pearson_decision'])
        target = output/subdir/(row['query_id']+'_similarity.png')
        target.parent.mkdir(parents=True, exist_ok=True)
        coverage, coverage_reasons = coverage_data(row, source)
        reasons.extend(coverage_reasons)
        render(row, panels, reasons, target, coverage)
        return dict(result, plot_status='PLOTTED_WITH_GAPS' if reasons else 'PLOTTED',
                    reason=' | '.join(reasons), plot_path=str(target.relative_to(output)), plot_sha256=sha(target),
                    null_sha256=null_sha, **{m+'_reference_n': len(panels[m]['values']) for m in METRICS})
    except Exception as error:
        return dict(result, plot_status='ERROR', reason=f'{type(error).__name__}: {error}')


def run(args):
    required = (*IDENTITY, 'status', 'detection_status',
                *(f'{m}_{f}' for m in METRICS for f in METRIC_FIELDS))
    rows = read_table(args.summary, required)
    ids = [r['query_id'] for r in rows]
    if not ids or len(set(ids)) != len(ids) or any(not QUERY_ID.fullmatch(q) for q in ids):
        raise ValueError('Summary must contain unique, valid query IDs')
    if args.jobs < 1:
        raise ValueError('--jobs must be positive')
    if args.query_id:
        if not set(args.query_id) <= set(ids):
            raise ValueError('Requested query ID absent from the summary')
        rows = [r for r in rows if r['query_id'] in set(args.query_id)]
    sources = null_sources(args, rows)
    args.outdir.mkdir(parents=True, exist_ok=True)
    methods = [name.strip() for name in getattr(args, 'methods', 'normal').split(',')]
    if len(set(methods)) != len(methods) or not set(methods) <= {'normal', 'minimum', 'empirical', 'bootstrap'}:
        raise ValueError('Invalid or duplicate diagnostic methods')
    method_path = getattr(args, 'decision_methods', None)
    methods_by_key = {}
    if method_path is not None:
        needed = {r['query_id']: r for r in rows}
        fields = (*IDENTITY, 'method', 'jaccard_score', 'reference_lower_limit', 'bootstrap_lower_ci',
                  'bootstrap_upper_ci', 'jaccard_decision', 'pearson_p', 'pearson_decision',
                  'detection_status', 'method_status')
        for entry in read_table(method_path, fields):
            q, method = entry['query_id'], entry['method']
            if q not in needed or method not in methods:
                continue
            if (q, method) in methods_by_key or any(entry[k] != needed[q][k] for k in IDENTITY):
                raise ValueError(f'Duplicate or mismatched method result: {q}/{method}')
            methods_by_key[q, method] = entry
        if len(methods_by_key) != len(rows)*len(methods):
            raise ValueError('Missing diagnostic method results')
    elif methods != ['normal']:
        raise ValueError('Non-normal diagnostics require --decision-methods')
    tasks = [(r, sources.get(r['query_id'], (None, None)), args.outdir, args.stub, method,
              methods_by_key.get((r['query_id'], method)), method if len(methods) > 1 else '')
             for r in rows for method in methods]
    if args.jobs > 1 and len(tasks) > 1:
        with ProcessPoolExecutor(max_workers=min(args.jobs, len(tasks))) as executor:
            results = list(executor.map(plot_one, tasks))
    else:
        results = [plot_one(task) for task in tasks]
    with NamedTemporaryFile(dir=args.outdir, mode='w', encoding='utf-8', newline='', delete=False) as handle:
        writer = csv.DictWriter(handle, PLOT_FIELDS, delimiter='\t', lineterminator='\n')
        writer.writeheader()
        writer.writerows(results)
    os.replace(handle.name, args.outdir/'plot_index.tsv')
    counts = dict(Counter(row['plot_status'] for row in results))
    print(json.dumps(counts, sort_keys=True))
    return 1 if counts.get('ERROR') else 0


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--version', action='version', version=VERSION)
    parser.add_argument('--summary', required=True, type=Path)
    parser.add_argument('--decision-methods', type=Path)
    parser.add_argument('--methods', default='normal')
    parser.add_argument('--outdir', required=True, type=Path)
    sources = parser.add_mutually_exclusive_group(required=True)
    sources.add_argument('--null-root', type=Path, help='Directory with <query_id>/null.tsv.gz')
    sources.add_argument('--work-dir', type=Path, help='Nextflow work directory, with --trace')
    sources.add_argument('--staged-nulls', type=Path, help=argparse.SUPPRESS)
    parser.add_argument('--staged-ids', type=Path, help=argparse.SUPPRESS)
    parser.add_argument('--trace', type=Path)
    parser.add_argument('--query-id', action='append')
    parser.add_argument('--jobs', type=int, default=1)
    parser.add_argument('--stub', action='store_true', help=argparse.SUPPRESS)
    args = parser.parse_args()
    try:
        return run(args)
    except (ValueError, OSError, KeyError, csv.Error) as error:
        parser.exit(2, f'ERROR: {error}\n')


if __name__ == '__main__':
    raise SystemExit(main())
