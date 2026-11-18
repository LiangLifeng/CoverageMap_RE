#!/usr/bin/env python3
"""Exploratory Jaccard decision methods on the frozen N-draw reference null.

The original normal-tail decision in summary.tsv remains authoritative.
Here all methods use the SAME query mean and N reference row means. Repeated
source IDs are deliberately retained, not treated as independent observations.
Prediction/Bootstrap labels are descriptive, not calibrated error guarantees.
"""
import argparse
from collections import defaultdict
import csv
import gzip
import hashlib
import json
import math
from pathlib import Path
from statistics import NormalDist

import numpy as np


METHODS = ('normal', 'empirical', 'minimum', 'bootstrap')
COLUMNS = ('query_id', 'Sample', 'GCFID', 'status', 'arm', 'method', 'alpha',
           'null_count', 'unique_references', 'jaccard_score', 'reference_lower_limit',
           'bootstrap_lower_ci', 'bootstrap_upper_ci', 'bootstrap_pass_fraction',
           'jaccard_decision', 'pearson_p', 'pearson_decision',
           'detection_status', 'method_status')


def parse_methods(spec):
    names = tuple(name.strip() for name in spec.split(','))
    if not names or any(name not in METHODS for name in names) or len(set(names)) != len(names):
        raise ValueError('Jaccard methods must be distinct names from: ' + ','.join(METHODS))
    return names


def validate_options(alpha, bootstrap_replicates, bootstrap_confidence):
    if not math.isfinite(alpha) or not 0 < alpha < 1:
        raise ValueError('alpha must be a finite probability in (0,1)')
    if (type(bootstrap_replicates) is not int or not 10 <= bootstrap_replicates <= 10000
            or not math.isfinite(bootstrap_confidence) or not 0 < bootstrap_confidence < 1):
        raise ValueError('Bootstrap requires 10..10000 replicates and confidence in (0,1)')


def _number(value):
    if value is None or value == '':
        return None
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    return number if math.isfinite(number) else None


def _decision(p, alpha):
    return 'NOT_EVALUATED' if p is None else 'PASS' if p >= alpha else 'FAIL'


def _combine(jaccard, pearson):
    if jaccard in ('NOT_EVALUATED', 'UNCERTAIN') or pearson == 'NOT_EVALUATED':
        return 'NOT_EVALUATED'
    return 'DETECTED' if jaccard == pearson == 'PASS' else 'NOT_DETECTED'


def empirical_lower(values, alpha):
    """Noninterpolated observed lower rank, ceil(alpha*n), 1-based."""
    if not len(values):
        raise ValueError('An empirical limit needs nonempty scores')
    return float(np.sort(values)[max(0, math.ceil(alpha * len(values)) - 1)])


def _null_inputs(summary, draws, nulls):
    """Reject incomplete or undefined draw/null sets, without shrinking denominators."""
    expected = int(summary.get('requested_draws') or 0)
    if expected < 2 or len(draws) != expected:
        return None, None, 'incomplete_draw_slots'
    try:
        draw_by_id = {int(row['draw_id']): row for row in draws}
    except (KeyError, TypeError, ValueError):
        return None, None, 'invalid_draw_identity'
    if set(draw_by_id) != set(range(expected)) or any(
            row['status'] != 'READY' or not row.get('reference_id')
            for row in draw_by_id.values()):
        return None, None, 'incomplete_or_invalid_draws'
    rows = [row for row in nulls if row.get('metric') == 'jaccard']
    if len(rows) != expected:
        return None, None, 'incomplete_reference_null'
    try:
        by_id = {int(row['left_draw']): row for row in rows}
    except (KeyError, TypeError, ValueError):
        return None, None, 'invalid_null_identity'
    if set(by_id) != set(range(expected)):
        return None, None, 'duplicate_or_missing_null_rows'
    values = [_number(by_id[i].get('value')) for i in range(expected)]
    if any(value is None for value in values):
        return None, None, 'nonestimable_reference_null'
    return np.asarray(values, dtype=np.float64), [
        str(draw_by_id[i]['reference_id']) for i in range(expected)], 'estimable'


def _bootstrap_limits(values, reference_ids, alpha, replicates, confidence, seed, query_id):
    """Source-cluster bootstrap conditional on the saved draw/null statistics.

    Resample distinct reference source groups, retaining all their draw slots.
    This is a stability diagnostic, NOT a fresh independent prediction interval.
    """
    clusters = defaultdict(list)
    for value, reference_id in zip(values, reference_ids):
        clusters[reference_id].append(float(value))
    if len(clusters) < 2:
        return None
    groups = [np.asarray(values, dtype=np.float64) for values in clusters.values()]
    token = hashlib.sha256(json.dumps(
        [seed, query_id, 'reference-cluster-bootstrap-v1', alpha, replicates, confidence],
        separators=(',', ':')).encode()).digest()
    rng = np.random.default_rng(np.frombuffer(token, dtype='<u4'))
    cutoffs = np.empty(replicates, dtype=np.float64)
    for i in range(replicates):
        sampled = rng.integers(0, len(groups), size=len(groups))
        cutoffs[i] = empirical_lower(np.concatenate([groups[j] for j in sampled]), alpha)
    tail = (1 - confidence) / 2
    return (float(np.quantile(cutoffs, tail)),
            float(np.quantile(cutoffs, 1 - tail)), cutoffs)


def evaluate(summary, draws, nulls, *, methods='normal', alpha=.05,
             bootstrap_replicates=1000, bootstrap_confidence=.95, seed=401):
    """One row per method. No resampling of alignments or references occurs."""
    names = parse_methods(methods)
    validate_options(alpha, bootstrap_replicates, bootstrap_confidence)
    score = _number(summary.get('jaccard'))
    pearson_p = _number(summary.get('pearson_p'))
    pearson = _decision(pearson_p, alpha)
    values, refs, null_state = _null_inputs(summary, draws, nulls)
    complete = (values is not None and score is not None
                and summary.get('status') not in ('ERROR', 'ZERO_QUERY_DEPTH'))
    if not complete and null_state == 'estimable':
        null_state = 'nonestimable_query'
    unique = len(set(refs)) if refs is not None else 0
    output = []
    for method in names:
        row = dict(query_id=summary.get('query_id',''), Sample=summary.get('Sample',''),
                   GCFID=summary.get('GCFID',''), status=summary.get('status',''),
                   arm=summary.get('arm',''), method=method, alpha=alpha,
                   null_count=0 if values is None else len(values),
                   unique_references=unique, jaccard_score=score,
                   reference_lower_limit=None, bootstrap_lower_ci=None,
                   bootstrap_upper_ci=None, bootstrap_pass_fraction=None,
                   jaccard_decision='NOT_EVALUATED', pearson_p=pearson_p,
                   pearson_decision=pearson, detection_status='NOT_EVALUATED',
                   method_status=null_state)
        if not complete:
            output.append(row)
            continue
        if method == 'normal':
            p = _number(summary.get('jaccard_p'))
            row['jaccard_decision'] = _decision(p, alpha)
            if p is None:
                row['method_status'] = str(summary.get('jaccard_status') or 'nonestimable_normal')
            else:
                mean, sd = _number(summary.get('jaccard_null_mean')), _number(summary.get('jaccard_null_sd'))
                if mean is not None and sd is not None:
                    row['reference_lower_limit'] = mean + NormalDist().inv_cdf(alpha) * sd
                row['method_status'] = 'estimable'
        elif method in ('empirical', 'minimum'):
            limit = empirical_lower(values, alpha) if method == 'empirical' else float(np.min(values))
            row['reference_lower_limit'] = limit
            row['jaccard_decision'] = 'PASS' if score >= limit else 'FAIL'
            row['method_status'] = 'exploratory_dependent_reference_draws'
        else:
            boot = _bootstrap_limits(values, refs, alpha, bootstrap_replicates,
                                     bootstrap_confidence, seed, row['query_id'])
            if boot is None:
                row['method_status'] = 'insufficient_unique_references_for_bootstrap'
            else:
                lower, upper, cutoffs = boot
                row.update(reference_lower_limit=empirical_lower(values, alpha),
                           bootstrap_lower_ci=lower, bootstrap_upper_ci=upper,
                           bootstrap_pass_fraction=float(np.mean(cutoffs <= score)))
                # Declare unresolved comparisons; do not silently force a call.
                row['jaccard_decision'] = ('PASS' if score >= upper else
                                           'FAIL' if score < lower else 'UNCERTAIN')
                row['method_status'] = ('exploratory_bootstrap_consensus'
                                        if row['jaccard_decision'] != 'UNCERTAIN'
                                        else 'bootstrap_boundary_uncertain')
        row['detection_status'] = _combine(row['jaccard_decision'], pearson)
        output.append(row)
    return output


def _read_table(path):
    opener = gzip.open if path.suffix == '.gz' else open
    with opener(path, 'rt', encoding='utf-8', newline='') as handle:
        return list(csv.DictReader(handle, delimiter='\t'))


def evaluate_saved(directory, summary, *, methods='normal', alpha=.05,
                   bootstrap_replicates=1000, bootstrap_confidence=.95, seed=401):
    """Score a finished query without decoding BAM/CRAM or changing its draws."""
    directory = Path(directory)
    draws_file, null_file = directory/'draws.tsv.gz', directory/'null.tsv.gz'
    draws, nulls = (_read_table(draws_file), _read_table(null_file)) if (
        draws_file.is_file() and null_file.is_file()) else ([], [])
    return evaluate(summary, draws, nulls, methods=methods, alpha=alpha,
                    bootstrap_replicates=bootstrap_replicates,
                    bootstrap_confidence=bootstrap_confidence, seed=seed)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--run-dir', type=Path, required=True,
                        help='Completed v2 query output root with queries/<query_id>/summary.json')
    parser.add_argument('--out', type=Path, required=True, help='New TSV output path')
    parser.add_argument('--methods', default='normal,empirical,minimum,bootstrap')
    parser.add_argument('--alpha', type=float, default=.05)
    parser.add_argument('--bootstrap-replicates', type=int, default=1000)
    parser.add_argument('--bootstrap-confidence', type=float, default=.95)
    parser.add_argument('--seed', type=int, default=401)
    args = parser.parse_args()
    parse_methods(args.methods)
    validate_options(args.alpha, args.bootstrap_replicates, args.bootstrap_confidence)
    root = args.run_dir/'queries'
    if not root.is_dir():
        parser.error('Missing query output directory: ' + str(root))
    if args.out.exists():
        parser.error('Output already exists; choose a new file')
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open('w', encoding='utf-8', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=COLUMNS, delimiter='\t',
                                lineterminator='\n', extrasaction='raise')
        writer.writeheader()
        for directory in sorted(root.iterdir()):
            if not directory.is_dir():
                continue
            summary_file = directory/'summary.json'
            if not summary_file.is_file():
                raise ValueError('Missing completed query summary: ' + str(directory))
            summary = json.loads(summary_file.read_text())
            writer.writerows(evaluate_saved(directory, summary, methods=args.methods,
                                            alpha=args.alpha, bootstrap_replicates=args.bootstrap_replicates,
                                            bootstrap_confidence=args.bootstrap_confidence, seed=args.seed))


if __name__ == '__main__':
    main()
