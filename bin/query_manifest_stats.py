"""Position-aligned Tools-style statistics with explicit non-estimable states."""
from concurrent.futures import ThreadPoolExecutor
import math
from statistics import NormalDist

import numpy as np
from bin.jaccard_breadth import covered_state, packed, pair


def similarity(x, y, *, covered_x=None, covered_y=None, lengths=None, threshold=0.):
    x, y = np.asarray(x, dtype=float), np.asarray(y, dtype=float)
    if x.ndim != 1 or x.shape != y.shape or not x.size:
        raise ValueError('Similarity requires matching nonempty vectors')
    if not np.isfinite(x).all() or not np.isfinite(y).all() or min(x.min(), y.min()) < 0:
        raise ValueError('Coverage must be finite and nonnegative')
    if covered_x is None or covered_y is None:
        if threshold != 0 or covered_x is not None or covered_y is not None:
            raise ValueError('Breadth thresholds require both covered-bases vectors and lengths')
        states = x > 0, y > 0
    else:
        states = covered_state(covered_x, lengths, threshold), covered_state(covered_y, lengths, threshold)
        if states[0].shape != x.shape or states[1].shape != y.shape:
            raise ValueError('Covered-bases shape differs from aligned bases')
    jaccard = pair(packed(states[0]), packed(states[1]))
    a, b = x-x.mean(), y-y.mean()
    norm = float(np.linalg.norm(a)*np.linalg.norm(b))
    return dict(**jaccard,
                pearson=float(a@b/(norm+1e-12)) if norm else None,
                pearson_status='estimable' if norm else 'constant_or_single_window')


def normal_test(observed, values, reason=None, *, stable_tail=False):
    """Original lower normal CDF, population SD; no invented degenerate p-values."""
    finite = [float(v) for v in values if v is not None and math.isfinite(v)]
    mean = float(np.mean(finite)) if finite else None
    # A constant repeated float may have nonzero np.std from mean rounding.
    # Exact equality is sufficient here; do not invent a near-zero threshold.
    sd = (0.0 if min(finite) == max(finite) else float(np.std(finite))) if finite else None
    state = reason or ('nonfinite_similarity' if observed is None or len(finite) != len(values)
                       else 'insufficient_null' if len(finite) < 2
                       else 'zero_variance_null' if sd == 0 else 'estimable')
    p = None
    if state == 'estimable':
        z = (observed-mean)/sd
        # erfc avoids cancellation deep in the left tail. Keep legacy replay exact.
        p = 0.5 * math.erfc(-z / math.sqrt(2)) if stable_tail else NormalDist().cdf(z)
    return dict(value=observed, null_count=len(finite), null_mean=mean, null_sd=sd, p=p, status=state)


def off_diagonal_means(size, pair_rows):
    """Reduce in upper-triangle order regardless of how pairs were calculated."""
    sums = {k: np.zeros(size) for k in ('jaccard', 'pearson')}
    valid = {k: np.ones(size, dtype=bool) for k in sums}
    for i, pairs in enumerate(pair_rows):
        for j, scores in enumerate(pairs, start=i+1):
            for key, value in zip(sums, scores):
                if value is None:
                    valid[key][[i, j]] = False
                else:
                    sums[key][i] += value
                    sums[key][j] += value
    return {key: [float(value/(size-1)) if valid[key][i] and size > 1 else None
                  for i, value in enumerate(sums[key])] for key in sums}


def low_null(matrix, *, covered=None, lengths=None, threshold=0.):
    """Off-diagonal row means, never subtract an assumed diagonal of 1."""
    size = len(matrix)

    def pair_row(i):
        for j in range(i+1, size):
            scores = similarity(matrix[i], matrix[j], covered_x=None if covered is None else covered[i],
                                covered_y=None if covered is None else covered[j], lengths=lengths, threshold=threshold)
            yield scores['jaccard'], scores['pearson']

    return off_diagonal_means(size, map(pair_row, range(size-1)))


def reference_row_means(matrix, threads, *, covered=None, lengths=None, threshold=0.):
    """Prepare each coverage map once, then score independent rows of pairs."""
    size = len(matrix)

    states = matrix > 0 if covered is None else covered_state(covered, lengths, threshold)
    if covered is None and threshold != 0:
        raise ValueError('Positive breadth threshold requires covered bases')

    def prepare(item):
        row, state = item
        covered = packed(state)
        # Allocate individual vectors as similarity() does for its dot products.
        centered = row.astype(np.float64)
        centered -= centered.mean()
        return covered, covered.bit_count(), centered, np.linalg.norm(centered)

    def pair_row(i):
        covered, count, centered, norm = profiles[i]
        scores = []
        for other in range(i+1, size):
            other_covered, other_count, other_centered, other_norm = profiles[other]
            intersection = (covered & other_covered).bit_count()
            union = count + other_count - intersection
            denominator = float(norm * other_norm)
            pearson = float(centered @ other_centered / (denominator+1e-12)) if denominator else None
            scores.append((intersection / union if union else None, pearson))
        return scores

    rows = range(size-1)
    workers = min(threads, len(rows))
    if workers <= 1:
        profiles = list(map(prepare, zip(matrix, states)))
        return off_diagonal_means(size, map(pair_row, rows))
    with ThreadPoolExecutor(max_workers=workers) as executor:
        profiles = list(executor.map(prepare, zip(matrix, states)))
        pair_rows = list(executor.map(pair_row, rows))
    return off_diagonal_means(size, pair_rows)


def aggregate(draws, nulls, arm, requested, *, stable_tail=False):
    result = {}
    for metric in ('jaccard', 'pearson'):
        values = [d.get(metric) for d in draws if d['status'] == 'READY']
        finite = [v for v in values if v is not None]
        complete = len(finite) == requested
        mean = float(np.mean(finite)) if finite else None
        maximum = float(max(finite)) if finite else None
        observed = mean if arm == 'low' else maximum
        reason = 'incomplete_or_nonestimable_draws' if not complete else None
        entry = normal_test(observed, nulls.get(metric, []), reason, stable_tail=stable_tail)
        entry.update(mean=mean, max=maximum, valid_draws=len(finite),
                     aggregation='mean' if arm == 'low' else 'max')
        result[metric] = entry
    return result
