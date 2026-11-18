"""Coverage-breadth states shared by query/reference and reference/reference."""
import math
import numpy as np

DEFAULT_GRID = (0., .01, .05, .10, .15, .20, .30)


def thresholds(single=0., grid=None):
    values = [single]
    if grid:
        values.extend(DEFAULT_GRID if grid == 'default' else map(float, grid.split(',')))
    if any(not math.isfinite(v) or not 0 <= v <= 1 for v in values):
        raise ValueError('Jaccard coverage thresholds must be finite in [0,1]')
    return list(dict.fromkeys(values))


def covered_state(covered, lengths, threshold=0.):
    covered, lengths = np.asarray(covered), np.asarray(lengths)
    if (covered.ndim not in (1, 2) or covered.shape[-1:] != lengths.shape or lengths.ndim != 1
            or not lengths.size or covered.dtype.kind not in 'iu' or lengths.dtype.kind not in 'iu'
            or np.any(lengths <= 0) or np.any(covered < 0) or np.any(covered > lengths)):
        raise ValueError('Covered bases must be integer counts within positive window lengths')
    thresholds(threshold)
    return covered > 0 if threshold == 0 else covered / lengths >= threshold


def packed(state):
    return int.from_bytes(np.packbits(state).tobytes(), 'big')


def pair(left, right):
    intersection = (left & right).bit_count()
    union = (left | right).bit_count()
    return dict(jaccard=intersection/union if union else None,
                intersection=intersection, union=union,
                jaccard_status='estimable' if union else 'zero_union')


def reference_means(states):
    """Ordered upper triangle, incomplete rows remain non-estimable."""
    profiles = [packed(row) for row in states]
    sums, valid = np.zeros(len(profiles)), np.ones(len(profiles), dtype=bool)
    for i, left in enumerate(profiles):
        for j in range(i+1, len(profiles)):
            score = pair(left, profiles[j])['jaccard']
            if score is None:
                valid[[i, j]] = False
            else:
                sums[i] += score
                sums[j] += score
    return [float(v/(len(profiles)-1)) if valid[i] and len(profiles)>1 else None for i,v in enumerate(sums)]
