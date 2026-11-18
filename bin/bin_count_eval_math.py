"""Exact windows and descriptive coverage-space scores; ordinations never set scores."""
from collections import Counter
from concurrent.futures import ThreadPoolExecutor
from contextlib import nullcontext
import math
import numpy as np

METRICS = ('euclidean', 'cosine', 'jensenshannon', 'correlation')
DISTANCE_TILE_ROWS = 8
DISTANCE_PANEL_ROWS = 256


class UndefinedCorrelationError(ValueError):
    """Undefined Pearson profiles; never remove rows or substitute a metric."""

    def __init__(self, issues):
        self.issues = issues
        invalid = sum(r['reason'] != 'NEAR_CONSTANT_PROFILE' for r in issues)
        super().__init__(f'Pearson correlation undefined for {invalid} row/center profiles; '
                         'no records removed or correlations imputed')


def correlation_diagnostics(values, kind='row'):
    """Bounded-memory numerical diagnostics; exact constants are undefined."""
    values = np.atleast_2d(values)
    if values.ndim != 2:
        raise ValueError('Finite 2D Pearson profiles required')
    issues = []
    batch = max(1, (16 * 1024**2) // (8 * max(1, values.shape[1])))
    for start in range(0, len(values), batch):
        block = np.asarray(values[start:start+batch], dtype=np.float64)
        if not np.isfinite(block).all():
            raise ValueError('Finite 2D Pearson profiles required')
        if block.shape[1] < 2:
            issues.extend(dict(kind=kind, index=start+i, reason='FEWER_THAN_TWO_BINS')
                          for i in range(len(block)))
            continue
        constant = np.ptp(block, axis=1) == 0
        means = block.mean(axis=1)
        norm = np.linalg.norm(block - means[:, None], axis=1)
        near = (~constant) & (norm < 1e-13 * np.abs(means))
        issues.extend(dict(kind=kind, index=start+int(i),
                           reason='CONSTANT_PROFILE' if constant[i] else 'NEAR_CONSTANT_PROFILE')
                      for i in np.flatnonzero(constant | near))
    return issues


def require_correlation_defined(values):
    issues = correlation_diagnostics(values)
    if any(r['reason'] != 'NEAR_CONSTANT_PROFILE' for r in issues):
        raise UndefinedCorrelationError(issues)


def exact_edges(length, count):
    if type(length) is not int or type(count) is not int or not 1 <= count <= length:
        raise ValueError('Require integer 1 <= bin count <= genome length')
    width = length // count
    edges = np.arange(count + 1, dtype=np.int64) * width
    edges[-1] = length
    return edges


def normalized(values, target=1e6):
    values = np.asarray(values, dtype=np.float64)
    if values.ndim != 1 or not np.isfinite(values).all() or (values < 0).any():
        raise ValueError('Nonnegative finite 1D coverage required')
    total = values.sum()
    if not math.isfinite(total) or total <= 0 or not math.isfinite(target) or target <= 0:
        raise ValueError('Cannot normalize zero/nonfinite coverage')
    return values * (target / total)


def distances(left, right, metric):
    from scipy.spatial.distance import cdist
    if metric not in METRICS:
        raise ValueError('Unsupported distance')
    if metric == 'correlation':
        require_correlation_defined(left)
        require_correlation_defined(right)
    result = cdist(np.atleast_2d(left), np.atleast_2d(right), metric=metric)
    if not np.isfinite(result).all():
        raise ValueError('Nonfinite distances (zero norms or invalid distribution)')
    return np.clip(result, 0, 2) if metric == 'correlation' else np.maximum(result, 0)


def score(matrix, labels, metric, diagnostics=None, *, jobs=1):
    """Mean-profile centers; equal-GCF within mean, unordered pair between mean."""
    labels = np.asarray(labels)
    groups = sorted(set(labels.tolist()))
    if matrix.ndim != 2 or len(matrix) != len(labels) or not groups:
        raise ValueError('Matrix/labels shape mismatch')
    centers = np.empty((len(groups), matrix.shape[1]), dtype=np.float64)
    per_row = np.empty(len(matrix), dtype=np.float64)
    batch = max(1, (16 * 1024**2) // (8 * matrix.shape[1]))

    def make_center(i):
        indices = np.flatnonzero(labels == groups[i])
        total = np.zeros(matrix.shape[1], dtype=np.float64)
        for k in range(0, len(indices), batch):
            block = np.asarray(matrix[indices[k:k+batch]])
            if not np.isfinite(block).all() or (block < 0).any() or (block.sum(axis=1) <= 0).any():
                raise ValueError('Invalid/zero normalized coverage')
            total += block.sum(axis=0)
        centers[i] = total / len(indices)

    def score_group(i):
        indices = np.flatnonzero(labels == groups[i])
        for k in range(0, len(indices), batch):
            ids = indices[k:k+batch]
            per_row[ids] = distances(matrix[ids], centers[i], metric)[:, 0]
        return dict(gcf=groups[i], samples=len(indices), within_mean=float(per_row[indices].mean()))

    workers = min(jobs, len(groups))
    with (ThreadPoolExecutor(workers) if workers > 1 else nullcontext()) as executor:
        map_groups = executor.map if executor else map
        for _ in map_groups(make_center, range(len(groups))):
            pass
        if metric == 'correlation':
            issues = correlation_diagnostics(matrix) + correlation_diagnostics(centers, 'center')
            if diagnostics is not None:
                diagnostics.extend(issues)
            if any(r['reason'] != 'NEAR_CONSTANT_PROFILE' for r in issues):
                raise UndefinedCorrelationError(issues)
        within = list(map_groups(score_group, range(len(groups))))
    center_distances = distances(centers, centers, metric)
    pairs = [(groups[i], groups[j], float(center_distances[i, j]))
             for i in range(len(groups)) for j in range(i+1, len(groups))]
    intra = float(np.mean([r['within_mean'] for r in within]))
    inter = float(np.mean([p[2] for p in pairs])) if pairs else None
    state = 'SINGLE_GCF' if not pairs else ('ZERO_WITHIN_DISTANCE' if intra == 0 else 'OK')
    result = dict(distance=metric, gcfs=len(groups), records=len(matrix), within_mean=intra,
                  between_mean=inter, ratio=inter/intra if state == 'OK' else None, status=state,
                  center_definition='arithmetic mean profile; not a minimizing metric barycenter',
                  within_weighting='equal GCFs', between_weighting='equal unordered GCF pairs')
    if metric == 'correlation':
        result.update(correlation_definition='1 - Pearson r; centered profiles',
                      within_pearson_r_mean=1-intra,
                      between_pearson_r_mean=1-inter if inter is not None else None,
                      near_constant_profiles=sum(r['reason']=='NEAR_CONSTANT_PROFILE' for r in issues))
    return result, groups, centers, within, per_row, pairs


def pairwise(matrix, metric, jobs=1):
    from scipy.spatial.distance import pdist, squareform
    if metric == 'correlation':
        require_correlation_defined(matrix)
        result = np.empty((len(matrix), len(matrix)), dtype=np.float64)

        def store(rows, columns, block):
            result[rows, columns] = block

        distance_tiles(matrix, metric, jobs, store)
    elif jobs > 1 and len(matrix) > 64:
        result = parallel_pairwise(matrix, metric, jobs)
    else:
        result = squareform(pdist(matrix, metric=metric))
    if not np.isfinite(result).all():
        raise ValueError('Nonfinite pairwise distances')
    np.maximum(result, 0, out=result)
    if metric == 'correlation':
        np.minimum(result, 2, out=result)
    np.fill_diagonal(result, 0)
    return result


def parallel_pairwise(matrix, metric, jobs):
    from scipy.spatial.distance import cdist

    matrix = np.ascontiguousarray(matrix, dtype=np.float64)
    count = len(matrix)
    block_rows = min(64, max(1, count // (4 * jobs)))
    result = np.empty((count, count), dtype=np.float64)

    def write_block(start):
        stop = min(start + block_rows, count)
        block = cdist(matrix[start:stop], matrix[start:], metric=metric)
        width = stop - start
        # Keep the same i < j direction as pdist, including within this tile.
        diagonal = block[:, :width]
        lower = np.tril_indices(width, -1)
        diagonal[lower] = diagonal.T[lower]
        # Each worker owns a disjoint upper strip and its reflected lower strip.
        result[start:stop, start:] = block
        result[stop:, start:stop] = block[:, width:].T

    with ThreadPoolExecutor(max_workers=jobs) as executor:
        for _ in executor.map(write_block, range(0, count, block_rows)):
            pass
    return result


def distance_profiles(matrix, rows, correlation):
    if not correlation:
        return np.ascontiguousarray(matrix[rows], dtype=np.float64)
    block = np.array(matrix[rows], dtype=np.float64, order='C', copy=True)
    block -= block.mean(axis=1, keepdims=True)
    return block


def distance_tiles(matrix, metric, jobs, consume):
    """Visit each pair in its original upper-triangle direction and reflect it."""
    from scipy.spatial.distance import cdist

    correlation = metric == 'correlation'
    kernel = 'cosine' if correlation else metric
    with ThreadPoolExecutor(max_workers=jobs) as executor:
        for panel_start in range(0, len(matrix), DISTANCE_PANEL_ROWS):
            end = min(panel_start + DISTANCE_PANEL_ROWS, len(matrix))
            right = distance_profiles(matrix, slice(panel_start, end), correlation)

            def calculate(start):
                stop = min(start + DISTANCE_TILE_ROWS, end)
                rows = slice(start, stop)
                first = max(start, panel_start)
                if start < panel_start:
                    left = distance_profiles(matrix, rows, correlation)
                else:
                    left = right[start-panel_start:stop-panel_start]
                block = cdist(left, right[first-panel_start:], metric=kernel)
                if first == start:
                    diagonal = block[:, :stop-start]
                    lower = np.tril_indices(stop-start, -1)
                    diagonal[lower] = diagonal.T[lower]
                    np.fill_diagonal(diagonal, 0)
                return rows, slice(first, end), block

            for rows, columns, block in executor.map(calculate, range(0, end, DISTANCE_TILE_ROWS)):
                consume(rows, columns, block)
                mirror = max(rows.stop, columns.start)
                if mirror < end:
                    consume(slice(mirror, end), rows, block[:, mirror-columns.start:].T)
            del right


def tsne_neighbor_count(n, perplexity):
    return min(n - 1, int(3 * perplexity + 1))


def distance_replay_rows(bins):
    return max(DISTANCE_TILE_ROWS, min(4096, 1024**2 // bins))


def nearest_partition(distances, size):
    """Match the installed sklearn dense precomputed selector, including ties."""
    selected = np.argpartition(distances, size-1, axis=1)[:, :size]
    values = np.take_along_axis(distances, selected, axis=1)
    order = np.argsort(values, axis=1)
    return np.take_along_axis(values, order, axis=1), np.take_along_axis(selected, order, axis=1)


def tsne_neighbor_graph(matrix, metric, perplexity, jobs):
    """Build the exact sparse precomputed input without keeping an n-by-n matrix."""
    from scipy.sparse import csr_array
    from scipy.spatial.distance import cdist

    n = len(matrix)
    correlation = metric == 'correlation'
    if correlation:
        require_correlation_defined(matrix)
    kernel = 'cosine' if correlation else metric
    # sklearn queries k+1 before dropping self. Keep one extra to detect cutoff ties.
    size = tsne_neighbor_count(n, perplexity) + 1
    keep = min(n, size + 1)
    distances = np.full((n, keep), np.inf)
    indices = np.empty((n, keep), dtype=np.int32)
    indices.fill(-1)

    def clamp(block):
        if not np.isfinite(block).all():
            raise ValueError('Nonfinite pairwise distances')
        np.maximum(block, 0, out=block)
        if correlation:
            np.minimum(block, 2, out=block)

    def update(rows, columns, block):
        clamp(block)
        combined = np.concatenate((distances[rows], block), axis=1)
        selected = np.argpartition(combined, keep-1, axis=1)[:, :keep]
        previous = np.take_along_axis(indices[rows], np.minimum(selected, keep-1), axis=1)
        indices[rows] = np.where(selected < keep, previous, selected - keep + columns.start)
        distances[rows] = np.take_along_axis(combined, selected, axis=1)

    distance_tiles(matrix, metric, jobs, update)
    order = np.argsort(distances, axis=1)
    distances = np.take_along_axis(distances, order, axis=1)
    indices = np.take_along_axis(indices, order, axis=1)
    # Full-row replay preserves NumPy's tie selection and ordering, without storing other rows.
    tied = np.flatnonzero(np.any(distances[:, 1:] == distances[:, :-1], axis=1))
    distances = np.ascontiguousarray(distances[:, :size])
    indices = np.ascontiguousarray(indices[:, :size])
    replay_rows = distance_replay_rows(matrix.shape[1])

    def replay(index):
        row = np.empty((1, n), dtype=np.float64)
        profile = distance_profiles(matrix, slice(index, index+1), correlation)
        for low, high in ((0, index), (index+1, n)):
            for start in range(low, high, replay_rows):
                stop = min(start + replay_rows, high)
                other = distance_profiles(matrix, slice(start, stop), correlation)
                block = cdist(other, profile, metric=kernel).T if low == 0 else cdist(profile, other, metric=kernel)
                row[:, start:stop] = block
                del other
        row[0, index] = 0
        clamp(row)
        return nearest_partition(row, size)

    if len(tied):
        with ThreadPoolExecutor(max_workers=jobs) as executor:
            for start in range(0, len(tied), jobs):
                batch = tied[start:start+jobs]
                for index, (values, selected) in zip(batch, executor.map(replay, batch)):
                    distances[index] = values[0]
                    indices[index] = selected[0]
    pointers = np.arange(0, (n+1)*size, size, dtype=np.int64)
    return csr_array((distances.ravel(), indices.ravel(), pointers), shape=(n, n))


def pcoa(distance_matrix):
    from scipy.linalg import eigh
    squared = np.square(distance_matrix)
    row_mean = squared.mean(axis=1)
    gram = -.5 * (squared - row_mean[:, None] - row_mean[None, :] + squared.mean())
    eigenvalues, vectors = eigh(gram, check_finite=False)
    order = np.argsort(eigenvalues)[::-1]
    eigenvalues, vectors = eigenvalues[order], vectors[:, order]
    tolerance = max(1., float(np.abs(eigenvalues).max())) * np.finfo(float).eps * max(1, len(eigenvalues))
    positive = np.flatnonzero(eigenvalues > tolerance)
    coords = np.zeros((len(gram), 2))
    for k, index in enumerate(positive[:2]):
        col = vectors[:, index]
        if col[np.argmax(np.abs(col))] < 0:
            col = -col
        coords[:, k] = col * np.sqrt(eigenvalues[index])
    negative = eigenvalues[eigenvalues < -tolerance]
    denom = float(eigenvalues[positive].sum())
    return coords, dict(eigenvalues=eigenvalues.tolist(),
                        negative_eigenvalue_count=len(negative),
                        negative_absolute_mass=float(np.abs(negative).sum()),
                        positive_variance_explained=[float(max(v, 0)/denom) if denom else 0
                                                     for v in eigenvalues[:2]],
                        warning='NEGATIVE_EIGENVALUES' if len(negative) else '')


def tsne(distance_matrix, seed, perplexity, iterations, learning_rate, jobs):
    from sklearn.manifold import TSNE
    if not 0 < perplexity < distance_matrix.shape[0]:
        raise ValueError('t-SNE requires 0 < perplexity < number of rows; no automatic adjustment')
    model = TSNE(n_components=2, metric='precomputed', init='random', random_state=seed,
                 perplexity=perplexity, max_iter=iterations, learning_rate=learning_rate,
                 n_jobs=jobs, method='barnes_hut')
    coords = model.fit_transform(distance_matrix)
    return coords, dict(kl_divergence=float(model.kl_divergence_), iterations=int(model.n_iter_),
                        seed=seed, perplexity=perplexity, learning_rate=learning_rate,
                        initialization='random', metric='precomputed; sklearn squares distances internally')
