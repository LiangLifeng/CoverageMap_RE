# Exact function source extracted from Tools.snapshot.py (user source).
# Original full-file SHA256: b94569f5ae821b58c46bd46f9149d5acb753f56feef2532bc1ea009a84023d3e
# Loaded by tests with NumPy, typing, and a standard normal CDF provider.

def bin_coverage(
    cov: np.ndarray,
    bin_count: Optional[int] = 1000,
    bin_size: Optional[int] = None,
    mode: str = "sum",
    drop_tail: bool = False,
) -> np.ndarray:
    """Bin per-base coverage along the last axis.

    Parameters
    ----------
    cov : np.ndarray
        Coverage array of shape ``(..., genome_len)``.
    bin_size : int, optional
        Number of bases per bin (mutually exclusive with ``bin_count``).
    bin_count : int, optional
        Total number of bins (mutually exclusive with ``bin_size``).
    mode : {"sum", "mean"}
        Aggregation within each bin.
    drop_tail : bool
        If ``False`` (default), keep the remainder as an extra bin.

    Returns
    -------
    np.ndarray
        Binned coverage of shape ``(..., n_bins)``.
    """
    cov = np.asarray(cov)

    if cov.ndim == 0:
        raise ValueError("cov must have at least one dimension.")
    length = cov.shape[-1]

    if length == 0:
        return cov.copy()

    if (bin_size is None) == (bin_count is None):
        raise ValueError("Must specify exactly one of bin_size or bin_count.")

    if bin_size is None:
        if bin_count <= 0:
            raise ValueError("bin_count must be a positive integer.")
        bin_size = int(np.ceil(length / bin_count))
    elif bin_size <= 0:
        raise ValueError("bin_size must be a positive integer.")

    if mode not in ("sum", "mean"):
        raise ValueError("mode must be 'sum' or 'mean'.")

    n_full = length // bin_size
    end = n_full * bin_size

    main = cov[..., :end]
    new_shape = (*main.shape[:-1], n_full, bin_size)
    main = main.reshape(new_shape)
    binned = main.sum(axis=-1) if mode == "sum" else main.mean(axis=-1)

    tail = length - end
    if tail > 0 and not drop_tail:
        tail_slice = cov[..., end:]
        if mode == "sum":
            tail_val = tail_slice.sum(axis=-1, keepdims=True)
        else:
            tail_val = tail_slice.mean(axis=-1, keepdims=True)
        binned = np.concatenate([binned, tail_val], axis=-1)

    return binned.astype(np.int32 if mode == "sum" else np.float32)

def z_test(
    test_value: float,
    ref_value_pool: np.ndarray,
    one_side: bool = True,
) -> float:
    """Z-test of ``test_value`` against ``ref_value_pool``.

    Returns the p-value. When the reference pool has zero variance, returns
    ``1.0`` if ``test_value`` equals the mean, otherwise ``0.0``.
    """
    ref = np.asarray(ref_value_pool, dtype=float)
    if ref.size == 0:
        raise ValueError("ref_value_pool must not be empty.")

    mu = ref.mean()
    sigma = ref.std(ddof=0)
    if sigma == 0:
        return 1.0 if test_value == mu else 0.0

    z = (test_value - mu) / sigma
    if one_side:
        return float(norm.cdf(z))
    return float(2.0 * (1.0 - norm.cdf(abs(z))))

def pearson_corr(
    x: np.ndarray,
    y: Optional[np.ndarray] = None,
    eps: float = 1e-12,
) -> Union[float, np.ndarray]:
    """Pearson correlation coefficient(s).

    Supports:

    - vector vs vector  -> scalar
    - vector vs matrix  -> ``(n,)``
    - matrix vs matrix  -> ``(n, k)``
    - matrix internal   -> ``(n, n)``
    """
    x = np.asarray(x, dtype=float)
    if x.ndim == 1:
        if y is None:
            raise ValueError("y must be provided when x is 1-D.")
        x = x[None, :]

    x_centered = x - x.mean(axis=1, keepdims=True)
    x_norm = np.linalg.norm(x_centered, axis=1, keepdims=True)

    if y is None:
        corr = (x_centered @ x_centered.T) / (x_norm @ x_norm.T + eps)
        return corr if corr.shape != (1, 1) else float(corr[0, 0])

    y = np.asarray(y, dtype=float)
    if y.ndim == 1:
        y = y[None, :]

    y_centered = y - y.mean(axis=1, keepdims=True)
    y_norm = np.linalg.norm(y_centered, axis=1, keepdims=True)

    corr = (x_centered @ y_centered.T) / (x_norm @ y_norm.T + eps)

    if corr.shape == (1, 1):
        return float(corr[0, 0])
    if corr.shape[0] == 1:
        return corr[0]
    if corr.shape[1] == 1:
        return corr[:, 0]
    return corr

def jaccard_similarity(
    x: np.ndarray,
    y: Optional[np.ndarray] = None,
    threshold: float = 0,
) -> Union[float, np.ndarray]:
    """Jaccard similarity on presence/absence (values > ``threshold``).

    Supports:

    - vector vs vector  -> scalar
    - vector vs matrix  -> ``(n,)``
    - matrix vs matrix  -> ``(n, k)``
    - matrix internal   -> ``(n, n)``
    """
    x = np.asarray(x)
    if x.ndim == 1:
        if y is None:
            raise ValueError("y must be provided when x is 1-D.")
        x = x[None, :]

    x_bin = x > threshold

    if y is None:
        n = x_bin.shape[0]
        sim = np.zeros((n, n), dtype=float)
        for i in range(n):
            inter = np.logical_and(x_bin[i], x_bin).sum(axis=1)
            union = np.logical_or(x_bin[i], x_bin).sum(axis=1)
            sim[i] = inter / np.maximum(union, 1)
        return sim if sim.shape != (1, 1) else float(sim[0, 0])

    y = np.asarray(y)
    if y.ndim == 1:
        y = y[None, :]
    y_bin = y > threshold

    sim = np.zeros((x_bin.shape[0], y_bin.shape[0]), dtype=float)
    for i in range(x_bin.shape[0]):
        inter = np.logical_and(x_bin[i], y_bin).sum(axis=1)
        union = np.logical_or(x_bin[i], y_bin).sum(axis=1)
        sim[i] = inter / np.maximum(union, 1)

    if sim.shape == (1, 1):
        return float(sim[0, 0])
    if sim.shape[0] == 1:
        return sim[0]
    if sim.shape[1] == 1:
        return sim[:, 0]
    return sim

def concordance_test_high_cov(
    test_cov: np.ndarray,
    ref_cov_pool: np.ndarray,
    ref_simi_jcd: np.ndarray,
    ref_simi_r: np.ndarray,
    alpha: float = 0.01,
    n_samples: int = 1000,
) -> Tuple[bool, float, bool, float]:
    """High-coverage concordance test against precomputed reference statistics.

    Returns ``(accept_jcd, p_jcd, accept_r, p_r)``.
    """
    rng = np.random.default_rng()

    ref_simi_jcd = np.asarray(ref_simi_jcd)
    ref_simi_r = np.asarray(ref_simi_r)

    if ref_simi_jcd.size > n_samples:
        ref_simi_jcd = ref_simi_jcd[
            rng.choice(ref_simi_jcd.size, size=n_samples, replace=False)
        ]
    if ref_simi_r.size > n_samples:
        ref_simi_r = ref_simi_r[
            rng.choice(ref_simi_r.size, size=n_samples, replace=False)
        ]

    test_cov = bin_coverage(test_cov, bin_count=1000)
    ref_cov_pool = bin_coverage(ref_cov_pool, bin_count=1000)

    max_sim_jcd = jaccard_similarity(test_cov, ref_cov_pool).max()
    max_sim_r = pearson_corr(test_cov, ref_cov_pool).max()

    p_value_jcd = z_test(max_sim_jcd, ref_simi_jcd)
    accept_jcd = p_value_jcd >= alpha

    p_value_r = z_test(max_sim_r, ref_simi_r)
    accept_r = p_value_r >= alpha

    return accept_jcd, p_value_jcd, accept_r, p_value_r
