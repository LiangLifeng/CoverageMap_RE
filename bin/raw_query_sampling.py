"""Deterministic whole-unit depth matching for Raw Query vNext.

Each draw is a random permutation prefix; select its nearest nonempty prefix.
Reads and fragments are never truncated or duplicated within a draw.
"""
from dataclasses import dataclass
import numpy as np


@dataclass(frozen=True)
class QueryArm:
    name: str
    mate_side: str
    query_ids: np.ndarray
    target_bases: int
    slots: int


@dataclass(frozen=True)
class UnitDraw:
    indices: np.ndarray
    actual_bases: int
    target_bases: int
    available_bases: int
    delta_bases: int
    within_tolerance: bool
    sampling_status: str


def query_arms(pool, sampling_unit, reference_draws):
    """Separate PE mate-stratified draws without dividing query depth by two."""
    if sampling_unit not in ('read', 'fragment'):
        raise ValueError('sampling_unit must be read or fragment')
    if reference_draws < 2:
        raise ValueError('reference_draws must be >=2')
    indices = np.arange(len(pool), dtype=np.intp)
    if sampling_unit == 'fragment':
        return [QueryArm('paired', 'both', indices, int(pool.total_bases), reference_draws)]
    roles = getattr(pool, 'unit_mate_sides', None)
    if roles is None:
        raise ValueError('Read-unit geometry lacks mate-side metadata')
    roles = np.asarray(roles)
    if len(roles) != len(pool):
        raise ValueError('Mate-side metadata length mismatch')
    if (roles == 1).any() and (roles == 2).any() and not np.any(roles == 0):
        left, right = reference_draws // 2, reference_draws - reference_draws // 2
        return [QueryArm('read1', 'read1', indices[roles == 1],
                         int(pool.unit_bases[roles == 1].sum()), left),
                QueryArm('read2', 'read2', indices[roles == 2],
                         int(pool.unit_bases[roles == 2].sum()), right)]
    # Mixed/partial mate information is handled as unstratified independent reads;
    # never drop valid reads silently just to force a two-arm split.
    return [QueryArm('single', 'all', indices, int(pool.total_bases), reference_draws)]


def compatible_reference_ids(pool, mate_side):
    if mate_side in ('both', 'unpaired', 'all'):
        return np.arange(len(pool), dtype=np.intp)
    roles = getattr(pool, 'unit_mate_sides', None)
    if roles is None:
        raise ValueError('Reference read geometry lacks mate-side metadata')
    side = 1 if mate_side == 'read1' else 2
    return np.flatnonzero(np.asarray(roles) == side).astype(np.intp)


def select_nearest_units(pool, eligible_ids, rng, target_bases, tolerance_bp):
    """Closest nonempty permutation prefix; undershoot is a reported outcome.

    If reference capacity is lower than requested, use every eligible unit and
    report the deficit. Out-of-tolerance draws retain their geometry for QC,
    but must not silently enter an accepted reference distribution.
    """
    if not isinstance(target_bases, (int, np.integer)) or target_bases < 0:
        raise ValueError('target_bases must be nonnegative integer')
    if not isinstance(tolerance_bp, (int, np.integer)) or tolerance_bp < 0:
        raise ValueError('tolerance_bp must be nonnegative integer')
    ids = np.asarray(eligible_ids, dtype=np.intp)
    if ids.ndim != 1 or np.any(ids < 0) or np.any(ids >= len(pool)) or len(np.unique(ids)) != len(ids):
        raise ValueError('eligible_ids must be unique in-range unit IDs')
    if target_bases == 0 or not len(ids):
        return UnitDraw(ids[:0], 0, int(target_bases), 0,
                        -int(target_bases), False, 'NO_ELIGIBLE_UNITS')
    available = int(pool.unit_bases[ids].sum(dtype=np.int64))
    order = rng.permutation(ids)
    running = np.cumsum(pool.unit_bases[order], dtype=np.int64)
    at = int(np.searchsorted(running, target_bases, side='left'))
    best = min(len(ids), at + 1)
    if at > 0:
        before = int(running[at - 1])
        after = int(running[at]) if at < len(ids) else before
        if abs(before - target_bases) <= abs(after - target_bases):
            best = at
    chosen = order[:best]
    actual = int(running[best - 1])
    delta = actual - int(target_bases)
    status = 'EXACT' if delta == 0 else 'OVERSHOOT' if delta > 0 else 'UNDERSHOOT'
    return UnitDraw(chosen, actual, int(target_bases), available, delta,
                    abs(delta) <= tolerance_bp, status)
