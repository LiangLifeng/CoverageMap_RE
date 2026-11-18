"""Bounded CPU-parallel simulation and original normal-approximation concordance."""
from contextlib import nullcontext
from concurrent.futures import ThreadPoolExecutor, ProcessPoolExecutor, wait, FIRST_COMPLETED
import hashlib
import multiprocessing
from pathlib import Path
import time
import numpy as np
from bin.coverage_core import BlockPool, select_units, window_sums, choose_depth_sampler
from Tools import z_test

RNG_POLICY = 'depth-value-role-replicate-v1'
_PROCESS_POOL = None
_PROCESS_EDGES = None


def task_rng(seed, depth, role, replicate):
    digest = hashlib.sha256(float(depth).hex().encode()).digest()[:16]
    words = np.frombuffer(digest, dtype='<u4').tolist()
    return np.random.default_rng(np.random.SeedSequence([seed, *words, role, replicate]))


def evaluate_windows(test, simulated, alpha=.01, block_size=64, backend="cpu"):
    """Same off-diagonal row means / normal z calculation; no retained R x R matrix."""
    test = np.asarray(test, dtype=np.float64)
    x = np.asarray(simulated, dtype=np.float64)
    if x.ndim != 2 or len(x) < 2 or test.shape != (x.shape[1],):
        raise ValueError('Require N>=2 simulations and matching test windows')
    if not 0 < alpha < 1 or block_size < 1:
        raise ValueError('Invalid alpha/block size')
    if np.any(~np.isfinite(x)) or np.any(~np.isfinite(test)) or np.any(x < 0) or np.any(test < 0):
        raise ValueError('Coverage must be finite and nonnegative')
    if test.sum() == 0 or np.any(x.sum(axis=1) == 0):
        raise ValueError('Zero coverage cannot provide concordance evidence')
    if backend not in ('cpu', 'gpu'):
        raise ValueError('Unknown statistics backend')
    xp = np
    if backend == 'gpu':
        from bin.gpu_core import require_cupy
        xp = require_cupy()
        test, x = xp.asarray(test), xp.asarray(x)
    n = len(x)
    positive = (x > 0).astype(xp.float64 if backend == "gpu" else xp.int64)
    totals = positive.sum(axis=1)
    test_positive = (test > 0).astype(xp.float64 if backend == "gpu" else xp.int64)
    test_intersections = positive @ test_positive
    test_unions = totals + test_positive.sum() - test_intersections
    test_j = float(xp.mean(test_intersections / test_unions))
    j_rows = xp.zeros(n)
    centered = x - x.mean(axis=1, keepdims=True)
    norms = xp.linalg.norm(centered, axis=1)
    tc = test - test.mean()
    tn = float(xp.linalg.norm(tc))
    pearson_ok = tn > 0 and xp.all(norms > 0)
    r_rows = xp.zeros(n) if pearson_ok else None
    if pearson_ok:
        test_r = float(xp.mean(xp.clip((centered @ tc) / (norms * tn), -1, 1)))
    for start in range(0, n, block_size):
        stop = min(start + block_size, n)
        intersections = positive[start:stop] @ positive.T
        union = totals[start:stop, None] + totals[None, :] - intersections
        values = intersections / union
        values[xp.arange(stop-start), xp.arange(start, stop)] = 0
        j_rows[start:stop] = values.sum(axis=1) / (n-1)
        if pearson_ok:
            values = xp.clip((centered[start:stop] @ centered.T) / (norms[start:stop, None] * norms[None, :]), -1, 1)
            values[xp.arange(stop-start), xp.arange(start, stop)] = 0
            r_rows[start:stop] = values.sum(axis=1) / (n-1)
    def result(value, distribution):
        if backend == "gpu":
            distribution = xp.asnumpy(distribution)
        p = z_test(value, distribution)
        return dict(status='estimable', similarity=value, reference_mean=float(distribution.mean()),
                    reference_sd=float(distribution.std()), p_value=p, accepted=bool(p >= alpha))
    return dict(jaccard=result(test_j, j_rows),
                pearson=result(test_r, r_rows) if pearson_ok else
                    dict(status='non_estimable', reason='constant_test_or_reference_windows',
                         similarity=None, p_value=None, accepted=None),
                alpha=alpha, model='existing-normal-approximation-not-calibrated',
                simulations=n, windows=x.shape[1])


def _batch(pool, edges, job, gpu=None):
    index, depth, begin, end, sampler, seed = job
    rows = np.empty((end-begin, len(edges)-1), dtype=np.int64)
    summaries = []
    audits = {}
    started = time.perf_counter()
    selected = []
    for offset, replicate in enumerate(range(begin, end)):
        ids = select_units(pool, task_rng(seed, depth, 0, replicate),
                           target_bases=depth*pool.genome_len, depth_sampler=sampler)
        if gpu is None:
            rows[offset] = window_sums(pool, ids, edges)
        else:
            selected.append(ids)
        mass = int(pool.unit_bases[ids].sum(dtype=np.int64))
        if gpu is None and int(rows[offset].sum()) != mass:
            raise AssertionError('Coverage mass mismatch')
        summaries.append(dict(replicate=replicate, units=len(ids), aligned_bases=mass,
                              actual_depth=mass/pool.genome_len,
                              selected_ids_sha256=hashlib.sha256(ids.astype('<i8').tobytes()).hexdigest()))
        if replicate in (0, end-1):
            audits[str(replicate)] = ids
    if gpu is not None:
        rows = gpu._count_selected_batch(selected)
        if not np.array_equal(rows.sum(axis=1), [r['aligned_bases'] for r in summaries]):
            raise AssertionError('GPU coverage mass mismatch')
    return index, begin, rows, summaries, audits, time.perf_counter()-started


def _process_init(directory, genome_len, edges):
    global _PROCESS_POOL, _PROCESS_EDGES
    directory = Path(directory)
    _PROCESS_POOL = BlockPool(*(np.load(directory/f'{name}.npy', mmap_mode='r')
                               for name in ('starts', 'ends', 'indptr')), genome_len)
    _PROCESS_EDGES = edges


def _process_batch(job):
    return _batch(_PROCESS_POOL, _PROCESS_EDGES, job)


def simulate_grid(pool, edges, depths, replicates, seed, workers, executor, batch_size,
                  sampler, output, complete, gpu=None, timing=None):
    """Bound pending jobs to 2*workers; complete callbacks consume one depth at a time."""
    if replicates < 2 or workers < 1 or batch_size < 1 or executor not in ('serial','thread','process'):
        raise ValueError('Invalid simulation/executor configuration')
    if len(set(depths)) != len(depths) or any(not np.isfinite(d) or d <= 0 or d*pool.genome_len > pool.total_bases for d in depths):
        raise ValueError('Depths must be positive, unique and supported by the pool')
    if gpu is not None and executor == 'process' and workers > 1:
        raise ValueError('GPU simulation supports serial/thread CPU sampling, not process workers')
    output = Path(output)
    calibration = {}
    for index, depth in enumerate(depths):
        tick = time.perf_counter()
        chosen, decision = choose_depth_sampler(pool, depth*pool.genome_len) if sampler == 'auto' else (sampler, None)
        if chosen not in ('lazy','permutation'):
            raise ValueError('Unknown sampler')
        calibration[index] = dict(sampler=chosen, decision=decision, calibration_seconds=time.perf_counter()-tick)
    def jobs():
        for index, depth in enumerate(depths):
            for begin in range(0, replicates, batch_size):
                yield index, depth, begin, min(replicates,begin+batch_size), calibration[index]['sampler'], seed
    active = {}
    def consume(answer):
        index, begin, rows, summaries, audits, seconds = answer
        if index not in active:
            directory=output/f'depth-{index:05d}'
            directory.mkdir()
            active[index]=dict(directory=directory,
                matrix=np.lib.format.open_memmap(directory/'simulations.npy',mode='w+',dtype='<i8',shape=(replicates,len(edges)-1)),
                summaries=[None]*replicates, received=0, seconds=0.)
        state=active[index];end=begin+len(rows)
        state['matrix'][begin:end]=rows
        state['summaries'][begin:end]=summaries
        state['received']+=len(rows);state['seconds']+=seconds
        with timing.measure('audit_ids_write') if timing else nullcontext():
            np.savez(state['directory']/f'ids-{begin:06d}.npz', **audits)
        if state['received']==replicates:
            state['matrix'].flush()
            complete(index,depths[index],state['directory'],state['matrix'],state['summaries'],calibration[index],state['seconds'])
            state['matrix']._mmap.close()
            del active[index]
    pool_executor = None
    try:
        if workers==1 or executor=='serial':
            for job in jobs():consume(_batch(pool,edges,job,gpu))
        else:
            if executor=='process':
                shared=output/'shared-pool';shared.mkdir()
                for name in ('starts','ends','indptr'):
                    np.save(shared/f'{name}.npy',getattr(pool,name),allow_pickle=False)
                pool_executor=ProcessPoolExecutor(max_workers=workers,
                    mp_context=multiprocessing.get_context('spawn'),initializer=_process_init,
                    initargs=(str(shared),pool.genome_len,edges))
                submit=lambda job:pool_executor.submit(_process_batch,job)
            else:
                pool_executor=ThreadPoolExecutor(max_workers=workers)
                submit=lambda job:pool_executor.submit(_batch,pool,edges,job,gpu)
            pending=set();iterator=iter(jobs());exhausted=False
            while pending or not exhausted:
                while len(pending)<2*workers and not exhausted:
                    try:pending.add(submit(next(iterator)))
                    except StopIteration:exhausted=True
                if pending:
                    done,pending=wait(pending,return_when=FIRST_COMPLETED)
                    for future in done:consume(future.result())
    finally:
        if pool_executor is not None:pool_executor.shutdown(wait=True,cancel_futures=True)
        for state in active.values():state['matrix']._mmap.close()
    return calibration
