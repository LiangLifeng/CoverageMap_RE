#!/usr/bin/env python3
"""Sample reference and test at matching depths, then evaluate concordance."""
import argparse
from decimal import Decimal
import hashlib
import json
import os
from pathlib import Path
import sys
import time


def parse_args():
    p=argparse.ArgumentParser(description=__doc__)
    source=p.add_mutually_exclusive_group(required=True)
    source.add_argument('--alignment',type=Path,help='Reference high-depth alignment pool')
    source.add_argument('--reference-bundle',type=Path,help='Saved reference directory or its storage.json')
    source.add_argument('--reference-depth-manifest',type=Path,help='Verified sample_manifest from READY_DEPTH_INDEX.tsv; reuse stored depth matrices')
    source.add_argument('--reference-sqlite',type=Path,help='Read-only SQLite reference index for batched cached comparisons')
    query=p.add_mutually_exclusive_group()
    query.add_argument('--test-alignment',type=Path)
    query.add_argument('--query-manifest',type=Path,help='CSV: Sample,CramFilePath,GCFID')
    p.add_argument('--query-sample',help='Explicit sample identity for a single SQLite query')
    p.add_argument('--gcf',help='Explicit GCF for a single SQLite query')
    p.add_argument('--query-jobs',type=int,default=1)
    p.add_argument('--query-scheduling',choices=['gcf','query'],default='gcf',
                   help='gcf: keep same-GCF queries in one worker; query: independent query tasks')
    p.add_argument('--reference-limit',type=int,default=1,help='Native index: nearest eligible references per target (default1)')
    p.add_argument('--comparison-policy',choices=['depth-band','nearest'],default='depth-band',
                   help='Native: depth-band uses the <5X simulation / >=5X reference-pool routes; nearest replays the previous workflow')
    p.add_argument('--depth-metric',choices=['mean'],default='mean',help='Full-FAI depth metric (currently mean only)')
    geometry = p.add_mutually_exclusive_group()
    geometry.add_argument('--window-size', '--bin-size', type=int,
                          help='Native depth-band fixed window length in bp (default 1000)')
    geometry.add_argument('--bin-count', type=int, help='Native depth-band exact GCF bin count')
    p.add_argument('--depth-tolerance',type=float,default=1.,help='Inclusive absolute reference depth difference in X for high-depth queries')
    p.add_argument('--low-depth-threshold',type=float,default=5.,help='Queries strictly below this depth use simulation')
    p.add_argument('--simulation-depths',help='Low-depth target layers in X, comma separated; default query mean; N repeats per layer via --replicates')
    p.add_argument('--simulation-references',default='all',help='Per-GCF source sample selection: all, positive integer count, or percentage such as 20%%')
    p.add_argument('--simulation-source-depth-range',help='Optional original source mean-depth bounds MIN:MAX (inclusive); default unrestricted')
    p.add_argument('--simulation-match-depth-band',action='store_true',help='Also restrict low-depth source samples to query depth +/- --depth-tolerance')
    p.add_argument('--reference-cache-mb',type=float,default=256)
    p.add_argument('--reference-batch-size',type=int,default=32)
    p.add_argument('--reference-samples',help='Optional comma-separated reference sample filter')
    p.add_argument('--database-root',type=Path,help='Relocated reference database root for SQLite relative paths')
    identity=p.add_mutually_exclusive_group()
    identity.add_argument('--include-self',action='store_true',default=None,
                          help='Include same-sample/source references (default in native depth-band mode)')
    identity.add_argument('--exclude-self',dest='include_self',action='store_false',default=None,
                          help='Explicitly exclude same-sample IDs and matching source-content SHA')
    ref=p.add_mutually_exclusive_group()
    ref.add_argument('--reference',type=Path);ref.add_argument('--reference-index',type=Path)
    grid=p.add_mutually_exclusive_group()
    grid.add_argument('--depth',type=float);grid.add_argument('--depths',help='Comma-separated positive depth values')
    grid.add_argument('--depth-grid',help='start:stop:step; include stop endpoint')
    p.add_argument('--replicates',type=int,default=100)
    p.add_argument('--bins',type=int,default=1000)
    p.add_argument('--seed',type=int,default=401)
    p.add_argument('--workers',type=int,default=1)
    p.add_argument('--executor',choices=['serial','thread','process'],default='process')
    p.add_argument('--extraction-python',type=Path,help='Explicit pysam-compatible Python for ABI-separated extraction')
    p.add_argument('--coverage-backend', choices=['cpu','gpu'], default='cpu')
    p.add_argument('--statistics-backend', choices=['cpu','gpu'], default='cpu')
    p.add_argument('--batch-size',type=int,default=16)
    p.add_argument('--depth-sampler',choices=['auto','permutation','lazy'],default='permutation')
    p.add_argument('--alpha',type=float,default=.01)
    p.add_argument('--similarity-block-size',type=int,default=64)
    p.add_argument('--mate-policy',choices=['both','read1','read2'])
    p.add_argument('--sampling-unit',choices=['read','fragment'])
    p.add_argument('--min-mapq',type=int)
    p.add_argument('--exclude-duplicates',action='store_true',default=None)
    p.add_argument('--exclude-qcfail',action='store_true',default=None)
    p.add_argument('--outdir',type=Path,required=True)
    args=p.parse_args()
    if args.reference_sqlite is None and (args.test_alignment is None or args.query_manifest is not None):
        p.error('Existing comparison modes require --test-alignment; --query-manifest requires --reference-sqlite')
    return args


def depth_values(args):
    if args.depth is not None:return [args.depth]
    if args.depths:return [float(x) for x in args.depths.split(',')]
    if args.depth_grid is None:return []
    start,stop,step=map(Decimal,args.depth_grid.split(':'))
    if not start.is_finite() or not stop.is_finite() or not step.is_finite() or start<=0 or stop<start or step<=0:
        raise ValueError('Invalid start:stop:step')
    count=int((stop-start)//step)+1
    if count>100000:raise ValueError('Depth grid exceeds 100000 points')
    values=[start+i*step for i in range(count)]
    if values[-1]!=stop:values.append(stop)
    return list(map(float,values))


def main():
    a=parse_args()
    if a.workers<1 or a.seed<0 or a.replicates<2 or a.batch_size<1 or not 0<a.alpha<1 or a.similarity_block_size<1:
        raise ValueError('Invalid workers/seed/replicates/batch/alpha')
    # Set native thread budgets before importing numerical libraries (also inherited by spawned workers).
    for key in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS']:
        os.environ[key]='1'
    sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
    import numpy as np
    from bin.coverage_core import window_edges,select_units,window_sums
    from bin.concordance_core import simulate_grid,task_rng,evaluate_windows,RNG_POLICY
    from bin.concordance_inputs import load_inputs
    depths=depth_values(a)
    if a.outdir.exists():raise FileExistsError(a.outdir)
    if getattr(a, 'reference_sqlite', None):
        from bin.reference_index import ReferenceIndex
        index=ReferenceIndex(a.reference_sqlite,a.database_root)
        native=index.metadata.get('storage_mode')=='native'
        index.close()
        if native:
            from bin.native_concordance import run_native
            return run_native(a,depths)
        if not depths:raise ValueError('Legacy reference index requires explicit depths')
        from bin.concordance_batch import run_indexed
        return run_indexed(a, depths)
    if not depths:raise ValueError('Explicit depths required outside native SQLite mode')
    started=time.perf_counter()
    cached = None
    if getattr(a, 'reference_depth_manifest', None):
        from bin.ref_depth_worker import load_depth
        first_matrix, first_edges, first_metadata = load_depth(a.reference_depth_manifest, depths[0])
        recorded = first_metadata['sample']
        del first_matrix, first_edges, first_metadata
        if any(not any(row['depth'] == d and row['status'] == 'READY' for row in recorded['depths']) for d in depths):
            raise ValueError('Exact stored depth is unavailable')
        cached = True
        a.reference_bundle = a.reference_depth_manifest.resolve().parent / 'geometry'
        # Database geometry and replicate dimensions are authoritative for this mode.
        a.bins, a.replicates = recorded['bins'], recorded['replicates']
    inputs=load_inputs(a)
    reference_pool,test_pool,layout=inputs.reference_pool,inputs.test_pool,inputs.layout
    if any(not np.isfinite(d) or d<=0 or d*layout.length>min(reference_pool.total_bases,test_pool.total_bases) for d in depths):
        raise ValueError('Requested depth exceeds reference or test aligned mass')
    if len(set(depths))!=len(depths):raise ValueError('Duplicate depth values')
    extraction_seconds=time.perf_counter()-started
    edges=window_edges(layout.length,a.bins)
    a.outdir.mkdir(parents=True);status=a.outdir/'status';status.mkdir();(status/'RUNNING').touch()
    results={}
    gpu = None
    test_gpu = None
    def complete(index,depth,directory,matrix,summaries,calibration,worker_seconds):
        tick=time.perf_counter()
        test_ids=select_units(test_pool,task_rng(a.seed,depth,1,0),target_bases=depth*layout.length,
                              depth_sampler=calibration['sampler'])
        test=(test_gpu._count_selected_batch([test_ids])[0] if test_gpu is not None else window_sums(test_pool,test_ids,edges))
        if int(test.sum())!=int(test_pool.unit_bases[test_ids].sum()):raise AssertionError('Test mass mismatch')
        np.save(directory/'test.npy',test);np.save(directory/'test_ids.npy',test_ids)
        test_seconds=time.perf_counter()-tick
        stat_tick=time.perf_counter()
        result=evaluate_windows(test,matrix,a.alpha,a.similarity_block_size,backend=a.statistics_backend)
        result.update(depth=depth,depth_sampler=calibration,reference_replicates=summaries,
                      test_units=len(test_ids),test_actual_depth=int(test.sum())/layout.length,
                      test_selected_ids_sha256=hashlib.sha256(test_ids.astype('<i8').tobytes()).hexdigest(),
                      worker_simulation_seconds_sum=worker_seconds,test_draw_and_write_seconds=test_seconds,statistics_seconds=time.perf_counter()-stat_tick)
        (directory/'concordance.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
        results[index]=dict(depth=depth,jaccard=result['jaccard'],pearson=result['pearson'],test_actual_depth=result['test_actual_depth'])
    try:
        gpu_started=time.perf_counter()
        if a.coverage_backend == 'gpu':
            from bin.gpu_core import GpuWindowAccumulator
            gpu=GpuWindowAccumulator(reference_pool,edges)
            test_gpu=gpu if test_pool is reference_pool else GpuWindowAccumulator(test_pool,edges)
        if a.statistics_backend == 'gpu':
            from bin.gpu_core import require_cupy
            require_cupy().cuda.Stream.null.synchronize()
        gpu_initialization_seconds=time.perf_counter()-gpu_started
        np.save(a.outdir/'edges.npy',edges)
        if cached is None:
            simulate_grid(reference_pool,edges,depths,a.replicates,a.seed,a.workers,a.executor,
                          a.batch_size,a.depth_sampler,a.outdir,complete,gpu=gpu)
        else:
            for index, depth in enumerate(depths):
                matrix, saved_edges, metadata = load_depth(a.reference_depth_manifest, depth)
                if not np.array_equal(edges, saved_edges):
                    raise ValueError('Stored depth window geometry differs')
                directory = a.outdir / f'depth-{index:05d}'
                directory.mkdir()
                reference = metadata['reference']
                complete(index, depth, directory, matrix, reference['reference_replicates'], reference['sampler'], 0.)
                del matrix
                load_depth(a.reference_depth_manifest, depth)  # Verify cached payloads remained unchanged.
        inputs.verify()
        provenance=dict(schema='coverage-concordance-1',rng_policy=RNG_POLICY,
                        inputs_sha256=inputs.identities,depths=depths,replicates=a.replicates,seed=a.seed,
                        workers=a.workers,executor=a.executor,native_threads_per_worker=1,
                        reference_units=len(reference_pool),test_units=len(test_pool),
                        reference_pool_numeric_bytes=reference_pool.nbytes,test_pool_numeric_bytes=test_pool.nbytes,
                        bins=a.bins,policies=inputs.policies,test_policy='independent-complete-unit-draw-at-matching-depth',
                        elapsed_seconds=time.perf_counter()-started,extraction_seconds=extraction_seconds,
                        python=sys.version,numpy=np.__version__,
                        actual_executor='serial' if a.workers==1 or a.executor=='serial' else a.executor,
                        comparison_scope=inputs.comparison_scope,
                        reference_input_mode='depth_database' if cached is not None else 'bundle' if inputs.bundle_info else 'alignment',
                        reference_depth_manifest=str(a.reference_depth_manifest.resolve()) if cached is not None else None,
                        cached_reference_seed=recorded['seed'] if cached is not None else None,
                        reference_simulations_reused=cached is not None,
                        reference_bundle=inputs.bundle_info,derived_test_layout=inputs.derived_layout,
                        layout_comparison=inputs.layout_comparison,
                        extraction_bridge=inputs.bridge_metadata,
                        coverage_backend=a.coverage_backend,statistics_backend=a.statistics_backend,
                        gpu_initialization_seconds=gpu_initialization_seconds,
                        gpu=gpu.metrics() if gpu is not None else None,
                        statistics_calibrated=False)
        (a.outdir/'provenance.json').write_text(json.dumps(provenance,indent=2,allow_nan=False)+'\n')
        (a.outdir/'summary.json').write_text(json.dumps([results[i] for i in range(len(depths))],indent=2,allow_nan=False)+'\n')
        (status/'RUNNING').unlink();(status/'SUCCESS').touch();(status/'exit_code.txt').write_text('0\n')
    except BaseException:
        (status/'RUNNING').unlink(missing_ok=True);(status/'FAILED').touch();(status/'exit_code.txt').write_text('1\n');raise


if __name__=='__main__':sys.exit(main())
