#!/usr/bin/env python3
"""Plan, build and evaluate stratified genome coverage Bin Counts (geometry only)."""
import argparse
import json
import math
import os
from pathlib import Path
import sys
from time import perf_counter

for _name in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS'):
    os.environ.setdefault(_name,'1')
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from bin import bin_count_eval_io as io
from bin import bin_count_eval_math as maths
import numpy as np


def write_pearson_diagnostics(dest, issues, rows, groups):
    records = []
    for issue in issues:
        row = rows[issue['index']] if issue['kind']=='row' else None
        records.append(dict(kind=issue['kind'], index=issue['index'],
                            row_id=row['row_id'] if row else '',
                            row_key=row['row_key'] if row else '',
                            sample=row['sample'] if row else '',
                            gcf=row['gcf'] if row else groups[issue['index']],
                            reason=issue['reason']))
    io.write_tsv(dest/'pearson_diagnostics.tsv',
                 ['kind','index','row_id','row_key','sample','gcf','reason'],records)


def score_candidate(matrix, labels, metric, rows, candidate, dest, root, jobs):
    print(f"SCORE {candidate['stratum']} B={candidate['bin_count']} {metric}",file=sys.stderr,flush=True)
    diagnostics=[]
    try:
        result,groups,centers,within,per_row,pairs=maths.score(matrix,labels,metric,diagnostics,jobs=jobs)
    except maths.UndefinedCorrelationError as error:
        write_pearson_diagnostics(dest,error.issues,rows,sorted(set(labels)))
        io.write_json(dest/'FAILURE.json',dict(status='UNDEFINED_PEARSON',
            stratum=candidate['stratum'],bin_count=candidate['bin_count'],distance=metric,
            matrices_dir=str(root),matrix_manifest_sha256=io.sha(root/'result.json'),
            code_sha256=io.code_hashes(),reason=str(error),
            diagnostics='pearson_diagnostics.tsv',records_removed=0))
        print(f'Pearson failure details: {dest}/pearson_diagnostics.tsv',
              file=sys.stderr,flush=True)
        raise
    if metric=='correlation':
        write_pearson_diagnostics(dest,diagnostics,rows,groups)
        if diagnostics:
            print(f'WARNING Pearson: {len(diagnostics)} near-constant profiles; '
                  f'see {dest}/pearson_diagnostics.tsv',file=sys.stderr,flush=True)
    result.update(stratum=candidate['stratum'],bin_count=candidate['bin_count'],bin_sizes_bp=candidate['bin_sizes_bp'])
    np.save(dest/'centers.npy',centers,allow_pickle=False)
    io.write_tsv(dest/'centers.tsv',['center_row','gcf'],[dict(center_row=i,gcf=g) for i,g in enumerate(groups)])
    io.write_tsv(dest/'within_gcf.tsv',['gcf','samples','within_mean'],within)
    extra=['pearson_r'] if metric=='correlation' else []
    io.write_tsv(dest/'within_rows.tsv',['row_id','row_key','gcf','distance']+extra,
                 [dict(row_id=r['row_id'],row_key=r['row_key'],gcf=r['gcf'],distance=float(d),
                       **({'pearson_r':1-float(d)} if extra else {})) for r,d in zip(rows,per_row)])
    io.write_tsv(dest/'between_centers.tsv',['left_gcf','right_gcf','distance']+extra,
                 [dict(left_gcf=a,right_gcf=b,distance=d,
                       **({'pearson_r':1-d} if extra else {})) for a,b,d in pairs])
    return result


def evaluation_bytes(n, b, g, jobs, ordination, correlation, *, plots_only=False, perplexity=30):
    base = 32 * 1024**2
    panel = min(n, maths.DISTANCE_PANEL_ROWS)
    tile = maths.DISTANCE_TILE_ROWS
    workers = min(jobs, n)
    profiles = panel + min(n-panel, workers*tile)
    if plots_only:
        if n <= perplexity:
            return base
        keep = min(n, maths.tsne_neighbor_count(n, perplexity) + 2)
        graph = 256 * n * keep
        distances = 8 * n * panel + 24 * n * workers
        # Tied rows may be replayed, with one query and one tile per worker.
        profiles = max(profiles, workers * (1 + min(n, maths.distance_replay_rows(b))))
        centered = 8 * b * profiles if correlation else 0
        return base + graph + distances + centered
    dense = 64 * n * n if ordination else 0
    scoring = (0 if ordination else 8 * g * g) + 8 * g * b
    scoring += max(64 * 1024**2, 32 * b) * max(0, min(jobs, g) - 1)
    centered = 8 * b * profiles if ordination and correlation else 0
    return base + dense + scoring + centered


def check_grid(preflight_dir, strata_config, metrics, ordination, jobs, matrix_gib, working_gib, *,
               index=None, normalized_only=False):
    _, summary, rows, _ = io.frozen_preflight(preflight_dir)
    if index is not None:
        path, fingerprint = io.validate_frozen_index(index, summary)
        print(f'Validated frozen SQLite: {path} (SHA-256 {fingerprint})', flush=True)
    layouts = {}
    for row in rows:
        length = int(row['genome_length_bp'])
        layout = layouts.setdefault(row['gcf'], {'lengths': [length]})
        if length <= 0 or layout['lengths'][0] != length:
            raise ValueError(f"Inconsistent or nonpositive genome length for {row['gcf']}")
    config = io.read_strata_config(strata_config)
    breaks = [item['min_genome_bp'] for item in config]
    strata = io.strata_plan(rows, layouts, breaks, None, config)
    candidates = [(layer, candidate) for layer in strata for candidate in layer['candidates']]
    if not candidates:
        raise ValueError('No matrix candidates')
    matrix_bytes = sum(candidate['raw_plus_normalized_bytes'] for _, candidate in candidates)
    matrix_bytes = io.matrix_payload_bytes(matrix_bytes, normalized_only)
    working_bytes = 0
    plots_only = ordination == 'tsne'
    for metric in metrics:
        minimum = 2 if plots_only or metric == 'correlation' else 1
        selected = [(layer, candidate) for layer, candidate in candidates
                    if candidate['bin_count'] >= minimum]
        if not selected:
            raise ValueError(f'No {metric} candidates with B >= {minimum}')
        working_bytes = max(working_bytes, *(evaluation_bytes(
            len(layer['row_indices']), candidate['bin_count'], len(layer['gcfs']),
            jobs, plots_only, metric == 'correlation', plots_only=plots_only)
            for layer, candidate in selected))
    gib = 1024**3
    print(f'Projected final matrices: {matrix_bytes/gib:.2f} GiB; '
          f'shards + merged: {2*matrix_bytes/gib:.2f} GiB; '
          f'peak evaluation per job: {working_bytes/gib:.2f} GiB', flush=True)
    io.check_sharded_matrix_budget(matrix_bytes, matrix_gib)
    if not math.isfinite(working_gib) or working_gib <= 0 or working_bytes > working_gib * gib:
        raise ValueError(f'Projected evaluation needs {working_bytes/gib:.2f} GiB within --max-working-gib')


def evaluate(matrices_dir, outdir, metrics, ordinations, seeds, perplexity, iterations,
             learning_rate, jobs, max_working_gib, *, tasks=None, min_bin_count=1, plots_only=False):
    import scipy
    import sklearn
    from threadpoolctl import threadpool_limits
    from bin.bin_count_eval_render import scatter
    metrics, ordinations, seeds = list(dict.fromkeys(metrics)), list(dict.fromkeys(ordinations)), list(dict.fromkeys(seeds))
    if plots_only and (ordinations != ['tsne'] or tasks is not None):
        raise ValueError('Plot-only evaluation requires t-SNE on one node')
    selected = None
    if tasks is not None:
        selected = {}
        for folder, metric in tasks:
            if metric in selected.setdefault(folder, []) or metric not in metrics:
                raise ValueError('Duplicate or unknown evaluation task')
            selected[folder].append(metric)
    root, manifest = io.load(matrices_dir, 'matrices', folders=selected)
    if min_bin_count < 1:
        raise ValueError('Minimum bin count must be positive')
    if plots_only:
        min_bin_count = max(2, min_bin_count)
    candidates = [c for c in manifest['candidates'] if c['bin_count'] >= min_bin_count]
    if not candidates and min_bin_count > 1:
        raise ValueError('No matrix candidates meet the minimum bin count')
    if selected is not None:
        candidates = [c for c in candidates if c['folder'] in selected]
        if len(candidates) != len(selected):
            raise ValueError('Unknown evaluation candidate')
    if (not metrics or not set(metrics)<=set(maths.METRICS) or not ordinations or
            not set(ordinations)<={'tsne','pcoa','none'} or ('none' in ordinations and len(ordinations)>1)):
        raise ValueError('Invalid distance/ordination options')
    if (not seeds or any(not 0<=s<2**32 for s in seeds) or not math.isfinite(perplexity) or perplexity<=0 or
            iterations<250 or jobs<1 or not math.isfinite(max_working_gib) or max_working_gib<=0 or
            (learning_rate!='auto' and (not math.isfinite(learning_rate) or learning_rate<=0))):
        raise ValueError('Invalid ordination/resource parameters')
    for candidate in candidates:
        n,b=candidate['shape']
        rows=io.read_tsv(io.child(root,candidate['folder'])/'rows.tsv')
        g=len({r['gcf'] for r in rows})
        required = evaluation_bytes(n, b, g, jobs, ordinations != ['none'],
                                    'correlation' in metrics, plots_only=plots_only,
                                    perplexity=perplexity)
        if required>max_working_gib*1024**3:
            raise ValueError(f'Projected dense distance/ordination working memory {required/1024**3:.2f}GiB exceeds --max-working-gib')
    out=io.new_directory(outdir)
    scores, embedded = [], []
    with threadpool_limits(limits=jobs):
        for candidate in candidates:
            folder=io.child(root,candidate['folder'])
            rows=io.read_tsv(folder/'rows.tsv')
            with io.mapped_matrix(folder/'normalized.npy') as matrix:
                if list(matrix.shape)!=candidate['shape'] or matrix.dtype!=np.dtype('float64') or len(rows)!=len(matrix):
                    raise ValueError('Matrix artifact shape/dtype mismatch')
                for start in range(0,len(matrix),16):
                    if not np.allclose(matrix[start:start+16].sum(axis=1),manifest['normalization_target'],rtol=1e-12,atol=1e-7):
                        raise ValueError('Normalized row sums differ from recorded target')
                labels=[r['gcf'] for r in rows]
                for metric in metrics if selected is None else selected[candidate['folder']]:
                    dest=out/candidate['folder']/metric
                    dest.mkdir(parents=True)
                    if not plots_only:
                        scores.append(score_candidate(matrix, labels, metric, rows, candidate,
                                                      dest, root, jobs))
                    if ordinations==['none']:
                        continue
                    started = perf_counter()
                    if plots_only:
                        distance_matrix = (maths.tsne_neighbor_graph(matrix, metric, perplexity, jobs)
                                           if len(rows) > perplexity else None)
                    else:
                        distance_matrix = maths.pairwise(matrix, metric, jobs=jobs)
                    calculated = perf_counter()
                    if not plots_only:
                        np.save(dest/'pairwise_distances.npy',distance_matrix,allow_pickle=False)
                    print(f"TIMING pairwise {candidate['stratum']} B={candidate['bin_count']} {metric} "
                          f"compute={calculated-started:.2f}s write={perf_counter()-calculated:.2f}s",
                          file=sys.stderr,flush=True)
                    for method in ordinations:
                        for seed in seeds if method=='tsne' else [None]:
                            entry=dict(stratum=candidate['stratum'],bin_count=candidate['bin_count'],distance=metric,
                                       ordination=method,seed=seed,status='COMPLETE')
                            if method=='tsne' and len(rows)<=perplexity:
                                entry.update(status='UNAVAILABLE_PERPLEXITY',reason='Require perplexity < record count; no auto adjustment')
                                embedded.append(entry);continue
                            view=dest/(method+(f'_seed{seed}' if seed is not None else ''))
                            view.mkdir()
                            if method=='pcoa':
                                coords,details=maths.pcoa(distance_matrix)
                            else:
                                coords,details=maths.tsne(distance_matrix,seed,perplexity,iterations,learning_rate,jobs)
                            if not np.isfinite(coords).all():
                                raise ValueError('Nonfinite ordination coordinates')
                            io.write_tsv(view/'coordinates.tsv',['row_id','row_key','sample','gcf','x','y'],
                                [dict(row_id=r['row_id'],row_key=r['row_key'],sample=r['sample'],gcf=r['gcf'],x=float(x),y=float(y))
                                 for r,(x,y) in zip(rows,coords)])
                            io.write_json(view/'ordination.json',details)
                            scatter(coords,rows,view,f"{candidate['stratum']} | B={candidate['bin_count']} | {metric} | {method}",method)
                            entry.update(folder=str(view.relative_to(out)),details=details)
                            embedded.append(entry)
                    del distance_matrix
    metadata = dict(kind='evaluation' if tasks is None else 'evaluation-shard', schema=1,
        matrices_dir=str(root), matrix_manifest_sha256=io.sha(root/'result.json'),
        distance=metrics,ordination=ordinations,seeds=seeds,perplexity=perplexity,iterations=iterations,
        learning_rate=learning_rate,jobs=jobs,max_working_gib=max_working_gib,min_bin_count=min_bin_count,
        embeddings=embedded,
        scores_space='original normalized coverage vectors',scoring_does_not_depend_on_ordination=True,
        versions=dict(scipy=scipy.__version__,sklearn=sklearn.__version__),
        status='COMPLETE_WITH_ORDINATION_GAPS' if any(e['status']!='COMPLETE' for e in embedded) else 'COMPLETE')
    if plots_only:
        metadata['plots_only'] = True
        metadata.pop('scores_space')
        metadata.pop('scoring_does_not_depend_on_ordination')
    elif tasks is None:
        write_score_outputs(out, scores, metrics)
    else:
        io.write_json(out/'scores.json', scores)
        metadata['tasks'] = tasks
    return io.seal(out, metadata)


def write_score_outputs(out, scores, metrics):
    from bin.bin_count_eval_render import curves

    io.write_json(out/'scores.json',scores)
    columns=['stratum','bin_count','bin_sizes_bp','distance','gcfs','records','within_mean','between_mean','ratio','status']
    io.write_tsv(out/'scores.tsv',columns,[{k:(','.join(map(str,r[k])) if k=='bin_sizes_bp' else r[k]) for k in columns} for r in scores])
    winners=[]
    for stratum in sorted({r['stratum'] for r in scores}):
        for metric in metrics:
            eligible=[r for r in scores if r['stratum']==stratum and r['distance']==metric and r['status']=='OK']
            maximum=max((r['ratio'] for r in eligible),default=None)
            winners.extend(dict(stratum=stratum,distance=metric,bin_count=r['bin_count'],ratio=r['ratio'])
                           for r in eligible if math.isclose(r['ratio'],maximum,rel_tol=1e-12,abs_tol=1e-12))
    io.write_json(out/'candidate_maxima.json',winners)
    curves(scores,out)


def add_plan_options(p):
    p.add_argument('--preflight-dir',required=True)
    p.add_argument('--index',required=True)
    p.add_argument('--strata-config',type=Path,help='JSON file defining each stratum and its nominal bin widths')
    p.add_argument('--genome-size-breaks-bp',type=int,nargs='+')
    p.add_argument('--bin-sizes-bp',type=int,nargs='+',default=None)
    p.add_argument('--bin-size-start-bp',type=int)
    p.add_argument('--bin-size-stop-bp',type=int)
    p.add_argument('--bin-size-step-bp',type=int)
    p.add_argument('--source-root-map',help='JSON file mapping known source IDs to absolute relocated roots')


def add_eval_options(p):
    p.add_argument('--distance',nargs='+',choices=maths.METRICS,default=['euclidean'],
                   help='correlation = 1 - Pearson r; all scores use normalized coverage space')
    p.add_argument('--ordination',nargs='+',choices=['tsne','pcoa','none'],default=['tsne'])
    p.add_argument('--seeds',type=int,nargs='+',default=[20261103])
    p.add_argument('--perplexity',type=float,default=30)
    p.add_argument('--iterations',type=int,default=1000)
    p.add_argument('--learning-rate',default='auto',help='auto or positive float')
    p.add_argument('--jobs',type=int,default=1,help='Worker limit for build (run), scoring, and pairwise distances')
    p.add_argument('--max-working-gib',type=float,default=8)
    p.add_argument('--min-bin-count',type=int,default=1,help='Score only matrix candidates at or above this bin count')


def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__)
    sub=p.add_subparsers(dest='command',required=True)
    plan=sub.add_parser('plan',help='Validate frozen metadata and expand candidates; no geometry extraction')
    add_plan_options(plan)
    check = sub.add_parser('check-grid', help='Check frozen SQLite and resource budgets before submitting')
    check.add_argument('--preflight-dir', required=True)
    check.add_argument('--index', help='Validate this SQLite fingerprint against the frozen preflight')
    check.add_argument('--strata-config', required=True)
    check.add_argument('--distance', nargs='+', choices=maths.METRICS, required=True)
    check.add_argument('--ordination', choices=['none', 'tsne'], required=True)
    check.add_argument('--jobs', type=int, required=True)
    check.add_argument('--max-matrix-gib', type=float, required=True)
    check.add_argument('--max-working-gib', type=float, required=True)
    build=sub.add_parser('build',help='Build raw+normalized matrices using verified geometry only')
    build.add_argument('--plan',required=True)
    build.add_argument('--max-matrix-gib',type=float,default=64)
    build.add_argument('--jobs',type=int,default=1,help='Concurrent sample-GCF records')
    shard = sub.add_parser('build-shard', help='Build one independent row shard')
    merge = sub.add_parser('merge-build', help='Validate all row shards and publish one matrix set')
    for stage in (shard, merge):
        stage.add_argument('--plan', required=True)
        stage.add_argument('--shards', type=int, required=True)
        stage.add_argument('--max-matrix-gib', type=float, default=64)
        stage.add_argument('--jobs', type=int, default=1, help='Worker count per build shard')
    shard.add_argument('--rank', type=int, required=True, help='Zero-based shard index')
    shard.add_argument('--shard-dir', required=True, help='Shared directory for independent shard outputs')
    merge.add_argument('--shard-dir', required=True)
    merge.add_argument('--outdir', required=True, help='New exclusive matrix output directory')
    ev=sub.add_parser('evaluate',help='Reuse saved matrices for distance scores and ordinations')
    ev.add_argument('--matrices',required=True)
    add_eval_options(ev)
    ev.add_argument('--plots-only', action='store_true', help='Write t-SNE views for B >= 2 without scores')
    run=sub.add_parser('run',help='Run plan, build and evaluate into three separate subdirectories')
    add_plan_options(run);add_eval_options(run)
    run.add_argument('--max-matrix-gib',type=float,default=64)
    for parser in (check, build, shard, merge, run):
        parser.add_argument('--normalized-only', action='store_true',
                            help='Store and budget only normalized matrices; use on all shard and merge stages')
    for parser in (build, ev, run):
        parser.add_argument('--nodes',type=int,default=1,help='Slurm nodes, one worker per node')
    for parser in (plan,build,ev,run):
        parser.add_argument('--outdir',required=True,help='New exclusive output directory')
    a=p.parse_args(argv)
    try:
        if a.command == 'check-grid':
            if a.jobs < 1:
                raise ValueError('Jobs must be positive')
            check_grid(a.preflight_dir, a.strata_config, a.distance, a.ordination,
                       a.jobs, a.max_matrix_gib, a.max_working_gib,
                       index=a.index, normalized_only=a.normalized_only)
            return
        if a.command == 'evaluate' and a.plots_only and (a.nodes != 1 or a.ordination != ['tsne']):
            raise ValueError('--plots-only requires --ordination tsne and one node')
        if a.command in ('build-shard', 'merge-build'):
            from bin.bin_count_multinode import run_build_shard_stage
            output = run_build_shard_stage(a)
            print(f'COMPLETE {output}', flush=True)
            return
        if a.command in ('build', 'evaluate', 'run'):
            if a.nodes < 1:
                raise ValueError('Nodes must be positive')
            if a.nodes > 1:
                from bin.bin_count_multinode import check_allocation
                check_allocation(a.nodes, a.jobs)
        if a.command in ('evaluate', 'run') and a.min_bin_count < 1:
            raise ValueError('Minimum bin count must be positive')
        if a.command in ('plan','run'):
            legacy_sizes = (a.bin_sizes_bp, a.bin_size_start_bp, a.bin_size_stop_bp, a.bin_size_step_bp)
            strata_config = None
            if a.strata_config is None:
                breaks = a.genome_size_breaks_bp or [1000,10000,100000]
                bin_sizes = io.sizes(*legacy_sizes)
            else:
                if any(value is not None for value in (a.genome_size_breaks_bp, *legacy_sizes)):
                    raise ValueError('--strata-config cannot be combined with legacy bin-width or stratum arguments')
                strata_config = io.read_strata_config(a.strata_config)
                breaks = [item['min_genome_bp'] for item in strata_config]
                bin_sizes = None
            mapping=json.loads(Path(a.source_root_map).read_text()) if a.source_root_map else None
            if mapping is not None and not isinstance(mapping,dict):
                raise ValueError('Root map must be a JSON object')
            output=io.new_directory(a.outdir) if a.command=='run' else Path(a.outdir)
            io.create_plan(a.preflight_dir,a.index,output/'plan' if a.command=='run' else output,
                           breaks,bin_sizes,mapping,strata_config=strata_config)
        if a.command == 'run' and a.nodes > 1:
            from bin.bin_count_multinode import run_multinode
            run_multinode(output, a)
            print(f'COMPLETE {output}', flush=True)
            return
        if a.command in ('build', 'evaluate') and a.nodes > 1:
            from bin.bin_count_multinode import run_stage
            run_stage(a)
            print(f'COMPLETE {Path(a.outdir).absolute()}', flush=True)
            return
        if a.command in ('build','run'):
            io.build(output/'plan' if a.command=='run' else a.plan,
                     output/'matrices' if a.command=='run' else a.outdir,a.max_matrix_gib,
                     jobs=a.jobs, normalized_only=a.normalized_only)
        if a.command in ('evaluate','run'):
            learning_rate=a.learning_rate if a.learning_rate=='auto' else float(a.learning_rate)
            evaluate(output/'matrices' if a.command=='run' else a.matrices,
                     output/'evaluation' if a.command=='run' else a.outdir,
                     a.distance,a.ordination,a.seeds,a.perplexity,a.iterations,learning_rate,a.jobs,
                     a.max_working_gib,min_bin_count=a.min_bin_count,
                     plots_only=a.command == 'evaluate' and a.plots_only)
        if a.command=='run':
            io.seal(output,dict(kind='run',schema=1,status=json.loads((output/'evaluation/result.json').read_text())['status'],components=['plan','matrices','evaluation']))
        print(f'COMPLETE {Path(a.outdir).absolute()}',flush=True)
    except (ValueError,OSError,KeyError,TypeError,ImportError) as error:
        p.exit(1,f'ERROR: {error}\n')


if __name__=='__main__':
    main()
