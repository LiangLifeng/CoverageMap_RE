#!/usr/bin/env python3
"""Slurm orchestration for exact row-sharded build and task-sharded evaluation."""
import argparse
from contextlib import ExitStack
import hashlib
from io import BytesIO
import json
import os
from pathlib import Path
import shutil
import socket
import subprocess
import sys
from time import perf_counter

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from bin import bin_count_eval_io as io
from bin.bin_count_evaluate import evaluate, write_score_outputs


def check_allocation(nodes, jobs):
    if not (nodes > 1 and 1 <= jobs <= int(os.environ.get('SLURM_CPUS_PER_TASK', '0'))
            and nodes == int(os.environ.get('SLURM_JOB_NUM_NODES', '0'))):
        raise ValueError('Multi-node run requires a Slurm allocation with exactly --nodes nodes '
                         'and --cpus-per-task at least --jobs')


def distribute(candidates, metrics, ordinations, nodes):
    """Assign longest estimated tasks first, with a stable tie-break."""
    costs = []
    for candidate in candidates:
        n, b = candidate['shape']
        for metric in metrics:
            weight = n * b
            if ordinations != ['none']:
                weight += n * n * b * (6 if metric == 'jensenshannon' else 1)
                if 'pcoa' in ordinations:
                    weight += n ** 3
            costs.append((-weight, candidate['folder'], metric))
    tasks, load = [[] for _ in range(nodes)], [0] * nodes
    for neg_weight, folder, metric in sorted(costs):
        rank = min(range(nodes), key=load.__getitem__)
        tasks[rank].append([folder, metric])
        load[rank] -= neg_weight
    return tasks


def rank_folder(work, phase, rank):
    return work / phase / f'rank-{rank:04d}'


def worker(work, dispatch, phase, rank):
    config = json.loads(dispatch.read_text())
    nodes = config['nodes']
    if not 0 <= rank < nodes:
        raise ValueError('Worker rank is outside the dispatch')
    target = rank_folder(work, phase, rank)
    started = perf_counter()
    print(f'RANK {phase} {rank}/{nodes} {socket.gethostname()}', flush=True)
    if phase == 'build':
        io.build(config['plan_dir'], target, config['max_matrix_gib'],
                 jobs=config['jobs'], row_shard=(rank, nodes),
                 normalized_only=config.get('normalized_only', False))
    else:
        evaluate(config['matrices_dir'], target, **config['evaluation'], tasks=config['tasks'][rank])
    io.write_json(work / f'{phase}-{rank:04d}.json',
                  dict(rank=rank, host=socket.gethostname(), seconds=perf_counter()-started))


def launch_slurm(work, dispatch, phase, nodes, jobs):
    command = ['/usr/bin/srun', f'--nodes={nodes}', f'--ntasks={nodes}',
               '--ntasks-per-node=1', f'--cpus-per-task={jobs}', '--cpu-bind=cores',
               '--kill-on-bad-exit=1', sys.executable, '-B', str(Path(__file__).resolve()),
               str(work), str(dispatch), phase]
    try:
        subprocess.run(command, check=True)
    except subprocess.CalledProcessError as error:
        raise ValueError(f'Slurm {phase} workers failed with exit {error.returncode}') from error
    hosts = {json.loads((work/f'{phase}-{rank:04d}.json').read_text())['host'] for rank in range(nodes)}
    if len(hosts) != nodes:
        raise ValueError(f'Expected {nodes} distinct worker hosts for {phase}, observed {len(hosts)}')


MERGE_BLOCK_BYTES = 32 * 1024**2
# a95bfe7 and eceb78b built compatible raw+normalized shards.
LEGACY_MERGE_CODE = {
    'bin_count_evaluate.py': '7ab3af216cd5609b5457f679493b0140ea9110c396d891c724a37352f6606947',
    'bin_count_eval_io.py': '4723d88304c442d108b5e27309a6db4cc064d147df81562cbe433b34da6291e9',
    'bin_count_multinode.py': '23922c3145275a797ff2d1308004bea3f6a48b0043b42538f3b381ed4389c588',
}
STREAMED_MERGE_CODE = {
    **LEGACY_MERGE_CODE,
    'bin_count_eval_io.py': '4ad2a22ab4ff04672b933327f9fafbe45da99ad36a4fb8f57de3a43e934310b5',
    'bin_count_multinode.py': 'fce3fa372580991a5fbf6fdf6d4e5de85e67251d39c6c9ea2d32f5e02ab0d658',
}
# Last code identity before gcf-concat-v1, for frozen legacy shards.
PRE_CONCAT_MERGE_CODE = {
    'coverage_core.py': 'f57e93a9bddfa27b2a72561f77aad56af1406f691a36ed81f085d065d2678ee9',
    'reference_store.py': '3c51b4e0d465d28fe5e5695b08d2e38310f5114582ed2441a945a1cb2f0a1595',
    'bin_count_eval_io.py': '7fd6647ce9459ec7284b8ce1d8d9f5716751dc49937734460fc96737a3357ee2',
    'bin_count_multinode.py': 'b32e6efa6ef98555f230200f3c7a374ed733385f88f048925bd8ccc2cc684435',
}


def copy_verified_edges(source, target, expected):
    fingerprint = hashlib.sha256()
    with source.open('rb') as src, target.open('wb') as dst:
        for block in iter(lambda: src.read(1024**2), b''):
            fingerprint.update(block)
            dst.write(block)
    actual = fingerprint.hexdigest()
    if actual != expected:
        raise ValueError(f'Changed artifact: {source}')
    return actual


def merge_array(target, inputs, shape):
    """Consume each shard in order and write canonical row-major NPY bytes once."""
    import numpy as np

    nrows, bins = shape
    dtype = np.dtype('float64')
    nodes = len(inputs)
    batch_rows = min(nrows, max(1, MERGE_BLOCK_BYTES // (bins * dtype.itemsize)))
    batch = np.empty((batch_rows, bins), dtype=dtype)
    rank_rows = (batch_rows + nodes - 1) // nodes
    scratch = np.empty((rank_rows, bins), dtype=dtype) if rank_rows > 1 and nodes > 1 else None

    header = BytesIO()
    np.lib.format.write_array_header_1_0(
        header, dict(descr=np.lib.format.dtype_to_descr(dtype), fortran_order=False, shape=tuple(shape)))
    header_bytes = header.getvalue()
    fingerprint = hashlib.sha256(header_bytes)
    with ExitStack() as stack:
        streams = {}
        for rank, entry in enumerate(inputs):
            if entry is None:
                continue
            path, expected, count = entry
            # Inspect the existing NPY header without reading the mapped payload.
            with io.mapped_matrix(path) as matrix:
                if (matrix.shape != (count, bins) or matrix.dtype != dtype
                        or not matrix.flags.c_contiguous):
                    raise ValueError(f'Matrix shard shape/dtype/layout mismatch: {path}')
                offset = matrix.offset
            stream = stack.enter_context(path.open('rb'))
            checksum = hashlib.sha256(stream.read(offset))
            streams[rank] = stream, checksum, path, expected

        with target.open('wb') as output:
            output.write(header_bytes)
            for first in range(0, nrows, batch_rows):
                block = batch[:min(batch_rows, nrows-first)]
                for rank, (source, checksum, path, _) in streams.items():
                    rows = block[(rank-first) % nodes::nodes]
                    if not len(rows):
                        continue
                    destination = rows if rows.flags.c_contiguous else scratch[:len(rows)]
                    raw = memoryview(destination).cast('B')
                    if source.readinto(raw) != len(raw):
                        raise ValueError(f'Truncated matrix shard: {path}')
                    checksum.update(raw)
                    if destination is not rows:
                        rows[:] = destination
                raw = memoryview(block).cast('B')
                output.write(raw)
                fingerprint.update(raw)

        for stream, checksum, path, expected in streams.values():
            for trailing in iter(lambda: stream.read(1024**2), b''):
                checksum.update(trailing)
            if checksum.hexdigest() != expected:
                raise ValueError(f'Changed artifact: {path}')
    return fingerprint.hexdigest()


def merge_matrix_payloads(out, plan, shards, available, nodes, normalized_only):
    remaining = [dict(manifest['artifact_sha256']) for _, manifest in shards]
    written_sha256 = {}
    merged = []
    columns = ['row_id', 'row_key', 'sample', 'gcf', 'reference_id']
    names = ('normalized.npy',) if normalized_only else ('raw_mean_depth.npy', 'normalized.npy')
    for layer in plan['strata']:
        rows = [plan['rows'][i] for i in layer['row_indices']]
        expected_rows = [{**{k: str(row[k]) for k in columns[:-1]}, 'reference_id': str(row['id'])}
                         for row in rows]
        for candidate in layer['candidates']:
            b = candidate['bin_count']
            relative = f"{layer['id']}/B{b}"
            target = out / relative
            target.mkdir(parents=True)
            expected = dict(stratum=layer['id'], bin_count=b, bin_sizes_bp=candidate['bin_sizes_bp'],
                            folder=relative, shape=[len(rows), b])
            io.write_tsv(target/'rows.tsv', columns, expected_rows)
            inputs = {name: [None] * nodes for name in names}
            edge_key = f'{relative}/edges.json'
            for rank, (root, _) in enumerate(shards):
                expected_slice = expected_rows[rank::nodes]
                if not expected_slice:
                    continue
                metadata = available[rank].pop(relative, None)
                if metadata != {**expected, 'shape': [len(expected_slice), b]}:
                    raise ValueError(f'Missing or unexpected matrix shard for {relative}, rank {rank}')
                row_key = f'{relative}/rows.tsv'
                row_file = io.child(root, row_key)
                if io.sha(row_file) != remaining[rank].pop(row_key):
                    raise ValueError(f'Changed artifact: {row_file}')
                if io.read_tsv(row_file) != expected_slice:
                    raise ValueError(f'Matrix shard row order mismatch for {relative}, rank {rank}')
                edge_file = io.child(root, edge_key)
                edge_expected = remaining[rank].pop(edge_key)
                if rank == 0:
                    written_sha256[edge_key] = copy_verified_edges(edge_file, target/'edges.json', edge_expected)
                else:
                    if io.sha(edge_file) != edge_expected:
                        raise ValueError(f'Changed artifact: {edge_file}')
                    if edge_expected != written_sha256[edge_key]:
                        raise ValueError(f'Matrix shard edges mismatch for {relative}, rank {rank}')
                for name in names:
                    key = f'{relative}/{name}'
                    inputs[name][rank] = (io.child(root, key), remaining[rank].pop(key), len(expected_slice))
            for name in names:
                written_sha256[f'{relative}/{name}'] = merge_array(target/name, inputs[name], expected['shape'])
            merged.append(expected)
    if any(available):
        raise ValueError('Unexpected matrix shard candidates')
    for (root, _), unused in zip(shards, remaining):
        for name, fingerprint in unused.items():
            if io.sha(io.child(root, name)) != fingerprint:
                raise ValueError(f'Changed artifact: {name}')
    return merged, written_sha256


def merge_matrices(plan_dir, work, output, nodes, jobs, *, normalized_only=False):
    import numpy as np

    plan_root, plan = io.load(plan_dir, 'plan')
    policy = io.matrix_policy(plan)
    plan_hash = io.sha(plan_root/'result.json')
    code_hashes = io.code_hashes()
    shards = [io.read_manifest(rank_folder(work, 'build', rank), 'matrix-shard') for rank in range(nodes)]
    for rank, (_, manifest) in enumerate(shards):
        if manifest.get('normalized_only', False) is not normalized_only:
            raise ValueError(f'Shard {rank} storage mode does not match merge; use --normalized-only on all stages')
    shard_code = shards[0][1]['code_sha256']
    compatible_code = [code_hashes]
    if policy == io.LEGACY_MATRIX_WINDOWS:
        legacy_base = dict(code_hashes, **PRE_CONCAT_MERGE_CODE)
        del legacy_base['gcf_windows.py']
        compatible_code.append(legacy_base)
        if not normalized_only:
            compatible_code.extend({**base, **previous}
                                   for base in (legacy_base, code_hashes)
                                   for previous in (LEGACY_MERGE_CODE, STREAMED_MERGE_CODE))
    if shard_code not in compatible_code:
        raise ValueError('Matrix shards were built with a different code version')
    available = []
    for rank, (root, manifest) in enumerate(shards):
        candidates = {c['folder']: c for c in manifest['candidates']}
        if (manifest['plan_manifest_sha256'] != plan_hash or manifest['row_shard'] != [rank, nodes]
                or manifest['build_jobs'] != jobs or manifest['normalization_target'] != 1e6
                or manifest['code_sha256'] != shard_code or manifest['numpy'] != np.__version__
                or manifest['window_policy'] != policy
                or manifest['source_geometry_verified'] is not True
                or len(candidates) != len(manifest['candidates'])):
            raise ValueError('Invalid or duplicate matrix shard')
        available.append(candidates)

    out = io.new_directory(output)
    try:
        merged, written_sha256 = merge_matrix_payloads(out, plan, shards, available, nodes, normalized_only)
        metadata = dict(kind='matrices', schema=1, plan_dir=str(plan_root),
            plan_manifest_sha256=plan_hash, candidates=merged, normalization_target=1e6, build_jobs=jobs,
            window_policy=policy, source_geometry_verified=True)
        if normalized_only:
            metadata['normalized_only'] = True
        return io.seal(out, metadata, written_sha256=written_sha256)
    except BaseException:
        shutil.rmtree(out, ignore_errors=True)
        raise


def merge_evaluation(matrices_dir, work, output, tasks, options):
    matrix_root, matrices = io.read_manifest(matrices_dir, 'matrices')
    matrix_hash = io.sha(matrix_root/'result.json')
    folders = {(c['stratum'], c['bin_count']): c['folder'] for c in matrices['candidates']}
    expected = [(c['folder'], metric) for c in matrices['candidates']
                if c['bin_count'] >= options['min_bin_count'] for metric in options['metrics']]
    if sorted(tuple(task) for rank_tasks in tasks for task in rank_tasks) != sorted(expected):
        raise ValueError('Evaluation dispatch has missing or duplicate tasks')

    scores, embeddings, artifacts = {}, {}, {}
    metadata = None
    option_fields = {'metrics': 'distance', 'ordinations': 'ordination'}
    for rank, rank_tasks in enumerate(tasks):
        root, manifest = io.load(rank_folder(work, 'evaluate', rank), 'evaluation-shard')
        if (manifest['matrix_manifest_sha256'] != matrix_hash or manifest['tasks'] != rank_tasks
                or any(manifest[option_fields.get(k, k)] != v for k, v in options.items())):
            raise ValueError(f'Evaluation shard provenance mismatch for rank {rank}')
        if metadata is not None and manifest['versions'] != metadata['versions']:
            raise ValueError('Evaluation ranks use different library versions')
        metadata = manifest
        assigned = {tuple(task) for task in rank_tasks}
        observed = set()
        for score in json.loads((root/'scores.json').read_text()):
            key = folders[score['stratum'], score['bin_count']], score['distance']
            if key not in assigned or key in observed:
                raise ValueError('Unexpected or duplicate score in evaluation shard')
            scores[key] = score
            artifacts[key] = io.child(root, f'{key[0]}/{key[1]}')
            observed.add(key)
        if observed != assigned:
            raise ValueError('Evaluation shard has missing scores')
        for entry in manifest['embeddings']:
            key = folders[entry['stratum'], entry['bin_count']], entry['distance']
            if key not in assigned:
                raise ValueError('Unexpected embedding in evaluation shard')
            embeddings.setdefault(key, []).append(entry)
        expected_views = [(method, seed) for method in options['ordinations'] if method != 'none'
                          for seed in (options['seeds'] if method == 'tsne' else [None])]
        for key in assigned:
            if [(e['ordination'], e['seed']) for e in embeddings.get(key, [])] != expected_views:
                raise ValueError(f'Missing or duplicate embeddings for {key}')

    out = io.new_directory(output)
    for key in expected:
        target = out / key[0] / key[1]
        target.parent.mkdir(parents=True, exist_ok=True)
        artifacts[key].rename(target)
    write_score_outputs(out, [scores[key] for key in expected], options['metrics'])
    embedded = [entry for key in expected for entry in embeddings.get(key, [])]
    metadata = {k: v for k, v in metadata.items() if k not in ('tasks', 'artifact_sha256')}
    metadata.update(kind='evaluation', embeddings=embedded,
                    status='COMPLETE_WITH_ORDINATION_GAPS' if any(e['status'] != 'COMPLETE' for e in embedded) else 'COMPLETE')
    return io.seal(out, metadata)


def check_build_disk(plan_dir, max_matrix_gib, *, normalized_only=False):
    plan_root, plan = io.read_manifest(plan_dir, 'plan')
    payload_bytes = io.matrix_payload_bytes(plan['planned_matrix_bytes'], normalized_only)
    io.check_sharded_matrix_budget(payload_bytes, max_matrix_gib)
    return plan_root


def run_build_shard_stage(args):
    if args.shards < 2 or args.jobs < 1:
        raise ValueError('Independent build requires at least two shards and positive jobs')
    normalized_only = getattr(args, 'normalized_only', False)
    plan = check_build_disk(args.plan, args.max_matrix_gib, normalized_only=normalized_only)
    work = Path(args.shard_dir).resolve()
    if args.command == 'build-shard':
        if not 0 <= args.rank < args.shards:
            raise ValueError('Shard rank must be between zero and --shards minus one')
        target = rank_folder(work, 'build', args.rank)
        io.build(plan, target, args.max_matrix_gib, jobs=args.jobs,
                 row_shard=(args.rank, args.shards), normalized_only=normalized_only)
    else:
        target = Path(args.outdir).resolve()
        if target == work or target.is_relative_to(work):
            raise ValueError('Merged output must be outside the shard directory')
        merge_matrices(plan, work, target, args.shards, args.jobs, normalized_only=normalized_only)
    return target


def build_phase(work, plan_dir, matrices_dir, args, launch):
    normalized_only = getattr(args, 'normalized_only', False)
    plan_root = check_build_disk(plan_dir, args.max_matrix_gib, normalized_only=normalized_only)
    dispatch = work/'dispatch.json'
    io.write_json(dispatch, dict(nodes=args.nodes, jobs=args.jobs, plan_dir=str(plan_root),
                                max_matrix_gib=args.max_matrix_gib, normalized_only=normalized_only))
    started = perf_counter()
    launch(work, dispatch, 'build', args.nodes, args.jobs)
    built = perf_counter()
    merge_matrices(plan_root, work, matrices_dir, args.nodes, args.jobs, normalized_only=normalized_only)
    merged = perf_counter()
    shutil.rmtree(work/'build')
    return dict(build_workers=built-started, matrix_merge=merged-built,
                build_shard_cleanup=perf_counter()-merged)


def evaluate_phase(work, matrices_dir, output, args, launch):
    matrix_root, matrices = io.read_manifest(matrices_dir, 'matrices')
    options = dict(metrics=list(dict.fromkeys(args.distance)),
                   ordinations=list(dict.fromkeys(args.ordination)), seeds=list(dict.fromkeys(args.seeds)),
                   perplexity=args.perplexity, iterations=args.iterations,
                   learning_rate=args.learning_rate if args.learning_rate == 'auto' else float(args.learning_rate),
                   jobs=args.jobs, max_working_gib=args.max_working_gib,
                   min_bin_count=args.min_bin_count)
    candidates = [c for c in matrices['candidates'] if c['bin_count'] >= args.min_bin_count]
    if not candidates:
        raise ValueError('No matrix candidates meet the minimum bin count')
    dispatch = work/'dispatch.json'
    tasks = distribute(candidates, options['metrics'], options['ordinations'], args.nodes)
    io.write_json(dispatch, dict(nodes=args.nodes, matrices_dir=str(matrix_root), evaluation=options, tasks=tasks))
    evaluate_started = perf_counter()
    launch(work, dispatch, 'evaluate', args.nodes, args.jobs)
    evaluated = perf_counter()
    result = merge_evaluation(matrix_root, work, output, tasks, options)
    merged_evaluation = perf_counter()
    shutil.rmtree(work/'evaluate')
    return result, dict(evaluation_workers=evaluated-evaluate_started,
                        evaluation_merge=merged_evaluation-evaluated)


def run_stage(args, *, launch=launch_slurm):
    output = Path(args.outdir).absolute()
    if output.exists():
        raise FileExistsError(f'Output directory already exists: {output}')
    work = io.new_directory(output.with_name(output.name + '.multinode'))
    if args.command == 'build':
        timings = build_phase(work, args.plan, output, args, launch)
    else:
        _, timings = evaluate_phase(work, args.matrices, output, args, launch)
    io.write_json(work/'timings.json', timings)


def run_multinode(output, args, *, launch=launch_slurm):
    """Build once per frozen row, then distribute each distance/candidate task."""
    work = io.new_directory(output/'_multinode')
    timings = build_phase(work, output/'plan', output/'matrices', args, launch)
    result, evaluation_timings = evaluate_phase(work, output/'matrices', output/'evaluation', args, launch)
    io.write_json(work/'timings.json', {**timings, **evaluation_timings})
    return io.seal(output, dict(kind='run', schema=1, status=result['status'],
        components=['plan', 'matrices', 'evaluation'], nodes=args.nodes))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('work', type=Path)
    parser.add_argument('dispatch', type=Path)
    parser.add_argument('phase', choices=('build', 'evaluate'))
    args = parser.parse_args()
    worker(args.work, args.dispatch, args.phase, int(os.environ['SLURM_PROCID']))
