"""Exercise shard seams and verify full results against the single-node path."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import json
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

from bin import bin_count_eval_io as io, bin_count_multinode as multinode
from bin import bin_count_eval_math as maths
from bin.bin_count_evaluate import evaluate, main
from bin.bin_count_multinode import run_multinode, run_stage, worker
from tests.test_bin_count_evaluate import fixture


def compare_runs(single, multiple):
    _, left = io.load(single/'matrices', 'matrices')
    _, right = io.load(multiple/'matrices', 'matrices')
    if left['candidates'] != right['candidates'] or left['artifact_sha256'] != right['artifact_sha256']:
        raise ValueError('Distributed matrices or row/edge metadata differ')
    _, left = io.load(single/'evaluation', 'evaluation')
    _, right = io.load(multiple/'evaluation', 'evaluation')
    def numeric_artifacts(manifest):
        return {name: sha for name, sha in manifest['artifact_sha256'].items()
                if Path(name).suffix not in ('.pdf', '.png')}

    if numeric_artifacts(left) != numeric_artifacts(right):
        raise ValueError('Distributed evaluation numeric artifacts differ')
    if left['embeddings'] != right['embeddings'] or left['status'] != right['status']:
        raise ValueError('Distributed ordination metadata differs')
    return len(numeric_artifacts(left))


class MultinodeTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        self.pre, self.db = fixture(self.root)
        self.args = argparse.Namespace(nodes=2, jobs=2, max_matrix_gib=1, max_working_gib=1,
                                      distance=list(maths.METRICS), ordination=['pcoa'],
                                      seeds=[11], perplexity=3, iterations=300, learning_rate='auto',
                                      min_bin_count=1)

    def plan(self, name):
        output = self.root / name
        io.create_plan(self.pre, self.db, output/'plan', [10, 30], [5, 7])
        return output

    def launch_local(self, work, dispatch, phase, nodes, jobs):
        with ThreadPoolExecutor(nodes) as executor:
            # Matplotlib runs in separate processes on Slurm; serialize it here.
            if phase == 'build':
                results = [executor.submit(worker, work, dispatch, phase, rank) for rank in range(nodes)]
                for result in results:
                    result.result()
            else:
                for rank in range(nodes):
                    worker(work, dispatch, phase, rank)

    def test_single_vs_two_shards_exact(self):
        single, multiple = self.plan('single'), self.plan('multiple')
        io.build(single/'plan', single/'matrices', jobs=self.args.jobs)
        evaluate(single/'matrices', single/'evaluation', self.args.distance, self.args.ordination,
                 self.args.seeds, 3, 300, 'auto', self.args.jobs, 1)
        run_multinode(multiple, self.args, launch=self.launch_local)
        io.load(multiple, 'run')
        self.assertGreater(compare_runs(single, multiple), 16)
        config = json.loads((multiple/'_multinode/dispatch.json').read_text())
        self.assertTrue(all(config['tasks']))

    def test_missing_rank_does_not_seal(self):
        output = self.plan('missing')

        def missing(work, dispatch, phase, nodes, jobs):
            self.launch_local(work, dispatch, phase, nodes, jobs)
            shutil.rmtree(work/'build/rank-0001')

        with self.assertRaises(FileNotFoundError):
            run_multinode(output, self.args, launch=missing)
        self.assertFalse((output/'SUCCESS').exists())
        self.assertFalse((output/'matrices/SUCCESS').exists())

    def test_independent_build_cli_merges_exactly(self):
        output = self.plan('array')
        script = Path(__file__).resolve().parents[1] / 'bin/bin_count_evaluate.py'
        shared = ['--plan', str(output/'plan'), '--shard-dir', str(output/'shards'),
                  '--shards', '3', '--jobs', '2']

        def cli(command, *options, check=True):
            return subprocess.run([sys.executable, '-B', str(script), command, *shared, *options],
                                  capture_output=True, text=True, check=check)

        combined = output/'matrices'
        failed = cli('merge-build', '--outdir', str(combined), check=False)
        self.assertNotEqual(failed.returncode, 0)
        self.assertFalse(combined.exists())
        with ThreadPoolExecutor(3) as executor:
            list(executor.map(lambda rank: cli('build-shard', '--rank', str(rank)), range(3)))
        failed = cli('build-shard', '--rank', '0', check=False)
        self.assertNotEqual(failed.returncode, 0)
        failed = cli('merge-build', '--outdir', str(combined), '--jobs', '1', check=False)
        self.assertIn('Invalid or duplicate matrix shard', failed.stderr)
        self.assertFalse(combined.exists())
        cli('merge-build', '--outdir', str(combined))
        baseline = output/'baseline'
        io.build(output/'plan', baseline, jobs=2)
        self.assertEqual(io.load(baseline, 'matrices')[1], io.load(combined, 'matrices')[1])

    def test_larger_benchmark_fixture_has_sorted_sample_inventory(self):
        root = self.root/'larger'
        root.mkdir()
        preflight, _ = fixture(root, samples_per_gcf=12)
        _, summary, rows, _ = io.frozen_preflight(preflight)
        self.assertEqual(summary['counts']['retained_unique_samples'], 12)
        self.assertEqual(len(rows), 48)

    def test_streamed_merge_verifies_bytes_and_retries_with_previous_shards(self):
        output = self.plan('streamed')
        plan = output/'plan'
        shards = output/'shards'
        combined = output/'matrices'
        baseline = output/'baseline'
        io.build(plan, baseline)
        for rank in range(3):
            shard = multinode.rank_folder(shards, 'build', rank)
            io.build(plan, shard, row_shard=(rank, 3))
            _, manifest = io.read_manifest(shard, 'matrix-shard')
            manifest['code_sha256'].update(multinode.LEGACY_MERGE_CODE)
            io.write_json(shard/'result.json', manifest)
            (shard/'SUCCESS').write_text(io.sha(shard/'result.json') + '\n')

        first = multinode.rank_folder(shards, 'build', 0)
        raw = first / manifest['candidates'][0]['folder'] / 'raw_mean_depth.npy'
        original = raw.read_bytes()
        raw.write_bytes(original[:-1] + bytes([original[-1] ^ 1]))
        with self.assertRaisesRegex(ValueError, 'Changed artifact'):
            multinode.merge_matrices(plan, shards, combined, 3, 1)
        self.assertFalse(combined.exists())
        raw.write_bytes(original)

        original_sha = io.sha

        def reject_full_array_rescan(path):
            self.assertNotEqual(Path(path).suffix, '.npy')
            self.assertFalse(Path(path).is_relative_to(combined) and Path(path).name == 'edges.json')
            return original_sha(path)

        # Force several output blocks and exercise unaligned starts across the three ranks.
        with patch.object(multinode, 'MERGE_BLOCK_BYTES', 128), patch.object(io, 'sha', reject_full_array_rescan):
            multinode.merge_matrices(plan, shards, combined, 3, 1)
        self.assertEqual(io.load(baseline, 'matrices')[1], io.load(combined, 'matrices')[1])

    def test_normalized_only_cli_budget_shards_merge_and_tsne_parity(self):
        config = self.root/'strata.json'
        io.write_json(config, {'strata': [dict(
            min_genome_bp=0, max_genome_bp_exclusive=None, nominal_bin_widths_bp=[2, 5, 23])]})
        plan = self.root/'plan'
        common = ['--preflight-dir', str(self.pre), '--index', str(self.db),
                  '--strata-config', str(config)]
        main(['plan', *common, '--outdir', str(plan)])
        _, planned = io.read_manifest(plan, 'plan')
        full_bytes = planned['planned_matrix_bytes']
        sharded_limit = str(full_bytes / 1024**3)
        check = ['check-grid', *common, '--distance', 'euclidean', '--ordination', 'tsne',
                 '--jobs', '1', '--max-matrix-gib', sharded_limit, '--max-working-gib', '1']
        with self.assertRaises(SystemExit):
            main(check)
        main([*check, '--normalized-only'])

        baseline, single, merged = (self.root/name for name in ('full', 'single', 'merged'))
        io.build(plan, baseline)
        main(['build', '--plan', str(plan), '--outdir', str(single), '--normalized-only',
              '--max-matrix-gib', str((full_bytes // 2) / 1024**3)])
        work = self.root/'shards'
        shared = ['--plan', str(plan), '--shard-dir', str(work), '--shards', '2',
                  '--max-matrix-gib', sharded_limit]
        for rank in range(2):
            main(['build-shard', *shared, '--rank', str(rank), '--normalized-only'])
            shard = multinode.rank_folder(work, 'build', rank)
            _, manifest = io.load(shard, 'matrix-shard')
            self.assertIs(manifest['normalized_only'], True)
            self.assertFalse(list(shard.rglob('raw_mean_depth.npy')))
        with self.assertRaisesRegex(ValueError, 'storage mode'):
            multinode.merge_matrices(plan, work, merged, 2, 1)
        self.assertFalse(merged.exists())
        main(['merge-build', *shared, '--outdir', str(merged), '--normalized-only'])

        _, full_manifest = io.load(baseline, 'matrices')
        retained = {name: digest for name, digest in full_manifest['artifact_sha256'].items()
                    if not name.endswith('/raw_mean_depth.npy')}
        for matrices in (single, merged):
            _, manifest = io.load(matrices, 'matrices')
            self.assertIs(manifest['normalized_only'], True)
            self.assertEqual(manifest['candidates'], full_manifest['candidates'])
            self.assertEqual(manifest['artifact_sha256'], retained)
            self.assertFalse(list(matrices.rglob('raw_mean_depth.npy')))

        options = ['--ordination', 'tsne', '--plots-only', '--distance',
                   'euclidean', 'cosine', 'correlation', '--perplexity', '3',
                   '--iterations', '300', '--seeds', '19']
        outputs = []
        for matrices in (baseline, merged):
            plots = self.root/f'plots-{matrices.name}'
            main(['evaluate', '--matrices', str(matrices), '--outdir', str(plots), *options])
            outputs.append(plots)
        old, new = outputs
        coords = sorted(path.relative_to(old) for path in old.rglob('coordinates.tsv'))
        self.assertEqual(len(coords), 6)
        for path in coords:
            self.assertEqual((old/path).read_bytes(), (new/path).read_bytes())
            metadata = path.with_name('ordination.json')
            self.assertEqual((old/metadata).read_bytes(), (new/metadata).read_bytes())

    def test_independent_distance_jobs_reuse_one_sealed_matrix(self):
        io.create_plan(self.pre, self.db, self.root/'plan', [10, 30], [5, 20])
        matrices = self.root/'matrices'
        args = argparse.Namespace(**vars(self.args))
        args.nodes = 3
        args.command, args.plan, args.outdir = 'build', self.root/'plan', matrices
        run_stage(args, launch=self.launch_local)
        _, manifest = io.load(matrices, 'matrices')
        self.assertEqual([c['bin_count'] for c in manifest['candidates']], [1, 4, 1, 6])
        digest = io.sha(matrices/'result.json')
        for metric in ('euclidean', 'cosine', 'correlation'):
            args.command, args.matrices = 'evaluate', matrices
            args.outdir = self.root/f'distributed-{metric}'
            args.distance, args.ordination = [metric], ['none']
            args.min_bin_count = 2 if metric == 'correlation' else 1
            run_stage(args, launch=self.launch_local)
            single = self.root/f'single-{metric}'
            evaluate(matrices, single, [metric], ['none'], [11], 3, 300, 'auto', 2, 1,
                     min_bin_count=args.min_bin_count)
            self.assertEqual(io.sha(single/'scores.tsv'), io.sha(args.outdir/'scores.tsv'))
            scores = io.read_tsv(args.outdir/'scores.tsv')
            self.assertEqual(len(scores), 2 if metric == 'correlation' else 4)
            self.assertEqual(io.sha(matrices/'result.json'), digest)
        io.load(matrices, 'matrices')
