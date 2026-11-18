"""Saved-reference integration without local real data or a pysam dependency."""
import json
from pathlib import Path
import sys
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from bin import coverage_concordance as cli
from bin.concordance_inputs import POLICY_DEFAULTS, load_inputs, resolve_bundle
from bin.coverage_core import BlockPool, Reference
from bin.reference_store import build_bundle


class ConcordanceBundleTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary.name)
        rng = np.random.default_rng(912)
        starts = rng.integers(0, 470, size=500)
        self.pool = BlockPool.from_intervals(np.column_stack([starts, starts + 15]), 500)
        self.layout = Reference(('a', 'b'), (300, 200), ('one', 'two'))
        self.alignment = self.root / 'test.cram'
        self.alignment.write_bytes(b'synthetic-only-test-alignment')
        self.index = self.root / 'ref.fai'
        self.index.write_text('a\t300\nb\t200\n')
        self.bundle = self.root / 'reference'
        build_bundle(self.pool, self.layout, self.bundle, 'source-sample', dict(POLICY_DEFAULTS),
                     [{'path': '/not-present/original-reference.cram', 'sha256': 'reference-content'}],
                     workers=2, chunk_units=100, storage_format='packed')

    def tearDown(self):
        self.temporary.cleanup()

    def args(self, **kwargs):
        values = dict(alignment=None, reference_bundle=self.bundle / 'storage.json',
            test_alignment=self.alignment, reference=None, reference_index=None,
            extraction_python=None, outdir=self.root / 'result', workers=3,
            seed=401, replicates=9, batch_size=4, alpha=.01, similarity_block_size=8,
            bins=23, depth=None, depths='.1,1,3', depth_grid=None, executor='thread',
            depth_sampler='permutation', coverage_backend='cpu', statistics_backend='cpu')
        values.update({key: None for key in POLICY_DEFAULTS})
        values.update(kwargs)
        return SimpleNamespace(**values)

    def test_index_loading_does_not_open_original_reference(self):
        tables = []
        def extract(alignment, reference, reference_index, **policies):
            self.assertEqual(alignment, self.alignment)
            self.assertIsNone(reference)
            self.assertEqual(policies, POLICY_DEFAULTS)
            tables.append(reference_index)
            self.assertEqual(reference_index.read_text(), 'a\t300\nb\t200\n')
            return self.pool, self.layout
        with patch('bin.concordance_inputs.extract_pool', side_effect=extract) as mocked:
            inputs = load_inputs(self.args())
        self.assertEqual(mocked.call_count, 1)
        self.assertFalse(tables[0].exists())
        self.assertTrue(inputs.derived_layout)
        self.assertEqual(inputs.bundle_info['sample_id'], 'source-sample')
        self.assertEqual(inputs.comparison_scope, 'separate-input-pools')
        inputs.verify()

    def test_cli_same_counts_ids_and_statistics_as_alignment(self):
        outputs = []
        for mode in ('alignment', 'bundle'):
            for sampler in ('permutation', 'lazy'):
                out = self.root / f'{mode}-{sampler}'
                args = self.args(outdir=out, depth_sampler=sampler)
                if mode == 'alignment':
                    args.alignment = self.alignment
                    args.reference_bundle = None
                    args.reference_index = self.index
                with patch.object(cli, 'parse_args', return_value=args), patch(
                        'bin.concordance_inputs.extract_pool', return_value=(self.pool, self.layout)):
                    cli.main()
                outputs.append(out)
                self.assertTrue((out / 'status/SUCCESS').exists())
        for sampler in ('permutation', 'lazy'):
            raw, saved = self.root / f'alignment-{sampler}', self.root / f'bundle-{sampler}'
            self.assertEqual(json.loads((raw / 'summary.json').read_text()),
                             json.loads((saved / 'summary.json').read_text()))
            for depth in range(3):
                left, right = raw / f'depth-{depth:05d}', saved / f'depth-{depth:05d}'
                for name in ('simulations.npy', 'test.npy', 'test_ids.npy'):
                    np.testing.assert_array_equal(np.load(left / name), np.load(right / name))
                a = json.loads((left / 'concordance.json').read_text())
                b = json.loads((right / 'concordance.json').read_text())
                self.assertEqual(a['reference_replicates'], b['reference_replicates'])
                self.assertEqual(a['test_selected_ids_sha256'], b['test_selected_ids_sha256'])

    def test_saved_nondefault_policies_are_inherited(self):
        policies = dict(POLICY_DEFAULTS, min_mapq=20, exclude_duplicates=True, mate_policy='read1')
        bundle = self.root / 'nondefault'
        build_bundle(self.pool, self.layout, bundle, 's', policies, [], storage_format='npz')
        with patch('bin.concordance_inputs.extract_pool', return_value=(self.pool, self.layout)) as mocked:
            inputs = load_inputs(self.args(reference_bundle=bundle))
        self.assertEqual(inputs.policies, policies)
        self.assertEqual(mocked.call_args.kwargs['min_mapq'], 20)
        self.assertTrue(mocked.call_args.kwargs['exclude_duplicates'])

    def test_policy_conflicts_rejected_before_extraction(self):
        with patch('bin.concordance_inputs.extract_pool') as mocked:
            for kwargs in ({'min_mapq': 20}, {'mate_policy': 'read1'}, {'exclude_qcfail': True}):
                with self.subTest(kwargs=kwargs), self.assertRaisesRegex(ValueError, 'Policy'):
                    load_inputs(self.args(**kwargs))
            mocked.assert_not_called()

    def test_layout_and_md5_conflicts_rejected(self):
        for layout in (Reference(('a', 'b'), (299, 201), ('one', 'two')),
                       Reference(('a', 'b'), (300, 200), ('wrong', 'two')),
                       Reference(('b', 'a'), (200, 300), ('two', 'one'))):
            with patch('bin.concordance_inputs.extract_pool', return_value=(self.pool, layout)):
                with self.assertRaisesRegex(ValueError, 'layouts differ'):
                    load_inputs(self.args())

    def test_corrupt_and_incomplete_bundle_fail(self):
        (self.bundle / 'status/SUCCESS').unlink()
        with self.assertRaises(ValueError):
            load_inputs(self.args())
        (self.bundle / 'status/SUCCESS').touch()
        chunk = next(self.bundle.glob('chunk-*.zlib'))
        chunk.write_bytes(b'corrupt')
        with self.assertRaisesRegex(ValueError, 'Checksum'):
            load_inputs(self.args())

    def test_input_changes_are_detected(self):
        with patch('bin.concordance_inputs.extract_pool', return_value=(self.pool, self.layout)):
            inputs = load_inputs(self.args())
        self.alignment.write_bytes(b'changed')
        with self.assertRaisesRegex(ValueError, 'Input changed'):
            inputs.verify()

    def test_bundle_status_change_detected(self):
        with patch('bin.concordance_inputs.extract_pool', return_value=(self.pool, self.layout)):
            inputs = load_inputs(self.args())
        (self.bundle / 'status/FAILED').touch()
        with self.assertRaisesRegex(ValueError, 'status changed'):
            inputs.verify()

    def test_explicit_index_and_legacy_defaults(self):
        args = self.args(alignment=self.alignment, reference_bundle=None, reference_index=self.index)
        with patch('bin.concordance_inputs.extract_pool', return_value=(self.pool, self.layout)) as mocked:
            inputs = load_inputs(args)
        self.assertEqual(mocked.call_count, 1)
        self.assertEqual(inputs.policies, POLICY_DEFAULTS)
        self.assertIs(inputs.reference_pool, inputs.test_pool)
        self.assertEqual(inputs.comparison_scope, 'same-source-self-consistency')
        with self.assertRaisesRegex(ValueError, 'requires'):
            load_inputs(self.args(alignment=self.alignment, reference_bundle=None))

    def test_cli_rejects_unsupported_depth_and_output_collision(self):
        with patch.object(cli, 'parse_args', return_value=self.args(depths='1000')), patch(
                'bin.concordance_inputs.extract_pool', return_value=(self.pool, self.layout)):
            with self.assertRaisesRegex(ValueError, 'exceeds'):
                cli.main()
        out = self.root / 'exists'
        out.mkdir()
        with patch.object(cli, 'parse_args', return_value=self.args(outdir=out)):
            with self.assertRaises(FileExistsError):
                cli.main()

    def test_index_must_belong_to_complete_bundle(self):
        with self.assertRaisesRegex(ValueError, 'storage.json'):
            resolve_bundle(self.bundle / 'manifest.json')
        isolated = self.root / 'isolated'
        isolated.mkdir()
        (isolated / 'storage.json').write_text((self.bundle / 'storage.json').read_text())
        with self.assertRaises(ValueError):
            load_inputs(self.args(reference_bundle=isolated / 'storage.json'))


if __name__ == '__main__':
    unittest.main()
