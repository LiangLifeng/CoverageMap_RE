"""Lossless storage and downstream sampling invariants; stdlib unittest runner."""
import json
from pathlib import Path
import sys
import tempfile
import unittest

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from bin.coverage_core import (BlockPool, Reference, dense_coverage, select_units,
                               window_edges, window_sums)
from bin.reference_store import (bounded_map, build_bundle, load_reference_bundle,
                                 read_storage, save_storage, sha256)


def fixture():
    # Fragment-like unit 0 has two distant ends; unit 1 contains a CIGAR gap.
    pool = BlockPool(np.array([0, 21, 3, 10, 2], dtype=np.uint32),
                     np.array([4, 26, 5, 15, 29], dtype=np.uint32),
                     np.array([0, 2, 4, 5], dtype=np.uint64), 31)
    return pool, Reference(('a', 'b'), (20, 11), ('abc', ''))


class ReferenceStoreTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary.name)

    def tearDown(self):
        self.temporary.cleanup()

    def build(self, pool=None, layout=None, name='bundle', **kwargs):
        if pool is None:
            pool, layout = fixture()
        path = self.root / name
        report = build_bundle(pool, layout, path, 'sample-A', {'sampling_unit': 'fragment'},
                              [{'path': '/immutable/input.cram', 'sha256': 'source-identity'}],
                              chunk_units=1, **kwargs)
        return path, report

    def test_all_formats_dense_and_preserved_metadata(self):
        original, layout = fixture()
        for codec in ('packed', 'npz', 'npy', 'auto'):
            with self.subTest(codec=codec):
                path, report = self.build(name=codec, storage_format=codec, workers=3)
                pool, actual_layout, manifest = load_reference_bundle(path, workers=2)
                self.assertEqual(actual_layout, layout)
                self.assertEqual(manifest['reference_m5_missing'], ['b'])
                self.assertEqual(manifest['sample_id'], 'sample-A')
                self.assertEqual(manifest['policies'], {'sampling_unit': 'fragment'})
                self.assertEqual(manifest['sources'][0]['sha256'], 'source-identity')
                self.assertEqual(report['selected_codec'], min(report['candidates'], key=lambda c: c['storage_bytes'])['codec'])
                self.assertEqual(len(report['candidates']), 5 if codec == 'auto' else 1)
                if codec == 'packed':
                    self.assertEqual(report['selected_level'], 9)

                for seed in range(5):
                    ids = select_units(pool, np.random.default_rng(seed), target_bases=17)
                    expected_ids = select_units(original, np.random.default_rng(seed), target_bases=17)
                    np.testing.assert_array_equal(ids, expected_ids)
                    dense = dense_coverage(original, ids)
                    for bins in (1, 7, 31):
                        edges = window_edges(31, bins)
                        expected = [dense[a:b].sum() for a, b in zip(edges[:-1], edges[1:])]
                        np.testing.assert_array_equal(window_sums(pool, ids, edges), expected)

    def test_production_default_is_single_verified_packed9(self):
        path, report = self.build()
        self.assertEqual((report['requested_format'], report['selected_codec'], report['selected_level']),
                         ('packed', 'packed', 9))
        self.assertEqual(len(report['candidates']), 1)
        self.assertTrue(report['replay_pass'])
        self.assertEqual(json.loads((path/'storage.json').read_text())['level'], 9)
        original, _ = fixture()
        loaded, _, _ = load_reference_bundle(path)
        np.testing.assert_array_equal(dense_coverage(loaded, np.arange(len(loaded))),
                                      dense_coverage(original, np.arange(len(original))))

    def test_thread_count_deterministic_payload(self):
        first, _ = self.build(name='serial', storage_format='packed', workers=1)
        second, _ = self.build(name='parallel', storage_format='packed', workers=4)
        first_files = {p.name: sha256(p) for p in first.glob('chunk-*.zlib')}
        self.assertEqual(first_files, {p.name: sha256(p) for p in second.glob('chunk-*.zlib')})
        self.assertEqual((first / 'storage.json').read_bytes(), (second / 'storage.json').read_bytes())

    def test_empty_pool_all_formats(self):
        pool = BlockPool.from_intervals([], 31)
        layout = Reference(('a',), (31,), ('',))
        for codec in ('packed', 'npz', 'npy', 'auto'):
            path, _ = self.build(pool, layout, name=codec, storage_format=codec)
            loaded, _, _ = load_reference_bundle(path)
            self.assertEqual(len(loaded), 0)
            self.assertEqual(loaded.total_bases, 0)

    def test_uint64_coordinates_and_counts(self):
        genome = 2**32 + 100
        pool = BlockPool.from_intervals(np.array([[2**32, 2**32 + 21], [0, 3]], dtype=np.uint64), genome)
        layout = Reference(('large',), (genome,), ('m5',))
        path, _ = self.build(pool, layout, storage_format='packed')
        loaded, _, _ = load_reference_bundle(path)
        np.testing.assert_array_equal(loaded.starts, pool.starts)
        np.testing.assert_array_equal(loaded.ends, pool.ends)

    def test_many_blocks_and_length_widths(self):
        starts = np.arange(300, dtype=np.uint32) * 2
        pool = BlockPool(starts, starts + 1, np.array([0, 300], dtype=np.uint64), 100000)
        layout = Reference(('one',), (100000,), ('',))
        path, _ = self.build(pool, layout, storage_format='packed')
        loaded, _, _ = load_reference_bundle(path)
        np.testing.assert_array_equal(loaded.indptr, [0, 300])
        long = BlockPool.from_intervals(np.array([[1, 70001]]), 100000)
        path, _ = self.build(long, layout, name='long', storage_format='packed')
        np.testing.assert_array_equal(load_reference_bundle(path)[0].ends, [70001])

    def test_corruption_rejected(self):
        path, _ = self.build(storage_format='packed')
        chunk = next(path.glob('chunk-*.zlib'))
        chunk.write_bytes(chunk.read_bytes() + b'bad')
        with self.assertRaisesRegex(ValueError, 'Checksum'):
            load_reference_bundle(path)

    def test_incomplete_and_failed_rejected(self):
        path, _ = self.build()
        (path / 'status/SUCCESS').unlink()
        with self.assertRaisesRegex(ValueError, 'incomplete'):
            load_reference_bundle(path)
        (path / 'status/SUCCESS').touch()
        (path / 'status/FAILED').touch()
        with self.assertRaises(ValueError):
            load_reference_bundle(path)

    def test_collision_invalid_arguments_and_layout(self):
        path, _ = self.build()
        with self.assertRaises(FileExistsError):
            self.build()
        for kwargs in ({'workers': 0}, {'storage_format': 'unknown'}):
            with self.assertRaises(ValueError):
                self.build(name='bad' + str(kwargs), **kwargs)
        pool, _ = fixture()
        with self.assertRaises(ValueError):
            self.build(pool, Reference(('x',), (5,), ('',)), name='mismatch')

    def test_source_changed_retains_failure(self):
        pool, layout = fixture()
        def changed():
            raise ValueError('Source changed')
        path = self.root / 'changed'
        with self.assertRaisesRegex(ValueError, 'Source changed'):
            build_bundle(pool, layout, path, 'x', {}, [], source_check=changed)
        self.assertTrue((path / 'status/FAILED').exists())
        self.assertFalse((path / 'status/SUCCESS').exists())

    def test_escape_and_unchecksummed_payload_rejected(self):
        path, _ = self.build(storage_format='npy')
        checks = json.loads((path / 'checksums.json').read_text())
        checks['../escape'] = 'bad'
        (path / 'checksums.json').write_text(json.dumps(checks))
        with self.assertRaisesRegex(ValueError, 'path'):
            load_reference_bundle(path)
        del checks['../escape']
        del checks['starts.npy']
        (path / 'checksums.json').write_text(json.dumps(checks))
        with self.assertRaisesRegex(ValueError, 'Unchecksummed'):
            load_reference_bundle(path)

    def test_truncated_chunk_rejected(self):
        pool, _ = fixture()
        directory = self.root / 'raw'
        save_storage(pool, directory, 'packed', 2, 1)
        chunk = next(directory.glob('chunk-*.zlib'))
        chunk.write_bytes(chunk.read_bytes()[:-1])
        with self.assertRaises(ValueError):
            read_storage(directory)

    def test_bounded_map_preserves_order_and_propagates_errors(self):
        self.assertEqual(list(bounded_map(lambda i: i * i, range(100), 4)), [i*i for i in range(100)])
        with self.assertRaises(ZeroDivisionError):
            list(bounded_map(lambda i: 1 / i, range(4), 2))


if __name__ == '__main__':
    unittest.main()
