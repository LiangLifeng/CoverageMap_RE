"""Parallel metadata preparation preserves immutable deterministic SQLite indexes."""
from contextlib import closing
from collections import deque
from concurrent.futures import ProcessPoolExecutor
import json
import multiprocessing
from pathlib import Path
import shutil
import sqlite3
import subprocess
import sys
import unittest
from unittest.mock import patch

import test_ref_depth_database as fixtures
from bin import reference_index as indexer
from bin.ref_depth_catalog import collect, create_plan, digest
from bin.ref_depth_worker import run_gcf


class ParallelIndexTests(unittest.TestCase):
    def setUp(self):
        self.f = fixtures.DepthDatabaseTests('runTest')
        self.f.setUp()
        self.root = self.f.root

    def tearDown(self):
        self.f.tearDown()

    def database(self, native=True):
        # More than two batches, two GCFs; exercise ordered cross-GCF boundaries.
        rows = [dict(self.f.rows[n % 2], Sample='sample-{:03d}'.format(n),
                     GCFID=fixtures.GCF1 if n < 34 else fixtures.GCF2) for n in range(36)]
        fixtures.write_csv(self.f.records, ['Sample', 'CramFilePath', 'GCFID'], rows)
        fixtures.write_csv(self.f.refs, ['GCF_ID', 'GCF_path'],
                           [dict(GCF_ID=g, GCF_path=str(self.f.fai)) for g in (fixtures.GCF1, fixtures.GCF2)])
        config = dict(self.f.config)
        if native:
            config.update(storage_mode='native', depths=[], bins=None, window_size=7)
        database = self.root / ('native' if native else 'legacy')
        plan = create_plan(self.f.records, self.f.refs, database, config, fixtures.COLUMNS)
        for gcf in plan['groups']:
            self.assertEqual(run_gcf(database, gcf, self.f.extractor), 0)
        self.assertTrue(collect(database)['complete'])
        return database, plan

    def assert_same_index(self, left, right):
        # Metadata, source IDs, row IDs, ordering, selection, integrity, and file bytes.
        with closing(sqlite3.connect(left.as_uri() + '?mode=ro', uri=True)) as a, a, closing(sqlite3.connect(right.as_uri() + '?mode=ro', uri=True)) as b, b:
            self.assertEqual(list(a.iterdump()), list(b.iterdump()))
            for connection in (a, b):
                self.assertEqual(connection.execute('PRAGMA integrity_check').fetchone()[0], 'ok')
                self.assertEqual(connection.execute('PRAGMA foreign_key_check').fetchall(), [])
        self.assertEqual(digest(left), digest(right))
        self.assertEqual(left.stat().st_mode & 0o777, 0o444)
        self.assertEqual(right.stat().st_mode & 0o777, 0o444)

    def test_native_legacy_serial_process_exact_and_sources_unchanged(self):
        pids = set()
        def executor_factory(*args, **kwargs):
            executor = ProcessPoolExecutor(*args, **kwargs)
            original = executor.submit
            def submit(*args, **kwargs):
                result = original(*args, **kwargs)
                pids.update(executor._processes)
                return result
            executor.submit = submit
            return executor
        for native in (True, False):
            database, plan = self.database(native)
            before = {p: digest(p) for p in database.rglob('*') if p.is_file()}
            for verify in (False, True):
                serial = self.root / '{}-{}-serial.sqlite'.format(native, verify)
                parallel = self.root / '{}-{}-parallel.sqlite'.format(native, verify)
                indexer.build_index(database, serial, verify, min_depth=2, min_samples=2)
                with patch.object(indexer, 'ProcessPoolExecutor', side_effect=executor_factory):
                    indexer.build_index(database, parallel, verify, min_depth=2, min_samples=2, jobs=2)
                self.assert_same_index(serial, parallel)
            self.assertEqual(before, {p: digest(p) for p in database.rglob('*') if p.is_file()})
        self.assertGreaterEqual(len(pids), 2)
        self.assertNotIn(multiprocessing.current_process().pid, pids)

    def test_corruption_strict_skip_and_multiroot_duplicate_equivalence(self):
        database, plan = self.database()
        copy = self.root / 'copy'
        shutil.copytree(database, copy)
        attempt = next(database.glob('*/samples/*/attempt-*'))
        (attempt / 'native_coverage.npz').write_bytes(b'corrupt payload')
        for jobs in (1, 2):
            output = self.root / 'strict-{}.sqlite'.format(jobs)
            with self.assertRaises(ValueError):
                indexer.build_index(database, output, True, jobs=jobs)
            self.assertFalse(output.exists())
            self.assertFalse(list(self.root.glob('.reference-registry-*')))
        serial, parallel = self.root / 'skip1.sqlite', self.root / 'skip2.sqlite'
        for output, jobs in ((serial, 1), (parallel, 2)):
            info = indexer.build_index([database, copy], output, True, skip_invalid_samples=True, jobs=jobs)
            self.assertEqual((info['skipped_samples'], info['counts']['references'], info['duplicate_copies']), (1, 36, 35))
            self.assertTrue(info['partial_index'])
        self.assert_same_index(serial, parallel)
        # All invalid stays inspectable, and empty valid-result batches are handled.
        for success in database.glob('*/samples/*/attempt-*/status/SUCCESS'):
            success.unlink()
        left, right = self.root / 'empty1.sqlite', self.root / 'empty2.sqlite'
        for out, jobs in ((left, 1), (right, 2)):
            info = indexer.build_index(database, out, True, skip_invalid_samples=True, jobs=jobs)
            self.assertEqual((info['counts']['references'], info['skipped_samples']), (0, 36))
        self.assert_same_index(left, right)

    def test_source_changes_remain_fatal_and_output_exclusive(self):
        database, plan = self.database()
        summary = database / 'SUMMARY.json'
        original = summary.read_bytes()
        prepare = indexer._prepared_samples
        def mutate(*args, **kwargs):
            yield from prepare(*args, **kwargs)
            summary.write_bytes(original + b'\n')
        for jobs in (1, 2):
            output = self.root / 'changed{}.sqlite'.format(jobs)
            with patch.object(indexer, '_prepared_samples', side_effect=mutate):
                with self.assertRaisesRegex(ValueError, 'changed while indexing'):
                    indexer.build_index(database, output, jobs=jobs)
            self.assertFalse(output.exists())
            summary.write_bytes(original)
        output = self.root / 'exists.sqlite'
        output.write_bytes(b'user content')
        with self.assertRaises(FileExistsError):
            indexer.build_index(database, output, jobs=2)
        self.assertEqual(output.read_bytes(), b'user content')

    def test_cli_jobs_and_external_audit(self):
        database, plan = self.database()
        next(database.glob('*/samples/*/attempt-*/status/SUCCESS')).unlink()
        entry = Path(indexer.__file__).resolve()
        output, audit = self.root / 'cli.sqlite', self.root / 'errors.tsv'
        result = subprocess.run([sys.executable, '-B', str(entry), 'build', '--database', str(database),
            '--out', str(output), '--jobs', '2', '--skip-invalid-samples', '--verify-payloads',
            '--skipped-samples-tsv', str(audit)], capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(json.loads(result.stdout)['skipped_samples'], 1)
        self.assertEqual(len(audit.read_text().splitlines()), 2)
        with closing(sqlite3.connect(output)) as connection, connection:
            self.assertIsNone(connection.execute("SELECT 1 FROM sqlite_master WHERE name='skipped_samples'").fetchone())
        for jobs in (0, -1, True, 1.5):
            with self.assertRaisesRegex(ValueError, 'positive integer'):
                indexer.build_index(database, self.root / 'invalid' / 'never.sqlite', jobs=jobs)
        self.assertFalse((self.root / 'invalid').exists())
        result = subprocess.run([sys.executable, '-B', str(entry), 'build', '--database', str(database),
            '--out', str(self.root / 'zero.sqlite'), '--jobs', '0'], capture_output=True, text=True)
        self.assertNotEqual(result.returncode, 0)
        self.assertFalse((self.root / 'zero.sqlite').exists())


class BoundedPreparationTests(unittest.TestCase):
    def test_order_bounded_prefetch_and_cleanup_on_writer_close(self):
        class Future:
            def __init__(self, owner, batch):
                self.owner, self.batch = owner, batch
                self.cancelled = False
            def result(self):
                self.owner.outstanding -= 1
                # Results represent work finishing out of order, but remain queued.
                return [(g, s, None, None) for g, layout, s in self.batch]
            def cancel(self):
                self.cancelled = True
        class Executor:
            def __init__(self, **kwargs):
                self.outstanding = self.peak = 0
                self.futures = []
                self.closed = False
            def submit(self, function, batch):
                self.outstanding += 1
                self.peak = max(self.peak, self.outstanding)
                future = Future(self, batch)
                self.futures.append(future)
                return future
            def shutdown(self, **kwargs):
                self.closed = kwargs == dict(wait=True, cancel_futures=True)
        plan = dict(config={}, groups={'g':dict(layout={}, samples=[dict(sample=str(n)) for n in range(1000)])})
        executor = Executor()
        with patch.object(indexer, 'ProcessPoolExecutor', return_value=executor):
            result = list(indexer._prepared_samples(Path('.'), plan, 'sha', False, 2))
        self.assertEqual([r[1]['sample'] for r in result], [str(n) for n in range(1000)])
        self.assertLessEqual(executor.peak, 4)
        self.assertTrue(executor.closed)
        executor = Executor()
        with patch.object(indexer, 'ProcessPoolExecutor', return_value=executor):
            iterator = indexer._prepared_samples(Path('.'), plan, 'sha', False, 2)
            next(iterator)
            iterator.close()
        self.assertTrue(executor.closed)
        self.assertTrue(all(f.cancelled for f in executor.futures[1:]))
        self.assertEqual(len(executor.futures), 4)
