"""Shared-filesystem visibility failures must not publish incomplete manifests."""
from contextlib import redirect_stderr
import errno
import io
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from bin import reference_store as store
from bin import ref_depth_worker as worker
from bin.ref_depth_catalog import completed_sample, digest
from bin.ref_depth_worker import run_gcf
import test_ref_depth_database as fixtures
from test_reference_store import fixture


class PayloadRetryTests(unittest.TestCase):
    def test_transient_access_retries_and_logs_exact_path(self):
        for code in store.TRANSIENT_PAYLOAD_ERRNOS:
            with self.subTest(errno=code):
                operation = unittest.mock.Mock(side_effect=[OSError(code, 'transient'),
                                                            OSError(code, 'transient'), 'ok'])
                capture = io.StringIO()
                with patch.object(store.time, 'sleep') as sleep, redirect_stderr(capture):
                    self.assertEqual(store.retry_payload_io(operation, '/payload/chunk.zlib'), 'ok')
                self.assertEqual(sleep.call_args_list, [unittest.mock.call(.2), unittest.mock.call(.5)])
                self.assertEqual(capture.getvalue().count('PAYLOAD RETRY'), 2)
                self.assertIn('/payload/chunk.zlib', capture.getvalue())

    def test_permanent_failure_bounded_and_permissions_or_validation_not_retried(self):
        with patch.object(store.time, 'sleep') as sleep, redirect_stderr(io.StringIO()):
            operation = unittest.mock.Mock(side_effect=FileNotFoundError(errno.ENOENT, 'missing'))
            with self.assertRaises(FileNotFoundError):
                store.retry_payload_io(operation, '/missing')
            self.assertEqual(operation.call_count, 5)
            self.assertEqual(sleep.call_count, 4)
        for error in [PermissionError(errno.EACCES, 'denied'), ValueError('Checksum mismatch')]:
            with patch.object(store.time, 'sleep') as sleep:
                operation = unittest.mock.Mock(side_effect=error)
                with self.assertRaises(type(error)):
                    store.retry_payload_io(operation, '/bad')
                self.assertEqual(operation.call_count, 1)
                sleep.assert_not_called()

    def test_native_both_manifests_complete_despite_omitted_directory_entries(self):
        f = fixtures.DepthDatabaseTests(methodName='runTest')
        f.setUp()
        self.addCleanup(f.tearDown)
        database, plan = f.plan(storage_mode='native', window_size=7, depths=[], geometry_format='packed')
        original_iterdir, original_rglob = Path.iterdir, Path.rglob

        def hidden_iterdir(path):
            return (p for p in original_iterdir(path) if p.suffix != '.zlib')

        def hidden_rglob(path, pattern):
            return (p for p in original_rglob(path, pattern) if p.suffix != '.zlib')

        with patch.object(Path, 'iterdir', hidden_iterdir), patch.object(Path, 'rglob', hidden_rglob):
            self.assertEqual(run_gcf(database, fixtures.GCF1, f.extractor), 0)
        for record in plan['groups'][fixtures.GCF1]['samples']:
            attempt, result = completed_sample(database / fixtures.GCF1 / 'samples' / record['folder'],
                                                digest(database / 'database.json'), verify=True)
            storage = json.loads((attempt / 'geometry/storage.json').read_text())
            inner = json.loads((attempt / 'geometry/checksums.json').read_text())
            outer = json.loads((attempt / 'checksums.json').read_text())
            self.assertTrue(storage['chunks'])
            for chunk in storage['chunks']:
                name = chunk['file']
                self.assertEqual(inner[name], outer['geometry/' + name])
            loaded, _, _ = store.load_reference_bundle(attempt / 'geometry')
            self.assertEqual(loaded.total_bases, result['aligned_bases'])

    def test_bundle_recovers_transient_chunk_read_then_roundtrips(self):
        with tempfile.TemporaryDirectory() as temporary:
            out = Path(temporary) / 'bundle'
            pool, layout = fixture()
            original_open = Path.open
            failures = []

            def transient_open(path, *args, **kwargs):
                mode = args[0] if args else kwargs.get('mode', 'r')
                if path.name == 'chunk-000000000000.zlib' and mode == 'rb' and len(failures) < 2:
                    failures.append(str(path))
                    raise OSError(errno.ESTALE, 'stale handle', str(path))
                return original_open(path, *args, **kwargs)

            with patch.object(Path, 'open', transient_open), patch.object(store.time, 'sleep') as sleep:
                store.build_bundle(pool, layout, out, 'sample', {}, [], workers=1, chunk_units=1)
            self.assertEqual(len(failures), 2)
            self.assertEqual(sleep.call_count, 2)
            loaded, _, _ = store.load_reference_bundle(out)
            self.assertEqual(loaded.total_bases, pool.total_bases)

    def test_missing_written_chunk_exhausts_retry_and_never_publishes_success(self):
        with tempfile.TemporaryDirectory() as temporary:
            out = Path(temporary) / 'bundle'
            pool, layout = fixture()
            original_encode = store.encode_chunk

            def missing_chunk(job):
                result = original_encode(job)
                (job[0] / result['file']).unlink()
                return result

            with (patch.object(store, 'encode_chunk', missing_chunk),
                  patch.object(store.time, 'sleep') as sleep, redirect_stderr(io.StringIO())):
                with self.assertRaises(FileNotFoundError):
                    store.build_bundle(pool, layout, out, 'sample', {}, [], workers=1, chunk_units=1)
            self.assertEqual(sleep.call_count, 4)
            self.assertTrue((out / 'status/FAILED').exists())
            self.assertFalse((out / 'status/SUCCESS').exists())
            self.assertFalse((out / 'checksums.json').exists())

    def test_outer_manifest_retries_transient_access_and_rejects_persistent_missing(self):
        for persistent in (False, True):
            with self.subTest(persistent=persistent):
                f = fixtures.DepthDatabaseTests(methodName='runTest')
                f.setUp()
                try:
                    database, plan = f.plan(storage_mode='native', window_size=7, depths=[], geometry_format='packed')
                    record = plan['groups'][fixtures.GCF1]['samples'][0]
                    original_required, original_open = worker.required_sample_payloads, Path.open
                    state = dict(outer=False, errors=0)

                    def required(*args):
                        names = original_required(*args)
                        state['outer'] = True
                        return names

                    def transient_open(path, *args, **kwargs):
                        mode = args[0] if args else kwargs.get('mode', 'r')
                        if (state['outer'] and path.suffix == '.zlib' and mode == 'rb'
                                and (persistent or state['errors'] < 2)):
                            state['errors'] += 1
                            raise FileNotFoundError(errno.ENOENT, 'temporarily invisible', str(path))
                        return original_open(path, *args, **kwargs)

                    with (patch.object(worker, 'required_sample_payloads', required),
                          patch.object(Path, 'open', transient_open), patch.object(store.time, 'sleep') as sleep,
                          redirect_stderr(io.StringIO())):
                        if persistent:
                            with self.assertRaises(FileNotFoundError):
                                worker.build_sample(database, plan, fixtures.GCF1, record, extractor=f.extractor)
                        else:
                            worker.build_sample(database, plan, fixtures.GCF1, record, extractor=f.extractor)
                    attempt = database / fixtures.GCF1 / 'samples' / record['folder'] / 'attempt-000001'
                    self.assertEqual(sleep.call_count, 4 if persistent else 2)
                    self.assertEqual((attempt / 'status/SUCCESS').exists(), not persistent)
                    self.assertEqual((attempt / 'status/FAILED').exists(), persistent)
                    if not persistent:
                        completed_sample(attempt.parent, digest(database / 'database.json'), verify=True)
                finally:
                    f.tearDown()


if __name__ == '__main__':
    unittest.main()
