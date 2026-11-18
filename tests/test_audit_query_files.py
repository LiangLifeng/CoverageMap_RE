"""File audit diagnosis without real CRAMs, scientific calculations or writes to inputs."""
from contextlib import closing
import argparse
import csv
import errno
import gzip
import importlib.util
import json
from pathlib import Path
import sqlite3
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

SCRIPT = Path(__file__).resolve().parents[1]/'bin/audit_query_files.py'
spec = importlib.util.spec_from_file_location('file_audit', SCRIPT)
audit = importlib.util.module_from_spec(spec)
spec.loader.exec_module(audit)


class FileAuditTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name).resolve()
        self.args = argparse.Namespace(rechecks=0, recheck_delay=0, hash_payloads=False)
        self.library = self.root/'library'
        self.attempt = self.library/'gcf/sample/attempt-000001'
        (self.attempt/'geometry').mkdir(parents=True)
        (self.attempt/'status').mkdir()
        self.write(self.library/'database.json', {})
        self.write(self.library/'SUMMARY.json', {})
        record = dict(storage_mode='native', sources=[dict(path=str(self.root/'original.cram'))])
        self.write(self.attempt/'result.json', record)
        self.write(self.attempt/'geometry/storage.json', dict(codec='packed', chunks=[dict(file='part.zlib')]))
        self.write(self.attempt/'geometry/checksums.json', {})
        for name in ('edges.npy', 'native_coverage.npz', 'geometry/part.zlib',
                     'geometry/manifest.json', 'geometry/report.json'):
            (self.attempt/name).write_bytes(b'payload')
        names = [str(p.relative_to(self.attempt)) for p in self.attempt.rglob('*') if p.is_file()]
        self.write(self.attempt/'checksums.json', {n: audit.sha256(self.attempt/n) for n in names})
        (self.attempt/'status/SUCCESS').touch()
        (self.attempt/'status/exit_code.txt').write_text('0\n')
        self.ref = dict(id=17, sample='sample', gcf='gcf', selected=0,
                        manifest_path='gcf/sample/attempt-000001/result.json',
                        manifest_sha=audit.sha256(self.attempt/'result.json'),
                        checksums_sha=audit.sha256(self.attempt/'checksums.json'))
        self.index = self.root/'index.sqlite'
        self.build_index()
        self.cram = self.root/'query.cram'; self.cram.write_bytes(b'CRAM fixture')
        self.manifest = self.root/'input.tsv'
        self.manifest.write_text(f'Sample\tGCFID\tDepth\tCramPath\nQ\tgcf\t2\t{self.cram}\n')

    def tearDown(self):
        self.temp.cleanup()

    def write(self, path, data):
        path.write_text(json.dumps(data))

    def build_index(self, federated=False):
        with closing(sqlite3.connect(self.index)) as c, c:
            c.execute('CREATE TABLE metadata (key TEXT,value TEXT)')
            metadata = dict(storage_mode='native')
            metadata.update(sources=[dict(id=5, root=str(self.library))]) if federated else metadata.update(database_root=str(self.library))
            c.executemany('INSERT INTO metadata VALUES (?,?)', [(k,json.dumps(v)) for k,v in metadata.items()])
            columns = ','.join(k+' '+('INTEGER' if k in ('id', 'selected') else 'TEXT') for k in self.ref)
            c.execute('CREATE TABLE references_table ('+columns+(',source_id INTEGER' if federated else '')+')')
            values = list(self.ref.values())+([5] if federated else [])
            c.execute('INSERT INTO references_table VALUES ('+','.join('?' for _ in values)+')',values)

    def rows(self, name, out):
        with gzip.open(out/name, 'rt') as f:
            return list(csv.DictReader(f, delimiter='\t'))

    def run_cli(self, name, *options):
        out = self.root/name
        command = [sys.executable, '-B', str(SCRIPT), '--outdir', str(out), *map(str, options)]
        result = subprocess.run(command, capture_output=True, text=True)
        self.assertTrue((out/'summary.json').is_file(), result.stderr)
        return result, out, json.loads((out/'summary.json').read_text())

    def test_single_and_federated_full_closure_unselected_provenance_and_immutable(self):
        for federated in (False, True):
            if federated:
                self.index.unlink(); self.build_index(True)
            before = audit.sha256(self.index)
            result, out, summary = self.run_cli('federated' if federated else 'single', '--index', self.index,
                                                 '--manifest', self.manifest, '--hash-payloads')
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertTrue(summary['enumeration_complete'])
            rows = self.rows('files.tsv.gz', out)
            chunks = [r for r in rows if r['path'].endswith('part.zlib')]
            self.assertEqual(len(chunks), 1)
            self.assertEqual(chunks[0]['selected'], '0')
            self.assertTrue(any(r['role']=='reference_original_provenance' and r['status']=='MISSING' for r in rows))
            self.assertTrue(any(r['path'].endswith('SUCCESS') and r['status']=='OK' for r in rows))
            self.assertEqual(before, audit.sha256(self.index))

    def test_missing_chunk_first_bad_component_and_hash_mismatch(self):
        (self.attempt/'geometry/part.zlib').unlink()
        result, out, summary = self.run_cli('missing', '--index', self.index)
        self.assertEqual(result.returncode, 2, result.stderr)
        self.assertTrue(summary['enumeration_complete'])
        rows = self.rows('issues.tsv.gz', out)
        chunk = next(r for r in rows if r['path'].endswith('part.zlib'))
        self.assertEqual(chunk['errno'], str(errno.ENOENT))
        self.assertEqual(chunk['first_bad_component'], str(self.attempt/'geometry/part.zlib'))
        (self.attempt/'geometry/part.zlib').write_bytes(b'changed')
        result, out, _ = self.run_cli('corrupt', '--index', self.index, '--hash-payloads')
        self.assertEqual(result.returncode, 2)
        self.assertTrue(any(r['status']=='SHA_MISMATCH' for r in self.rows('issues.tsv.gz', out)))

    def test_missing_manifest_and_malformed_ledger_incomplete(self):
        (self.attempt/'result.json').unlink()
        result, out, summary = self.run_cli('missing_manifest', '--index', self.index)
        self.assertEqual(result.returncode, 1)
        self.assertFalse(summary['enumeration_complete'])
        self.assertTrue(any(r['path'].endswith('part.zlib') for r in self.rows('files.tsv.gz', out)))
        (self.attempt/'checksums.json').write_text('{broken')
        result, _, summary = self.run_cli('bad_ledger', '--index', self.index)
        self.assertEqual(result.returncode, 1)
        self.assertFalse(summary['enumeration_complete'])

    def test_permission_errno_is_not_missing(self):
        with patch.object(Path, 'stat', side_effect=PermissionError(errno.EACCES, 'denied')):
            row = audit.probe(audit.context('/denied/input.cram', 'test'))
        self.assertEqual(row['status'], 'PERMISSION_DENIED')
        self.assertEqual(row['errno'], errno.EACCES)

    def test_broken_symlink_loop_directory_empty_and_relative(self):
        link = self.root/'broken'; link.symlink_to(self.root/'absent')
        row = audit.probe(audit.context(link, 'test'))
        self.assertEqual(row['status'], 'MISSING')
        self.assertIn('broken', row['symlinks'])
        loop = self.root/'loop'; loop.symlink_to(loop)
        self.assertIn(audit.probe(audit.context(loop, 'test'))['status'], ('SYMLINK_LOOP', 'IO_ERROR'))
        self.assertEqual(audit.probe(audit.context(self.root, 'test'))['status'], 'NOT_FILE')
        empty = self.root/'empty'; empty.touch()
        self.assertEqual(audit.probe(audit.context(empty, 'test'))['status'], 'EMPTY')
        self.assertEqual(audit.probe(audit.context('relative', 'test'))['status'], 'INVALID_PATH')

    def test_initial_failure_is_preserved_after_recovery(self):
        args = argparse.Namespace(rechecks=1, recheck_delay=0)
        with patch.object(audit, 'probe', side_effect=[dict(status='MISSING'),dict(status='OK')]):
            rows = audit.observe(audit.context(self.cram, 'test'), args)
        self.assertEqual([r['status'] for r in rows], ['MISSING', 'OK'])
        self.assertEqual([r['observation'] for r in rows], [0,1])

    def test_unsafe_ledger_and_escaping_symlink(self):
        with self.assertRaises(ValueError):
            audit.stored_path(self.library, '../outside')
        escape = self.library/'escape'; escape.symlink_to(self.root)
        with self.assertRaises(ValueError):
            audit.stored_path(self.library, 'escape/query.cram')
        self.write(self.attempt/'checksums.json', {'../outside':'a'})
        result, out, _ = self.run_cli('unsafe', '--index', self.index)
        self.assertEqual(result.returncode, 2)
        self.assertTrue(any(r['status']=='UNSAFE_PATH' for r in self.rows('issues.tsv.gz', out)))

    def test_serial_thread_equivalence(self):
        results = []
        for jobs in (1, 3):
            result, out, _ = self.run_cli('jobs'+str(jobs), '--index', self.index, '--jobs', jobs)
            self.assertEqual(result.returncode, 0)
            results.append([{k:v for k,v in r.items() if k!='utc'} for r in self.rows('files.tsv.gz', out)])
        self.assertEqual(*results)

    def test_frozen_run_history_and_early_return(self):
        run = self.root/'run'; (run/'queries/q00000001').mkdir(parents=True)
        frozen = {k:dict(path=str(p), sha256=audit.sha256(p)) for k,p in
                  [('manifest',self.manifest),('database',self.index)]}
        self.write(run/'run.json', frozen)
        (run/'summary.tsv').write_text('query_id\nq00000001\n')
        query = run/'queries/q00000001'
        self.write(query/'audit.json', dict(query=dict(CramPath=str(self.cram))))
        self.write(query/'checksums.json', {'draws.tsv.gz':'not-used-for-file-access'})
        with gzip.open(query/'draws.tsv.gz','wt') as f:
            f.write('error\treference_id\tdraw_id\nMissing file: '+str(self.attempt/'result.json')+'\t17\t0\n')
        result, out, summary = self.run_cli('history', '--query-run', run)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(summary['historical_current_status'], {'OK':1})
        history = self.rows('historical_missing.tsv.gz', out)
        self.assertEqual(history[0]['reference_id'], '17')
        self.write(query/'checksums.json', {})
        result, _, summary = self.run_cli('early', '--query-run', run)
        self.assertEqual(result.returncode, 0)
        self.assertEqual(summary['historical_observations'], 0)

    def test_invalid_manifest_and_exclusive_output(self):
        self.manifest.write_text('Wrong\ninput\n')
        result, out, summary = self.run_cli('invalid', '--manifest', self.manifest)
        self.assertEqual(result.returncode, 1)
        self.assertFalse(summary['enumeration_complete'])
        before = (out/'summary.json').read_bytes()
        second = subprocess.run([sys.executable, '-B', str(SCRIPT), '--manifest', str(self.manifest),
                                 '--outdir',str(out)], capture_output=True)
        self.assertNotEqual(second.returncode, 0)
        self.assertEqual(before, (out/'summary.json').read_bytes())


if __name__ == '__main__':
    unittest.main()
