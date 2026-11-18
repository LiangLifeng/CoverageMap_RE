"""Opt-in sample omission preserves validation, atomic rows and source provenance."""
import csv
import json
from pathlib import Path
import sqlite3
import unittest
from unittest.mock import patch

import test_ref_depth_database as fixtures
from bin.ref_depth_catalog import collect, create_plan, digest
from bin.ref_depth_worker import run_gcf
from bin.reference_index import ReferenceIndex, build_index
from bin.reference_registry import export_report, filter_registry


class SkipInvalidTests(unittest.TestCase):
    def setUp(self):
        self.f=fixtures.DepthDatabaseTests('runTest');self.f.setUp();self.root=self.f.root
        self.config=dict(self.f.config,storage_mode='native',depths=[],bins=None,window_size=7)
        self.database=self.root/'native'
        self.plan=create_plan(self.f.records,self.f.refs,self.database,self.config,fixtures.COLUMNS)
        self.assertEqual(run_gcf(self.database,fixtures.GCF1,self.f.extractor),0)
        self.assertTrue(collect(self.database)['complete'])
        self.path=self.root/'index.sqlite'

    def tearDown(self):self.f.tearDown()

    def attempt(self,sample='b'):
        row=next(s for s in self.plan['groups'][fixtures.GCF1]['samples'] if s['sample']==sample)
        return next((self.database/fixtures.GCF1/'samples'/row['folder']).glob('attempt-*'))

    def snapshot(self):return {str(p):digest(p) for p in self.database.rglob('*') if p.is_file()}

    def test_incomplete_source_skips_failed_only_and_retains_audit_after_refilter(self):
        attempt=self.attempt();(attempt/'status/SUCCESS').unlink();(attempt/'status/FAILED').touch()
        (attempt/'status/exit_code.txt').write_text('1\n')
        summary_path=self.database/'SUMMARY.json';summary=json.loads(summary_path.read_text())
        summary.update(complete=False,ready_samples=1);summary_path.write_text(json.dumps(summary))
        before=self.snapshot()
        with self.assertRaisesRegex(ValueError,'completed reference database'):
            build_index(self.database,self.root/'strict.sqlite',True)
        info=build_index(self.database,self.path,True,skip_invalid_samples=True)
        self.assertEqual(info['counts']['references'],1)
        self.assertEqual((info['skipped_samples'],info['planned_input_references'],info['partial_index']),(1,2,True))
        self.assertFalse(info['sources'][0]['source_summary_complete'])
        index=ReferenceIndex(self.path)
        try:
            self.assertEqual([r['sample'] for r in index.nearest(fixtures.GCF1,1,10)],['a'])
            row=index.connection.execute('SELECT * FROM skipped_samples').fetchone()
            self.assertEqual(row['sample'],'b');self.assertIn('Missing successful',row['error'])
            self.assertIsNone(index.connection.execute('PRAGMA foreign_key_check').fetchone())
        finally:index.close()
        filtered=self.root/'filtered.sqlite';filter_registry(self.path,filtered,0,2)
        report=self.root/'report';export_report(filtered,report)
        with (report/'skipped_samples.tsv').open() as handle:rows=list(csv.DictReader(handle,delimiter='\t'))
        self.assertEqual(len(rows),1);self.assertEqual(rows[0]['sample'],'b')
        with (report/'gcf_depths.tsv').open() as handle:groups=list(csv.DictReader(handle,delimiter='\t'))
        self.assertEqual((groups[0]['total_samples'],groups[0]['selected_samples'],groups[0]['skipped_copies']),('1','0','1'))
        self.assertTrue(json.loads((report/'selection.json').read_text())['partial_index'])
        self.assertEqual(before,self.snapshot())

    def test_complete_source_corrupt_sample_is_omitted_without_bypassing_checksum(self):
        attempt=self.attempt();checks=attempt/'checksums.json';data=json.loads(checks.read_text())
        del data['edges.npy'];checks.write_text(json.dumps(data));before=self.snapshot()
        with self.assertRaisesRegex(ValueError,'checksum'):
            build_index(self.database,self.root/'strict.sqlite',True)
        info=build_index(self.database,self.path,True,skip_invalid_samples=True)
        self.assertEqual((info['counts']['references'],info['skipped_samples']),(1,1))
        index=ReferenceIndex(self.path)
        try:
            self.assertIn('checksum',index.connection.execute('SELECT error FROM skipped_samples').fetchone()[0])
            self.assertEqual(index.metadata['selected_counts']['references'],1)
        finally:index.close()
        self.assertEqual(before,self.snapshot())

    def test_sample_transaction_rolls_back_references_depths_and_counts(self):
        # A late legacy metadata error occurs after the reference INSERT and first depth INSERT.
        legacy,plan=self.f.plan('legacy');self.assertEqual(run_gcf(legacy,fixtures.GCF1,self.f.extractor),0);collect(legacy)
        sample=next(s for s in plan['groups'][fixtures.GCF1]['samples'] if s['sample']=='b')
        attempt=next((legacy/fixtures.GCF1/'samples'/sample['folder']).glob('attempt-*'))
        record=json.loads((attempt/'result.json').read_text());bad=attempt/record['depths'][1]['reference_metadata']
        from bin.reference_index import relative_file
        def fail(root,path):
            if Path(path)==bad:raise ValueError('injected late metadata failure')
            return relative_file(root,path)
        with patch('bin.reference_index.relative_file',side_effect=fail):
            info=build_index(legacy,self.path,True,skip_invalid_samples=True)
        self.assertEqual(info['counts'],dict(gcfs=1,references=1,depths=3,ready=2))
        index=ReferenceIndex(self.path)
        try:
            self.assertEqual(index.connection.execute('SELECT COUNT(*) FROM depths').fetchone()[0],3)
            self.assertEqual(index.connection.execute('SELECT COUNT(*) FROM references_table').fetchone()[0],1)
            self.assertIsNone(index.connection.execute('PRAGMA foreign_key_check').fetchone())
        finally:index.close()

    def test_all_invalid_group_retained_as_empty_and_systemic_errors_still_fail(self):
        for name in ('a','b'):(self.attempt(name)/'status/SUCCESS').unlink()
        info=build_index(self.database,self.path,True,skip_invalid_samples=True)
        self.assertEqual((info['counts']['references'],info['skipped_samples']),(0,2))
        self.assertEqual(info['selected_counts'],dict(gcfs=0,references=0))
        index=ReferenceIndex(self.path)
        try:
            self.assertIsNone(index.group(fixtures.GCF1,True))
            self.assertEqual(index.connection.execute('SELECT total_samples FROM gcf_statistics').fetchone()[0],0)
        finally:index.close()
        for error in (KeyboardInterrupt(),sqlite3.OperationalError('disk full')):
            with patch('bin.reference_index.completed_sample',side_effect=error):
                with self.assertRaises(type(error)):
                    build_index(self.database,self.root/'system-error.sqlite',skip_invalid_samples=True)
            self.assertFalse((self.root/'system-error.sqlite').exists())

    def test_three_sources_keep_healthy_copy_of_invalid_pair_and_skipped_source_audit(self):
        import shutil
        roots=[self.database]
        for name in ('zzcopy-a','zzcopy-b'):
            copy=self.root/name;shutil.copytree(self.database,copy);roots.append(copy)
        (self.attempt()/'status/SUCCESS').unlink()
        info=build_index(roots,self.path,True,skip_invalid_samples=True)
        self.assertEqual((info['counts']['references'],info['skipped_samples'],info['duplicate_copies']),(2,1,3))
        self.assertEqual(info['planned_input_references'],6)
        index=ReferenceIndex(self.path)
        try:
            self.assertEqual(len(index.nearest(fixtures.GCF1,1,10)),2)
            self.assertEqual(index.connection.execute('SELECT COUNT(*) FROM reference_sources').fetchone()[0],5)
            self.assertEqual(index.connection.execute('SELECT database_root FROM skipped_samples').fetchone()[0],str(self.database))
        finally:index.close()

    def test_healthy_source_with_flag_is_not_marked_partial(self):
        info=build_index(self.database,self.path,True,skip_invalid_samples=True)
        self.assertEqual((info['skipped_samples'],info['partial_index']),(0,False))
        report=self.root/'report';export_report(self.path,report)
        self.assertEqual(len((report/'skipped_samples.tsv').read_text().splitlines()),1)

    def test_existing_schema3_without_skip_table_still_reports(self):
        build_index(self.database,self.path,True)
        current=ReferenceIndex(self.path);older_path=self.root/'old-schema3.sqlite'
        connection=sqlite3.connect(older_path)
        try:
            current.connection.backup(connection)
            connection.execute('DROP TABLE skipped_samples');connection.commit()
        finally:current.close();connection.close()
        report=self.root/'old-report';export_report(older_path,report)
        self.assertEqual(len((report/'skipped_samples.tsv').read_text().splitlines()),1)

    def test_cli_flag_publishes_partial_index_and_exports_errors(self):
        import subprocess
        import sys
        (self.attempt()/'status/SUCCESS').unlink()
        summary_path=self.database/'SUMMARY.json';summary=json.loads(summary_path.read_text())
        summary.update(complete=False,ready_samples=1);summary_path.write_text(json.dumps(summary))
        entry=Path(__file__).resolve().parents[1]/'bin/reference_index.py'
        command=[sys.executable,'-B',str(entry),'build','--database',str(self.database),
                 '--out',str(self.path),'--verify-payloads']
        strict=subprocess.run(command,capture_output=True,text=True)
        self.assertNotEqual(strict.returncode,0);self.assertFalse(self.path.exists())
        allowed=subprocess.run([*command,'--skip-invalid-samples'],capture_output=True,text=True)
        self.assertEqual(allowed.returncode,0,allowed.stderr)
        result=json.loads(allowed.stdout)
        self.assertEqual((result['skipped_samples'],result['partial_index']),(1,True))
        self.assertIn('INDEX SKIP '+fixtures.GCF1+'/b',allowed.stderr)
        self.assertIn('INDEX PARTIAL',allowed.stderr)
        report=self.root/'cli-report'
        subprocess.run([sys.executable,'-B',str(entry),'report','--index',str(self.path),'--outdir',str(report)],
                       check=True,capture_output=True,text=True)
        self.assertIn('\tb\t',(report/'skipped_samples.tsv').read_text())

    def test_external_only_audit_removes_table_and_deleted_error_bytes(self):
        (self.attempt()/'status/SUCCESS').unlink()
        audit=self.root/'skipped_samples.tsv'
        marker='PRIVATE_ERROR_TEXT_NOT_TO_BE_STORED_IN_SQLITE'
        from bin.reference_index import completed_sample
        def inspect(path,*args,**kwargs):
            if path==self.attempt().parent:raise ValueError(marker)
            return completed_sample(path,*args,**kwargs)
        with patch('bin.reference_index.completed_sample',side_effect=inspect):
            info=build_index(self.database,self.path,True,skip_invalid_samples=True,skipped_samples_tsv=audit)
        self.assertIn(marker,audit.read_text());self.assertNotIn(marker.encode(),self.path.read_bytes())
        self.assertEqual(info['skipped_samples_audit']['sha256'],digest(audit))
        index=ReferenceIndex(self.path)
        try:
            self.assertIsNone(index.connection.execute("SELECT 1 FROM sqlite_master WHERE name='skipped_samples'").fetchone())
            self.assertEqual(index.connection.execute('SELECT copies FROM skipped_sample_counts').fetchone()[0],1)
            self.assertEqual(index.metadata['skipped_samples'],1)
        finally:index.close()
        filtered=self.root/'external-filtered.sqlite';filter_registry(self.path,filtered,0,2)
        index=ReferenceIndex(filtered)
        try:
            self.assertIsNone(index.connection.execute("SELECT 1 FROM sqlite_master WHERE name='skipped_samples'").fetchone())
        finally:index.close()
        report=self.root/'external-report';export_report(filtered,report)
        self.assertIn(marker,(report/'skipped_samples.tsv').read_text())
        with (report/'gcf_depths.tsv').open() as handle:row=next(csv.DictReader(handle,delimiter='\t'))
        self.assertEqual(row['skipped_copies'],'1')

    def test_existing_index_can_externalize_during_filter_and_relocate_audit(self):
        (self.attempt()/'status/SUCCESS').unlink()
        build_index(self.database,self.path,True,skip_invalid_samples=True)
        before=digest(self.path);out=self.root/'only-valid.sqlite';audit=self.root/'errors.tsv'
        filter_registry(self.path,out,0,1,skipped_samples_tsv=audit)
        self.assertEqual(digest(self.path),before)
        index=ReferenceIndex(out)
        try:self.assertIsNone(index.connection.execute("SELECT 1 FROM sqlite_master WHERE name='skipped_samples'").fetchone())
        finally:index.close()
        copied=self.root/'copied.sqlite';new_audit=self.root/'copied-errors.tsv'
        filter_registry(out,copied,0,1,skipped_samples_tsv=new_audit)
        self.assertEqual(new_audit.read_bytes(),audit.read_bytes())
        # Missing or mutated sidecars cannot silently produce a blank error report.
        audit.chmod(0o644);audit.write_text('corrupted')
        with self.assertRaisesRegex(ValueError,'checksum'):
            export_report(out,self.root/'bad-report')
        self.assertFalse((self.root/'bad-report').exists())
        export_report(copied,self.root/'good-report')

    def test_audit_targets_are_exclusive_and_failed_index_publish_cleans_own_tsv(self):
        audit=self.root/'exists.tsv';audit.write_text('user content')
        with self.assertRaises(FileExistsError):
            build_index(self.database,self.path,True,skipped_samples_tsv=audit)
        self.assertEqual(audit.read_text(),'user content');self.assertFalse(self.path.exists())
        with self.assertRaises(ValueError):
            build_index(self.database,self.path,True,skipped_samples_tsv=self.path)
        with self.assertRaises(ValueError):
            build_index(self.database,self.path,True,skipped_samples_tsv=self.database/'bad.tsv')
        import os
        original_link=os.link
        def fail_index(source,target):
            if Path(target)==self.path:raise FileExistsError('injected publication race')
            return original_link(source,target)
        new_audit=self.root/'new-errors.tsv'
        with patch('bin.reference_registry.os.link',side_effect=fail_index):
            with self.assertRaises(FileExistsError):
                build_index(self.database,self.path,True,skipped_samples_tsv=new_audit)
        self.assertFalse(self.path.exists());self.assertFalse(new_audit.exists())

    def test_external_audit_cli_build_and_filter_parameters(self):
        import subprocess
        import sys
        (self.attempt()/'status/SUCCESS').unlink()
        entry=Path(__file__).resolve().parents[1]/'bin/reference_index.py'
        audit=self.root/'external-cli.tsv'
        result=subprocess.run([sys.executable,'-B',str(entry),'build','--database',str(self.database),
            '--out',str(self.path),'--skip-invalid-samples','--verify-payloads','--skipped-samples-tsv',str(audit)],
            check=True,capture_output=True,text=True)
        self.assertEqual(json.loads(result.stdout)['skipped_samples_audit']['path'],str(audit))
        filtered=self.root/'external-cli-filtered.sqlite';new_audit=self.root/'external-cli-filtered.tsv'
        subprocess.run([sys.executable,'-B',str(entry),'filter','--index',str(self.path),'--out',str(filtered),
            '--min-depth','0','--min-samples','2','--skipped-samples-tsv',str(new_audit)],
            check=True,capture_output=True,text=True)
        self.assertEqual(new_audit.read_bytes(),audit.read_bytes())
        index=ReferenceIndex(filtered)
        try:
            self.assertEqual(index.metadata['selected_counts']['references'],0)
            self.assertIsNone(index.connection.execute("SELECT 1 FROM sqlite_master WHERE name='skipped_samples'").fetchone())
        finally:index.close()
