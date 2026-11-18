"""Concurrency bounds, native integrity and deterministic single-writer catalogs."""
from contextlib import redirect_stdout
import io
import json
import multiprocessing
import os
from pathlib import Path
import sys
import threading
import time
import unittest
from unittest.mock import patch

sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import test_ref_depth_database as fixtures
from bin import native_execution as execution, native_reference as native
from bin.ref_depth_catalog import collect, create_plan, digest, exclusive
from bin.ref_depth_database import ProgressLog
from bin.ref_depth_worker import run_sample
from bin.coverage_core import Reference


def delayed_task(gcf,sample):
    """Picklable synthetic scheduling probe; no genomic computation."""
    from bin.native_execution import ROOT,SAMPLES
    parent=ROOT/gcf/'samples'/SAMPLES[gcf,sample]['folder']
    start=time.monotonic();time.sleep(.15);end=time.monotonic()
    (parent/'interval.json').write_text(json.dumps(dict(start=start,end=end,pid=os.getpid(),gcf=gcf)))
    return dict(gcf=gcf,sample=sample,state='COMPLETE',exit_code=0,reused=False,
                worker_pid=os.getpid(),log=str(parent/'worker.log'))


def interrupted_task(gcf,sample):
    raise InterruptedError('synthetic interrupted worker')


class SchedulingTests(unittest.TestCase):
    def setUp(self):
        self.fixture=fixtures.DepthDatabaseTests('runTest');self.fixture.setUp()
        self.root=self.fixture.root
        rows=[dict(Sample=f's{i:03d}',CramFilePath=str(self.fixture.a),GCFID=fixtures.GCF1) for i in range(8)]
        fixtures.write_csv(self.fixture.records,['Sample','CramFilePath','GCFID'],rows)
        config=dict(self.fixture.config,storage_mode='native',depths=[],window_size=7,bins=None,geometry_format='packed')
        self.database=self.root/'db'
        self.plan=create_plan(self.fixture.records,self.fixture.refs,self.database,config,fixtures.COLUMNS)
        self.group=self.plan['groups'][fixtures.GCF1]

    def tearDown(self):self.fixture.tearDown()

    def build(self):
        with redirect_stdout(io.StringIO()):
            for sample in self.group['samples']:
                self.assertEqual(run_sample(self.database,self.plan,fixtures.GCF1,sample,self.fixture.extractor)['exit_code'],0)

    def test_single_gcf_uses_bounded_persistent_workers(self):
        with patch.object(execution,'perform',delayed_task):
            report=execution.execute(self.database,self.plan,2,ProgressLog(self.database,self.plan,persist=False))
        intervals=[json.loads(p.read_text()) for p in self.database.glob('*/samples/*/interval.json')]
        events=sorted([(r['start'],1) for r in intervals]+[(r['end'],-1) for r in intervals])
        active=maximum=0
        for _,delta in events:active+=delta;maximum=max(maximum,active)
        self.assertEqual(maximum,2)
        self.assertEqual(len(report['worker_pids']),2)
        self.assertEqual(report['sample_task_count'],8)
        self.assertEqual(report['local_exit_codes'],{fixtures.GCF1:0})
        self.assertEqual(json.loads((self.database/fixtures.GCF1/'progress.json').read_text())['samples_finished'],8)

    def test_worker_interruption_propagates_and_pool_is_reaped(self):
        before={p.pid for p in multiprocessing.active_children()}
        with patch.object(execution,'perform',interrupted_task):
            with self.assertRaisesRegex(InterruptedError,'synthetic interrupted'):
                execution.execute(self.database,self.plan,2,ProgressLog(self.database,self.plan,persist=False))
        self.assertFalse({p.pid for p in multiprocessing.active_children()}-before)
        self.assertEqual(json.loads((self.database/fixtures.GCF1/'progress.json').read_text())['state'],'INTERRUPTED')
        self.assertFalse((self.database/'local_sample_exit_codes.json').exists())

    def test_parallel_collection_bounded_deterministic_and_fai_cached(self):
        self.build();self.assertTrue(collect(self.database,jobs=1)['complete'])
        names=['GCF_INDEX.tsv','SAMPLE_INDEX.tsv','NATIVE_INDEX.tsv','RESOURCE_INDEX.tsv']
        before={n:(self.database/n).read_bytes() for n in names}
        original=native.verify_native_sample;guard=threading.Lock();active=maximum=0
        def observe(item):
            nonlocal active,maximum
            with guard:active+=1;maximum=max(maximum,active)
            try:time.sleep(.03);return original(item)
            finally:
                with guard:active-=1
        original_digest=native.digest;fai_calls=[]
        def counted(path):
            if str(path)==str(self.fixture.fai):fai_calls.append(path)
            return original_digest(path)
        with patch.object(native,'verify_native_sample',observe),patch.object(native,'digest',counted):
            self.assertTrue(collect(self.database,jobs=3)['complete'])
        self.assertEqual(maximum,3);self.assertEqual(len(fai_calls),2)
        self.assertEqual(before,{n:(self.database/n).read_bytes() for n in names})
        chunk=next(self.database.glob('*/samples/*/attempt-*/geometry/chunk-*.zlib'))
        chunk.write_bytes(chunk.read_bytes()+b'broken')
        result=collect(self.database,jobs=3)
        self.assertFalse(result['complete']);self.assertEqual(result['failed_samples'],1)

    def test_sample_lock_reuse_and_corrupt_attempt_retry(self):
        sample=self.group['samples'][0];parent=self.database/fixtures.GCF1/'samples'/sample['folder']
        with exclusive(parent/'.sample.lock'):
            with self.assertRaisesRegex(RuntimeError,'Locked'):
                run_sample(self.database,self.plan,fixtures.GCF1,sample,self.fixture.extractor)
        self.assertFalse((parent/'progress.json').exists())
        with redirect_stdout(io.StringIO()):
            first=run_sample(self.database,self.plan,fixtures.GCF1,sample,self.fixture.extractor)
            self.assertFalse(first['reused'])
            self.assertTrue(run_sample(self.database,self.plan,fixtures.GCF1,sample,self.fixture.extractor)['reused'])
            chunk=next(parent.glob('attempt-*/geometry/chunk-*.zlib'));chunk.write_bytes(b'corrupted')
            repaired=run_sample(self.database,self.plan,fixtures.GCF1,sample,self.fixture.extractor)
        self.assertEqual(repaired['exit_code'],0);self.assertFalse(repaired['reused'])
        self.assertEqual(len(list(parent.glob('attempt-*'))),2)
        self.assertEqual(chunk.read_bytes(),b'corrupted')

    def test_same_gcf_conflicting_known_m5_is_not_published_complete(self):
        def different(path,reference_index,**policies):
            pool,layout=self.fixture.extractor(path,reference_index,**policies)
            return pool,Reference(layout.names,layout.lengths,('different-known-m5',layout.md5[1]))
        with redirect_stdout(io.StringIO()):
            for sample in self.group['samples'][:-1]:
                self.assertEqual(run_sample(self.database,self.plan,fixtures.GCF1,sample,self.fixture.extractor)['exit_code'],0)
            self.assertEqual(run_sample(self.database,self.plan,fixtures.GCF1,self.group['samples'][-1],different)['exit_code'],0)
        result=collect(self.database,jobs=3)
        self.assertFalse(result['complete']);self.assertEqual(result['failed_samples'],1)

    def test_invalid_collection_limit_rejected_before_writes(self):
        before=(self.database/'PLAN.sha256').read_bytes()
        for jobs in (0,-1,True,1.5):
            with self.assertRaises(ValueError):collect(self.database,jobs=jobs)
        self.assertEqual((self.database/'PLAN.sha256').read_bytes(),before)

    def test_fai_change_during_parallel_collection_rejects_entire_group(self):
        self.build();original=native.verify_native_sample
        def changed(item):
            result=original(item)
            if item[3]['sample']==self.group['samples'][-1]['sample']:
                self.fixture.fai.write_text(self.fixture.fai.read_text()+'changed\t1\n')
            return result
        observed=[]
        with patch.object(native,'verify_native_sample',changed):
            result=collect(self.database,jobs=3,progress=lambda *values:observed.append(values))
        self.assertFalse(result['complete']);self.assertEqual(result['failed_samples'],8)
        self.assertEqual(observed[-1],(8,8,8))
