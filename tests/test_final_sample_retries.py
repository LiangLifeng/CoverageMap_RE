"""Final native-local retries use verified failures and preserve healthy bytes."""
from contextlib import redirect_stderr, redirect_stdout
import io
import json
from pathlib import Path
import sys
import unittest
from unittest.mock import patch

sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from bin import native_execution as execution, ref_depth_database as cli
from bin.coverage_core import BlockPool, Reference
from bin.ref_depth_catalog import collect, digest
from bin.ref_depth_worker import run_sample
import test_ref_depth_database as fixtures


def synthetic_task(gcf, sample):
    """Picklable probe runs real build/validation; only extraction is synthetic."""
    from bin.native_execution import ROOT, PLAN, SAMPLES
    record=SAMPLES[gcf,sample]
    parent=ROOT/gcf/'samples'/record['folder']
    def extractor(path, reference_index, **policies):
        fail_attempts=json.loads((ROOT/'fault.json').read_text())['fail_attempts']
        if sample=='a' and len(list(parent.glob('attempt-*'))) <= fail_attempts:
            raise OSError('injected sample failure')
        group=PLAN['groups'][gcf]
        layout=Reference(tuple(group['layout']['names']),tuple(group['layout']['lengths']),
                         tuple('' for _ in group['layout']['names']))
        return BlockPool.from_intervals([[0,3],[5,9]],layout.length),layout
    result=run_sample(ROOT,PLAN,gcf,record,extractor)
    result['log']=str(parent/'worker.log')
    return result


class FinalSampleRetryTests(unittest.TestCase):
    def setUp(self):
        self.fixture=fixtures.DepthDatabaseTests('runTest');self.fixture.setUp()
        self.database,self.plan=self.fixture.plan(storage_mode='native',window_size=7,depths=[],geometry_format='packed')
        self.samples={s['sample']:self.database/fixtures.GCF1/'samples'/s['folder']
                      for s in self.plan['groups'][fixtures.GCF1]['samples']}

    def tearDown(self):
        self.fixture.tearDown()

    def run_case(self, fail_attempts, corrupt_after_worker=False):
        (self.database/'fault.json').write_text(json.dumps(dict(fail_attempts=fail_attempts)))
        pins={};calls=[];old_corrupt=[]
        def observed_collect(*args,**kwargs):
            if not calls and corrupt_after_worker:
                chunk=next(self.samples['a'].glob('attempt-*/geometry/chunk-*.zlib'))
                chunk.write_bytes(b'injected corruption');old_corrupt.append(chunk)
            summary=collect(*args,**kwargs)
            calls.append(summary)
            current={str(p):digest(p) for p in self.samples['b'].rglob('*') if p.is_file()}
            if pins:self.assertEqual(pins,current)
            else:pins.update(current)
            return summary
        with (patch.object(execution,'perform',synthetic_task),patch.object(cli,'collect',observed_collect),
              redirect_stdout(io.StringIO()),redirect_stderr(io.StringIO())):
            result=cli.execute_local(self.database,jobs=2,collect_jobs=2)
        history=json.loads((self.database/'sample_retries.json').read_text())['invocations'][-1]
        self.assertEqual(len(list(self.samples['b'].glob('attempt-*'))),1)
        if old_corrupt:self.assertEqual(old_corrupt[0].read_bytes(),b'injected corruption')
        return result,history,calls

    def test_second_retry_succeeds_with_three_attempts_and_healthy_bytes_preserved(self):
        result,history,calls=self.run_case(2)
        self.assertTrue(result['complete'])
        self.assertEqual(result['sample_retry_rounds'],2)
        self.assertEqual(result['sample_retry_tasks'],2)
        self.assertEqual(result['local_exit_codes'],{fixtures.GCF1:0})
        self.assertEqual(len(calls),3)
        self.assertEqual(len(list(self.samples['a'].glob('attempt-*'))),3)
        self.assertTrue((self.samples['a']/'attempt-000001/status/FAILED').exists())
        self.assertTrue((self.samples['a']/'attempt-000002/status/FAILED').exists())
        self.assertTrue((self.samples['a']/'attempt-000003/status/SUCCESS').exists())
        self.assertEqual([r['round'] for r in history['rounds']],[1,2])
        self.assertTrue(all(r['samples']==[dict(gcf=fixtures.GCF1,sample='a')] for r in history['rounds']))
        self.assertTrue(history['complete'])

    def test_persistent_failure_stops_after_two_additional_attempts(self):
        result,history,calls=self.run_case(99)
        self.assertFalse(result['complete'])
        self.assertEqual(result['failed_samples'],1)
        self.assertEqual(result['sample_retry_rounds'],2)
        self.assertEqual(len(list(self.samples['a'].glob('attempt-*'))),3)
        self.assertFalse(list(self.samples['a'].glob('attempt-*/status/SUCCESS')))
        self.assertEqual(result['local_exit_codes'],{fixtures.GCF1:1})
        self.assertEqual(history['remaining_failed_samples'],1)

    def test_collection_only_failure_is_rebuilt_and_recollected(self):
        result,history,calls=self.run_case(0,corrupt_after_worker=True)
        self.assertTrue(result['complete'])
        self.assertEqual(result['sample_retry_rounds'],1)
        self.assertEqual(len(list(self.samples['a'].glob('attempt-*'))),2)
        self.assertEqual(len(calls),2)

    def test_successful_first_pass_does_not_retry(self):
        result,history,calls=self.run_case(0)
        self.assertTrue(result['complete'])
        self.assertEqual(result['sample_retry_rounds'],0)
        self.assertEqual(history['rounds'],[])
        self.assertEqual(len(calls),1)

    def test_interruption_never_starts_final_retries(self):
        with patch.object(execution,'execute',side_effect=InterruptedError('stop')) as execute:
            with self.assertRaises(InterruptedError):cli.execute_local(self.database,jobs=2)
        self.assertEqual(execute.call_count,1)
        self.assertFalse((self.database/'sample_retries.json').exists())
