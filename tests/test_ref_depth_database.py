"""Scientific invariants, restart/catalog behavior and local/Slurm interfaces."""
import csv
from contextlib import redirect_stderr, redirect_stdout
import io
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from bin import ref_depth_database as cli
from bin.ref_depth_catalog import (collect, completed_sample, create_plan, digest, exclusive,
                                   load_plan, sample_folder)
from bin.ref_depth_worker import load_depth, run_gcf
from bin.coverage_core import BlockPool, Reference, canonical_reference, dense_coverage
from bin.concordance_core import task_rng
from bin.concordance_core import evaluate_windows

GCF1, GCF2 = 'GCF_000000001.1', 'GCF_000000002.1'
COLUMNS = dict(sample='Sample', cram='CramFilePath', gcf='GCFID', ref_gcf='GCF_ID', fai='GCF_path')


def write_csv(path, columns, rows):
    with path.open('w', newline='') as handle:
        writer = csv.DictWriter(handle, columns)
        writer.writeheader(); writer.writerows(rows)


class DepthDatabaseTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary.name).resolve()
        self.records, self.refs = self.root / 'records.csv', self.root / 'refs.csv'
        self.fai = self.root / 'one.fai'
        self.fai.write_text('one\t20\ntwo\t11\n')
        self.a, self.b = self.root / 'a.cram', self.root / 'b.cram'
        self.a.write_bytes(b'toy a'); self.b.write_bytes(b'toy b')
        for path in (self.a,self.b):
            Path(str(path)+'.crai').write_bytes(b'fixture sidecar for provenance checks')
        self.rows = [dict(Sample='a', CramFilePath=str(self.a), GCFID=GCF1),
                     dict(Sample='b', CramFilePath=str(self.b), GCFID=GCF1)]
        write_csv(self.records, ['Sample', 'CramFilePath', 'GCFID'], self.rows)
        write_csv(self.refs, ['GCF_ID', 'GCF_path'], [dict(GCF_ID=GCF1, GCF_path=str(self.fai))])
        self.config = dict(depths=[.1, 1., 10.], cpus=2, memory='1G', python=sys.executable,
            replicates=5, bins=7, seed=401, batch_size=2, sampler='permutation',
            geometry_format='npz', matrix_format='npz', insufficient_depth='record',
            policies=dict(mate_policy='both', sampling_unit='read', min_mapq=0,
                          exclude_duplicates=False, exclude_qcfail=False))
        self.layout = Reference(('one', 'two'), (20, 11), ('aaa', 'bbb'))
        self.small = BlockPool(np.array([0, 21, 3, 10, 2], dtype=np.uint32),
            np.array([4, 26, 5, 15, 29], dtype=np.uint32), np.array([0, 2, 4, 5], dtype=np.uint64), 31)
        self.big = BlockPool.from_intervals(np.tile(np.array([[1, 11]]), (50, 1)), 31)

    def tearDown(self):
        self.temporary.cleanup()

    def plan(self, name='db', **config):
        output = self.root / name
        plan = create_plan(self.records, self.refs, output, dict(self.config, **config), COLUMNS)
        return output, plan

    def extractor(self, path, reference_index, sort_contigs=False, **policies):
        self.assertEqual(policies, self.config['policies'])
        pool = self.small if Path(path) == self.a else self.big
        if sort_contigs:
            layout = canonical_reference(self.layout)
            self.assertEqual(layout.names, self.layout.names)
            return pool, layout
        return pool, self.layout

    def run_args(self, mode='local'):
        return cli.parser().parse_args(['run', '--records', str(self.records), '--references', str(self.refs),
            '--outdir', str(self.root / 'one-command'), '--mode', mode,
            '--window-size', '7', '--cpus', '2', '--interval', '.01', '--alignment-capture', 'legacy'])

    def test_progress_counts_samples_not_equal_weight_gcfs(self):
        plan={'groups':{'a':{'samples':[1,2,3]}, 'b':{'samples':[1]}}}
        rows=[dict(gcf='a',state='RUNNING',samples_finished=1,failed_samples=1,
                   samples_reused=1,current_sample='third'),dict(gcf='b',state='COMPLETE',samples_finished=1)]
        value=cli.progress_counts(plan,rows)
        self.assertEqual((value['samples_processed'],value['samples_total']),(3,4))
        self.assertEqual((value['samples_complete'],value['samples_failed'],value['samples_reused']),(2,1,1))
        self.assertEqual((value['gcfs_terminal'],value['gcfs_total']),(1,2))
        self.assertEqual(value['current_samples'],['a/third'])

    def test_progress_heartbeat_is_throttled_and_persisted(self):
        database,plan=self.plan('progress-log')
        capture=io.StringIO()
        with patch.object(cli.time,'monotonic',side_effect=[0,0,1,10]),redirect_stderr(capture):
            logger=cli.ProgressLog(database,plan)
            logger.update();logger.update();logger.update()
        lines=capture.getvalue().splitlines()
        self.assertEqual(len(lines),2)
        self.assertIn('samples=0/2',lines[0]);self.assertIn('elapsed=10.0s',lines[1])
        self.assertEqual((database/'build.log').read_text(),capture.getvalue())

    def test_read_only_watch_terminal_does_not_collect_or_write(self):
        database,plan=self.plan('read-only-watch')
        with redirect_stdout(io.StringIO()):
            run_gcf(database,GCF1,self.extractor)
        before={str(p):digest(p) for p in database.rglob('*') if p.is_file()}
        args=['ref_depth_database.py','status','--database',str(database),'--watch','--no-collect']
        stderr=io.StringIO()
        with patch.object(sys,'argv',args),patch.object(cli,'collect') as collect_mock,redirect_stderr(stderr),redirect_stdout(io.StringIO()):
            self.assertEqual(cli.main(),0)
        collect_mock.assert_not_called()
        self.assertIn('samples=2/2',stderr.getvalue())
        self.assertEqual(before,{str(p):digest(p) for p in database.rglob('*') if p.is_file()})

    def test_worker_logs_failure_and_final_progress(self):
        database,plan=self.plan('failure-progress')
        def extract(path,reference_index,**policies):
            if Path(path)==self.a:
                raise ValueError('test input rejected')
            return self.extractor(path,reference_index,**policies)
        output=io.StringIO()
        with redirect_stdout(output):
            self.assertEqual(run_gcf(database,GCF1,extract),1)
        self.assertIn('SAMPLE FAILED a: ValueError: test input rejected',output.getvalue())
        self.assertIn('SAMPLE OK b',output.getvalue())
        p=json.loads((database/GCF1/'progress.json').read_text())
        self.assertEqual((p['state'],p['samples_finished'],p['failed_samples']),('FAILED',1,1))

    def test_crai_optional_before_plan_and_existing_archive_preserved(self):
        crai=Path(str(self.a)+'.crai'); saved=crai.read_bytes()
        crai.unlink()
        self.a.with_suffix('.crai').write_bytes(saved)  # Alternate name is not required or implicitly selected.
        missing_db,missing_plan=self.plan('missing-optional')
        self.assertIsNone(missing_plan['groups'][GCF1]['samples'][0]['cram_index'])
        self.assertEqual(missing_plan['cram_index_policy'],'optional-source-cram-crai-v1')
        self.assertEqual(run_gcf(missing_db,GCF1,self.extractor),0)
        self.assertTrue(collect(missing_db)['complete'])
        parent=missing_db/GCF1/'samples'/missing_plan['groups'][GCF1]['samples'][0]['folder']
        attempt,result=completed_sample(parent,digest(missing_db/'database.json'))
        self.assertIsNone(result['cram_index'])
        self.assertFalse((attempt/'source.cram.crai').exists())
        crai.write_bytes(b'')
        with self.assertRaisesRegex(ValueError,'empty'):
            self.plan('empty')
        crai.write_bytes(saved)
        database,plan=self.plan('valid-index')
        self.assertEqual(run_gcf(database,GCF1,self.extractor),0)
        sample=plan['groups'][GCF1]['samples'][0]
        self.assertEqual(sample['cram_index']['sha256'],digest(crai))
        attempt,result=completed_sample(database/GCF1/'samples'/sample['folder'],digest(database/'database.json'))
        self.assertEqual((attempt/'source.cram.crai').read_bytes(),saved)
        self.assertEqual(result['cram_index'],sample['cram_index'])
        self.assertEqual(crai.read_bytes(),saved)
        self.assertTrue(collect(database)['complete'])
        # Same-length changed contents and restored mtime must still be rejected by SHA.
        stat=crai.stat();crai.write_bytes(b'X'*len(saved));os.utime(crai,ns=(stat.st_atime_ns,stat.st_mtime_ns))
        self.assertEqual(run_gcf(database,GCF1,self.extractor),1)
        self.assertFalse(collect(database)['complete'])

    def test_empty_cram_is_still_rejected(self):
        self.a.write_bytes(b'')
        with self.assertRaisesRegex(ValueError,'Input file is empty'):
            self.plan('empty-source')
        self.assertFalse((self.root/'empty-source').exists())

    def test_crai_deletion_and_mutation_during_extraction_fail(self):
        database,plan=self.plan('deleted-index')
        Path(str(self.a)+'.crai').unlink()
        self.assertEqual(run_gcf(database,GCF1,self.extractor),1)
        self.assertFalse(collect(database)['complete'])
        Path(str(self.a)+'.crai').write_bytes(b'fixture sidecar')
        database,plan=self.plan('mutate-index')
        def changed(path,**kwargs):
            answer=self.extractor(path,**kwargs)
            Path(str(path)+'.crai').write_bytes(b'changed during extraction')
            return answer
        self.assertEqual(run_gcf(database,GCF1,changed),1)
        self.assertFalse(collect(database)['complete'])

    def test_required_payload_empty_or_unregistered_is_rejected(self):
        database,plan=self.plan('payload-validation')
        self.assertEqual(run_gcf(database,GCF1,self.extractor),0)
        parent=database/GCF1/'samples'/plan['groups'][GCF1]['samples'][0]['folder']
        fingerprint=digest(database/'database.json')
        attempt,_=completed_sample(parent,fingerprint)
        archive=attempt/'source.cram.crai';original=archive.read_bytes()
        archive.write_bytes(b'')
        with self.assertRaises(ValueError):completed_sample(parent,fingerprint,verify=False)
        archive.write_bytes(original)
        registry=attempt/'checksums.json';checks=json.loads(registry.read_text())
        removed=checks.pop('geometry/storage.json');registry.write_text(json.dumps(checks))
        with self.assertRaisesRegex(ValueError,'Missing required payload checksum'):
            completed_sample(parent,fingerprint,verify=False)
        checks['geometry/storage.json']=removed;registry.write_text(json.dumps(checks))
        reference=next(attempt.glob('coverage/*/reference.json'));reference.write_bytes(b'')
        with self.assertRaisesRegex(ValueError,'empty'):
            completed_sample(parent,fingerprint,verify=False)

    def test_slurm_dry_run_omits_memory_and_never_starts_jobs(self):
        args=self.run_args('slurm');args.dry_run=True;args.partition='cpu';args.cpus=10
        with patch('bin.ref_depth_database.subprocess.run') as called, patch('bin.ref_depth_database.subprocess.Popen') as spawned:
            result=cli.run(args)
        called.assert_not_called();spawned.assert_not_called()
        self.assertEqual((result['state'],result['submitted'],result['scheduler_checked']),('DRY_RUN',False,False))
        self.assertIsNone(result['resources']['memory'])
        self.assertEqual(result['resources']['cpus'],10)
        script=(args.outdir/GCF1/'job.sbatch.sh').read_text()
        self.assertNotIn('--mem',script);self.assertIn('--partition=cpu',script)
        self.assertIn('--cpus-per-task=10',script)
        self.assertFalse((args.outdir/'submission.json').exists())
        self.assertFalse(list(args.outdir.glob('*/samples/*/attempt-*')))
        self.assertEqual(json.loads((args.outdir/'pipeline.json').read_text())['state'],'DRY_RUN')
        command=result['jobs'][0]['command']
        self.assertIn('--cpus-per-task=10',command);self.assertEqual(command[command.index('-p')+1],'cpu')
        self.assertFalse(any(c.startswith('--mem') for c in command))

    def test_existing_plan_preview_and_scheduler_default_environment(self):
        database,plan=self.plan('default-memory',memory=None)
        with patch('bin.ref_depth_database.subprocess.run') as called:
            result=cli.submit(database,preview=True)
        called.assert_not_called();self.assertEqual(result['state'],'DRY_RUN')
        self.assertFalse((database/'submission.json').exists())
        with patch.dict(os.environ,dict(SBATCH_MEM_PER_NODE='999G',SBATCH_MEM_PER_CPU='999G',OTHER_EXAMPLE='keep')):
            with patch('bin.ref_depth_database.subprocess.run',return_value=SimpleNamespace(returncode=0,stdout='123',stderr='')) as called:
                cli.submit(database)
                environment=called.call_args.kwargs['env']
                self.assertNotIn('SBATCH_MEM_PER_NODE',environment);self.assertNotIn('SBATCH_MEM_PER_CPU',environment)
                self.assertEqual(environment['OTHER_EXAMPLE'],'keep')
        before=(database/'submission.json').read_bytes()
        cli.submit(database,preview=True)
        self.assertEqual((database/'submission.json').read_bytes(),before)

    def two_group_inputs(self):
        rows=[dict(self.rows[0]),dict(self.rows[1],GCFID=GCF2)]
        write_csv(self.records,['Sample','CramFilePath','GCFID'],rows)
        write_csv(self.refs,['GCF_ID','GCF_path'],[dict(GCF_ID=g,GCF_path=str(self.fai)) for g in (GCF1,GCF2)])

    def test_submission_cap_and_finished_slot_refill(self):
        self.two_group_inputs();database,plan=self.plan('cap',memory=None)
        with patch('bin.ref_depth_database.subprocess.run',return_value=SimpleNamespace(returncode=0,stdout='123',stderr='')) as called:
            records=cli.submit(database,jobs=1)
            self.assertEqual(len(records),1);self.assertEqual(called.call_count,1)
            command=called.call_args.args[0]
            self.assertEqual(command[command.index('-p')+1],'cpu')
            self.assertFalse(any(c.startswith('--mem') for c in command))
            cli.submit(database,jobs=1);self.assertEqual(called.call_count,1)
            self.assertEqual(run_gcf(database,GCF1,self.extractor),0)
            records=cli.submit(database,jobs=1)
            self.assertEqual(len(records),2)
            self.assertEqual(sum(c.args[0][0]=='sbatch' for c in called.call_args_list),2)
        limits=json.loads((database/'submission_limits.json').read_text())
        self.assertEqual((limits['tracked_active'],limits['submitted_this_call']),(1,1))
        with self.assertRaisesRegex(ValueError,'positive'):cli.submit(database,jobs=0)

    def test_failed_submit_burst_is_capped_and_not_automatically_retried(self):
        self.two_group_inputs();database,plan=self.plan('failed-burst')
        with patch('bin.ref_depth_database.subprocess.run',return_value=SimpleNamespace(returncode=1,stdout='',stderr='rejected')) as called:
            cli.submit(database,jobs=1);self.assertEqual(called.call_count,1)
            cli.submit(database,jobs=1);self.assertEqual(called.call_count,2)
            cli.submit(database,jobs=1);self.assertEqual(called.call_count,2)

    def test_run_slurm_refills_until_all_groups_finish(self):
        self.two_group_inputs();args=self.run_args('slurm');args.jobs=1
        actual_status=cli.status;observations=[]
        def poll(database,scheduler):
            records=json.loads((database/'submission.json').read_text())
            pending=[g for g,r in records.items() if not (database/g/'progress.json').exists()]
            self.assertLessEqual(len(pending),1)
            for gcf in pending:self.assertEqual(run_gcf(database,gcf,self.extractor),0)
            observations.append(len(records));return actual_status(database)
        with patch('bin.ref_depth_database.subprocess.run',return_value=SimpleNamespace(returncode=0,stdout='123',stderr='')), patch.object(cli,'status',side_effect=poll), patch.object(cli.time,'sleep'):
            result=cli.run(args)
        self.assertTrue(result['complete']);self.assertEqual(observations,[1,2])

    def test_status_refill_is_explicit(self):
        database,plan=self.plan('watch-refill')
        for refill in (False,True):
            arguments=['ref_depth_database.py','status','--database',str(database),'--watch','--scheduler','--jobs','1']
            if refill:arguments.append('--refill')
            with patch.object(sys,'argv',arguments),patch.object(cli,'submit') as submitted,patch.object(cli,'status',return_value=dict(terminal=True,rows=[dict(gcf=GCF1,state='COMPLETE',samples_finished=2)])),patch.object(cli,'collect',return_value=dict(complete=True)):
                self.assertEqual(cli.main(),0)
                if refill:submitted.assert_called_once_with(database,jobs=1)
                else:submitted.assert_not_called()

    def test_single_jobs_option_has_no_historical_alias(self):
        for command in ('submit','execute','status'):
            base=[command,'--database',str(self.root/'unused')]
            if command=='execute':base+=['--mode','slurm']
            self.assertEqual(cli.parser().parse_args(base+['--jobs','4']).jobs,4)
            with patch('sys.stderr'),self.assertRaises(SystemExit):
                cli.parser().parse_args(base+['--max-inflight','4'])
        args=self.run_args('slurm');args.jobs=4;args.dry_run=True
        self.assertEqual(cli.run(args)['concurrency_limit'],4)

    def test_explicit_memory_is_preserved(self):
        database,plan=self.plan('explicit-memory',memory='8G',partition='alternate')
        preview=cli.dry_run(database)
        self.assertEqual(preview['resources']['memory'],'8G')
        command=preview['jobs'][0]['command'];self.assertEqual(command[command.index('-p')+1],'alternate')
        self.assertIn('--mem=8G',preview['jobs'][0]['command'])
        self.assertIn('#SBATCH --mem=8G',(database/GCF1/'job.sbatch.sh').read_text())
        args=self.run_args('slurm');args.memory='default';args.dry_run=True
        self.assertIsNone(cli.run(args)['resources']['memory'])

    def test_one_command_local_prepares_builds_and_collects(self):
        args = self.run_args()
        def execute(database, jobs, collect_jobs=4):
            self.assertTrue((database / 'database.json').is_file())
            run_gcf(database, GCF1, self.extractor)
            return collect(database)
        with patch.object(cli, 'execute_local', side_effect=execute) as worker:
            result = cli.run(args)
        self.assertTrue(result['complete'])
        self.assertEqual(result['exit_code'], 0)
        worker.assert_called_once_with(args.outdir.resolve(), 1, collect_jobs=4)
        self.assertEqual(json.loads((args.outdir / 'pipeline.json').read_text())['state'], 'COMPLETE')
        self.assertTrue((args.outdir / 'NATIVE_INDEX.tsv').is_file())

    def test_one_command_slurm_waits_and_collects(self):
        args = self.run_args('slurm')
        actual_status = cli.status
        polls = []
        def observed(database, scheduler):
            polls.append(scheduler)
            if len(polls) == 1:
                return dict(actual_status(database), terminal=False, needs_attention=[], scheduler_returncode=0)
            run_gcf(database, GCF1, self.extractor)
            return actual_status(database)
        with patch.object(cli, 'submit') as submit, patch.object(cli, 'status', side_effect=observed), patch.object(cli.time, 'sleep') as sleep:
            result = cli.run(args)
        self.assertEqual(result['exit_code'], 0)
        self.assertTrue(result['complete'])
        self.assertEqual(submit.call_count,2)
        submit.assert_called_with(args.outdir.resolve(),jobs=1)
        sleep.assert_called_once_with(.01)
        self.assertEqual(polls, [True, True])

    def test_one_command_lost_job_is_not_success_or_infinite_wait(self):
        args = self.run_args('slurm')
        with patch.object(cli, 'submit'), patch.object(cli, 'status', return_value=dict(
                rows=[dict(gcf=GCF1,state='PENDING')], terminal=False, needs_attention=[GCF1], scheduler_returncode=0)), patch.object(cli.time, 'sleep'):
            result = cli.run(args)
        self.assertEqual(result['exit_code'], 1)
        self.assertFalse(result['complete'])
        self.assertEqual(json.loads((args.outdir / 'pipeline.json').read_text())['state'], 'NEEDS_ATTENTION')

    def test_one_command_failures_propagate_and_arguments_fail_before_prepare(self):
        args = self.run_args()
        args.interval = float('nan')
        with self.assertRaisesRegex(ValueError, 'interval'):
            cli.run(args)
        self.assertFalse(args.outdir.exists())
        args.interval = .01
        with patch.object(cli, 'execute_local', return_value=dict(complete=False, state='INCOMPLETE')):
            self.assertEqual(cli.run(args)['exit_code'], 1)
        self.assertEqual(json.loads((args.outdir / 'pipeline.json').read_text())['state'], 'FAILED')

    def test_manifest_counts_and_expected_slurm_resources(self):
        root, plan = self.plan()
        self.assertEqual(plan['counts'], dict(gcfs=1, samples=2, sample_gcf_pairs=2))
        script = root / GCF1 / 'job.sbatch.sh'
        directives = [line for line in script.read_text().splitlines() if line.startswith('#SBATCH')]
        self.assertEqual(directives, ['#SBATCH --cpus-per-task=2', '#SBATCH --mem=1G', '#SBATCH --partition=cpu'])
        subprocess.run(['bash', '-n', str(script)], check=True)
        result = collect(root)
        self.assertFalse(result['complete'])
        self.assertEqual(result['pending_samples'], 2)
        with (root / 'DEPTH_INDEX.tsv').open() as handle:
            self.assertEqual(len(list(csv.DictReader(handle, delimiter='\t'))), 6)

    def test_builder_rejects_removed_parameters_before_output(self):
        base=['run','--mode','local','--records',str(self.records),'--references',str(self.refs),'--outdir',str(self.root/'unused')]
        removed={'--depths':'1','--depth-grid':'1:2:1','--replicates':'5','--bins':'1000',
                 '--seed':'1','--batch-size':'16','--sampler':'auto','--matrix-format':'npz',
                 '--insufficient-depth':'record','--max-inflight':'4'}
        for flag,value in removed.items():
            with self.subTest(flag=flag),patch('sys.stderr'),self.assertRaises(SystemExit):
                cli.parser().parse_args(base+[flag,value])
        args=cli.parser().parse_args(base+['--dry-run'])
        with patch.object(cli.subprocess,'Popen') as spawned,patch.object(cli.subprocess,'run') as called:
            cli.run(args)
        spawned.assert_not_called();called.assert_not_called()
        config=json.loads((args.outdir/'database.json').read_text())['config']
        self.assertEqual(config['window_size'],1000)
        for key in ('depths','depth_grid','replicates','bins','seed','batch_size','sampler','matrix_format','insufficient_depth'):
            self.assertNotIn(key,config)

    def test_missing_duplicate_ambiguous_inputs_fail_before_output(self):
        write_csv(self.records, ['Sample', 'CramFilePath', 'GCFID'], self.rows + [self.rows[0]])
        with self.assertRaisesRegex(ValueError, 'Duplicate'):
            self.plan()
        write_csv(self.records, ['Sample', 'CramFilePath', 'GCFID'], self.rows)
        write_csv(self.refs, ['GCF_ID', 'GCF_path'], [dict(GCF_ID=GCF2, GCF_path=str(self.fai))])
        with self.assertRaisesRegex(ValueError, 'Missing FAI'):
            self.plan()
        write_csv(self.refs, ['GCF_ID', 'GCF_path'], [dict(GCF_ID=GCF1, GCF_path=str(self.fai))] * 2)
        with self.assertRaisesRegex(ValueError, 'Duplicate reference'):
            self.plan()
        self.assertFalse((self.root / 'db').exists())

    def test_success_nested_samples_depths_dense_and_catalog(self):
        root, plan = self.plan()
        self.assertEqual(run_gcf(root, GCF1, self.extractor), 0)
        summary = collect(root)
        self.assertTrue(summary['complete'])
        self.assertEqual(summary['state'], 'COMPLETE_WITH_GAPS')
        self.assertEqual(summary['unsupported_sample_depths'], 1)
        self.assertEqual(summary['ready_sample_depths'], 5)
        fingerprint = (root / 'PLAN.sha256').read_text().strip()
        for sample in plan['groups'][GCF1]['samples']:
            attempt, result = completed_sample(root / GCF1 / 'samples' / sample['folder'], fingerprint)
            pool = self.small if sample['sample'] == 'a' else self.big
            for row in result['depths']:
                if row['status'] != 'READY':
                    with self.assertRaises(ValueError):
                        load_depth(attempt / 'result.json', row['depth'])
                    continue
                matrix, edges, meta = load_depth(attempt / 'result.json', row['depth'])
                for archive in (attempt / row['matrix']).parent.glob('ids-*.npz'):
                    with np.load(archive, allow_pickle=False) as ids:
                        for replicate in ids:
                            dense = dense_coverage(pool, ids[replicate])
                            expected = [dense[a:b].sum() for a, b in zip(edges[:-1], edges[1:])]
                            np.testing.assert_array_equal(matrix[int(replicate)], expected)
                self.assertEqual(meta['sample']['sample'], sample['sample'])
                # Existing downstream statistics consume the stored reference rows directly.
                self.assertEqual(evaluate_windows(matrix[0], matrix)['simulations'], 5)
            with self.assertRaises(ValueError):
                load_depth(attempt / 'result.json', .2)
        with (root / 'READY_DEPTH_INDEX.tsv').open() as handle:
            ready = list(csv.DictReader(handle, delimiter='\t'))
        self.assertEqual(len(ready), 5)
        self.assertTrue(all(Path(row['matrix_path']).is_file() for row in ready))

    def test_npy_and_npz_have_identical_counts(self):
        arrays = []
        for codec in ('npy', 'npz'):
            root, plan = self.plan(name=codec, matrix_format=codec)
            run_gcf(root, GCF1, self.extractor)
            attempt = root / GCF1 / 'samples' / plan['groups'][GCF1]['samples'][0]['folder'] / 'attempt-000001'
            arrays.append(np.array(load_depth(attempt / 'result.json', 1)[0]))
        np.testing.assert_array_equal(*arrays)

    def test_resume_and_corruption_rebuild_new_attempt(self):
        root, plan = self.plan()
        run_gcf(root, GCF1, self.extractor)
        with patch('bin.ref_depth_worker.build_sample', side_effect=AssertionError('should reuse')):
            self.assertEqual(run_gcf(root, GCF1, self.extractor), 0)
        parent = root / GCF1 / 'samples' / plan['groups'][GCF1]['samples'][0]['folder']
        attempt = parent / 'attempt-000001'
        result = json.loads((attempt / 'result.json').read_text())
        payload = attempt / next(row['matrix'] for row in result['depths'] if row['depth'] == .1)
        payload.write_bytes(b'corrupt')
        self.assertFalse(collect(root)['complete'])
        with self.assertRaisesRegex(ValueError, 'checksum'):
            load_depth(attempt / 'result.json', .1)
        self.assertEqual(run_gcf(root, GCF1, self.extractor), 0)
        self.assertTrue((parent / 'attempt-000002/status/SUCCESS').exists())
        self.assertEqual(payload.read_bytes(), b'corrupt')
        self.assertTrue(collect(root)['complete'])

    def test_failure_continues_other_samples_and_retry(self):
        root, plan = self.plan()
        def failing(path, **kwargs):
            if Path(path) == self.a:
                raise RuntimeError('controlled extraction failure')
            return self.extractor(path, **kwargs)
        self.assertEqual(run_gcf(root, GCF1, failing), 1)
        summary = collect(root)
        self.assertEqual(summary['failed_samples'], 1)
        self.assertEqual(summary['ready_sample_depths'], 3)
        self.assertEqual(run_gcf(root, GCF1, self.extractor), 0)
        self.assertTrue(collect(root)['complete'])

    def test_source_and_fai_mutation_block_resume(self):
        root, _ = self.plan()
        self.a.write_bytes(b'modified original')
        self.assertEqual(run_gcf(root, GCF1, self.extractor), 1)
        self.assertEqual(collect(root)['failed_samples'], 1)
        self.fai.write_text('one\t21\ntwo\t11\n')
        self.assertEqual(run_gcf(root, GCF1, self.extractor), 1)
        self.assertEqual(collect(root)['failed_samples'], 2)

    def test_error_policy_and_empty_pool(self):
        root, _ = self.plan(insufficient_depth='error')
        self.assertEqual(run_gcf(root, GCF1, self.extractor), 1)
        self.assertEqual(collect(root)['failed_samples'], 1)
        root, _ = self.plan(name='empty')
        self.assertEqual(run_gcf(root, GCF1, lambda *a, **k: (BlockPool.from_intervals([], 31), self.layout)), 0)
        summary = collect(root)
        self.assertEqual(summary['ready_sample_depths'], 0)
        self.assertEqual(summary['unsupported_sample_depths'], 6)

    def test_gcf_layout_conflict_not_pooled(self):
        root, _ = self.plan()
        def extract(path, **kwargs):
            pool, layout = self.extractor(path, **kwargs)
            if Path(path) == self.b:
                layout = Reference(layout.names, layout.lengths, ('different', 'bbb'))
            return pool, layout
        self.assertEqual(run_gcf(root, GCF1, extract), 1)
        self.assertEqual(collect(root)['failed_samples'], 1)

    def test_plan_mutation_lock_and_safe_ids(self):
        root, _ = self.plan()
        with exclusive(root / GCF1 / '.worker.lock'):
            with self.assertRaisesRegex(RuntimeError, 'Locked'):
                run_gcf(root, GCF1, self.extractor)
        self.assertFalse((root / GCF1 / '.worker.lock').exists())
        self.assertNotIn('/', sample_folder('../odd sample'))
        (root / 'database.json').write_text('{}')
        with self.assertRaisesRegex(ValueError, 'checksum'):
            load_plan(root)

    def test_submit_is_explicit_and_queue_absence_is_not_success(self):
        root, _ = self.plan()
        with patch('bin.ref_depth_database.subprocess.run', return_value=SimpleNamespace(returncode=0, stdout='123;cluster\n', stderr='')) as mocked:
            records = cli.submit(root)
            self.assertEqual(records[GCF1]['job_id'], '123')
            self.assertEqual(mocked.call_args.args[0][0:2], ['sbatch', '--parsable'])
            cli.submit(root)
            self.assertEqual(mocked.call_count, 1)
        with patch('bin.ref_depth_database.subprocess.run', return_value=SimpleNamespace(returncode=0, stdout='', stderr='')):
            observed = cli.status(root, scheduler=True)
        self.assertFalse(observed['terminal'])
        self.assertIn('do not establish', observed['scheduler_note'])
        with self.assertRaisesRegex(ValueError, 'Slurm submission'):
            cli.execute_local(root, jobs=1)

    def test_ambiguous_submission_never_retried_automatically(self):
        root, _ = self.plan()
        with patch('bin.ref_depth_database.subprocess.run', side_effect=subprocess.TimeoutExpired('sbatch', 60)):
            with self.assertRaises(subprocess.TimeoutExpired):
                cli.submit(root)
        with self.assertRaisesRegex(ValueError, 'Unknown submission'):
            cli.submit(root, retry=True)

    def test_submit_failure_is_recorded_and_can_retry(self):
        root, _ = self.plan()
        with patch('bin.ref_depth_database.subprocess.run', return_value=SimpleNamespace(returncode=1, stdout='', stderr='no permission')):
            self.assertEqual(cli.submit(root)[GCF1]['state'], 'SUBMIT_FAILED')
        self.assertTrue(cli.status(root)['terminal'])
        self.assertGreater(collect(root)['failed_samples'], 0)
        with patch('bin.ref_depth_database.subprocess.run', return_value=SimpleNamespace(returncode=0, stdout='234\n', stderr='')):
            record = cli.submit(root, retry=True)[GCF1]
            self.assertEqual(record['job_id'], '234')
            self.assertEqual(len(record['history']), 1)


try:
    import pysam
except ImportError:
    pysam = None


@unittest.skipIf(pysam is None, 'Synthetic CRAM/BAM integration requires pysam on NH')
class ExecutionIntegrationTests(unittest.TestCase):
    def test_real_synthetic_cram_crai_pair_is_preserved(self):
        with tempfile.TemporaryDirectory(prefix='indexed-cram-build-') as temporary:
            root=Path(temporary);fasta=root/'reference.fa';fasta.write_text('>synthetic\n'+'A'*100+'\n')
            pysam.faidx(str(fasta))
            cram=root/'reference.cram'
            with pysam.AlignmentFile(str(cram),'wc',reference_filename=str(fasta),
                    header=dict(HD=dict(VN='1.6',SO='coordinate'),SQ=[dict(SN='synthetic',LN=100)])) as handle:
                for i in range(20):
                    read=pysam.AlignedSegment();read.query_name=f'read-{i}';read.query_sequence='A'*10
                    read.flag=0;read.reference_id=0;read.reference_start=i;read.mapping_quality=40;read.cigar=[(0,10)]
                    handle.write(read)
            pysam.index(str(cram));crai=Path(str(cram)+'.crai');before=(digest(cram),digest(crai))
            records,refs=root/'records.csv',root/'refs.csv'
            write_csv(records,['Sample','CramFilePath','GCFID'],[dict(Sample='sample',CramFilePath=str(cram),GCFID=GCF1)])
            write_csv(refs,['GCF_ID','GCF_path'],[dict(GCF_ID=GCF1,GCF_path=str(fasta)+'.fai')])
            database=root/'database';entry=Path(__file__).resolve().parents[1]/'bin/ref_depth_database.py'
            subprocess.run([sys.executable,'-B',str(entry),'run','--mode','local','--records',str(records),
                '--references',str(refs),'--outdir',str(database),'--window-size','10',
                '--cpus','2','--geometry-format','packed'],check=True,stdout=subprocess.DEVNULL)
            self.assertTrue(json.loads((database/'SUMMARY.json').read_text())['complete'])
            archive=next(database.glob('*/samples/*/attempt-*/source.cram.crai'))
            self.assertEqual(digest(archive),before[1]);self.assertEqual((digest(cram),digest(crai)),before)
            # The archived index is an actual functional CRAI against this unchanged CRAM.
            with pysam.AlignmentFile(str(cram),'rc',reference_filename=str(fasta),index_filename=str(archive),require_index=True) as handle:
                self.assertEqual(sum(1 for _ in handle.fetch('synthetic',0,100)),20)
            crai.unlink()
            without=root/'database-without-crai'
            subprocess.run([sys.executable,'-B',str(entry),'run','--mode','local','--records',str(records),
                '--references',str(refs),'--outdir',str(without),'--window-size','10',
                '--cpus','2','--geometry-format','packed'],check=True,stdout=subprocess.DEVNULL)
            self.assertTrue(json.loads((without/'SUMMARY.json').read_text())['complete'])
            self.assertFalse(list(without.glob('*/samples/*/attempt-*/source.cram.crai')))
            self.assertEqual(digest(cram),before[0])
            with (without/'SAMPLE_INDEX.tsv').open() as handle:
                row=next(csv.DictReader(handle,delimiter='\t'))
            sample=json.loads(Path(row['sample_manifest']).read_text())
            self.assertIsNone(sample['cram_index'])
            original_manifest=next(database.glob('*/samples/*/attempt-*/result.json'))
            with np.load(original_manifest.parent/'native_coverage.npz') as a,np.load(Path(row['sample_manifest']).parent/'native_coverage.npz') as b:
                for key in a.files:np.testing.assert_array_equal(a[key],b[key])



    def test_local_and_generated_slurm_job_have_identical_science(self):
        with tempfile.TemporaryDirectory(prefix='depth-db-fixture-') as temporary:
            root = Path(temporary)
            fai = root / 'synthetic.fai'
            fai.write_text('synthetic\t100\n')
            bam = root / 'synthetic.bam'
            with pysam.AlignmentFile(str(bam), 'wb', header=dict(HD=dict(VN='1.6'), SQ=[dict(SN='synthetic', LN=100)])) as handle:
                for i in range(80):
                    read = pysam.AlignedSegment()
                    read.query_name = f'toy-{i}'; read.query_sequence = 'A' * 10
                    read.flag = 0; read.reference_id = 0; read.reference_start = (i * 7) % 90
                    read.mapping_quality = 40; read.cigar = [(0, 10)]
                    handle.write(read)
            records, refs = root / 'records.csv', root / 'refs.csv'
            write_csv(records, ['Sample', 'CramFilePath', 'GCFID'], [dict(Sample=s, CramFilePath=str(bam), GCFID=g) for g, s in ((GCF1, 'one'), (GCF1, 'two'), (GCF2, 'three'))])
            write_csv(refs, ['GCF_ID', 'GCF_path'], [dict(GCF_ID=g, GCF_path=str(fai)) for g in (GCF1, GCF2)])
            config = dict(storage_mode='native',window_size=10,cpus=2,memory='1G',python=sys.executable,geometry_format='packed',
                policies=dict(mate_policy='both',sampling_unit='read',min_mapq=0,exclude_duplicates=False,exclude_qcfail=False))
            outputs = []
            for mode in ('local', 'slurm-script'):
                database = root / mode
                if mode == 'local':
                    entry = Path(__file__).resolve().parents[1] / 'bin/ref_depth_database.py'
                    subprocess.run([sys.executable, '-B', str(entry), 'run', '--mode', 'local', '--jobs', '2',
                        '--records', str(records), '--references', str(refs), '--outdir', str(database),
                        '--window-size', '10', '--cpus', '2',
                        '--memory', '1G', '--geometry-format', 'packed'], check=True,
                        stdout=subprocess.DEVNULL)
                    self.assertTrue(json.loads((database / 'SUMMARY.json').read_text())['complete'])
                else:
                    create_plan(records, refs, database, config, COLUMNS)
                    for script in sorted(database.glob('*/job.sbatch.sh')):
                        subprocess.run(['bash', str(script)], check=True)
                    self.assertTrue(collect(database)['complete'])
                with (database / 'NATIVE_INDEX.tsv').open() as handle:
                    rows = list(csv.DictReader(handle, delimiter='\t'))
                arrays = {}
                for row in rows:
                    with np.load(Path(row['sample_manifest']).parent/'native_coverage.npz') as data:
                        arrays[(row['gcf'],row['sample'])]=data['aligned_bases']
                self.assertEqual(len(arrays),3)
                outputs.append(arrays)
                if mode == 'local':
                    from bin.reference_index import build_index
                    index=root/'reference.sqlite';build_index(database,index,True)
                    query=root/'query'
                    subprocess.run([sys.executable,'-B',str(Path(__file__).resolve().parents[1]/'bin/coverage_concordance.py'),
                        '--comparison-policy','nearest','--reference-sqlite',str(index),'--test-alignment',str(bam),'--query-sample','query','--gcf',GCF1,
                        '--include-self','--depths','0.1,1','--replicates','5','--outdir',str(query)],check=True)
                    report=json.loads((query/'provenance.json').read_text())
                    self.assertEqual(report['failed_query_count'],0)
                    self.assertEqual(report['ready_query_count'],1)
            self.assertEqual(outputs[0].keys(), outputs[1].keys())
            for key in outputs[0]:
                np.testing.assert_array_equal(outputs[0][key], outputs[1][key])


if __name__ == '__main__':
    unittest.main()
