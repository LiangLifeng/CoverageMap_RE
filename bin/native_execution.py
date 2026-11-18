"""Bounded persistent native sample workers; one controller owns GCF progress."""
from concurrent.futures import ProcessPoolExecutor, wait, FIRST_COMPLETED
from contextlib import contextmanager
import csv
import json
import multiprocessing
import os
from pathlib import Path
import signal
import sys
import time
import traceback

from bin.ref_depth_catalog import atomic_json, load_plan

ROOT = PLAN = SAMPLES = None


def initialize(database):
    global ROOT, PLAN, SAMPLES
    for name in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS'):
        os.environ[name]='1'
    def interrupted(signum, frame):
        raise InterruptedError('Sample worker interrupted')
    signal.signal(signal.SIGTERM,interrupted)
    ROOT,PLAN=load_plan(database,check_code=True)
    SAMPLES={(g,s['sample']):s for g,group in PLAN['groups'].items() for s in group['samples']}
    from bin import ref_depth_worker  # Import numpy/pysam once, after thread limits.


@contextmanager
def sample_log(path):
    """Capture Python and native HTSlib output without sharing a GCF log."""
    sys.stdout.flush();sys.stderr.flush()
    saved=[os.dup(1),os.dup(2)]
    try:
        with path.open('ab',buffering=0) as handle:
            os.dup2(handle.fileno(),1);os.dup2(handle.fileno(),2)
            yield
    finally:
        sys.stdout.flush();sys.stderr.flush()
        for descriptor,original in zip((1,2),saved):
            os.dup2(original,descriptor);os.close(original)


def perform(gcf, sample):
    from bin.ref_depth_worker import run_sample
    record=SAMPLES[gcf,sample]
    path=ROOT/gcf/'samples'/record['folder']/'worker.log'
    with sample_log(path):
        try:
            result=run_sample(ROOT,PLAN,gcf,record)
        except Exception as error:
            traceback.print_exc()
            result=dict(gcf=gcf,sample=sample,state='FAILED',exit_code=1,
                        error=f'{type(error).__name__}: {error}',worker_pid=os.getpid())
    result['log']=str(path)
    return result


def sample_tasks(plan):
    tasks=[(g,s) for g,group in plan['groups'].items() for s in group['samples']]
    return sorted(tasks,key=lambda item:(-item[1]['source_stat']['bytes'],item[0],item[1]['sample']))


def execute(database, plan, jobs, progress, selected=None):
    tasks=sample_tasks(plan)
    states={g:dict(state='PENDING',samples_finished=0,failed_samples=0,samples_reused=0,
                   samples_total=len(group['samples']),failures=[],active_samples=[])
            for g,group in plan['groups'].items()}
    started=time.monotonic();outcomes={};group_exit={};pids=set()
    if selected is not None:
        selected=set(selected)
        known={(g,s['sample']) for g,s in tasks}
        if not selected or not selected <= known:
            raise ValueError('Retry sample selection is empty or outside plan')
        with (database/'SAMPLE_INDEX.tsv').open() as handle:
            rows=list(csv.DictReader(handle,delimiter='\t'))
        if len(rows)!=len(known) or {(r['gcf'],r['sample']) for r in rows}!=known:
            raise ValueError('Retry catalog does not cover the immutable plan')
        for row in rows:
            key=(row['gcf'],row['sample'])
            if key in selected and row['state']!='FAILED':
                raise ValueError('Retry selection contains a sample without final failure')
            if key not in selected:
                if row['state']!='COMPLETE':
                    raise ValueError('Retry would omit an unfinished sample')
                states[row['gcf']]['samples_finished']+=1
        outcomes=json.loads((database/'local_sample_exit_codes.json').read_text())
        tasks=[(g,s) for g,s in tasks if (g,s['sample']) in selected]
    def publish(gcf):
        row=states[gcf];processed=row['samples_finished']+row['failed_samples']
        row['state']=('FAILED' if row['failed_samples'] else 'COMPLETE') if processed==row['samples_total'] else 'RUNNING' if row['active_samples'] or processed else 'PENDING'
        row['current_sample']=row['active_samples'][0] if row['active_samples'] else ''
        atomic_json(database/gcf/'progress.json',dict(row,updated=time.time(),elapsed_seconds=time.monotonic()-started))
        if processed==row['samples_total']:
            group_exit[gcf]=int(bool(row['failed_samples']))
            atomic_json(database/'local_exit_codes.json',group_exit)
    for gcf in states:publish(gcf)
    baseline={p.pid for p in multiprocessing.active_children()}
    executor=ProcessPoolExecutor(max_workers=min(jobs,len(tasks)),mp_context=multiprocessing.get_context('spawn'),
                                 initializer=initialize,initargs=(str(database),))
    pending={};iterator=iter(tasks);exhausted=False
    progress.event(f'LOCAL SAMPLE START jobs={jobs} tasks={len(tasks)} persistent_workers; logs: <database>/<GCF>/samples/<sample-folder>/worker.log')
    progress.update(force=True)
    try:
        while pending or not exhausted:
            while not exhausted and len(pending)<jobs:
                try:gcf,sample=next(iterator)
                except StopIteration:
                    exhausted=True;break
                future=executor.submit(perform,gcf,sample['sample'])
                pending[future]=(gcf,sample['sample'])
                states[gcf]['active_samples'].append(sample['sample']);publish(gcf)
            if pending:
                done,_=wait(pending,timeout=.1,return_when=FIRST_COMPLETED)
                for future in done:
                    gcf,sample=pending.pop(future);result=future.result()
                    if (result['gcf'],result['sample'])!=(gcf,sample):
                        raise RuntimeError('Sample worker returned wrong identity')
                    if result['state']=='INTERRUPTED':
                        raise InterruptedError(f'Sample worker interrupted: {gcf}/{sample}')
                    states[gcf]['active_samples'].remove(sample)
                    if result['exit_code']:
                        states[gcf]['failed_samples']+=1
                        states[gcf]['failures'].append(dict(sample=sample,error=result.get('error',result['state'])))
                        progress.event(f"SAMPLE FAILED {gcf}/{sample} log={result['log']}")
                    else:
                        states[gcf]['samples_finished']+=1
                        states[gcf]['samples_reused']+=int(result['reused'])
                    pids.add(result['worker_pid']);outcomes[gcf+'/'+sample]=result['exit_code'];publish(gcf)
            progress.update()
    except BaseException:
        children=[p for p in multiprocessing.active_children() if p.pid not in baseline]
        for child in children:child.terminate()
        for child in children:
            child.join(timeout=10)
            if child.is_alive():child.kill();child.join()
        for gcf,row in states.items():
            if row['active_samples']:
                atomic_json(database/gcf/'progress.json',dict(row,state='INTERRUPTED',updated=time.time()))
        raise
    finally:
        executor.shutdown(wait=True,cancel_futures=True)
        progress.update(force=True)
    atomic_json(database/'local_sample_exit_codes.json',outcomes)
    return dict(local_exit_codes=group_exit,local_scheduling='sample',sample_task_count=len(tasks),
                worker_pids=sorted(pids),requested_cpu_budget=jobs*plan['config']['cpus'])
