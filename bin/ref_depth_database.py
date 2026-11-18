#!/usr/bin/env python3
"""Build and collect original GCF/sample references with virtual GCF windows.
CRAM-only inputs are supported; optional CRAI and required output payloads are verified."""
import argparse
import csv
import json
import math
import os
from pathlib import Path
import re
import signal
import subprocess
import sys
import time

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from bin.ref_depth_catalog import (atomic_json, collect, create_plan, database_status,
    exclusive, generate_scripts, load_plan, worker_command)


def progress_counts(plan, rows):
    """Read worker counters only; never scan/rehash sample payloads for display."""
    states = ('COMPLETE', 'FAILED', 'INTERRUPTED', 'SUBMIT_FAILED')
    done = sum(int(r.get('samples_finished', 0)) for r in rows)
    failed = sum(int(r.get('failed_samples', 0)) for r in rows)
    total = sum(len(g['samples']) for g in plan['groups'].values())
    return dict(gcfs_total=len(plan['groups']), gcfs_complete=sum(r['state']=='COMPLETE' for r in rows),
        gcfs_terminal=sum(r['state'] in states for r in rows),
        gcfs_failed=sum(r['state'] in states[1:] for r in rows),
        gcfs_running=sum(r['state']=='RUNNING' for r in rows),
        gcfs_pending=sum(r['state']=='PENDING' for r in rows),
        samples_total=total, samples_complete=done, samples_failed=failed,
        samples_processed=done+failed,
        samples_reused=sum(int(r.get('samples_reused', 0)) for r in rows),
        current_samples=[r['gcf']+'/'+r['current_sample'] for r in rows
                         if r['state']=='RUNNING' and r.get('current_sample')])


class ProgressLog:
    """Line-oriented, flushed progress for terminals, nohup and Slurm logs."""
    def __init__(self, database, plan, persist=True, interval=10):
        self.database, self.plan = Path(database), plan
        self.path = self.database/'build.log' if persist else None
        self.started = time.monotonic()
        self.last = None
        self.interval = interval

    def event(self, message):
        line = time.strftime('[%Y-%m-%d %H:%M:%S] ') + message
        print(line, file=sys.stderr, flush=True)
        if self.path:
            with self.path.open('a') as handle:
                handle.write(line+'\n')

    def update(self, rows=None, force=False):
        now = time.monotonic()
        if not force and self.last is not None and now-self.last < self.interval:
            return
        self.last = now
        p = progress_counts(self.plan, rows if rows is not None else database_status(self.database, self.plan))
        fraction = p['samples_processed']/p['samples_total'] if p['samples_total'] else 0
        filled = min(20, int(20*fraction))
        self.event(f"[{'#'*filled}{'-'*(20-filled)}] {fraction:6.1%} "
            f"samples={p['samples_processed']}/{p['samples_total']} "
            f"ok={p['samples_complete']} failed={p['samples_failed']} reused={p['samples_reused']} | "
            f"GCF={p['gcfs_terminal']}/{p['gcfs_total']} ok={p['gcfs_complete']} "
            f"failed={p['gcfs_failed']} running={p['gcfs_running']} pending={p['gcfs_pending']} | "
            f"elapsed={now-self.started:.1f}s" +
            (' current='+','.join(p['current_samples'][:5]) if p['current_samples'] else ''))


def collection_progress(log):
    started=time.monotonic();last=None;signature=None
    def report(done,total,failed):
        nonlocal last,signature
        now=time.monotonic()
        if signature==(done,total,failed):return
        if done not in (0,total) and last is not None and now-last<10:
            return
        last=now
        signature=(done,total,failed)
        fraction=done/total if total else 0
        filled=min(20,int(20*fraction))
        log.event(f"COLLECT [{'#'*filled}{'-'*(20-filled)}] {fraction:6.1%} verified={done}/{total} failed={failed} elapsed={now-started:.1f}s")
    return report


def prepare(args):
    from bin.gcf_windows import GCF_CONCAT

    bins = getattr(args, 'bin_count', None)
    size = args.window_size if args.window_size is not None or bins is not None else 1000
    if args.cpus < 1 or (size if bins is None else bins) < 1 or args.min_mapq < 0 or args.collect_jobs < 1:
        raise ValueError('CPU/window-size/collect-jobs must be positive; min-mapq must be nonnegative')
    memory = args.memory.upper() if args.memory and args.memory.lower() != 'default' else None
    if memory and not re.fullmatch(r'[1-9]\d*[KMGTP]?', memory):
        raise ValueError('Memory must be positive, for example 4G')
    partition = getattr(args, 'partition', None)
    if partition and not re.fullmatch(r'[A-Za-z0-9_.-]+(?:,[A-Za-z0-9_.-]+)*', partition):
        raise ValueError('Invalid Slurm partition')
    python = Path(args.python).resolve()
    if not python.is_file() or not os.access(python, os.X_OK):
        raise ValueError('Python executable is not accessible')
    if args.sampling_unit == 'fragment' and args.mate_policy != 'both':
        raise ValueError('Fragment sampling requires both mates')
    config = dict(storage_mode='native', window_size=size, bin_count=bins,
        window_policy=GCF_CONCAT, cpus=args.cpus,
        memory=memory, partition=partition, python=str(python), geometry_format=args.geometry_format,
        policies=dict(mate_policy=args.mate_policy, sampling_unit=args.sampling_unit,
                      min_mapq=args.min_mapq, exclude_duplicates=args.exclude_duplicates,
                      exclude_qcfail=args.exclude_qcfail))
    if getattr(args, 'alignment_capture', 'legacy') == 'complete':
        if getattr(args, 'sampling_unit', None) not in (None, 'read'):
            raise ValueError('Complete geometry chooses sampling unit at query time')
        config['capture_alignment'] = True
        config['policies'].pop('sampling_unit', None)
        config['policies']['coverage_alignment_policy'] = 'primary_only'
    reuse = getattr(args, 'reuse_geometry_from', None)
    if reuse is not None:
        source, plan = load_plan(reuse)
        if plan['config']['policies'] != config['policies']:
            raise ValueError('Reusable geometry read policies differ from the new plan')
        config['reuse_geometry'] = dict(database=str(source),
                                        plan_sha256=(source/'PLAN.sha256').read_text().strip())
    columns = {key: getattr(args, key + '_column') for key in ('sample', 'cram', 'gcf', 'ref_gcf', 'fai')}
    plan = create_plan(args.records, args.references, args.outdir, config, columns)
    collect(args.outdir,jobs=args.collect_jobs)
    return dict(database=str(args.outdir.resolve()), counts=plan['counts'], state='PLANNED',
                window_size=size, bin_count=bins, window_policy=GCF_CONCAT,
                note='No worker or Slurm job has been started')


def execute_local(database, jobs, collect_jobs=4):
    if jobs < 1 or collect_jobs < 1:
        raise ValueError('jobs must be positive')
    database, plan = load_plan(database, check_code=True)
    progress = ProgressLog(database, plan)
    with exclusive(database / '.execution.lock'):
        submission = database / 'submission.json'
        if submission.exists():
            records = json.loads(submission.read_text())
            states = {r['gcf']: r['state'] for r in database_status(database, plan)}
            if any(r['state'] in ('SUBMITTED', 'UNKNOWN_SUBMISSION') and not states[gcf].startswith('COMPLETE')
                   for gcf, r in records.items()):
                raise ValueError('Existing Slurm submission; inspect jobs before using local mode')
        if plan['config'].get('storage_mode')=='native':
            from bin.native_execution import execute
            execution=execute(database,plan,jobs,progress)
            progress.event(f'COLLECTING: verifying payloads with collect-jobs={collect_jobs}')
            summary=collect(database,jobs=collect_jobs,progress=collection_progress(progress))
            summary,execution=finish_native_retries(database,plan,jobs,collect_jobs,progress,summary,execution)
            summary.update(execution)
            progress.event(f"LOCAL END state={summary['state']} complete={summary['complete']}")
            return summary
        queue, active, outcomes = iter(plan['groups']), {}, {}
        exhausted = False
        try:
            progress.event(f'LOCAL START jobs={jobs} database={database}; worker logs: <database>/<GCF>/worker.log')
            progress.update(force=True)
            while active or not exhausted:
                while not exhausted and len(active) < jobs:
                    try:
                        gcf = next(queue)
                    except StopIteration:
                        exhausted = True
                        break
                    log = (database / gcf / 'worker.log').open('ab')
                    child = subprocess.Popen(worker_command(database, gcf, plan), cwd=database / gcf,
                                             stdout=log, stderr=log)
                    active[child] = (gcf, log)
                for child in list(active):
                    status = child.poll()
                    if status is not None:
                        gcf, log = active.pop(child)
                        log.close(); outcomes[gcf] = status
                        atomic_json(database / 'local_exit_codes.json', outcomes)
                        if status:
                            progress.event(f'GCF FAILED {gcf} exit_code={status} log={database/gcf/"worker.log"}')
                progress.update()
                if active:
                    time.sleep(.1)
        finally:
            for child, (_, log) in active.items():
                child.terminate()
                try:
                    child.wait(timeout=10)
                except subprocess.TimeoutExpired:
                    child.kill(); child.wait()
                log.close()
            progress.update(force=True)
        progress.event('COLLECTING: verifying payloads and writing final catalogs')
        summary = collect(database,jobs=collect_jobs,progress=collection_progress(progress))
        progress.event(f"LOCAL END state={summary['state']} complete={summary['complete']}")
        return dict(summary, local_exit_codes=outcomes, requested_cpu_budget=jobs * plan['config']['cpus'])


def failed_sample_keys(database, plan):
    """Use final verification results, including failures missed by worker exits."""
    with (database/'SAMPLE_INDEX.tsv').open() as handle:
        rows=list(csv.DictReader(handle,delimiter='\t'))
    known={(g,s['sample']) for g,group in plan['groups'].items() for s in group['samples']}
    if len(rows)!=len(known) or {(r['gcf'],r['sample']) for r in rows}!=known:
        raise ValueError('Final sample catalog does not cover the immutable plan')
    return sorted((r['gcf'],r['sample']) for r in rows if r['state']=='FAILED')


def finish_native_retries(database, plan, jobs, collect_jobs, progress, summary, execution):
    """Two additional failed-sample rounds; keep successful samples untouched."""
    from bin.native_execution import execute
    path=database/'sample_retries.json'
    history=json.loads(path.read_text()) if path.exists() else dict(invocations=[])
    invocation=dict(started=time.time(),max_additional_attempts=2,rounds=[])
    history['invocations'].append(invocation)
    atomic_json(path,history)
    pids=set(execution['worker_pids'])
    retry_tasks=0
    for number in (1,2):
        if not summary['failed_samples']:
            break
        selected=failed_sample_keys(database,plan)
        if len(selected)!=summary['failed_samples']:
            raise ValueError('Final failed-sample count differs from catalog')
        record=dict(round=number,state='RUNNING',samples=[dict(gcf=g,sample=s) for g,s in selected])
        invocation['rounds'].append(record)
        atomic_json(path,history)
        progress.event(f'FINAL SAMPLE RETRY round={number}/2 failed_samples={len(selected)}')
        try:
            observed=execute(database,plan,jobs,progress,selected=selected)
            pids.update(observed['worker_pids']);retry_tasks+=len(selected)
            execution['local_exit_codes']=observed['local_exit_codes']
            summary=collect(database,jobs=collect_jobs,progress=collection_progress(progress))
        except BaseException as error:
            record.update(state='INTERRUPTED' if isinstance(error,(KeyboardInterrupt,InterruptedError)) else 'FAILED',
                          error=f'{type(error).__name__}: {error}',ended=time.time())
            atomic_json(path,history)
            raise
        record.update(state='FINISHED',ended=time.time(),failed_samples=summary['failed_samples'],
                      ready_samples=summary['ready_samples'],complete=summary['complete'])
        atomic_json(path,history)
    invocation.update(ended=time.time(),complete=summary['complete'],remaining_failed_samples=summary['failed_samples'])
    atomic_json(path,history)
    execution.update(worker_pids=sorted(pids),sample_retry_rounds=len(invocation['rounds']),
                     sample_retry_tasks=retry_tasks,max_additional_sample_attempts=2,
                     sample_retry_audit=str(path))
    return summary,execution


def submission_command(database, gcf, plan):
    config = plan['config']
    command = ['sbatch', '--parsable', f"--cpus-per-task={config['cpus']}"]
    command.extend(['-p', config.get('partition') or 'cpu'])
    if config.get('memory'):
        command.append('--mem=' + config['memory'])
    return command + [str(database / gcf / 'job.sbatch.sh')]


def dry_run(database, mode='slurm', jobs=1):
    if jobs < 1:
        raise ValueError('jobs must be positive')
    database, plan = load_plan(database, check_code=True)
    value = dict(state='DRY_RUN', submitted=False, workers_started=False, exit_code=0,
        mode=mode, database=str(database), resources=dict(cpus=plan['config']['cpus'], memory=plan['config'].get('memory'), partition=plan['config'].get('partition') or 'cpu'),
        concurrency_limit=jobs, first_round_max_submissions=min(jobs,len(plan['groups'])),
        scheduler_checked=False,
        jobs=[dict(gcf=gcf, script=str(database/gcf/'job.sbatch.sh'),
                   command=submission_command(database,gcf,plan) if mode=='slurm' else worker_command(database,gcf,plan))
              for gcf in plan['groups']])
    if mode=='local' and plan['config'].get('storage_mode')=='native':
        value.update(local_scheduling='sample',sample_task_count=plan['counts']['sample_gcf_pairs'],
            max_additional_sample_attempts=2,
            local_worker_limit=jobs,first_round_max_submissions=min(jobs,plan['counts']['sample_gcf_pairs']),
            local_note='Local mode uses persistent sample workers; GCF commands above are standalone worker/script references')
    atomic_json(database/'dry_run.json', value)
    return value


def submit(database, retry=False, preview=False, jobs=1):
    if jobs < 1:
        raise ValueError('jobs must be positive')
    if preview:
        return dry_run(database, jobs=jobs)
    database, plan = load_plan(database, check_code=True)
    with exclusive(database / '.execution.lock'):
        path = database / 'submission.json'
        records = json.loads(path.read_text()) if path.is_file() else {}
        states = {row['gcf']: row['state'] for row in database_status(database, plan)}
        terminal = {'COMPLETE', 'FAILED', 'INTERRUPTED', 'SUBMIT_FAILED'}
        unresolved = [gcf for gcf,r in records.items() if r['state']=='UNKNOWN_SUBMISSION' and states[gcf] not in terminal]
        if unresolved:
            raise ValueError('Unknown submission; reconcile with operator before continuing: '+','.join(unresolved))
        active = sum(r['state'] in ('SUBMITTED','UNKNOWN_SUBMISSION') and states[gcf] not in terminal for gcf,r in records.items())
        attempted = submitted_now = 0
        for gcf in plan['groups']:
            if active >= jobs or attempted >= jobs:
                break
            if states[gcf].startswith('COMPLETE'):
                continue
            previous = records.get(gcf)
            if previous:
                if previous['state'] == 'UNKNOWN_SUBMISSION':
                    raise ValueError(f'Unknown submission for {gcf}; reconcile with operator, do not duplicate it')
                if not retry or states[gcf] in ('RUNNING', 'PENDING') and previous['state'] == 'SUBMITTED':
                    continue
            attempted += 1
            command = submission_command(database, gcf, plan)
            # Persist uncertainty before contacting scheduler; a killed submitter cannot silently duplicate jobs.
            history = list(previous.get('history', [])) if previous else []
            if previous:
                history.append({k: v for k, v in previous.items() if k != 'history'})
            records[gcf] = dict(state='UNKNOWN_SUBMISSION', command=command, started=time.time(), history=history)
            atomic_json(path, records)
            try:
                environment = os.environ.copy()
                if not plan['config'].get('memory'):
                    for name in ('SBATCH_MEM_PER_NODE', 'SBATCH_MEM_PER_CPU', 'SBATCH_MEM_PER_GPU'):
                        environment.pop(name, None)
                result = subprocess.run(command, cwd=database / gcf, capture_output=True, text=True, timeout=60, env=environment)
            except OSError as error:
                records[gcf].update(state='SUBMIT_FAILED', error=str(error))
                atomic_json(path, records)
                continue
            except subprocess.TimeoutExpired as error:
                records[gcf]['error'] = str(error)
                atomic_json(path, records)
                raise
            if result.returncode:
                records[gcf].update(state='SUBMIT_FAILED', stderr=result.stderr, returncode=result.returncode)
                atomic_json(path, records)
                continue
            job_id = result.stdout.strip().split(';', 1)[0]
            if not re.fullmatch(r'\d+', job_id):
                records[gcf]['stdout'] = result.stdout
                atomic_json(path, records)
                raise ValueError('Ambiguous sbatch output; submission left UNKNOWN')
            records[gcf].update(state='SUBMITTED', job_id=job_id, stdout=result.stdout)
            active += 1; submitted_now += 1
            atomic_json(path, records)
        atomic_json(database/'submission_limits.json',dict(concurrency_limit=jobs, tracked_active=active, attempted_this_call=attempted, submitted_this_call=submitted_now, deferred_unsubmitted=sum(gcf not in records and not states[gcf].startswith('COMPLETE') for gcf in plan['groups'])))
        return records


def status(database, scheduler=False):
    database, plan = load_plan(database)
    rows = database_status(database, plan)
    result = dict(rows=rows, terminal=all(r['state'] in ('COMPLETE', 'FAILED', 'INTERRUPTED', 'SUBMIT_FAILED') for r in rows))
    result['progress'] = progress_counts(plan, rows)
    if scheduler:
        path = database / 'submission.json'
        records = json.loads(path.read_text()) if path.exists() else {}
        ids = [r['job_id'] for r in records.values() if r.get('job_id')]
        if ids:
            query = subprocess.run(['squeue', '--noheader', '--jobs', ','.join(ids), '--format', '%i|%T'],
                                   capture_output=True, text=True, timeout=30)
            result['scheduler_returncode'] = query.returncode
            result['scheduler_stdout'] = query.stdout
            result['scheduler_stderr'] = query.stderr
            queued = dict(line.split('|', 1) for line in query.stdout.splitlines() if '|' in line)
            result['scheduler_states'] = {gcf: queued.get(record['job_id'], 'NOT_IN_QUEUE')
                for gcf, record in records.items() if record.get('job_id')}
            result['needs_attention'] = [row['gcf'] for row in rows if row['state'] in ('PENDING', 'RUNNING')
                and result['scheduler_states'].get(row['gcf']) == 'NOT_IN_QUEUE'] if not query.returncode else []
            result['scheduler_note'] = 'Absent queue entries do not establish successful completion'
    return result


def run(args):
    """Build from manifests through final verified catalogs in one invocation."""
    if args.jobs < 1 or not math.isfinite(args.interval) or args.interval <= 0:
        raise ValueError('jobs and interval must be positive')
    print(time.strftime('[%Y-%m-%d %H:%M:%S] ')+
          'PREPARING: validating manifests, references and input identities', file=sys.stderr, flush=True)
    prepare(args)
    database = args.outdir.resolve()
    _, plan = load_plan(database)
    progress = ProgressLog(database, plan, persist=not getattr(args, 'dry_run', False))
    def phase(state, **details):
        value = dict(state=state, mode=args.mode, updated=time.time(), **details)
        atomic_json(database / 'pipeline.json', value)
        print(json.dumps(value, sort_keys=True), flush=True)
        progress.event('PHASE '+state)
    try:
        if getattr(args, 'dry_run', False):
            value = dry_run(database, args.mode, args.jobs)
            phase('DRY_RUN', exit_code=0, submitted=False, workers_started=False)
            return value
        if args.mode == 'local':
            phase('EXECUTING')
            summary = execute_local(database, args.jobs, collect_jobs=args.collect_jobs)
        else:
            phase('SUBMITTING')
            phase('MONITORING', concurrency_limit=args.jobs)
            absent = {}
            while True:
                submit(database, jobs=args.jobs)
                observed = status(database, scheduler=True)
                progress.update(observed['rows'], force=True)
                print(json.dumps(observed, sort_keys=True), flush=True)
                if observed['terminal']:
                    break
                if observed.get('scheduler_returncode', 0):
                    raise RuntimeError('squeue failed; inspect pipeline.json and submitted jobs')
                missing = set(observed.get('needs_attention', []))
                absent = {gcf: absent.get(gcf, 0) + 1 for gcf in missing}
                if any(count >= 3 for count in absent.values()):
                    summary = collect(database,jobs=args.collect_jobs,progress=collection_progress(progress))
                    summary.update(exit_code=1, missing_jobs=sorted(absent))
                    phase('NEEDS_ATTENTION', summary=summary, exit_code=1)
                    return summary
                time.sleep(args.interval)
            phase('COLLECTING')
            summary = collect(database,jobs=args.collect_jobs,progress=collection_progress(progress))
        failed = not summary['complete'] or any(summary.get('local_exit_codes', {}).values())
        summary['exit_code'] = int(failed)
        progress.update(force=True)
        phase('FAILED' if failed else summary['state'], summary=summary, exit_code=int(failed))
        return summary
    except BaseException as error:
        phase('INTERRUPTED' if isinstance(error, (KeyboardInterrupt, InterruptedError)) else 'FAILED',
              error=str(error), exit_code=1)
        raise


class ConcurrencyOption(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):
        if values < 1:
            raise argparse.ArgumentError(self, 'jobs must be positive')
        previous = getattr(namespace, '_jobs_seen', None)
        if previous is not None and previous[1] != values:
            raise argparse.ArgumentError(self, f'{option_string} conflicts with {previous[0]}; use one --jobs value')
        namespace.jobs = values
        namespace._jobs_seen = (option_string, values)


def add_jobs_option(target):
    target.add_argument('--jobs', type=int, default=1, action=ConcurrencyOption,
        help='Maximum local native sample tasks (legacy: GCFs), or Slurm GCF jobs (default1)')


def parser():
    root = argparse.ArgumentParser(description=__doc__)
    sub = root.add_subparsers(dest='command', required=True)
    prepare_parser = sub.add_parser('prepare', help='Validate manifests; write immutable plan and job scripts only')
    prepare_parser.add_argument('--records', type=Path, required=True)
    prepare_parser.add_argument('--references', type=Path, required=True)
    prepare_parser.add_argument('--outdir', type=Path, required=True)
    for flag, default in [('sample', 'Sample'), ('cram', 'CramFilePath'), ('gcf', 'GCFID'),
                          ('ref-gcf', 'GCF_ID'), ('fai', 'GCF_path')]:
        prepare_parser.add_argument('--' + flag + '-column', default=default)
    prepare_parser.add_argument('--cpus', type=int, default=4)
    prepare_parser.add_argument('--collect-jobs', type=int, default=4, help='Maximum native final-verification threads (default4)')
    prepare_parser.add_argument('--memory', default=None, help='Optional per-job memory, e.g. 8G; omit or use default for scheduler defaults')
    prepare_parser.add_argument('--partition', default='cpu', help='Slurm child-job partition (default cpu)')
    prepare_parser.add_argument('--python', default=sys.executable)
    geometry = prepare_parser.add_mutually_exclusive_group()
    geometry.add_argument('--window-size', '--bin-size', type=int,
                          help='Fixed window length in bp over the whole GCF (default 1000)')
    geometry.add_argument('--bin-count', type=int, help='Exact number of balanced windows per GCF')
    prepare_parser.add_argument('--reuse-geometry-from', type=Path,
                                help='Reproject verified saved geometry from an existing database into the new output')
    prepare_parser.add_argument('--geometry-format', choices=['auto', 'packed', 'npy', 'npz'], default='packed',
        help='Default packed: fixed zlib9, one verified candidate; auto is explicit format benchmarking')
    prepare_parser.add_argument('--mate-policy', choices=['both', 'read1', 'read2'], default='both')
    prepare_parser.add_argument('--alignment-capture', choices=['complete', 'legacy'], default='complete')
    prepare_parser.add_argument('--sampling-unit', choices=['read', 'fragment'], default='read',
        help='Legacy geometry only; complete geometry always supports runtime selection')
    prepare_parser.add_argument('--min-mapq', type=int, default=0)
    prepare_parser.add_argument('--exclude-duplicates', action='store_true')
    prepare_parser.add_argument('--exclude-qcfail', action='store_true')
    run_parser = sub.add_parser('run', parents=[prepare_parser], add_help=False,
        help='Prepare, execute/submit, wait and collect the database in one command')
    run_parser.add_argument('--mode', choices=['local', 'slurm'], required=True)
    run_parser.add_argument('--dry-run', action='store_true', help='Write plan and preview only; never start workers or contact Slurm')
    add_jobs_option(run_parser)
    run_parser.add_argument('--interval', type=float, default=10)
    for name in ('execute', 'submit', 'status', 'collect', 'worker'):
        p = sub.add_parser(name)
        p.add_argument('--database', type=Path, required=True)
        if name in ('execute','submit','status'):
            add_jobs_option(p)
        if name in ('execute','collect','status'):
            p.add_argument('--collect-jobs',type=int,default=4)
        if name in ('execute', 'submit'):
            p.add_argument('--dry-run', action='store_true')
        if name == 'execute':
            p.add_argument('--mode', choices=['local', 'slurm'], required=True)
        if name == 'submit':
            p.add_argument('--retry', action='store_true')
        if name == 'worker':
            p.add_argument('--gcf', required=True)
        if name == 'status':
            p.add_argument('--no-collect', action='store_true', help='Read-only watch: never collect or rewrite existing catalogs')
            p.add_argument('--refill',action='store_true',help='Explicitly submit deferred jobs while watching scheduler')
            p.add_argument('--scheduler', action='store_true')
            p.add_argument('--watch', action='store_true')
            p.add_argument('--interval', type=float, default=10)
    return root


def main():
    args = parser().parse_args()
    if args.command == 'prepare':
        value = prepare(args)
    elif args.command == 'run':
        value = run(args)
    elif args.command == 'worker':
        for key in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
            os.environ[key] = '1'
        def interrupted(signum, frame):
            raise InterruptedError('Worker interrupted')
        signal.signal(signal.SIGTERM, interrupted)
        from bin.ref_depth_worker import run_gcf
        return run_gcf(args.database, args.gcf)
    elif args.command == 'execute':
        if args.dry_run:
            value = dry_run(args.database, args.mode,args.jobs)
        elif args.mode == 'local':
            value = execute_local(args.database, args.jobs, collect_jobs=args.collect_jobs)
        else:
            database, plan = load_plan(args.database, check_code=True)
            generate_scripts(database, plan)
            value = dict(state='SCRIPTS_READY', submitted=False, job_count=len(plan['groups']))
    elif args.command == 'submit':
        value = submit(args.database, args.retry, preview=args.dry_run, jobs=args.jobs)
    elif args.command == 'collect':
        database,plan=load_plan(args.database)
        value = collect(database,jobs=args.collect_jobs,progress=collection_progress(ProgressLog(database,plan)))
    else:
        if not math.isfinite(args.interval) or args.interval <= 0:
            raise ValueError('interval must be positive')
        if args.refill and not (args.watch and args.scheduler):
            raise ValueError('--refill requires --watch --scheduler')
        if args.no_collect and args.refill:
            raise ValueError('--no-collect is read-only and cannot be combined with --refill')
        absent = {}
        database, plan = load_plan(args.database)
        progress = ProgressLog(database, plan, persist=False)
        while True:
            if args.refill:
                submit(args.database,jobs=args.jobs)
            value = status(args.database, args.scheduler)
            progress.update(value['rows'], force=True)
            if args.refill:
                if value.get('scheduler_returncode',0):
                    raise RuntimeError('squeue failed; inspect submitted jobs before refilling')
                missing = set(value.get('needs_attention',[]))
                absent = {gcf:absent.get(gcf,0)+1 for gcf in missing}
                if any(count>=3 for count in absent.values()):
                    raise RuntimeError('Submitted jobs absent from scheduler without terminal worker files; reconcile before continuing')
            if not args.watch or value['terminal']:
                break
            print(json.dumps(value), flush=True)
            time.sleep(args.interval)
        if args.watch and value['terminal'] and not args.no_collect:
            value['collection'] = collect(args.database,jobs=args.collect_jobs,progress=collection_progress(progress))
    print(json.dumps(value, sort_keys=True, indent=2))
    if args.command == 'run' or getattr(args,'dry_run',False):
        return value.get('exit_code',0)
    return 1 if (args.command in ('execute', 'collect') and value.get('complete') is False
                 or args.command == 'execute' and any(value.get('local_exit_codes', {}).values())
                 or args.command == 'status' and args.watch and value.get('collection', {}).get('complete') is False
                 or args.command == 'submit' and any(r['state'] != 'SUBMITTED' for r in value.values())) else 0


if __name__ == '__main__':
    try:
        sys.exit(main())
    except (ValueError, RuntimeError, OSError) as error:
        print(f'ERROR: {error}', file=sys.stderr)
        sys.exit(2)
