#!/usr/bin/env python3
"""Global planning, deterministic chunks and exclusive reference publication."""
import argparse
from concurrent.futures import ProcessPoolExecutor, wait, FIRST_COMPLETED
import hashlib
import json
import os
from pathlib import Path
import resource
import shutil
import sys
import time
import traceback

for _key in ('OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ[_key] = '1'
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from bin.ref_depth_catalog import (atomic_json, completed_sample, create_plan, digest,
    exclusive, file_stat, load_plan, validate_cram_index)
from bin.reference_store import bounded_map

VERSION = '3.1.2'


def fingerprint(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(',', ':'), allow_nan=False).encode()).hexdigest()


def record_plan(plan, gcf, sample):
    group = plan['groups'][gcf]
    group = {k: v for k, v in group.items() if k != 'samples'}
    group['samples'] = [sample]
    return dict(schema=plan['schema'], config=plan['config'], groups={gcf: group},
        counts=dict(gcfs=1, samples=1, sample_gcf_pairs=1))


def source_identity(item):
    gcf, group, sample = item
    before = file_stat(sample['cram'])
    if before != sample['source_stat']:
        raise ValueError('CRAM changed during plan validation')
    sha = digest(sample['cram'])
    if sample.get('source_sha256') is not None and sha != sample['source_sha256']:
        raise ValueError('Source content differs from immutable plan')
    if file_stat(sample['cram']) != before:
        raise ValueError('CRAM changed while hashing')
    validate_cram_index(sample)
    return gcf, dict(sample, source_sha256=sha)


def plan_stage(args, *, eager_sources=False):
    if min(args.chunk_size, args.cpus) < 1 or args.smoke_records < 0:
        raise ValueError('Require positive chunk/CPU budgets and nonnegative smoke record count')
    started = time.monotonic()
    print(f'PLAN VALIDATING records={args.records} references={args.references}',flush=True)
    config = dict(storage_mode='native', execution_model='nextflow-chunks', capture_alignment=True,
        window_policy='gcf-concat-v1', window_size=None if args.bin_count is not None else args.window_size, bin_count=args.bin_count,
        cpus=1, geometry_format='packed', insufficient_depth='error',
        policies=dict(mate_policy='both', min_mapq=0, exclude_duplicates=False, exclude_qcfail=False,
                      coverage_alignment_policy='primary_only'))
    from bin.gcf_windows import GCF_CONCAT
    config['window_policy'] = GCF_CONCAT
    columns = dict(sample='Sample', cram='CramFilePath', gcf='GCFID', ref_gcf='GCF_ID', fai='GCF_path')
    # Production preprocessing needs source byte sizes, not a global content
    # scan. Full content identity is bound by each independent batch worker.
    # eager_sources retains the frozen single-database path for acceptance.
    plan = create_plan(args.records, args.references, args.outdir, config, columns,
        validation_workers=args.cpus,
        source_validator=source_identity if eager_sources and not args.smoke_records else None,
        publish_plan=False)
    total = plan['counts']['sample_gcf_pairs']
    if total == 0:
        raise ValueError('Empty reference records')
    roster = [(gcf, sample) for gcf, group in plan['groups'].items() for sample in group['samples']]
    if args.smoke_records:
        roster = roster[:args.smoke_records]
        chosen = {(gcf, s['sample']) for gcf, s in roster}
        plan['groups'] = {g: dict(group, samples=[s for s in group['samples'] if (g,s['sample']) in chosen])
                          for g, group in plan['groups'].items() if any((g,s['sample']) in chosen for s in group['samples'])}
        plan['scope'] = dict(mode='explicit-smoke-subset', full_records=total, selected_records=len(roster),
                            excluded_records=total-len(roster), full_production_complete=False)
    else:
        plan['scope'] = dict(mode='full', full_records=total, selected_records=total, excluded_records=0)
    if args.smoke_records:
        jobs = ((g, plan['groups'][g], s) for g,s in roster)
        for completed, (gcf, sample) in enumerate(bounded_map(source_identity, jobs, args.cpus), 1):
            roster[completed - 1][1].update(sample)
            if completed % 1000 == 0 or completed == len(roster):
                print(f'PLAN HASHED {completed}/{len(roster)} elapsed={time.monotonic()-started:.1f}s', flush=True)
    plan['counts'] = dict(gcfs=len(plan['groups']), samples=len({s['sample'] for _,s in roster}), sample_gcf_pairs=len(roster))
    # No source revision or resource/concurrency knobs enter chunk scientific identity.
    plan['config'] = config
    if not eager_sources:
        from bin.ref_build_batches import MODEL
        plan.update(batch_model=MODEL, batching=dict(target_records=args.chunk_size,
            strategy='descending-gcf-bytes-to-lightest-batch', intact_gcfs=True), source_validation='per-batch-sha256')
    atomic_json(args.outdir/'database.json', plan)
    (args.outdir/'PLAN.sha256').write_text(digest(args.outdir/'database.json')+'\n')
    chunks = args.outdir/'chunks'
    chunks.mkdir()
    if not eager_sources:
        from bin.ref_build_batches import batches
        from bin.ref_depth_catalog import atomic_table
        summaries = []
        for chunk in batches(plan, args.chunk_size):
            atomic_json(chunks/(chunk['batch_id']+'.json'), chunk)
            summaries.append(dict(batch_id=chunk['batch_id'], records=chunk['counts']['sample_gcf_pairs'],
                gcfs=chunk['counts']['gcfs'], source_bytes=chunk['source_bytes']))
        atomic_table(args.outdir/'batches.tsv', ('batch_id','records','gcfs','source_bytes'), summaries)
        atomic_json(args.outdir/'planning.json', dict(records=len(roster), full_records=total, chunks=len(summaries),
            gcfs=len(plan['groups']), scope=plan['scope'], elapsed_seconds=time.monotonic()-started,
            validation_workers=args.cpus, source_validation='per-batch-sha256',
            source_bytes_total=sum(row['source_bytes'] for row in summaries),
            source_bytes_hashed=sum(s['source_stat']['bytes'] for _,s in roster if 'source_sha256' in s)))
        print(f'PLAN BATCHED batches={len(summaries)} records={len(roster)} full_source_hashes=deferred-to-batches',flush=True)
        return
    for i, start in enumerate(range(0, len(roster), args.chunk_size)):
        selected = roster[start:start+args.chunk_size]
        chunk_groups = {}
        for gcf, sample in selected:
            chunk_groups.setdefault(gcf, dict(plan['groups'][gcf], samples=[]))['samples'].append(sample)
        chunk = dict(schema=plan['schema'], config=config, groups=chunk_groups,
            counts=dict(gcfs=len(chunk_groups), samples=len({s['sample'] for _,s in selected}), sample_gcf_pairs=len(selected)))
        atomic_json(chunks/f'chunk-{i:06d}.json', chunk)
    atomic_json(args.outdir/'planning.json', dict(records=len(roster), full_records=total, chunks=i+1,
        gcfs=len(plan['groups']), scope=plan['scope'], elapsed_seconds=time.monotonic()-started,
        validation_workers=args.cpus,
        source_bytes_hashed=sum(sample['source_stat']['bytes'] for _, sample in roster)))


def build_record(task):
    plan, gcf, sample, checkpoints = task
    started = time.monotonic()
    mini = record_plan(plan, gcf, sample)
    key = fingerprint(mini)
    root = Path(checkpoints)/key
    receipt = dict(key=key, gcf=gcf, sample=sample['sample'], checkpoint=str(root), pid=os.getpid(), reused=False)
    receipt['started_unix'] = time.time()
    try:
        hashed_here = 'source_sha256' not in sample
        if hashed_here:
            _, sample = source_identity((gcf, plan['groups'][gcf], sample))
            mini = record_plan(plan, gcf, sample)
            key = fingerprint(mini)
            root = Path(checkpoints)/key
            receipt.update(key=key, checkpoint=str(root))
        root.mkdir(parents=True, exist_ok=True)
        receipt['source_sha256'] = sample['source_sha256']
        with exclusive(root/'.record.lock'):
            if file_stat(sample['cram']) != sample['source_stat']:
                raise ValueError('Source differs from immutable plan')
            if digest(plan['groups'][gcf]['fai']) != plan['groups'][gcf]['fai_sha256']:
                raise ValueError('FAI differs from immutable plan')
            if not (root/'database.json').exists():
                atomic_json(root/'database.json', mini)
                (root/'PLAN.sha256').write_text(digest(root/'database.json')+'\n')
            elif json.loads((root/'database.json').read_text()) != mini:
                raise ValueError('Record checkpoint identity collision')
            sha = digest(root/'database.json')
            parent = root/gcf/'samples'/sample['folder']
            parent.mkdir(parents=True, exist_ok=True)
            previous = completed_sample(parent, sha, verify=True)
            receipt['reused'] = previous is not None
            if previous is not None and not hashed_here:
                source_identity((gcf, plan['groups'][gcf], sample))
            if previous is None:
                from bin.ref_depth_worker import execute_sample
                previous = execute_sample(root, mini, gcf, sample)
            attempt, result = previous[:2]
            receipt.update(state='SUCCESS', attempt=str(attempt), record_plan_sha256=sha,
                           result_sha256=digest(attempt/'result.json'))
    except Exception as error:
        receipt.update(state='FAILED', error=f'{type(error).__name__}: {error}', traceback=traceback.format_exc())
    receipt.update(elapsed_seconds=time.monotonic()-started,
                   finished_unix=time.time(),
                   process_peak_rss_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*(1 if sys.platform=='darwin' else 1024))
    root.mkdir(parents=True, exist_ok=True)
    atomic_json(root/f'record_receipt-{os.getpid()}-{time.time_ns()}.json', receipt)
    return receipt


def chunk_stage(args):
    if args.cpus < 1:
        raise ValueError('cpus must be positive')
    started = time.monotonic()
    plan = json.loads(args.chunk.read_text())
    if plan.get('batch_model') and getattr(args, 'database', None) is None:
        raise ValueError('Whole-GCF batches require --database for independent publication')
    args.checkpoints = args.checkpoints.resolve()
    args.checkpoints.mkdir(parents=True, exist_ok=True)
    # Serialize one record, not the entire (possibly very large) GCF batch,
    # for every process-pool submission.
    tasks = ((record_plan(plan,g,s),g,s,str(args.checkpoints))
             for g, group in plan['groups'].items() for s in group['samples'])
    workers = min(args.cpus, plan['counts']['sample_gcf_pairs'])
    if not workers:
        raise ValueError('Empty chunk')
    # Executor.map is unbounded on older Python; keep at most 2*workers futures.
    results = []
    tasks = iter(tasks)
    exhausted = False
    print(f'CHUNK START {args.chunk.stem} records={plan["counts"]["sample_gcf_pairs"]} workers={workers}',flush=True)
    with ProcessPoolExecutor(max_workers=workers) as executor:
        pending = set()
        while pending or not exhausted:
            while not exhausted and len(pending)<2*workers:
                try:
                    task = next(tasks)
                except StopIteration:
                    exhausted = True
                    break
                pending.add(executor.submit(build_record,task))
            if pending:
                done,pending = wait(pending,return_when=FIRST_COMPLETED)
                for future in done:
                    row = future.result()
                    results.append(row)
                    print(f"RECORD {row['state']} {row['key']} reused={row['reused']} pid={row['pid']} seconds={row['elapsed_seconds']:.3f} progress={len(results)}/{plan['counts']['sample_gcf_pairs']}",flush=True)
    results.sort(key=lambda row:(row['gcf'],row['sample']))
    failed = sum(r['state']!='SUCCESS' for r in results)
    receipt = dict(chunk_id=args.chunk.stem, chunk_sha256=digest(args.chunk),
        records=results, completed=len(results)-failed, failed=failed, skipped=0,
        workers=workers, worker_pids=sorted({r['pid'] for r in results}),
        state='FAILED' if failed else 'SUCCESS', elapsed_seconds=time.monotonic()-started)
    if not failed and plan.get('batch_model'):
        from bin.ref_build_batches import publish_batch
        receipt.update(publish_batch(plan, results, args.database, args.cpus))
    atomic_json(args.receipt, receipt)
    if failed:
        raise ValueError(f'{failed} records failed; see {args.receipt}; verified successes retained in {args.checkpoints}')


def collect_stage(args):
    started = time.monotonic()
    _, plan = load_plan(args.plan)
    if plan.get('batch_model'):
        from bin.ref_build_batches import collect_batches
        collect_batches(args, plan)
        return
    expected = {(g,s['sample']): s for g,group in plan['groups'].items() for s in group['samples']}
    chunks = {p.stem: p for p in (args.plan/'chunks').glob('*.json')}
    receipts = {}
    records = {}
    for path in sorted(args.receipts.glob('*.json')):
        receipt = json.loads(path.read_text())
        cid = receipt['chunk_id']
        if cid in receipts or cid not in chunks or receipt['chunk_sha256'] != digest(chunks[cid]):
            raise ValueError('Duplicate, unknown or changed chunk receipt')
        if receipt['state'] != 'SUCCESS' or receipt['failed'] or receipt['skipped']:
            raise ValueError('Incomplete/failed chunk cannot be collected')
        receipts[cid] = receipt
        roster = json.loads(chunks[cid].read_text())['groups']
        roster = {(g,s['sample']) for g,v in roster.items() for s in v['samples']}
        observed = set()
        for row in receipt['records']:
            key = row['gcf'],row['sample']
            if key in records or key not in roster or row['state'] != 'SUCCESS':
                raise ValueError('Duplicate, unexpected or failed record receipt')
            records[key] = row
            observed.add(key)
        if observed != roster or receipt['completed'] != len(roster):
            raise ValueError('Missing record receipt')
    if set(receipts) != set(chunks) or set(records) != set(expected):
        raise ValueError('Missing chunk or record receipts')
    database = args.database.resolve()
    database.parent.mkdir(parents=True, exist_ok=True)
    with exclusive(database.parent/('.'+database.name+'.collector.lock')):
        if database.exists():
            if not (database/'BUILD_RECEIPT.json').exists():
                raise FileExistsError('Existing incomplete output retained; inspect before retry: '+str(database))
            old = json.loads((database/'BUILD_RECEIPT.json').read_text())
            if old['plan_sha256'] != digest(args.plan/'database.json'):
                raise ValueError('Published database belongs to another plan; use new outdir')
            summary = json.loads((database/'SUMMARY.json').read_text())
            if not summary['complete'] or summary['records_completed'] != len(expected):
                raise ValueError('Published database verification failed')
            sha = digest(database/'database.json')
            for (gcf,sample_id), sample in expected.items():
                value = completed_sample(database/gcf/'samples'/sample['folder'],sha,verify=True)
                if value is None or value[1]['sources'][0]['sha256'] != sample['source_sha256']:
                    raise ValueError('Published reference payload verification failed')
        else:
            stage = database.parent/('.'+database.name+'.collecting')
            stage.mkdir(exist_ok=False)
            atomic_json(stage/'database.json', plan)
            shutil.copytree(args.plan/'inputs',stage/'inputs')
            sha = digest(stage/'database.json')
            (stage/'PLAN.sha256').write_text(sha+'\n')
            for (gcf,sample_id), sample in expected.items():
                row = records[(gcf,sample_id)]
                if file_stat(sample['cram']) != sample['source_stat']:
                    raise ValueError('Source stat changed before collection')
                checkpoint = Path(row['checkpoint'])
                _, mini = load_plan(checkpoint)
                if mini != record_plan(plan,gcf,sample):
                    raise ValueError('Record checkpoint differs from expected roster')
                value = completed_sample(checkpoint/gcf/'samples'/sample['folder'], row['record_plan_sha256'], verify=True)
                if value is None or digest(value[0]/'result.json') != row['result_sha256']:
                    raise ValueError('Record payload differs from receipt')
                attempt,result = value
                if result['sources'][0]['sha256'] != sample['source_sha256']:
                    raise ValueError('Built source content differs from immutable plan')
                target = stage/gcf/'samples'/sample['folder']/'attempt-000001'
                target.parent.mkdir(parents=True,exist_ok=True)
                atomic_json(stage/gcf/'reference.json',{k:v for k,v in plan['groups'][gcf].items() if k!='samples'})
                shutil.copytree(attempt,target)
                result = dict(result,plan_sha256=sha)
                atomic_json(target/'result.json',result)
                checks = json.loads((target/'checksums.json').read_text())
                checks['result.json'] = digest(target/'result.json')
                atomic_json(target/'checksums.json',checks)
            from bin.native_reference import collect_native
            summary = collect_native(stage,plan,verify=True,jobs=args.cpus,publish_root=database)
            if not summary['complete']:
                raise ValueError('Incomplete database retained at '+str(stage))
            summary.update(build_scope=plan['scope'],records_total=len(expected),records_completed=len(records),
                           records_failed=0,records_skipped=0,chunks=len(chunks))
            atomic_json(stage/'SUMMARY.json',summary)
            atomic_json(stage/'BUILD_RECEIPT.json',dict(plan_sha256=sha,records=len(records),chunks=len(chunks),
                        elapsed_seconds=time.monotonic()-started,scope=plan['scope']))
            stage.rename(database)
    atomic_json(args.receipt,dict(database=str(database),state='COLLECTED',records=len(records),chunks=len(chunks)))


def index_stage(args):
    from bin.reference_index import build_index, ReferenceIndex
    database = args.database.resolve()
    batch_sha = None
    sources = database
    if (database/'BATCHES.json').exists():
        from bin.ref_build_batches import index_sources
        sources = index_sources(database)
        batch_sha = digest(database/'BATCHES.json')
    with exclusive(database/'.index.lock'):
        output = database.parent/'reference.sqlite'
        if output.exists():
            old = json.loads((database/'INDEX_RECEIPT.json').read_text())
            if old['sha256'] != digest(output) or old['build_plan_sha256'] != digest(database/'database.json'):
                raise ValueError('Existing index differs from published receipt')
            if batch_sha is not None and old.get('batches_sha256') != batch_sha:
                raise ValueError('Existing index belongs to different batches')
            index = ReferenceIndex(output)
            index.close()
            result = old
        else:
            result = build_index(sources,output,verify_payloads=True,jobs=max(1,args.cpus-1))
            result['build_plan_sha256'] = digest(database/'database.json')
            if batch_sha is not None:
                result['batches_sha256'] = batch_sha
            atomic_json(database/'INDEX_RECEIPT.json',result)
        (database/'SUCCESS').touch()
    atomic_json(args.receipt,dict(index=str(output),sha256=digest(output),state='COMPLETE',
        scope=json.loads((database/'database.json').read_text())['scope']))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--version',action='version',version=VERSION)
    sub = parser.add_subparsers(dest='stage',required=True)
    plan = sub.add_parser('plan')
    for name in ('records','references','outdir'):
        plan.add_argument('--'+name,type=Path,required=True)
    plan.add_argument('--chunk-size',type=int,required=True)
    plan.add_argument('--smoke-records',type=int,default=0)
    plan.add_argument('--window-size',type=int,default=1000)
    plan.add_argument('--bin-count',type=int)
    plan.set_defaults(action=plan_stage)
    chunk = sub.add_parser('chunk')
    for name in ('chunk','checkpoints','receipt'):
        chunk.add_argument('--'+name,type=Path,required=True)
    chunk.add_argument('--database',type=Path,help='Independent batch database output (required for whole-GCF batches)')
    chunk.set_defaults(action=chunk_stage)
    collect = sub.add_parser('collect')
    for name in ('plan','receipts','database','receipt'):
        collect.add_argument('--'+name,type=Path,required=True)
    collect.set_defaults(action=collect_stage)
    index = sub.add_parser('index')
    for name in ('database','receipt'):
        index.add_argument('--'+name,type=Path,required=True)
    index.set_defaults(action=index_stage)
    for stage in (plan,chunk,collect,index):
        stage.add_argument('--cpus',type=int,default=1)
    args = parser.parse_args()
    if args.cpus < 1:
        parser.error('--cpus must be positive')
    try:
        args.action(args)
    except Exception as error:
        print(f'{args.stage.upper()} FAILED: {type(error).__name__}: {error}',file=sys.stderr,flush=True)
        raise


if __name__ == '__main__':
    main()
