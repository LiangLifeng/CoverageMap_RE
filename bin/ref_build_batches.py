"""Whole-GCF byte balancing and verified multi-root build publication."""
import csv
import heapq
import io
import json
from pathlib import Path
import re
import shutil
from types import SimpleNamespace

from bin.ref_depth_catalog import atomic_json, digest, exclusive, load_plan

MODEL = 'whole-gcf-byte-balanced-v1'


def batches(plan, target_records):
    """Deterministic longest-processing-time assignment using measured file bytes."""
    if target_records < 1 or not plan['groups']:
        raise ValueError('Require a positive batch target and nonempty groups')
    count = min(len(plan['groups']), (plan['counts']['sample_gcf_pairs'] + target_records - 1) // target_records)
    heap = [(0, 0, i) for i in range(count)]
    heapq.heapify(heap)
    assigned = [[] for _ in range(count)]
    sizes = {g: sum(s['source_stat']['bytes'] for s in group['samples']) for g, group in plan['groups'].items()}
    for gcf in sorted(plan['groups'], key=lambda g: (-sizes[g], g)):
        size, records, index = heapq.heappop(heap)
        assigned[index].append(gcf)
        heapq.heappush(heap, (size + sizes[gcf], records + len(plan['groups'][gcf]['samples']), index))
    for index, gcfs in enumerate(sorted(assigned, key=lambda gs: min(gs))):
        groups = {g: plan['groups'][g] for g in sorted(gcfs)}
        samples = [s for group in groups.values() for s in group['samples']]
        value = {key: plan[key] for key in ('schema', 'config', 'columns', 'code_sha256', 'input_manifests') if key in plan}
        value.update(groups=groups, batch_model=MODEL, batch_id=f'chunk-{index:06d}',
            source_bytes=sum(sizes[g] for g in gcfs),
            counts=dict(gcfs=len(groups), samples=len({s['sample'] for s in samples}), sample_gcf_pairs=len(samples)),
            scope=dict(mode='whole-gcf-batch', parent_scope=plan['scope'],
                       parent_records=plan['counts']['sample_gcf_pairs'], selected_records=len(samples)))
        yield value


def resolved_batch(plan, records):
    """Bind each planned source stat/path to its verified content identity."""
    from bin.ref_build_pipeline import fingerprint, record_plan
    expected = {(g, s['sample']): s for g, group in plan['groups'].items() for s in group['samples']}
    observed = {}
    for row in records:
        key = row['gcf'], row['sample']
        if key not in expected or key in observed or row['state'] != 'SUCCESS':
            raise ValueError('Duplicate, foreign or failed batch record')
        sha = row.get('source_sha256', '')
        if not isinstance(sha, str) or not re.fullmatch('[0-9a-f]{64}', sha):
            raise ValueError('Missing/invalid verified source SHA256')
        sample = dict(expected[key], source_sha256=sha)
        if expected[key].get('source_sha256', sha) != sha or fingerprint(record_plan(plan, key[0], sample)) != row['key']:
            raise ValueError('Batch source/checkpoint identity differs from planned record')
        observed[key] = sample
    if set(observed) != set(expected):
        raise ValueError('Missing batch record receipt')
    groups = {g: dict(group, samples=[observed[g, s['sample']] for s in group['samples']])
              for g, group in plan['groups'].items()}
    return dict(plan, groups=groups)


def frozen_json(path, value):
    path = Path(path)
    if path.exists():
        if json.loads(path.read_text()) != value:
            raise ValueError('Existing batch identity differs; use a fresh outdir: ' + str(path))
    else:
        atomic_json(path, value)


def frozen_csv(path, fields, rows):
    buffer = io.StringIO(newline='')
    writer = csv.DictWriter(buffer, fields, lineterminator='\n')
    writer.writeheader()
    writer.writerows(rows)
    value = buffer.getvalue()
    if path.exists() and path.read_text() != value:
        raise ValueError('Existing batch input snapshot differs: ' + str(path))
    if not path.exists():
        path.write_text(value)


def publish_batch(chunk, records, database, cpus):
    from bin.ref_build_pipeline import collect_stage
    plan = resolved_batch(chunk, records)
    # Reuse the established strict single-database collector within this batch.
    # Scientific record identities exclude these operational partition markers.
    plan.pop('batch_model', None)
    plan.pop('batch_id', None)
    database = Path(database).resolve()
    root = database.parent
    plan_root, receipts = root / 'plan', root / 'receipts'
    for path in (plan_root / 'chunks', plan_root / 'inputs', receipts):
        path.mkdir(parents=True, exist_ok=True)
    with exclusive(root / '.batch-publish.lock'):
        frozen_json(plan_root / 'database.json', plan)
        (plan_root / 'PLAN.sha256').write_text(digest(plan_root / 'database.json') + '\n')
        frozen_csv(plan_root / 'inputs/sample_manifest.csv', ['Sample', 'CramFilePath', 'GCFID'],
            [dict(Sample=s['sample'], CramFilePath=s['cram'], GCFID=g)
             for g, group in plan['groups'].items() for s in group['samples']])
        frozen_csv(plan_root / 'inputs/reference_manifest.csv', ['GCF_ID', 'GCF_path'],
            [dict(GCF_ID=g, GCF_path=group['fai']) for g, group in plan['groups'].items()])
        chunk_path = plan_root / 'chunks/chunk-000000.json'
        frozen_json(chunk_path, plan)
        atomic_json(receipts / 'chunk-000000.json', dict(chunk_id=chunk_path.stem,
            chunk_sha256=digest(chunk_path), records=records, completed=len(records), failed=0, skipped=0, state='SUCCESS'))
        collect_stage(SimpleNamespace(plan=plan_root, receipts=receipts, database=database,
            cpus=cpus, receipt=root / 'collected.json'))
    return dict(database=str(database), database_plan_sha256=digest(database / 'database.json'),
                database_build_receipt_sha256=digest(database / 'BUILD_RECEIPT.json'))


def collect_batches(args, plan):
    expected = {chunk['batch_id']: chunk for chunk in batches(plan, plan['batching']['target_records'])}
    files = {path.stem: path for path in (args.plan / 'chunks').glob('*.json')}
    if set(files) != set(expected):
        raise ValueError('Missing or unexpected planned batches')
    observed = {}
    entries = []
    for path in sorted(args.receipts.glob('*.json')):
        receipt = json.loads(path.read_text())
        cid = receipt['chunk_id']
        if cid not in expected or cid in observed:
            raise ValueError('Duplicate or foreign batch receipt')
        chunk = expected[cid]
        if json.loads(files[cid].read_text()) != chunk or receipt['chunk_sha256'] != digest(files[cid]):
            raise ValueError('Changed batch plan or receipt')
        if receipt['state'] != 'SUCCESS' or receipt['failed'] or receipt['skipped']:
            raise ValueError('Incomplete or failed batch')
        resolved = resolved_batch(chunk, receipt['records'])
        resolved.pop('batch_model', None)
        resolved.pop('batch_id', None)
        if receipt['completed'] != resolved['counts']['sample_gcf_pairs']:
            raise ValueError('Batch record count differs')
        database = (args.database.parent / 'batches' / cid / 'database').resolve()
        if Path(receipt['database']).resolve() != database:
            raise ValueError('Foreign batch database path')
        _, actual = load_plan(database)
        if actual != resolved or digest(database / 'database.json') != receipt['database_plan_sha256']:
            raise ValueError('Published batch plan differs')
        if digest(database / 'BUILD_RECEIPT.json') != receipt['database_build_receipt_sha256']:
            raise ValueError('Published batch receipt differs')
        summary = json.loads((database / 'SUMMARY.json').read_text())
        if not summary['complete'] or summary['records_completed'] != receipt['completed']:
            raise ValueError('Published batch is incomplete')
        entries.append(dict(batch_id=cid, database=str(database),
            plan_sha256=receipt['database_plan_sha256'], build_receipt_sha256=receipt['database_build_receipt_sha256'],
            records=receipt['completed'], gcfs=sorted(chunk['groups']), source_bytes=chunk['source_bytes']))
        observed[cid] = receipt
    if set(observed) != set(expected):
        raise ValueError('Missing batch receipts')
    entries.sort(key=lambda value: value['batch_id'])
    summary = dict(complete=True, batch_model=MODEL, batches=len(entries), build_scope=plan['scope'],
        counts=plan['counts'],
        records_total=plan['counts']['sample_gcf_pairs'], records_completed=sum(e['records'] for e in entries),
        records_failed=0, records_skipped=0, gcfs=sum(len(e['gcfs']) for e in entries),
        source_bytes=sum(e['source_bytes'] for e in entries))
    if summary['records_completed'] != summary['records_total'] or summary['gcfs'] != len(plan['groups']):
        raise ValueError('Global batch conservation failed')
    args.database.parent.mkdir(parents=True, exist_ok=True)
    with exclusive(args.database.parent / ('.' + args.database.name + '.collector.lock')):
        if args.database.exists():
            if (json.loads((args.database / 'database.json').read_text()) != plan
                    or json.loads((args.database / 'BATCHES.json').read_text()) != entries
                    or json.loads((args.database / 'SUMMARY.json').read_text()) != summary):
                raise ValueError('Existing merged build differs; use a fresh outdir')
        else:
            stage = args.database.parent / ('.' + args.database.name + '.collecting')
            stage.mkdir(exist_ok=False)
            atomic_json(stage / 'database.json', plan)
            (stage / 'PLAN.sha256').write_text(digest(stage / 'database.json') + '\n')
            shutil.copytree(args.plan / 'inputs', stage / 'inputs')
            shutil.copy2(args.plan / 'batches.tsv', stage / 'batches.tsv')
            atomic_json(stage / 'BATCHES.json', entries)
            atomic_json(stage / 'SUMMARY.json', summary)
            atomic_json(stage / 'BUILD_RECEIPT.json', dict(plan_sha256=digest(stage / 'database.json'),
                batches_sha256=digest(stage / 'BATCHES.json'), records=summary['records_completed'], batches=len(entries)))
            stage.rename(args.database)
    atomic_json(args.receipt, dict(database=str(args.database), state='COLLECTED',
        records=summary['records_completed'], chunks=len(entries), batch_model=MODEL))


def index_sources(database):
    """Return immutable complete roots; actual payload validation stays in the registry."""
    manifest = database / 'BATCHES.json'
    entries = json.loads(manifest.read_text())
    if not entries or len({entry['batch_id'] for entry in entries}) != len(entries):
        raise ValueError('Empty or duplicate batch index inputs')
    receipt = json.loads((database / 'BUILD_RECEIPT.json').read_text())
    if receipt['batches_sha256'] != digest(manifest):
        raise ValueError('Batch aggregation changed before indexing')
    roots = []
    for entry in entries:
        root = Path(entry['database'])
        if (digest(root / 'database.json') != entry['plan_sha256']
                or digest(root / 'BUILD_RECEIPT.json') != entry['build_receipt_sha256']):
            raise ValueError('Batch database changed before indexing')
        roots.append(root)
    return roots
