"""Immutable plans, safe paths, job generation and catalogs for the depth DB."""
import csv
from concurrent.futures import FIRST_COMPLETED, ThreadPoolExecutor, wait
from contextlib import contextmanager
import hashlib
import io
import json
import os
from pathlib import Path
import re
import shlex
import socket
import stat as stat_module
import sys
import time
from urllib.parse import quote

SCHEMA = 'gcf-depth-database-1'
CODE_ROOT = Path(__file__).resolve().parents[1]
CODE_FILES = ('bin/ref_depth_database.py', 'bin/ref_depth_catalog.py',
              'bin/ref_depth_worker.py', 'bin/coverage_core.py',
              'bin/concordance_core.py', 'bin/reference_store.py', 'bin/native_reference.py',
              'bin/gcf_windows.py', 'bin/alignment_geometry.py', 'bin/native_execution.py', 'bin/stage_timing.py', 'Tools.py')


def digest(path):
    value = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b''):
            value.update(block)
    return value.hexdigest()


def atomic_json(path, value):
    path = Path(path)
    temporary = path.with_name(path.name + f'.tmp-{os.getpid()}')
    # Immutable cohort plans can contain hundreds of thousands of rows. Stream
    # the identical encoding instead of allocating a second complete JSON string.
    with temporary.open('w') as handle:
        json.dump(value, handle, sort_keys=True, indent=2, allow_nan=False)
        handle.write('\n')
    temporary.replace(path)


def atomic_table(path, columns, rows):
    path = Path(path)
    temporary = path.with_name(path.name + f'.tmp-{os.getpid()}')
    with temporary.open('w', newline='') as handle:
        writer = csv.DictWriter(handle, columns, delimiter='\t', extrasaction='ignore')
        writer.writeheader()
        writer.writerows(rows)
    temporary.replace(path)


@contextmanager
def exclusive(path):
    path = Path(path)
    try:
        path.mkdir()
    except FileExistsError:
        raise RuntimeError(f'Locked: {path}; check the worker/job before recovering a stale lock')
    try:
        atomic_json(path / 'owner.json', dict(pid=os.getpid(), host=socket.gethostname(), started=time.time()))
        yield
    finally:
        (path / 'owner.json').unlink(missing_ok=True)
        path.rmdir()


def sample_folder(sample):
    return 'sample-' + quote(sample, safe='._-')[:80] + '-' + hashlib.sha256(sample.encode()).hexdigest()[:12]


def safe_file(root, relative):
    root, relative = Path(root), Path(relative)
    path = root / relative
    if relative.is_absolute() or '..' in relative.parts or not path.resolve().is_relative_to(root.resolve()):
        raise ValueError('Unsafe stored path')
    if not path.is_file():
        raise ValueError(f'Missing file: {path}')
    return path


def read_rows(path, columns, content=None):
    with io.StringIO((Path(path).read_bytes() if content is None else content).decode('utf-8-sig'), newline='') as handle:
        reader = csv.DictReader(handle)
        if not reader.fieldnames or not set(columns) <= set(reader.fieldnames):
            raise ValueError(f'Missing required columns in {path}: {columns}')
        rows = list(reader)
    if not rows:
        raise ValueError(f'Empty manifest: {path}')
    for row in rows:
        for column in columns:
            value = row.get(column)
            if value is None or not value.strip() or value != value.strip() or any(c in value for c in '\r\n\t\x00'):
                raise ValueError(f'Empty/whitespace/control value in {column}')
    return rows


def input_file_info(value, manifest):
    """Resolve an input and reuse its one regular-file stat for planning."""
    path = Path(value)
    path = path if path.is_absolute() else Path(manifest).parent / path
    path = path.resolve()
    try:
        info = path.stat()
    except (FileNotFoundError, NotADirectoryError):
        raise FileNotFoundError(f'Input not readable: {path}') from None
    if not stat_module.S_ISREG(info.st_mode) or not os.access(path, os.R_OK):
        raise FileNotFoundError(f'Input not readable: {path}')
    if info.st_size == 0:
        raise ValueError(f'Input file is empty: {path}')
    return path, dict(bytes=info.st_size, mtime_ns=info.st_mtime_ns)


def input_file(value, manifest):
    return input_file_info(value, manifest)[0]


def file_stat(path):
    stat = Path(path).stat()
    return dict(bytes=stat.st_size, mtime_ns=stat.st_mtime_ns)


def source_cram_index(cram):
    """Record an optional adjacent CRAI; CRAM-only inputs remain supported."""
    cram = Path(cram)
    # A .cram suffix already selects this branch in the legacy policy. Only
    # unknown extensions need a magic lookup; full source hashing still reads
    # every byte when the immutable Nextflow plan is constructed.
    if cram.suffix.lower() != '.cram':
        with cram.open('rb') as handle:
            if handle.read(4) != b'CRAM':
                return None
    index = Path(str(cram) + '.crai')
    try:
        info = index.lstat()
    except FileNotFoundError:
        return None
    if stat_module.S_ISLNK(info.st_mode):
        try:
            info = index.stat()
        except FileNotFoundError:
            raise FileNotFoundError(f'Optional CRAM index exists but is unreadable: {index}') from None
    if not stat_module.S_ISREG(info.st_mode) or not os.access(index, os.R_OK):
        raise FileNotFoundError(f'Optional CRAM index exists but is unreadable: {index}')
    before = dict(bytes=info.st_size, mtime_ns=info.st_mtime_ns)
    if before['bytes'] == 0:
        raise ValueError(f'Optional CRAM index is empty: {index}')
    sha = digest(index)
    if file_stat(index) != before:
        raise ValueError(f'CRAM index changed while planning: {index}')
    return dict(path=str(index), sha256=sha, stat=before)


def validate_cram_index(sample):
    index = sample.get('cram_index')
    if index is None:  # Legacy frozen plans/BAM inputs remain readable.
        return
    path = Path(index['path'])
    if file_stat(path) != index['stat'] or digest(path) != index['sha256']:
        raise ValueError(f'CRAM index changed since plan: {path}; restore original pair or create a new plan')


def read_layout(path):
    names, lengths = [], []
    for line in Path(path).read_text().splitlines():
        fields = line.split('\t')
        if len(fields) < 2 or not fields[0] or any(c in fields[0] for c in '\r\n\t\x00'):
            raise ValueError('Invalid FAI row')
        length = int(fields[1])
        if length < 1:
            raise ValueError('Nonpositive FAI length')
        names.append(fields[0]); lengths.append(length)
    if not names or len(set(names)) != len(names):
        raise ValueError('Empty/duplicate FAI contigs')
    return dict(names=names, lengths=lengths)


def bounded_plan_map(function, jobs, workers):
    """Work-conserving planning I/O; at most 2*workers futures, no row reordering contract."""
    if workers < 1:
        raise ValueError('validation_workers must be positive')
    if workers == 1:
        for job in jobs:
            yield function(job)
        return
    jobs = iter(jobs)
    executor = ThreadPoolExecutor(max_workers=workers)
    pending = set()
    exhausted = False
    try:
        while pending or not exhausted:
            while not exhausted and len(pending) < 2 * workers:
                try:
                    job = next(jobs)
                except StopIteration:
                    exhausted = True
                    break
                pending.add(executor.submit(function, job))
            if pending:
                done, pending = wait(pending, return_when=FIRST_COMPLETED)
                for future in done:
                    yield future.result()
    finally:
        executor.shutdown(wait=True, cancel_futures=True)


def create_plan(records, references, outdir, config, columns, *,
                validation_workers=1, source_validator=None, publish_plan=True):
    from bin.coverage_core import canonical_reference, layout_metadata, reference_from_layout
    from bin.gcf_windows import GCF_CONCAT, window_policy

    records, references, outdir = Path(records).resolve(), Path(references).resolve(), Path(outdir).resolve()
    if validation_workers < 1:
        raise ValueError('validation_workers must be positive')
    if outdir.exists():
        raise FileExistsError(f'Plan output already exists: {outdir}')
    started = time.monotonic()
    snapshots = {p: p.read_bytes() for p in (records, references)}
    rows = read_rows(records, [columns[k] for k in ('sample', 'cram', 'gcf')], snapshots[records])
    refs = read_rows(references, [columns[k] for k in ('ref_gcf', 'fai')], snapshots[references])
    lookup = {}
    for row in refs:
        gcf = row[columns['ref_gcf']]
        if gcf in lookup:
            raise ValueError(f'Duplicate reference ID: {gcf}')
        lookup[gcf] = row[columns['fai']]
    # Validate cohort-wide keys cheaply before dispatching any per-file I/O.
    pairs, required_gcfs = set(), set()
    for row in rows:
        sample, gcf = row[columns['sample']], row[columns['gcf']]
        if not re.fullmatch(r'(?:GCF|GCA)_\d+\.\d+', gcf):
            raise ValueError(f'Invalid assembly ID: {gcf}')
        if (sample, gcf) in pairs:
            raise ValueError(f'Duplicate sample/GCF: {sample}/{gcf}')
        pairs.add((sample, gcf))
        if gcf not in lookup:
            raise ValueError(f'Missing FAI mapping: {gcf}')
        required_gcfs.add(gcf)
    sample_count = len({sample for sample, _ in pairs})
    del pairs
    progress = config.get('execution_model') == 'nextflow-chunks'
    if progress:
        mode = 'metadata+sha256' if source_validator is not None else 'metadata'
        print(f'PLAN PREFLIGHT records={len(rows)} gcfs={len(required_gcfs)} '
              f'workers={validation_workers} mode={mode}', file=sys.stderr, flush=True)

    def inspect_reference(gcf):
        fai, before = input_file_info(lookup[gcf], references)
        layout = read_layout(fai)
        if config.get('storage_mode') == 'native' and window_policy(config) == GCF_CONCAT:
            layout = layout_metadata(canonical_reference(reference_from_layout(layout)))
        bins = config.get('bin_count') if config.get('bin_count') is not None else config.get('bins')
        if bins is not None and not 1 <= bins <= sum(layout['lengths']):
            raise ValueError(f'Too many bins for {gcf}')
        sha = digest(fai)
        if file_stat(fai) != before:
            raise ValueError(f'FAI changed while planning: {fai}')
        return gcf, dict(gcf=gcf, fai=str(fai), fai_sha256=sha, layout=layout, samples=[])

    groups = {}
    for number, (gcf, group) in enumerate(bounded_plan_map(
            inspect_reference, sorted(required_gcfs), validation_workers), 1):
        groups[gcf] = group
        if progress and (number % 1000 == 0 or number == len(required_gcfs)):
            print(f'PLAN REFERENCES {number}/{len(required_gcfs)} '
                  f'elapsed={time.monotonic()-started:.1f}s', file=sys.stderr, flush=True)

    def inspect_sample(row):
        sample, gcf = row[columns['sample']], row[columns['gcf']]
        try:
            cram, source_stat = input_file_info(row[columns['cram']], records)
            folder = sample_folder(sample)
            value = dict(sample=sample, cram=str(cram), source_stat=source_stat,
                         cram_index=source_cram_index(cram), folder=folder)
            return (gcf, value) if source_validator is None else source_validator((gcf, groups[gcf], value))
        except Exception as error:
            error.add_note(f'PLAN input sample={sample} gcf={gcf} path={row[columns["cram"]]}')
            raise

    seen_paths = set()
    hashed_bytes = 0
    for row_number, (gcf, value) in enumerate(bounded_plan_map(
            inspect_sample, rows, validation_workers), 1):
        folder = value['folder']
        if (gcf, folder) in seen_paths:
            raise ValueError('Sample output path collision')
        seen_paths.add((gcf, folder))
        groups[gcf]['samples'].append(value)
        if source_validator is not None:
            hashed_bytes += value['source_stat']['bytes']
        if progress and (row_number % 1000 == 0 or row_number == len(rows)):
            print(f'PLAN VALIDATED {row_number}/{len(rows)} gcfs={len(groups)} '
                  f'workers={validation_workers} hashed_bytes={hashed_bytes} '
                  f'elapsed={time.monotonic()-started:.1f}s', file=sys.stderr, flush=True)
    del seen_paths
    groups = {gcf: dict(group, samples=sorted(group['samples'], key=lambda r: r['sample']))
              for gcf, group in sorted(groups.items())}
    code = {name: digest(CODE_ROOT / name) for name in CODE_FILES}
    plan = dict(schema='gcf-native-database-2' if config.get('storage_mode') == 'native' else SCHEMA, cram_index_policy='optional-source-cram-crai-v1', config=config, groups=groups, code_sha256=code,
                input_manifests=[dict(path=str(p), sha256=hashlib.sha256(snapshots[p]).hexdigest()) for p in (records, references)],
                counts=dict(gcfs=len(groups), samples=sample_count, sample_gcf_pairs=len(rows)),
                columns=columns)
    outdir.mkdir(parents=True, exist_ok=False)
    (outdir / 'inputs').mkdir()
    for original, name in ((records, 'sample_manifest.csv'), (references, 'reference_manifest.csv')):
        (outdir / 'inputs' / name).write_bytes(snapshots[original])
    if not publish_plan:
        return plan
    atomic_json(outdir / 'database.json', plan)
    (outdir / 'PLAN.sha256').write_text(digest(outdir / 'database.json') + '\n')
    if config.get('execution_model') == 'nextflow-chunks':
        return plan
    for gcf, group in groups.items():
        directory = outdir / gcf
        (directory / 'samples').mkdir(parents=True)
        atomic_json(directory / 'reference.json', {k: v for k, v in group.items() if k != 'samples'})
        for sample in group['samples']:
            (directory / 'samples' / sample['folder']).mkdir()
    generate_scripts(outdir, plan)
    return plan


def load_plan(database, check_code=False):
    database = Path(database).resolve()
    if digest(database / 'database.json') != (database / 'PLAN.sha256').read_text().strip():
        raise ValueError('Database plan checksum mismatch')
    plan = json.loads((database / 'database.json').read_text())
    if plan['schema'] not in (SCHEMA, 'gcf-native-database-2'):
        raise ValueError('Unsupported database schema')
    if check_code and any(digest(CODE_ROOT / name) != expected for name, expected in plan['code_sha256'].items()):
        raise ValueError('Source code changed; use a new database plan')
    return database, plan


def worker_command(database, gcf, plan):
    return [plan['config']['python'], '-B', str(CODE_ROOT / 'bin/ref_depth_database.py'),
            'worker', '--database', str(Path(database).resolve()), '--gcf', gcf]


def generate_scripts(database, plan):
    for gcf in plan['groups']:
        directory = Path(database) / gcf
        command = shlex.join(worker_command(database, gcf, plan))
        log = shlex.quote(str(directory / 'worker.log'))
        memory = plan['config'].get('memory')
        partition = plan['config'].get('partition') or 'cpu'
        directives = (f'#SBATCH --mem={memory}\n' if memory else '') + (f'#SBATCH --partition={partition}\n' if partition else '')
        text = (f"#!/usr/bin/env bash\n#SBATCH --cpus-per-task={plan['config']['cpus']}\n"
                f"{directives}set -euo pipefail\n"
                f"exec >> {log} 2>&1\nexec {command}\n")
        (directory / 'job.sbatch.sh').write_text(text)


def required_sample_payloads(attempt, result):
    """Declared content closure, independent of directory enumeration."""
    from bin.reference_store import payload_file, storage_payload_names
    attempt = Path(attempt)
    required = {'result.json', 'edges.npy', 'geometry/manifest.json', 'geometry/storage.json',
                'geometry/report.json', 'geometry/checksums.json'}
    storage = json.loads(payload_file(attempt, 'geometry/storage.json').read_text())
    required.update('geometry/' + name for name in storage_payload_names(storage))
    if result.get('storage_mode') == 'native':
        required.add('native_coverage.npz')
    for row in result.get('depths', []):
        if row['status'] != 'READY':
            continue
        required.update((row['matrix'], row['reference_metadata']))
        audit_dir = safe_file(attempt, row['matrix']).parent
        audits = list(audit_dir.glob('ids-*.npz'))
        if not audits:
            raise ValueError(f'Missing sampled-ID audit files: {audit_dir}')
        required.update(p.relative_to(attempt).as_posix() for p in audits)
    if 'cram_index' in result:
        required.add('timing.json')
    if result.get('cram_index'):
        required.add('source.cram.crai')
    return required


def validate_sample_payloads(attempt, result, checks, verify=True):
    """Check required content files independently of empty status-marker semantics."""
    from bin.reference_store import payload_file
    attempt = Path(attempt)
    required = required_sample_payloads(attempt, result)
    if result.get('storage_mode') == 'native':
        from bin.native_reference import validate_native
        validate_native(attempt, result)
    if result.get('cram_index'):
        if result.get('cram_index_archive') != 'source.cram.crai':
            raise ValueError('Missing CRAM index archive declaration')
        if digest(safe_file(attempt,'source.cram.crai')) != result['cram_index']['sha256']:
            raise ValueError('CRAM index archive differs from source index record')
    if not required <= checks.keys():
        raise ValueError('Missing required payload checksum: ' + ','.join(sorted(required-checks.keys())))
    for name, expected in checks.items():
        path = payload_file(attempt, name)
        if not os.access(path, os.R_OK) or path.stat().st_size == 0:
            raise ValueError(f'Required payload is unreadable or empty: {path}')
        if verify and digest(path) != expected:
            raise ValueError(f'Corrupted sample payload: {name}')


def completed_sample(root, fingerprint, verify=True):
    for attempt in sorted(Path(root).glob('attempt-*'), reverse=True):
        if not (attempt / 'status/SUCCESS').is_file():
            continue
        if (attempt / 'status/FAILED').exists() or (attempt / 'status/RUNNING').exists():
            raise ValueError('Contradictory successful sample state')
        if (attempt / 'status/exit_code.txt').read_text().strip() != '0':
            raise ValueError('Contradictory successful sample exit status')
        checks = json.loads((attempt / 'checksums.json').read_text())
        if 'result.json' not in checks or digest(attempt / 'result.json') != checks['result.json']:
            raise ValueError('Unverified sample result')
        result = json.loads((attempt / 'result.json').read_text())
        if result['plan_sha256'] != fingerprint:
            raise ValueError('Sample belongs to another plan')
        validate_sample_payloads(attempt, result, checks, verify)
        return attempt, result
    return None


def database_status(database, plan):
    rows = []
    submission = Path(database) / 'submission.json'
    submitted = json.loads(submission.read_text()) if submission.is_file() else {}
    for gcf in plan['groups']:
        path = Path(database) / gcf / 'progress.json'
        progress = json.loads(path.read_text()) if path.is_file() else dict(state='PENDING')
        if gcf in submitted:
            progress.update(submission_state=submitted[gcf]['state'], job_id=submitted[gcf].get('job_id'))
            if progress['state'] == 'PENDING' and submitted[gcf]['state'] == 'SUBMIT_FAILED':
                progress['state'] = 'SUBMIT_FAILED'
        rows.append(dict(gcf=gcf, **progress))
    return rows


def collect(database, verify=True, jobs=4, progress=None):
    if not isinstance(jobs,int) or isinstance(jobs,bool) or jobs < 1:
        raise ValueError('collect-jobs must be a positive integer')
    database, plan = load_plan(database)
    if plan['config'].get('storage_mode') == 'native':
        from bin.native_reference import collect_native
        return collect_native(database, plan, verify, jobs=jobs, progress=progress)
    fingerprint = (database / 'PLAN.sha256').read_text().strip()
    submission = database / 'submission.json'
    submitted = json.loads(submission.read_text()) if submission.is_file() else {}
    samples, depths, gcfs, resource_rows = [], [], [], []
    with exclusive(database / '.catalog.lock'):
        for gcf, group in plan['groups'].items():
            ready, failed, pending, gaps = 0, 0, 0, 0
            gcf_depths = []
            progress_file = database / gcf / 'progress.json'
            progress = json.loads(progress_file.read_text()) if progress_file.is_file() else {}
            reported_failures = {r['sample']: r['error'] for r in progress.get('failures', [])}
            try:
                reference_error = '' if digest(group['fai']) == group['fai_sha256'] else 'Reference FAI changed since plan'
            except OSError as error:
                reference_error = str(error)
            for sample in group['samples']:
                sample_root = database / gcf / 'samples' / sample['folder']
                failure = ''
                try:
                    if reference_error:
                        raise ValueError(reference_error)
                    if file_stat(sample['cram']) != sample['source_stat']:
                        raise ValueError('CRAM stat changed since plan')
                    validate_cram_index(sample)
                    value = completed_sample(sample_root, fingerprint, verify)
                except (ValueError, OSError, KeyError) as error:
                    value, failure = None, str(error)
                if value:
                    attempt, result = value
                    state = 'COMPLETE'
                    manifest = str(attempt / 'result.json')
                    rows = result['depths']
                    resource_rows.append(dict(gcf=gcf, sample=sample['sample'], sample_manifest=manifest,
                        elapsed_seconds=result['elapsed_seconds'], user_cpu_seconds=result['user_cpu_seconds'],
                        system_cpu_seconds=result['system_cpu_seconds'],
                        process_peak_rss_bytes=result['process_peak_rss_bytes'], rss_scope=result['rss_scope'],
                        geometry_bytes=sum(p.stat().st_size for p in (attempt / 'geometry').rglob('*') if p.is_file()),
                        matrix_bytes=sum((attempt / r['matrix']).stat().st_size for r in rows if r['status'] == 'READY'),
                        sample_output_bytes=sum(p.stat().st_size for p in attempt.rglob('*') if p.is_file())))
                    if any(r['status'] != 'READY' for r in rows):
                        gaps += 1
                else:
                    errors = sorted(sample_root.glob('attempt-*/status/FAILED'))
                    failure = failure or reported_failures.get(sample['sample'], '')
                    if submitted.get(gcf, {}).get('state') == 'SUBMIT_FAILED':
                        failure = failure or 'Scheduler submission failed; see submission.json'
                    state = 'FAILED' if failure or errors else 'PENDING'
                    failed += state == 'FAILED'; pending += state == 'PENDING'
                    manifest = ''
                    result = dict(full_mean_depth='', genome_len='')
                    rows = [dict(depth=d, status=state, matrix='') for d in plan['config']['depths']]
                samples.append(dict(gcf=gcf, sample=sample['sample'], state=state, sample_manifest=manifest,
                                    cram=sample['cram'], full_mean_depth=result['full_mean_depth'], error=failure))
                for row in rows:
                    matrix = str(attempt / row['matrix']) if row['status'] == 'READY' else ''
                    if row['status'] == 'READY':
                        ready += 1
                    entry = dict(gcf=gcf, sample=sample['sample'], depth=row['depth'], status=row['status'],
                        matrix_path=matrix, sample_manifest=manifest, genome_len=result['genome_len'],
                        full_mean_depth=result['full_mean_depth'], replicates=plan['config']['replicates'],
                        bins=plan['config']['bins'])
                    depths.append(entry); gcf_depths.append(entry)
            state = 'FAILED' if failed else 'PENDING' if pending else 'COMPLETE_WITH_GAPS' if gaps else 'COMPLETE'
            index_path = database / gcf / 'DEPTH_INDEX.tsv'
            columns = ['gcf', 'sample', 'depth', 'status', 'matrix_path', 'sample_manifest', 'genome_len', 'full_mean_depth', 'replicates', 'bins']
            atomic_table(index_path, columns, gcf_depths)
            gcfs.append(dict(gcf=gcf, state=state, sample_count=len(group['samples']),
                             ready_sample_depths=ready, failed_samples=failed, pending_samples=pending,
                             samples_with_depth_gaps=gaps, gcf_index=str(index_path)))
        atomic_table(database / 'GCF_INDEX.tsv', list(gcfs[0]), gcfs)
        atomic_table(database / 'SAMPLE_INDEX.tsv', list(samples[0]), samples)
        atomic_table(database / 'DEPTH_INDEX.tsv', columns, depths)
        atomic_table(database / 'READY_DEPTH_INDEX.tsv', columns, [r for r in depths if r['status'] == 'READY'])
        atomic_table(database / 'RESOURCE_INDEX.tsv', ['gcf', 'sample', 'sample_manifest', 'elapsed_seconds',
            'user_cpu_seconds', 'system_cpu_seconds', 'process_peak_rss_bytes', 'rss_scope',
            'geometry_bytes', 'matrix_bytes', 'sample_output_bytes'], resource_rows)
        complete = all(r['state'].startswith('COMPLETE') for r in gcfs)
        summary = dict(schema=SCHEMA, counts=plan['counts'], complete=complete,
                       state='COMPLETE_WITH_GAPS' if complete and any(r['samples_with_depth_gaps'] for r in gcfs)
                             else 'COMPLETE' if complete else 'INCOMPLETE',
                       ready_sample_depths=sum(r['ready_sample_depths'] for r in gcfs),
                       failed_samples=sum(r['failed_samples'] for r in gcfs),
                       pending_samples=sum(r['pending_samples'] for r in gcfs),
                       unsupported_sample_depths=sum(r['status'] == 'UNSUPPORTED_DEPTH' for r in depths),
                       payloads_verified=verify)
        atomic_json(database / 'SUMMARY.json', summary)
    return summary
