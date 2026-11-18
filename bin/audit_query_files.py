#!/usr/bin/env python3
"""Read-only query/SQLite dependency audit; run inside the failing environment."""
import argparse
from collections import Counter, deque
from concurrent.futures import ThreadPoolExecutor
import csv
from datetime import datetime, timezone
import errno
import gzip
import hashlib
import json
import os
from pathlib import Path
import platform
import socket
import sqlite3
import stat
import sys
import time


COLUMNS = ('origin', 'query_id', 'draw_id', 'reference_id', 'sample', 'gcf',
           'selected', 'role', 'query_required', 'path', 'observation', 'utc',
           'status', 'errno', 'error', 'is_file', 'bytes', 'mtime_ns', 'device', 'inode', 'uid', 'gid', 'mode',
           'first_bad_component', 'symlinks', 'component_details', 'resolved_path', 'expected_sha256',
           'actual_sha256', 'historical_error')


def utc():
    return datetime.now(timezone.utc).isoformat()


def sha256(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for block in iter(lambda: f.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def context(path, role, **values):
    return dict(path=str(path), role=role, query_required=1, **values)


def first_bad(path):
    links, details = [], []
    current = Path(path.anchor)
    for part in path.parts[1:]:
        current /= part
        try:
            info = current.lstat()
            details.append(dict(path=str(current), uid=info.st_uid, gid=info.st_gid,
                                mode=oct(stat.S_IMODE(info.st_mode))))
            if stat.S_ISLNK(info.st_mode):
                links.append([str(current), os.readlink(current)])
            current.stat()
        except OSError as error:
            details.append(dict(path=str(current), errno=error.errno, error=str(error)))
            return str(current), json.dumps(links, ensure_ascii=False), json.dumps(details, ensure_ascii=False)
    return '', json.dumps(links, ensure_ascii=False), json.dumps(details, ensure_ascii=False)


def probe(item, expected='', allow_empty=False):
    """Preserve errno: Path.is_file alone masks access errors in some Pythons."""
    row = dict(item, utc=utc(), status='OK', errno='', error='', is_file='',
               expected_sha256=expected, actual_sha256='')
    path = Path(item['path'])
    if not item['path'] or not path.is_absolute():
        return dict(row, status='INVALID_PATH', error='Expected explicit absolute path')
    try:
        row['resolved_path'] = str(path.resolve())
        info = path.stat()
        row.update(bytes=info.st_size, mtime_ns=info.st_mtime_ns, device=info.st_dev, inode=info.st_ino,
                   uid=info.st_uid, gid=info.st_gid,
                   mode=oct(stat.S_IMODE(info.st_mode)), is_file=stat.S_ISREG(info.st_mode))
        if not stat.S_ISREG(info.st_mode):
            return dict(row, status='NOT_FILE')
        with path.open('rb') as f:
            f.read(1)
        if not allow_empty and info.st_size == 0:
            row['status'] = 'EMPTY'
        if expected:
            row['actual_sha256'] = sha256(path)
            if row['actual_sha256'] != expected:
                row['status'] = 'SHA_MISMATCH'
        # Record actual runtime API behaviour separately from the explicit stat.
        row['is_file'] = path.is_file()
        if not row['is_file']:
            row['status'] = 'CHANGED_DURING_PROBE'
    except (OSError, ValueError, RuntimeError) as error:
        code = getattr(error, 'errno', None)
        row.update(status={errno.ENOENT: 'MISSING', errno.EACCES: 'PERMISSION_DENIED',
                           errno.EPERM: 'PERMISSION_DENIED', errno.ENOTDIR: 'NOT_DIRECTORY',
                           errno.ELOOP: 'SYMLINK_LOOP'}.get(code, 'IO_ERROR'),
                   errno=code if code is not None else '', error=str(error))
        try:
            row['is_file'] = path.is_file()
        except (OSError, ValueError, RuntimeError) as other:
            row['is_file'] = type(other).__name__
    if row['status'] != 'OK':
        try:
            row['first_bad_component'], row['symlinks'], row['component_details'] = first_bad(path)
        except (OSError, ValueError, RuntimeError) as error:
            row['first_bad_component'] = str(error)
    return row


def observe(item, args, expected='', allow_empty=False):
    rows = [dict(probe(item, expected, allow_empty), observation=0)]
    if rows[0]['status'] != 'OK':
        for number in range(1, args.rechecks + 1):
            time.sleep(args.recheck_delay)
            rows.append(dict(probe(item, expected, allow_empty), observation=number))
    return rows


def problem(item, status, error):
    return dict(item, utc=utc(), status=status, error=str(error), observation=0)


def read_json(item, args, rows, expected=''):
    observations = observe(item, args, expected)
    rows.extend(observations)
    if observations[-1]['status'] not in ('OK', 'SHA_MISMATCH'):
        rows.append(problem(item, 'UNENUMERATED', 'Metadata unavailable; descendants unknown'))
        return None
    try:
        data = json.loads(Path(item['path']).read_text())
        if not isinstance(data, dict):
            raise ValueError('Expected JSON object')
        return data
    except (OSError, ValueError) as error:
        rows.append(problem(item, 'UNENUMERATED', error))
        return None


def stored_path(root, name):
    root, relative = Path(root), Path(name)
    if not root.is_absolute() or relative.is_absolute() or '..' in relative.parts or not relative.parts:
        raise ValueError('Unsafe stored path: ' + str(name))
    path = root / relative
    if not path.resolve().is_relative_to(root.resolve()):
        raise ValueError('Stored symlink escapes root: ' + str(path))
    return path


def audit_reference(job, args):
    try:
        return reference_files(job, args)
    except (OSError, ValueError, KeyError, TypeError, AttributeError, RuntimeError) as error:
        root, ref, origin = job
        return [problem(context(root, 'reference_enumeration', origin=origin,
                                reference_id=ref.get('id')), 'UNENUMERATED', error)]


def reference_files(job, args):
    root, ref, origin = job
    rows = []
    base = dict(origin=origin, reference_id=ref['id'], sample=ref['sample'],
                gcf=ref['gcf'], selected=ref.get('selected', 1))
    try:
        manifest = stored_path(root, ref['manifest_path'])
        attempt = manifest.parent
    except (ValueError, OSError, RuntimeError) as error:
        return [problem(context(str(root) + '/' + ref['manifest_path'], 'reference_manifest', **base),
                        'UNENUMERATED', error)]
    record = read_json(context(manifest, 'reference_manifest', **base), args, rows, ref['manifest_sha'])
    checks = read_json(context(attempt/'checksums.json', 'reference_checksums', **base),
                       args, rows, ref['checksums_sha'])
    for name in ('status/SUCCESS', 'status/exit_code.txt'):
        rows.extend(observe(context(attempt/name, 'reference_status', **base), args,
                            allow_empty=name.endswith('SUCCESS')))
    try:
        if (attempt/'status/exit_code.txt').read_text().strip() != '0':
            rows.append(problem(context(attempt/'status/exit_code.txt', 'reference_status', **base),
                                'INVALID_STATUS', 'Exit code is not zero'))
        for name in ('FAILED', 'RUNNING'):
            if (attempt/'status'/name).exists():
                rows.append(problem(context(attempt/'status'/name, 'reference_status', **base),
                                    'INVALID_STATUS', 'Conflicting status marker'))
    except OSError:
        pass  # Already recorded by observe.
    required = {'edges.npy', 'geometry/manifest.json', 'geometry/storage.json',
                'geometry/report.json', 'geometry/checksums.json'}
    if record:
        if record.get('storage_mode') == 'native':
            required.add('native_coverage.npz')
        if 'cram_index' in record:
            required.add('timing.json')
        if record.get('cram_index'):
            required.add('source.cram.crai')
        for depth in record.get('depths', []):
            if depth['status'] == 'READY':
                required.update((depth['matrix'], depth['reference_metadata']))
        for source in record.get('sources', []):
            if source.get('path'):
                item = context(source['path'], 'reference_original_provenance', **base)
                item['query_required'] = 0
                rows.extend(observe(item, args))
    storage = read_json(context(attempt/'geometry/storage.json', 'geometry_storage', **base), args, rows)
    if storage:
        try:
            if storage['codec'] == 'packed':
                required.update('geometry/' + c['file'] for c in storage['chunks'])
            elif storage['codec'] == 'npz':
                required.add('geometry/pool.npz')
            elif storage['codec'] == 'npy':
                required.update('geometry/'+n+'.npy' for n in ('starts', 'ends', 'indptr'))
            else:
                raise ValueError('Unknown storage codec')
        except (KeyError, TypeError, ValueError) as error:
            rows.append(problem(context(attempt/'geometry/storage.json', 'geometry_storage', **base),
                                'UNENUMERATED', error))
    geometry_checks = read_json(context(attempt/'geometry/checksums.json', 'geometry_checksums', **base),
                               args, rows)
    names = dict(checks or {})
    if geometry_checks:
        names.update({'geometry/'+n: v for n, v in geometry_checks.items()})
    for name in required:
        names.setdefault(name, '')
    for name in sorted(names):
        item = context(attempt/name, 'reference_payload', **base)
        try:
            item['path'] = str(stored_path(attempt, name))
            expected = names[name] if args.hash_payloads else ''
            rows.extend(observe(item, args, expected))
        except (ValueError, OSError, RuntimeError) as error:
            rows.append(problem(item, 'UNSAFE_PATH', error))
    return rows


def ordered_map(function, jobs, workers):
    with ThreadPoolExecutor(max_workers=workers) as pool:
        pending = deque()
        for job in jobs:
            pending.append(pool.submit(function, job))
            if len(pending) >= workers * 2:
                yield pending.popleft().result()
        while pending:
            yield pending.popleft().result()


def environment(args):
    data = dict(utc=utc(), hostname=socket.gethostname(), platform=platform.platform(),
                uid=os.getuid(), euid=os.geteuid(), gid=os.getgid(), groups=os.getgroups(),
                python=sys.version, executable=sys.executable, cwd=os.getcwd(), argv=sys.argv,
                options=vars(args), pid=os.getpid())
    for name in ('/proc/self/mountinfo', '/proc/self/status'):
        try:
            data[name] = Path(name).read_text()
        except OSError as error:
            data[name] = str(error)
    try:
        data['mount_namespace'] = os.readlink('/proc/self/ns/mnt')
    except OSError as error:
        data['mount_namespace'] = str(error)
    return data


class Report:
    def __init__(self, outdir):
        self.handles = [gzip.open(outdir/name, 'wt', encoding='utf-8', newline='') for name in
                        ('files.tsv.gz', 'issues.tsv.gz', 'historical_missing.tsv.gz')]
        self.writers = [csv.DictWriter(h, COLUMNS, delimiter='\t', lineterminator='\n',
                                      extrasaction='ignore') for h in self.handles]
        for writer in self.writers:
            writer.writeheader()
        self.counts, self.bad_components, self.history_counts = Counter(), Counter(), Counter()
        self.incomplete = self.required_issues = self.historical = 0

    def emit(self, rows, historical=False):
        for row in rows:
            self.writers[2 if historical else 0].writerow(row)
            if historical:
                self.historical += 1
                self.history_counts[row['status']] += 1
                continue
            self.counts[row['status']] += 1
            if row['status'] != 'OK':
                self.writers[1].writerow(row)
                self.required_issues += int(bool(row.get('query_required', 1)))
                self.incomplete += int(row['status'] == 'UNENUMERATED')
                if row.get('first_bad_component'):
                    self.bad_components[row['first_bad_component']] += 1

    def close(self):
        for handle in self.handles:
            handle.close()


def audit_manifest(path, args, report):
    report.emit(observe(context(path, 'query_manifest', origin=path), args))
    count = 0
    with open(path, encoding='utf-8-sig', newline='') as f:
        reader = csv.DictReader(f, delimiter='\t')
        fields = reader.fieldnames or []
        if len(set(fields)) != len(fields) or not {'Sample', 'GCFID', 'Depth', 'CramPath'} <= set(fields):
            raise ValueError('Query manifest requires unique TSV Sample,GCFID,Depth,CramPath columns')
        for count, row in enumerate(reader, 1):
            base = dict(origin=path, query_id=f'q{count:08d}', sample=row.get('Sample'), gcf=row.get('GCFID'))
            if None in row or any(v is None for v in row.values()):
                report.emit([problem(context(path, 'query_manifest_row', **base), 'UNENUMERATED',
                                     'Wrong number of TSV fields')])
                continue
            for name in ('CramPath', 'CraiPath', 'BedGraphPath'):
                if name == 'CramPath' or row.get(name):
                    report.emit(observe(context(row.get(name, ''), name, **base), args))
    if not count:
        raise ValueError('Empty query manifest')
    return count


def audit_index(path, args, report):
    report.emit(observe(context(path, 'sqlite', origin=path), args))
    connection = sqlite3.connect(Path(path).as_uri()+'?mode=ro', uri=True)
    try:
        connection.execute('PRAGMA query_only=ON')
        connection.row_factory = sqlite3.Row
        metadata = {r['key']: json.loads(r['value']) for r in connection.execute('SELECT * FROM metadata')}
        if metadata.get('storage_mode') != 'native':
            raise ValueError('This dependency auditor supports native reference indexes only')
        sources = metadata.get('sources') or [dict(id=0, root=metadata['database_root'])]
        roots = {s['id']: s['root'] for s in sources}
        for root in roots.values():
            for name in ('database.json', 'SUMMARY.json'):
                report.emit(observe(context(Path(root)/name, 'database_metadata', origin=path), args))
        count = 0
        def jobs():
            for ref in connection.execute('SELECT * FROM references_table ORDER BY id'):
                row = dict(ref)
                yield roots[row.get('source_id', 0)], row, path
        for rows in ordered_map(lambda job: audit_reference(job, args), jobs(), args.jobs):
            report.emit(rows)
            count += 1
            if count % 1000 == 0:
                print(f'INDEX {path} references={count} observations={sum(report.counts.values())}', flush=True)
        return count
    finally:
        connection.close()


def audit_history(run, args, report):
    summary = Path(run)/'summary.tsv'
    report.emit(observe(context(summary, 'historical_summary', origin=run), args))
    with summary.open(newline='') as f:
        reader = csv.DictReader(f, delimiter='\t')
        if 'query_id' not in (reader.fieldnames or []):
            raise ValueError('summary.tsv lacks query_id')
        for query in reader:
            qid = query['query_id']
            try:
                if not isinstance(qid, str) or not qid.startswith('q') or not qid[1:].isdecimal():
                    raise ValueError('Invalid query_id')
                directory = Path(run)/'queries'/qid
                rows = []
                audit = read_json(context(directory/'audit.json', 'historical_audit', origin=run,
                                          query_id=qid), args, rows)
                checks = read_json(context(directory/'checksums.json', 'historical_checksums', origin=run,
                                           query_id=qid), args, rows)
                report.emit(rows)
                if audit:
                    saved = audit.get('query', {})
                    for key in ('CramPath', 'CraiPath', 'BedGraphPath'):
                        if saved.get(key):
                            report.emit(observe(context(saved[key], 'saved_'+key, origin=run, query_id=qid,
                                                        sample=saved.get('Sample'), gcf=saved.get('GCFID')), args))
                if checks is None or 'draws.tsv.gz' not in checks:
                    continue  # Early-return queries legitimately have no draw table.
                path = directory/'draws.tsv.gz'
                with gzip.open(path, 'rt', newline='') as handle:
                    draws = csv.DictReader(handle, delimiter='\t')
                    if not {'error', 'reference_id', 'draw_id'} <= set(draws.fieldnames or []):
                        raise ValueError('Draw log lacks required columns')
                    for draw in draws:
                        if None in draw or any(v is None for v in draw.values()):
                            raise ValueError('Malformed draw TSV row')
                        error = draw.get('error', '')
                        if error.startswith('Missing file: '):
                            item = context(error[len('Missing file: '):], 'historical_missing', origin=run,
                                           query_id=qid, draw_id=draw['draw_id'], reference_id=draw['reference_id'],
                                           sample=draw.get('reference_sample'), gcf=draw.get('GCFID'),
                                           historical_error=error)
                            report.emit(observe(item, args), historical=True)
            except (OSError, ValueError, EOFError) as error:
                report.emit([problem(context(Path(run)/'queries'/qid, 'historical_draws', origin=run),
                                     'UNENUMERATED', error)])


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    for name in ('query-run', 'manifest', 'index'):
        p.add_argument('--'+name, action='append', default=[], help='Repeatable absolute path')
    p.add_argument('--outdir', required=True)
    p.add_argument('--jobs', type=int, default=8)
    p.add_argument('--rechecks', type=int, default=0)
    p.add_argument('--recheck-delay', type=float, default=.2)
    p.add_argument('--hash-payloads', action='store_true')
    args = p.parse_args()
    if not (args.query_run or args.manifest or args.index):
        p.error('At least one --query-run, --manifest or --index is required')
    if args.jobs < 1 or args.rechecks < 0 or not 0 <= args.recheck_delay <= 60:
        p.error('jobs >=1, rechecks >=0 and delay in [0,60] required')
    if any(not Path(path).is_absolute() for path in [args.outdir, *args.query_run, *args.manifest, *args.index]):
        p.error('All command-line paths must be absolute')
    return args


def main():
    args = parse_args()
    out = Path(args.outdir)
    out.mkdir(parents=True, exist_ok=False)
    (out/'environment.json').write_text(json.dumps(environment(args), indent=2, sort_keys=True)+'\n')
    report = Report(out)
    manifests, indexes = list(args.manifest), list(args.index)
    totals, errors = {}, []
    try:
        for run in args.query_run:
            rows = []
            value = read_json(context(Path(run)/'run.json', 'frozen_run', origin=run), args, rows)
            report.emit(rows)
            if value:
                try:
                    manifests.append(value['manifest']['path'])
                    indexes.append(value['database']['path'])
                    for key in ('manifest', 'database'):
                        report.emit(observe(context(value[key]['path'], 'frozen_'+key, origin=run), args,
                                            value[key]['sha256']))
                except KeyError as error:
                    errors.append(f'{run}: invalid frozen run: {error}')
        for kind, paths, function in (('manifest', manifests, audit_manifest), ('index', indexes, audit_index)):
            for path in dict.fromkeys(paths):
                try:
                    totals[kind+':'+path] = function(path, args, report)
                except (OSError, ValueError, KeyError, TypeError, sqlite3.Error) as error:
                    report.emit([problem(context(path, kind, origin=path), 'UNENUMERATED', error)])
                    errors.append(f'{path}: {type(error).__name__}: {error}')
        for run in args.query_run:
            try:
                audit_history(run, args, report)
            except (OSError, ValueError) as error:
                report.emit([problem(context(run, 'historical_run', origin=run), 'UNENUMERATED', error)])
                errors.append(f'{run}: {error}')
    finally:
        report.close()
    complete = not (report.incomplete or errors)
    code = 1 if not complete else 2 if report.required_issues else 0
    summary = dict(schema='query-file-audit-1', finished_utc=utc(), enumeration_complete=complete,
                   exit_code=code, counts=dict(report.counts), input_counts=totals,
                   required_issue_observations=report.required_issues, errors=errors,
                   first_bad_components=dict(sorted(report.bad_components.items())),
                   historical_observations=report.historical,
                   historical_current_status=dict(report.history_counts),
                   limitation='Current access only; historical Missing file does not prove build failure. '
                   'Compare audits from the original execution account/node/container/mount namespace.')
    (out/'summary.json').write_text(json.dumps(summary, indent=2, sort_keys=True)+'\n')
    print(json.dumps(summary, sort_keys=True), flush=True)
    return code


if __name__ == '__main__':
    sys.exit(main())
