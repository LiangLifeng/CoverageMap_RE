#!/usr/bin/env python3
"""Publish a read-only SQLite metadata index for an existing depth database."""
import argparse
from collections import deque
from concurrent.futures import ProcessPoolExecutor
from itertools import islice
import multiprocessing
import json
import math
import os
from pathlib import Path
import sqlite3
import sys
import tempfile
import zipfile
import zlib

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from bin.coverage_core import reference_from_layout
from bin.gcf_windows import require_layout_policy, window_policy
from bin.ref_depth_catalog import completed_sample, digest, load_plan, safe_file

SCHEMA = 'coverage-reference-sqlite-1'


def depth_key(depth):
    value = float(depth)
    if not math.isfinite(value) or value <= 0:
        raise ValueError('Depth must be finite and positive')
    return value.hex()


def compact(value):
    return json.dumps(value, sort_keys=True, separators=(',', ':'), allow_nan=False)


def relative_file(root, path):
    path = Path(path).resolve()
    if not path.is_relative_to(root):
        raise ValueError('Index payload escapes database root')
    safe_file(root, path.relative_to(root))
    return path.relative_to(root).as_posix()


INDEX_SAMPLE_ERRORS = (ValueError, OSError, KeyError, TypeError, EOFError, zipfile.BadZipFile, zlib.error)
INDEX_BATCH_SIZE = 16
_INDEX_CONTEXT = None


def validate_jobs(jobs):
    if isinstance(jobs, bool) or not isinstance(jobs, int) or jobs < 1:
        raise ValueError('jobs must be a positive integer')


def _prepare_sample(root, config, fingerprint, verify_payloads, gcf, layout, sample):
    """Read/validate one reference; return SQL-ready values, without SQLite access."""
    native = config.get('storage_mode') == 'native'
    value = completed_sample(root / gcf / 'samples' / sample['folder'], fingerprint,
                             verify=verify_payloads)
    if value is None:
        raise ValueError(f'Missing successful reference: {gcf}/{sample["sample"]}')
    attempt, record = value
    if record['gcf'] != gcf or record['sample'] != sample['sample']:
        raise ValueError('Reference identity differs from plan')
    if (any(record['layout'].get(k) != layout.get(k) for k in ('names','lengths','sort_policy','offsets'))
            or record['policies'] != config['policies']
            or (not native and (record['bins'] != config['bins']
            or record['replicates'] != config['replicates']))
            or (native and (record.get('storage_mode') != 'native'
            or record['window_size'] != config['window_size']
            or record.get('bin_count') != config.get('bin_count')
            or window_policy(record['native']) != window_policy(config)))):
        raise ValueError('Reference geometry/configuration differs from plan')
    if native:
        require_layout_policy(reference_from_layout(record['layout']), window_policy(config))
    if sorted(row['depth'] for row in record.get('depths', [])) != config.get('depths', []):
        raise ValueError('Reference depth grid differs from plan')
    if not math.isfinite(record['full_mean_depth']) or record['full_mean_depth'] < 0:
        raise ValueError('Invalid full_mean_depth')
    checks = json.loads((attempt / 'checksums.json').read_text())
    manifest = attempt / 'result.json'
    reference = (gcf, sample['sample'], relative_file(root, manifest), digest(manifest),
                 digest(attempt / 'checksums.json'), record['sources'][0]['sha256'], compact(record))
    native_fields = None
    if native:
        native_fields = (record['full_mean_depth'], record['genome_len'], record['aligned_bases'], record['units'],
                         record['native']['covered_fraction'], relative_file(root, attempt / 'native_coverage.npz'),
                         (attempt / 'geometry').relative_to(root).as_posix())
    depths = []
    ready = int(native)
    for row in record.get('depths', []):
        matrix, metadata, matrix_sha, metadata_sha, matrix_bytes = '', '', '', '', 0
        if row['status'] == 'READY':
            matrix_file = safe_file(attempt, row['matrix'])
            metadata_file = safe_file(attempt, row['reference_metadata'])
            matrix = relative_file(root, matrix_file)
            metadata = relative_file(root, metadata_file)
            matrix_sha, metadata_sha = checks[row['matrix']], checks[row['reference_metadata']]
            if digest(metadata_file) != metadata_sha or digest(attempt / 'edges.npy') != checks['edges.npy']:
                raise ValueError('Depth metadata checksum mismatch')
            matrix_bytes = matrix_file.stat().st_size
            ready += 1
        elif row['status'] != 'UNSUPPORTED_DEPTH':
            raise ValueError('Unexpected depth state in successful sample')
        depths.append((gcf, depth_key(row['depth']), row['depth'], row['status'],
                       matrix, metadata, matrix_sha, checks['edges.npy'], metadata_sha, matrix_bytes))
    return reference, record['full_mean_depth'], native_fields, depths, ready


def _initialize_index_worker(root, config, fingerprint, verify_payloads):
    global _INDEX_CONTEXT
    _INDEX_CONTEXT = (Path(root), config, fingerprint, verify_payloads)


def _prepare_batch(context, tasks):
    result = []
    for gcf, layout, sample in tasks:
        try:
            prepared = _prepare_sample(*context, gcf, layout, sample)
        except INDEX_SAMPLE_ERRORS as error:
            result.append((gcf, sample, None, error))
        else:
            result.append((gcf, sample, prepared, None))
    return result


def _index_worker(tasks):
    return _prepare_batch(_INDEX_CONTEXT, tasks)


def _prepared_samples(root, plan, fingerprint, verify_payloads, jobs):
    """At most 2*jobs batches in flight; consume in plan order, never completion order."""
    context = (root, plan['config'], fingerprint, verify_payloads)
    tasks = ((gcf, group['layout'], sample) for gcf, group in plan['groups'].items() for sample in group['samples'])
    if jobs == 1:
        for task in tasks:
            yield from _prepare_batch(context, [task])
        return
    batches = iter(lambda: list(islice(tasks, INDEX_BATCH_SIZE)), [])
    pending = deque()
    executor = ProcessPoolExecutor(max_workers=jobs, mp_context=multiprocessing.get_context('spawn'),
                                   initializer=_initialize_index_worker, initargs=context)
    try:
        for batch in islice(batches, 2 * jobs):
            pending.append(executor.submit(_index_worker, batch))
        while pending:
            prepared = pending.popleft().result()
            yield from prepared
            batch = next(batches, None)
            if batch is not None:
                pending.append(executor.submit(_index_worker, batch))
    finally:
        for future in pending:
            future.cancel()
        executor.shutdown(wait=True, cancel_futures=True)


def _build_single(database, output, verify_payloads=False, skip_invalid_samples=False, jobs=1):
    validate_jobs(jobs)
    root, plan = load_plan(database)
    native = plan['config'].get('storage_mode') == 'native'
    fingerprint = digest(root / 'database.json')
    summary_fingerprint = digest(root / 'SUMMARY.json')
    summary = json.loads((root / 'SUMMARY.json').read_text())
    if not summary.get('complete') and not skip_invalid_samples:
        raise ValueError('Index requires a completed reference database; use --skip-invalid-samples to index only valid samples')
    output = Path(output).resolve()
    if output.exists():
        raise FileExistsError(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    handle, temporary = tempfile.mkstemp(prefix='.reference-index-', dir=output.parent)
    os.close(handle)
    temporary = Path(temporary)
    connection = None
    try:
        connection = sqlite3.connect(temporary)
        connection.execute('PRAGMA foreign_keys=ON')
        connection.executescript('''
            CREATE TABLE metadata(key TEXT PRIMARY KEY, value TEXT NOT NULL);
            CREATE TABLE gcfs(gcf TEXT PRIMARY KEY, layout_json TEXT NOT NULL);
            CREATE TABLE skipped_samples(
                database_root TEXT NOT NULL,gcf TEXT NOT NULL,sample TEXT NOT NULL,
                sample_path TEXT NOT NULL,error_type TEXT NOT NULL,error TEXT NOT NULL,
                PRIMARY KEY(database_root,gcf,sample));
            CREATE INDEX skipped_gcf_lookup ON skipped_samples(gcf);
            CREATE TABLE references_table(
                id INTEGER PRIMARY KEY, gcf TEXT NOT NULL REFERENCES gcfs(gcf),
                sample TEXT NOT NULL, manifest_path TEXT NOT NULL UNIQUE,
                manifest_sha TEXT NOT NULL, checksums_sha TEXT NOT NULL,
                source_sha TEXT NOT NULL, record_json TEXT NOT NULL,
                full_mean_depth REAL NOT NULL DEFAULT 0, genome_len INTEGER,
                aligned_bases INTEGER, units INTEGER, covered_fraction REAL,
                coverage_path TEXT, geometry_path TEXT,
                UNIQUE(gcf,sample));
            CREATE INDEX references_lookup ON references_table(gcf,sample);
            CREATE INDEX native_depth_lookup ON references_table(gcf,full_mean_depth,sample);
            CREATE TABLE depths(
                reference_id INTEGER NOT NULL REFERENCES references_table(id),
                gcf TEXT NOT NULL, depth_key TEXT NOT NULL, depth REAL NOT NULL,
                status TEXT NOT NULL, matrix_path TEXT NOT NULL,
                metadata_path TEXT NOT NULL, matrix_sha TEXT NOT NULL,
                edges_sha TEXT NOT NULL, metadata_sha TEXT NOT NULL,
                matrix_bytes INTEGER NOT NULL,
                PRIMARY KEY(reference_id,depth_key));
            CREATE INDEX depth_lookup ON depths(gcf,depth_key,status,reference_id);
        ''')
        counts = dict(gcfs=len(plan['groups']), references=0, depths=0, ready=0)
        processed = skipped = 0
        for gcf, group in plan['groups'].items():
            connection.execute('INSERT INTO gcfs VALUES (?,?)', (gcf, compact(group['layout'])))
        print(f'INDEX PREPARE {root.name} jobs={jobs} backend={"process" if jobs > 1 else "serial"}',
              file=sys.stderr, flush=True)
        prepared_iterator = _prepared_samples(root, plan, fingerprint, verify_payloads, jobs)
        try:
            for gcf, sample, prepared, error in prepared_iterator:
                before = counts.copy()
                connection.execute('SAVEPOINT index_sample')
                try:
                    if error is not None:
                        raise error
                    reference, full_mean_depth, native_fields, depths, ready = prepared
                    connection.execute('''INSERT INTO references_table
                        (gcf,sample,manifest_path,manifest_sha,checksums_sha,source_sha,record_json)
                        VALUES (?,?,?,?,?,?,?)''', reference)
                    reference_id = connection.execute('SELECT last_insert_rowid()').fetchone()[0]
                    connection.execute('UPDATE references_table SET full_mean_depth=? WHERE id=?',
                                       (full_mean_depth, reference_id))
                    if native_fields is not None:
                        connection.execute('''UPDATE references_table SET full_mean_depth=?,genome_len=?,
                            aligned_bases=?,units=?,covered_fraction=?,coverage_path=?,geometry_path=? WHERE id=?''',
                            (*native_fields, reference_id))
                    for row in depths:
                        connection.execute('INSERT INTO depths VALUES (?,?,?,?,?,?,?,?,?,?,?)', (reference_id, *row))
                    counts['depths'] += len(depths)
                    counts['ready'] += ready
                    counts['references'] += 1
                except INDEX_SAMPLE_ERRORS as error:
                    connection.execute('ROLLBACK TO index_sample')
                    connection.execute('RELEASE index_sample')
                    counts = before
                    if not skip_invalid_samples:
                        raise
                    sample_path = root / gcf / 'samples' / sample['folder']
                    connection.execute('INSERT INTO skipped_samples VALUES (?,?,?,?,?,?)',
                        (str(root),gcf,sample['sample'],str(sample_path),type(error).__name__,str(error)))
                    skipped += 1
                    print(f'INDEX SKIP {gcf}/{sample["sample"]}: {type(error).__name__}: {error}',
                          file=sys.stderr, flush=True)
                else:
                    connection.execute('RELEASE index_sample')
                processed += 1
                if processed % 100 == 0:
                    print(f'INDEX {root.name} samples={processed}/{plan["counts"]["sample_gcf_pairs"]} '
                          f'indexed={counts["references"]} skipped={skipped}', file=sys.stderr, flush=True)
        finally:
            prepared_iterator.close()
        if processed != plan['counts']['sample_gcf_pairs'] or counts['references'] + skipped != processed:
            raise ValueError('Plan sample count does not match processed references')
        summary_matches = (counts['references'] == plan['counts']['sample_gcf_pairs']
            and counts['ready'] == summary['ready_samples' if native else 'ready_sample_depths'])
        if not skip_invalid_samples and not summary_matches:
            raise ValueError('Database summary does not match indexed products')
        metadata = dict(schema='coverage-reference-sqlite-2' if native else SCHEMA, storage_mode='native' if native else 'legacy-depths', database_root=str(root), plan_sha256=fingerprint,
                        config=plan['config'], counts=counts, verify_payloads=verify_payloads, summary_sha256=summary_fingerprint,
                        skip_invalid_samples=skip_invalid_samples, skipped_samples=skipped,
                        planned_references=processed, source_summary_complete=bool(summary.get('complete')),
                        source_summary_matches_index=summary_matches,
                        partial_index=bool(skipped or not summary.get('complete') or not summary_matches))
        connection.executemany('INSERT INTO metadata VALUES (?,?)', [(k, compact(v)) for k, v in metadata.items()])
        connection.commit()
        connection.execute('ANALYZE')
        connection.commit()
        if connection.execute('PRAGMA integrity_check').fetchone()[0] != 'ok':
            raise ValueError('SQLite integrity check failed')
        connection.close(); connection = None
        if (digest(root / 'database.json') != fingerprint or digest(root / 'SUMMARY.json') != summary_fingerprint
                or json.loads((root / 'SUMMARY.json').read_text()) != summary):
            raise ValueError('Reference database changed while indexing')
        os.link(temporary, output)  # Exclusive atomic publication; never replace an existing index.
        output.chmod(0o444)
        return dict(index=str(output), sha256=digest(output), **metadata)
    finally:
        if connection is not None:
            connection.close()
        temporary.unlink(missing_ok=True)


class ReferenceIndex:
    def __init__(self, path, database_root=None):
        self.path = Path(path).resolve()
        if not self.path.is_file():
            raise FileNotFoundError(self.path)
        self.connection = sqlite3.connect(self.path.as_uri() + '?mode=ro', uri=True)
        self.connection.row_factory = sqlite3.Row
        try:
            self.metadata = {r['key']: json.loads(r['value']) for r in self.connection.execute('SELECT * FROM metadata')}
            if self.metadata['schema'] not in (SCHEMA, 'coverage-reference-sqlite-2', 'coverage-reference-sqlite-3'):
                raise ValueError('Unsupported reference SQLite schema')
            self.federated = self.metadata['schema'] == 'coverage-reference-sqlite-3'
            self.window_policy = window_policy(self.metadata['config'])
            roots = self.metadata.get('sources', [])
            if database_root is not None and len(roots) > 1:
                raise ValueError('--database-root cannot relocate a multi-database index')
            root = database_root or self.metadata.get('database_root')
            self.root = Path(root).resolve() if root else None
            self.roots = {s['id']: Path(database_root or s['root']).resolve() for s in roots}
        except BaseException:
            self.close(); raise

    def root_for(self, row):
        return self.roots[row['source_id']] if self.federated else self.root

    def close(self):
        self.connection.close()

    def group(self, gcf, missing_ok=False):
        clause = ' AND selected=1' if self.federated else ''
        row = self.connection.execute('SELECT layout_json FROM gcfs WHERE gcf=?' + clause, (gcf,)).fetchone()
        if row is None:
            if missing_ok:
                return None
            raise ValueError(f'No indexed GCF: {gcf}')
        layout = json.loads(row[0])
        if self.metadata.get('storage_mode') == 'native':
            require_layout_policy(reference_from_layout(layout), self.window_policy)
        return layout

    def nearest(self, gcf, depth, limit=1, exclude_sample=None, exclude_source=None, samples=None):
        """Two indexed range seeks; exclusions precede LIMIT and ranking."""
        depth_key(depth)
        if limit < 1:
            raise ValueError('reference-limit must be positive')
        filters = ['gcf=?', 'full_mean_depth>0']
        if self.federated:
            filters.append('selected=1')
        params = [gcf]
        for key, value in (('sample',exclude_sample),('source_sha',exclude_source)):
            if value is not None:
                filters.append(key+'<>?');params.append(value)
        if samples:
            filters.append('sample IN ('+','.join('?' for _ in samples)+')');params.extend(sorted(samples))
        found=[]
        for operator, order in (('>=','ASC'),('<','DESC')):
            lookup = 'selected_depth_lookup' if self.federated else 'native_depth_lookup'
            sql = ('SELECT * FROM references_table INDEXED BY ' + lookup + ' WHERE '
                +' AND '.join(filters)+' AND full_mean_depth'+operator+'? ORDER BY full_mean_depth '+order+',sample ASC LIMIT ?')
            found.extend(dict(row) for row in self.connection.execute(sql, [*params,depth,limit]))
        return sorted(found,key=lambda r:(abs(r['full_mean_depth']-depth),r['sample']))[:limit]

    def candidates(self, gcf, minimum=0., maximum=None, exclude_sample=None, exclude_source=None, samples=None):
        """Positive native originals in an inclusive depth band; exclusions precede sampling."""
        filters = ['gcf=?', 'full_mean_depth>0', 'full_mean_depth>=?']
        params = [gcf, minimum]
        if maximum is not None:
            filters.append('full_mean_depth<=?'); params.append(maximum)
        if self.federated:
            filters.append('selected=1')
        for key, value in (('sample', exclude_sample), ('source_sha', exclude_source)):
            if value is not None:
                filters.append(key+'<>?'); params.append(value)
        if samples:
            filters.append('sample IN ('+','.join('?' for _ in samples)+')'); params.extend(sorted(samples))
        lookup = 'selected_depth_lookup' if self.federated else 'native_depth_lookup'
        return [dict(row) for row in self.connection.execute(
            'SELECT * FROM references_table INDEXED BY '+lookup+' WHERE '+' AND '.join(filters)+' ORDER BY sample', params)]

    def references(self, gcf, depth, batch_size=32):
        # GCF lookup followed by exact reference/depth PK lookup; never scan all depth records.
        selected = ' AND r.selected=1' if self.federated else ''
        cursor = self.connection.execute('''SELECT r.*, d.status AS depth_status,
            d.matrix_path,d.metadata_path,d.matrix_sha,d.edges_sha,d.metadata_sha,d.matrix_bytes
            FROM references_table r INDEXED BY references_lookup LEFT JOIN depths d
            ON d.reference_id=r.id AND d.depth_key=?
            WHERE r.gcf=?''' + selected + ' ORDER BY r.sample', (depth_key(depth), gcf))
        while True:
            rows = cursor.fetchmany(batch_size)
            if not rows:
                return
            for row in rows:
                yield dict(row)


def build_index(database, output, verify_payloads=False, min_depth=0, min_samples=1, skip_invalid_samples=False, skipped_samples_tsv=None, jobs=1):
    validate_jobs(jobs)
    from bin.reference_registry import build_registry
    return build_registry(database, output, verify_payloads, min_depth, min_samples, skip_invalid_samples, skipped_samples_tsv, jobs=jobs)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest='command', required=True)
    build = sub.add_parser('build')
    build.add_argument('--database', type=Path, action='append', required=True,
                       help='Repeat for multiple completed databases')
    build.add_argument('--out', type=Path, required=True)
    build.add_argument('--jobs', type=int, default=1, help='Parallel read/validation worker processes (default: 1); SQLite writes remain ordered')
    build.add_argument('--verify-payloads', action='store_true')
    build.add_argument('--skip-invalid-samples', action='store_true',
                       help='Allow incomplete source summaries; omit invalid samples with an audit record')
    selection = sub.add_parser('filter', help='New immutable selection using existing indexed metadata')
    selection.add_argument('--index', type=Path, required=True)
    selection.add_argument('--out', type=Path, required=True)
    for command in (build, selection):
        command.add_argument('--skipped-samples-tsv', type=Path,
                             help='Store omitted sample details only in this new TSV, outside SQLite')
        command.add_argument('--min-depth', type=float, default=0,
                             help='Minimum original full-genome mean depth, inclusive')
        command.add_argument('--min-samples', type=int, default=1,
                             help='Minimum distinct samples AFTER depth filtering per GCF')
    report = sub.add_parser('report')
    report.add_argument('--index', type=Path, required=True)
    report.add_argument('--outdir', type=Path, required=True)
    inspect = sub.add_parser('inspect')
    inspect.add_argument('--index', type=Path, required=True)
    args = parser.parse_args()
    if args.command == 'build':
        result = build_index(args.database, args.out, args.verify_payloads, args.min_depth, args.min_samples, args.skip_invalid_samples, args.skipped_samples_tsv, jobs=args.jobs)
    elif args.command == 'filter':
        from bin.reference_registry import filter_registry
        result = filter_registry(args.index, args.out, args.min_depth, args.min_samples, args.skipped_samples_tsv)
    elif args.command == 'report':
        from bin.reference_registry import export_report
        result = export_report(args.index, args.outdir)
    else:
        index = ReferenceIndex(args.index)
        try:
            result = index.metadata
        finally:
            index.close()
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == '__main__':
    main()
