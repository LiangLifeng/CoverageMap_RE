"""Immutable multi-root reference registry and metadata-only quality selections."""
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import sqlite3
import tempfile
import sys

from bin.coverage_core import layout_metadata, merge_references, reference_from_layout
from bin.ref_depth_catalog import digest
from bin.reference_index import ReferenceIndex, _build_single, compact, validate_jobs
from bin.reference_audit import (prepare_external, preflight_destination, iter_skipped,
                                 has_table, validate_external, COLUMNS as SKIPPED_COLUMNS)

SCHEMA = 'coverage-reference-sqlite-3'


def validate_thresholds(min_depth, min_samples):
    if not math.isfinite(min_depth) or min_depth < 0:
        raise ValueError('min-depth must be finite and nonnegative')
    if isinstance(min_samples, bool) or not isinstance(min_samples, int) or min_samples < 1:
        raise ValueError('min-samples must be a positive integer')


def metadata(connection):
    return {k: json.loads(v) for k, v in connection.execute('SELECT key,value FROM metadata')}


def set_metadata(connection, values):
    connection.executemany('INSERT OR REPLACE INTO metadata VALUES (?,?)',
                           [(k, compact(v)) for k, v in values.items()])


def scientific_config(config):
    # Operational resource/path/codec changes do not change comparison semantics.
    keys = ('storage_mode', 'policies', 'window_policy', 'window_size', 'bin_count',
            'bins', 'depths', 'replicates', 'seed', 'sampler')
    return {k: config.get(k) for k in keys}


def duplicate_signature(root, manifest_path):
    """Only collapse copies with matching origin, semantics AND numerical file hashes.

    Different geometry encodings are conservatively a conflict; never infer
    equality from depth alone. Operational timing/path fields are excluded.
    """
    manifest = Path(root)/manifest_path
    record = json.loads(manifest.read_text())
    checks = json.loads((manifest.parent/'checksums.json').read_text())
    fields = ('gcf','sample','storage_mode','layout','genome_len','units','blocks','aligned_bases',
              'full_mean_depth','policies','bins','window_size','native','seed','replicates','rng_policy','depths')
    values = {key:record.get(key) for key in fields}
    # Optional CRAI archival differs across otherwise identical builds. Its SHA
    # remains in each original manifest, but is not a different biological sample.
    values['source_sha'] = record['sources'][0]['sha256']
    values['numerical_payloads'] = {key:value for key,value in checks.items()
                                  if Path(key).suffix in ('.npy','.npz','.bin','.zlib')}
    values['geometry_storage_sha'] = checks.get('geometry/storage.json')
    if not values['numerical_payloads']:
        raise ValueError('Duplicate has no numerical payload identities')
    return values


def depth_stats(connection, gcf, eligible=False):
    clause = ' AND eligible_depth=1' if eligible else ''
    where = ' FROM references_table WHERE gcf=?' + clause
    n, lo, mean, hi = connection.execute(
        'SELECT COUNT(*),MIN(full_mean_depth),AVG(full_mean_depth),MAX(full_mean_depth)' + where, (gcf,)).fetchone()
    def quantile(q):
        if not n:
            return None
        position = (n - 1) * q
        lower, upper = math.floor(position), math.ceil(position)
        rows = [r[0] for r in connection.execute('SELECT full_mean_depth' + where +
                ' ORDER BY full_mean_depth,sample LIMIT ? OFFSET ?', (gcf, upper-lower+1, lower))]
        return rows[0] + (rows[-1]-rows[0]) * (position-lower)
    return dict(n=n, min=lo, q25=quantile(.25), median=quantile(.5), q75=quantile(.75), max=hi, mean=mean)


def apply_selection(connection, min_depth, min_samples):
    validate_thresholds(min_depth, min_samples)
    connection.execute('UPDATE references_table SET eligible_depth=(full_mean_depth>=?),selected=0', (min_depth,))
    connection.execute('DELETE FROM gcf_statistics')
    for (gcf,) in connection.execute('SELECT gcf FROM gcfs ORDER BY gcf').fetchall():
        total = depth_stats(connection, gcf)
        eligible = depth_stats(connection, gcf, True)
        selected = eligible['n'] >= min_samples
        reason = '' if selected else 'INSUFFICIENT_SAMPLES_AFTER_DEPTH_FILTER'
        connection.execute('INSERT INTO gcf_statistics VALUES (?,?,?,?,?,?)',
                           (gcf, total['n'], eligible['n'], int(selected), compact(total), compact(eligible)))
        connection.execute('UPDATE gcfs SET selected=? WHERE gcf=?', (int(selected), gcf))
        connection.execute('''UPDATE references_table SET selected=(eligible_depth AND ?),
            exclusion_reason=CASE WHEN eligible_depth=0 THEN 'BELOW_MIN_DEPTH' ELSE ? END WHERE gcf=?''',
            (int(selected), reason, gcf))
    counts = dict(gcfs=connection.execute('SELECT COUNT(*) FROM gcfs WHERE selected=1').fetchone()[0],
                  references=connection.execute('SELECT COUNT(*) FROM references_table WHERE selected=1').fetchone()[0])
    set_metadata(connection, dict(selection=dict(min_depth=min_depth, min_samples=min_samples,
        depth_metric='full_mean_depth', comparison='>=', sample_count='distinct sample IDs after depth filter'),
        selected_counts=counts))


def publish(connection, temporary, output, skipped_tsv=None):
    staged_audit = None
    linked_audit = False
    try:
        if skipped_tsv is not None:
            staged_audit = prepare_external(connection, skipped_tsv)
        connection.commit()
        if staged_audit is not None:
            connection.execute('VACUUM')  # Purge deleted per-sample error text, not just its SQL table.
        connection.execute('ANALYZE')
        connection.commit()
        if connection.execute('PRAGMA integrity_check').fetchone()[0] != 'ok':
            raise ValueError('SQLite integrity check failed')
        if connection.execute('PRAGMA foreign_key_check').fetchone():
            raise ValueError('SQLite foreign key check failed')
        values = metadata(connection)
        connection.close()
        index_sha = digest(temporary)
        temporary.chmod(0o444)
        if staged_audit is not None:
            staged_audit.chmod(0o444)
            os.link(staged_audit, skipped_tsv)  # Exclusive; never replace the user's audit.
            linked_audit = True
        try:
            os.link(temporary, output)
        except BaseException:
            if linked_audit:
                try:
                    if os.path.samefile(staged_audit, skipped_tsv):
                        skipped_tsv.unlink()
                except FileNotFoundError:
                    pass
            raise
        return dict(index=str(output), sha256=index_sha, **values)
    finally:
        if staged_audit is not None:
            staged_audit.unlink(missing_ok=True)


def new_output(output):
    if Path(output).is_symlink():
        raise FileExistsError(output)
    output = Path(output).resolve()
    if output.exists() or output.is_symlink():
        raise FileExistsError(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    return output


def build_registry(databases, output, verify_payloads=False, min_depth=0, min_samples=1, skip_invalid_samples=False, skipped_samples_tsv=None, jobs=1):
    validate_jobs(jobs)
    validate_thresholds(min_depth, min_samples)
    databases = [databases] if isinstance(databases, (str, Path)) else list(databases)
    roots = sorted({str(Path(p).resolve()) for p in databases})
    if not roots or len(roots) != len(databases):
        raise ValueError('Supply one or more distinct database roots (no aliases/duplicates)')
    output = new_output(output)
    skipped_tsv = preflight_destination(output, skipped_samples_tsv, roots)
    if any(output.is_relative_to(Path(root)) for root in roots):
        raise ValueError('Index output must be outside source database directories')
    with tempfile.TemporaryDirectory(prefix='.reference-registry-', dir=output.parent) as directory:
        temporary = Path(directory) / 'registry.sqlite'
        connection = sqlite3.connect(temporary)
        try:
            sources, total_counts, config, mode = [], dict(gcfs=0, references=0, depths=0, ready=0), None, None
            known_layouts = {}
            input_references = planned_input_references = 0
            for number, root in enumerate(roots, 1):
                print(f'INDEX database={number}/{len(roots)} root={root}', file=sys.stderr, flush=True)
                single = Path(directory) / f'source-{number}.sqlite'
                info = _build_single(root, single, verify_payloads, skip_invalid_samples, jobs=jobs)
                incoming = sqlite3.connect(single.as_uri() + '?mode=ro', uri=True)
                try:
                    if config is None:
                        config, mode = info['config'], info['storage_mode']
                        incoming.backup(connection)
                        connection.executescript('''
                            ALTER TABLE references_table ADD COLUMN source_id INTEGER NOT NULL DEFAULT 1;
                            ALTER TABLE references_table ADD COLUMN eligible_depth INTEGER NOT NULL DEFAULT 1;
                            ALTER TABLE references_table ADD COLUMN selected INTEGER NOT NULL DEFAULT 1;
                            ALTER TABLE references_table ADD COLUMN exclusion_reason TEXT NOT NULL DEFAULT '';
                            ALTER TABLE gcfs ADD COLUMN selected INTEGER NOT NULL DEFAULT 1;
                            CREATE TABLE sources(id INTEGER PRIMARY KEY, root TEXT UNIQUE NOT NULL,
                                plan_sha256 TEXT NOT NULL, summary_sha256 TEXT NOT NULL, config_json TEXT NOT NULL);
                            CREATE TABLE reference_sources(reference_id INTEGER NOT NULL REFERENCES references_table(id),
                                source_id INTEGER NOT NULL REFERENCES sources(id), manifest_path TEXT NOT NULL,
                                manifest_sha TEXT NOT NULL, checksums_sha TEXT NOT NULL,
                                PRIMARY KEY(reference_id,source_id));
                            CREATE TABLE gcf_statistics(gcf TEXT PRIMARY KEY REFERENCES gcfs(gcf),
                                total_samples INTEGER, eligible_samples INTEGER, selected INTEGER,
                                all_depth_json TEXT, eligible_depth_json TEXT);
                            CREATE INDEX selected_depth_lookup ON references_table(gcf,selected,full_mean_depth,sample);
                        ''')
                        connection.execute('''INSERT INTO reference_sources SELECT id,1,manifest_path,manifest_sha,checksums_sha FROM references_table''')
                    else:
                        if mode != info['storage_mode'] or scientific_config(config) != scientific_config(info['config']):
                            raise ValueError('Incompatible database scientific configurations; index separately')
                        offset = connection.execute('SELECT COALESCE(MAX(id),0) FROM references_table').fetchone()[0]
                        for gcf, layout in incoming.execute('SELECT gcf,layout_json FROM gcfs ORDER BY gcf'):
                            old = connection.execute('SELECT layout_json FROM gcfs WHERE gcf=?', (gcf,)).fetchone()
                            if old and json.loads(old[0]) != json.loads(layout):
                                raise ValueError(f'Conflicting GCF layout: {gcf}')
                            connection.execute('INSERT OR IGNORE INTO gcfs(gcf,layout_json) VALUES (?,?)', (gcf, layout))
                        for row in incoming.execute('SELECT * FROM skipped_samples ORDER BY gcf,sample'):
                            connection.execute('INSERT INTO skipped_samples VALUES (?,?,?,?,?,?)',row)
                        added = set()
                        for row in incoming.execute('SELECT * FROM references_table ORDER BY gcf,sample'):
                            old = connection.execute('SELECT id,source_id,manifest_path FROM references_table WHERE gcf=? AND sample=?', (row[1],row[2])).fetchone()
                            if old:
                                original_root = sources[old[1]-1]['root']
                                if duplicate_signature(original_root, old[2]) != duplicate_signature(root, row[3]):
                                    raise ValueError(f'Conflicting duplicate GCF/sample: {row[1]}/{row[2]}; numerical payload or provenance differs')
                                reference_id = old[0]
                            else:
                                reference_id = row[0]+offset
                                values = (reference_id, *row[1:], number, 1, 1, '')
                                connection.execute('INSERT INTO references_table VALUES (' + ','.join('?' for _ in values) + ')', values)
                                added.add(row[0])
                            connection.execute('INSERT INTO reference_sources VALUES (?,?,?,?,?)',
                                (reference_id,number,row[3],row[4],row[5]))
                        for row in incoming.execute('SELECT * FROM depths ORDER BY reference_id,depth_key'):
                            if row[0] in added:
                                connection.execute('INSERT INTO depths VALUES (?,?,?,?,?,?,?,?,?,?,?)', (row[0]+offset, *row[1:]))
                    # Check known M5 across ALL references, including ones later excluded.
                    for gcf, encoded in incoming.execute('SELECT gcf,record_json FROM references_table ORDER BY gcf,sample'):
                        record = json.loads(encoded)
                        if not math.isfinite(record['full_mean_depth']) or record['full_mean_depth'] < 0:
                            raise ValueError('Invalid full_mean_depth')
                        layout = record['layout']
                        ref = reference_from_layout(layout)
                        known_layouts[gcf] = merge_references(known_layouts[gcf], ref) if gcf in known_layouts else ref
                    source = dict(id=number, root=root, plan_sha256=info['plan_sha256'], summary_sha256=info['summary_sha256'],
                        source_summary_complete=info['source_summary_complete'],
                        source_summary_matches_index=info['source_summary_matches_index'],
                        planned_references=info['planned_references'], skipped_samples=info['skipped_samples'],
                        partial_index=info['partial_index'])
                    sources.append(source)
                    connection.execute('INSERT INTO sources VALUES (?,?,?,?,?)',
                        (number, root, source['plan_sha256'], source['summary_sha256'], compact(info['config'])))
                    input_references += info['counts']['references']
                    planned_input_references += info['planned_references']
                finally:
                    incoming.close()
                    single.unlink()
            total_counts['gcfs'] = connection.execute('SELECT COUNT(*) FROM gcfs').fetchone()[0]
            total_counts['references'] = connection.execute('SELECT COUNT(*) FROM references_table').fetchone()[0]
            total_counts['depths'] = connection.execute('SELECT COUNT(*) FROM depths').fetchone()[0]
            total_counts['ready'] = total_counts['references'] if mode=='native' else connection.execute("SELECT COUNT(*) FROM depths WHERE status='READY'").fetchone()[0]
            for gcf, layout in known_layouts.items():
                connection.execute('UPDATE gcfs SET layout_json=? WHERE gcf=?',
                    (compact(layout_metadata(layout)), gcf))
            connection.execute("DELETE FROM metadata WHERE key IN ('source_summary_complete','source_summary_matches_index','planned_references','summary_sha256')")
            set_metadata(connection, dict(schema=SCHEMA, storage_mode=mode, config=config, counts=total_counts,
                database_root=roots[0] if len(roots)==1 else None, sources=sources,
                duplicate_copies=input_references-total_counts['references'], input_references=input_references,
                planned_input_references=planned_input_references, skip_invalid_samples=skip_invalid_samples,
                skipped_samples=connection.execute('SELECT COUNT(*) FROM skipped_samples').fetchone()[0],
                partial_index=any(source['partial_index'] for source in sources),
                plan_sha256=hashlib.sha256(compact(sources).encode()).hexdigest(), verify_payloads=verify_payloads))
            apply_selection(connection, min_depth, min_samples)
            set_metadata(connection, dict(geometry_capabilities=(dict(schema='coverage-geometry-2',
                sampling_units=['read','fragment'], supplementary_available=True,
                coverage_alignment_policy='primary_only', mate_relationships_available=True)
                if config.get('capture_alignment') else dict(schema='legacy-block-pool',
                    sampling_units=[config.get('policies',{}).get('sampling_unit','read')],
                    supplementary_available=False))))
            for source in sources:
                root = Path(source['root'])
                if digest(root/'database.json') != source['plan_sha256'] or digest(root/'SUMMARY.json') != source['summary_sha256']:
                    raise ValueError('Source database changed while indexing')
            if any(source['partial_index'] for source in sources):
                print(f'INDEX PARTIAL indexed={total_counts["references"]} '
                      f'skipped={planned_input_references-input_references}; see skipped_samples report',
                      file=sys.stderr, flush=True)
            return publish(connection, temporary, output, skipped_tsv)
        finally:
            connection.close()


def filter_registry(index_path, output, min_depth=0, min_samples=1, skipped_samples_tsv=None):
    validate_thresholds(min_depth, min_samples)
    output = new_output(output)
    source_hash = digest(index_path)
    index = ReferenceIndex(index_path)
    try:
        if not index.federated:
            raise ValueError('Rebuild old schema index once with the current builder before filtering')
        skipped_tsv = preflight_destination(output, skipped_samples_tsv, [s['root'] for s in index.metadata.get('sources',[])])
        with tempfile.TemporaryDirectory(prefix='.reference-filter-', dir=output.parent) as directory:
            temporary = Path(directory)/'registry.sqlite'
            connection = sqlite3.connect(temporary)
            try:
                index.connection.backup(connection)
                apply_selection(connection, min_depth, min_samples)
                set_metadata(connection, dict(selection_parent=dict(path=str(Path(index_path).resolve()), sha256=source_hash),
                    selection_validation='indexed-metadata-only; payload verification inherited from build'))
                if digest(index_path) != source_hash:
                    raise ValueError('Input index changed while filtering')
                return publish(connection, temporary, output, skipped_tsv)
            finally:
                connection.close()
    finally:
        index.close()


def export_report(index_path, outdir):
    index = ReferenceIndex(index_path)
    try:
        if not index.federated:
            raise ValueError('Rebuild old schema index once to obtain GCF depth reports')
        validate_external(index.connection)
        outdir = Path(outdir).resolve()
        outdir.mkdir(parents=True, exist_ok=False)
        has_skipped = index.connection.execute("SELECT 1 FROM sqlite_master WHERE type='table' AND name='skipped_samples'").fetchone() is not None
        with (outdir/'gcf_depths.tsv').open('w', newline='') as handle:
            fields = ['gcf','total_samples','eligible_samples','selected_samples','selected','exclusion_reason','skipped_copies']
            fields += [f'{scope}_depth_{key}' for scope in ('all','eligible') for key in ('min','q25','median','q75','max','mean')]
            writer = csv.DictWriter(handle, fieldnames=fields, delimiter='\t'); writer.writeheader()
            for row in index.connection.execute('SELECT * FROM gcf_statistics ORDER BY gcf'):
                value = dict(gcf=row['gcf'], total_samples=row['total_samples'], eligible_samples=row['eligible_samples'],
                    selected_samples=row['eligible_samples'] if row['selected'] else 0, selected=row['selected'],
                    exclusion_reason='' if row['selected'] else 'INSUFFICIENT_SAMPLES_AFTER_DEPTH_FILTER')
                if has_skipped:
                    copies = index.connection.execute('SELECT COUNT(*) FROM skipped_samples WHERE gcf=?', (row['gcf'],)).fetchone()
                elif has_table(index.connection, 'skipped_sample_counts'):
                    copies = index.connection.execute('SELECT copies FROM skipped_sample_counts WHERE gcf=?', (row['gcf'],)).fetchone()
                else:
                    copies = None
                value['skipped_copies'] = copies[0] if copies else 0
                for scope, column in (('all','all_depth_json'), ('eligible','eligible_depth_json')):
                    value.update({f'{scope}_depth_{k}': v for k,v in json.loads(row[column]).items() if k!='n'})
                writer.writerow(value)
        with (outdir/'sample_depths.tsv').open('w', newline='') as handle:
            cursor = index.connection.execute('''SELECT r.gcf,r.sample,s.root AS database_root,r.full_mean_depth,
                r.eligible_depth,r.selected,r.exclusion_reason,r.manifest_path,r.source_sha,
                (SELECT COUNT(*) FROM reference_sources rs WHERE rs.reference_id=r.id) AS source_copies
                FROM references_table r JOIN sources s ON s.id=r.source_id ORDER BY r.gcf,r.sample''')
            writer = csv.writer(handle, delimiter='\t'); writer.writerow([c[0] for c in cursor.description]); writer.writerows(cursor)
        with (outdir/'reference_sources.tsv').open('w', newline='') as handle:
            cursor = index.connection.execute('''SELECT r.gcf,r.sample,s.root AS database_root,rs.manifest_path,
                rs.manifest_sha,rs.checksums_sha,(r.source_id=rs.source_id) AS canonical
                FROM reference_sources rs JOIN references_table r ON r.id=rs.reference_id
                JOIN sources s ON s.id=rs.source_id ORDER BY r.gcf,r.sample,s.root''')
            writer = csv.writer(handle, delimiter='\t'); writer.writerow([c[0] for c in cursor.description]); writer.writerows(cursor)
        with (outdir/'skipped_samples.tsv').open('w', newline='') as handle:
            writer = csv.writer(handle, delimiter='\t')
            writer.writerow(SKIPPED_COLUMNS)
            writer.writerows(iter_skipped(index.connection))
        (outdir/'selection.json').write_text(json.dumps(index.metadata, indent=2, sort_keys=True)+'\n')
        return dict(outdir=str(outdir), **index.metadata['selected_counts'])
    finally:
        index.close()
