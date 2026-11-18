"""Partition a native schema-3 index using an explicit taxonomy manifest."""
import argparse
import csv
from contextlib import closing, contextmanager
import hashlib
import itertools
import json
import os
from pathlib import Path
import sqlite3
import tempfile
import sys
import time
from datetime import datetime


TABLES = {
    'metadata': ('1', 'key'),
    'sources': ('1', 'id'),
    'gcfs': ('gcf IN (SELECT gcf FROM keep)', 'gcf'),
    'references_table': ('gcf IN (SELECT gcf FROM keep)', 'id'),
    'depths': ('gcf IN (SELECT gcf FROM keep)', 'reference_id,depth_key'),
    'reference_sources': ('reference_id IN (SELECT id FROM parent.references_table WHERE gcf IN (SELECT gcf FROM keep))', 'reference_id,source_id'),
    'gcf_statistics': ('gcf IN (SELECT gcf FROM keep)', 'gcf'),
    'skipped_sample_counts': ('gcf IN (SELECT gcf FROM keep)', 'gcf'),
}


def progress(message):
    print(datetime.now().astimezone().isoformat(timespec='seconds') + ' ' + message,
          file=sys.stderr, flush=True)


@contextmanager
def phase(label):
    started = time.perf_counter()
    progress('START ' + label)
    try:
        yield
    except BaseException:
        progress(f'FAIL {label} elapsed_s={time.perf_counter()-started:.3f}')
        raise
    progress(f'DONE {label} elapsed_s={time.perf_counter()-started:.3f}')


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def read_taxonomy(path):
    result = {}
    with Path(path).open(encoding='utf-8-sig', newline='') as handle:
        reader = csv.DictReader(handle)
        if not {'assembly_accession', 'domain'}.issubset(reader.fieldnames or []):
            raise ValueError('Taxonomy requires assembly_accession and domain')
        for row in reader:
            key, domain = row['assembly_accession'], row['domain']
            if not key or not domain or key != key.strip() or domain != domain.strip():
                raise ValueError('Empty or whitespace-padded classification')
            if key in result and result[key] != domain:
                raise ValueError('Conflicting classification: ' + key)
            result[key] = domain
    return result


def write_json(path, value):
    Path(path).write_text(json.dumps(value, sort_keys=True, indent=2) + '\n')


def write_mapping(path, rows):
    with Path(path).open('w', newline='') as handle:
        writer = csv.writer(handle, delimiter='\t', lineterminator='\n')
        writer.writerow(['gcf', 'domain', 'references', 'selected_references'])
        writer.writerows(rows)


def partition(parent, destination, domain, mapping, metadata, provenance):
    connection = sqlite3.connect(destination, uri=True)
    try:
        connection.execute('ATTACH DATABASE ? AS parent', (parent.as_uri() + '?mode=ro',))
        connection.execute('CREATE TEMP TABLE keep(gcf TEXT PRIMARY KEY)')
        connection.executemany('INSERT INTO keep VALUES (?)', [(r[0],) for r in mapping if r[1] == domain])
        for name, (where, order) in TABLES.items():
            definition = connection.execute('SELECT sql FROM parent.sqlite_master WHERE type=\'table\' AND name=?', (name,)).fetchone()
            if definition is None:
                if name == 'skipped_sample_counts':
                    continue
                raise ValueError('Missing table: ' + name)
            with phase(f'{domain} COPY_TABLE {name}'):
                connection.execute(definition[0])
                connection.execute(f'INSERT INTO main.{name} SELECT * FROM parent.{name} WHERE {where} ORDER BY {order}')
        with phase(f'{domain} CREATE_INDEXES'):
            for (definition,) in connection.execute("SELECT sql FROM parent.sqlite_master WHERE type='index' AND sql IS NOT NULL ORDER BY name").fetchall():
                connection.execute(definition)
        # Independent ordered readback compares every retained cell, including IDs,
        # record JSON, hashes, source roots, statistics, selection flags and paths.
        verified = {}
        for name, (where, order) in TABLES.items():
            if name == 'metadata' or not connection.execute("SELECT 1 FROM main.sqlite_master WHERE name=?", (name,)).fetchone():
                continue
            with phase(f'{domain} CHECK_ROWS {name}'):
                left = connection.execute(f'SELECT * FROM main.{name} ORDER BY {order}')
                right = connection.execute(f'SELECT * FROM parent.{name} WHERE {where} ORDER BY {order}')
                n = 0
                for a, b in itertools.zip_longest(left, right):
                    if a != b:
                        raise ValueError('Row mismatch: ' + name)
                    n += 1
                    if n % 50000 == 0:
                        progress(f'{domain} CHECK_ROWS {name} rows={n}')
                verified[name] = n
        count = lambda sql: connection.execute(sql).fetchone()[0]
        copies = count('SELECT COUNT(*) FROM reference_sources')
        references = verified['references_table']
        counts = dict(gcfs=verified['gcfs'], references=references, depths=verified['depths'], ready=references)
        selected = dict(gcfs=count('SELECT COUNT(*) FROM gcfs WHERE selected=1'), references=count('SELECT COUNT(*) FROM references_table WHERE selected=1'))
        updates = dict(counts=counts, selected_counts=selected, input_references=copies,
                       planned_input_references=copies, duplicate_copies=copies-references,
                       taxonomy_partition=dict(domain=domain, **provenance,
                           source_metadata_scope='sources retain original full-build provenance',
                           parent_aggregates={k: metadata.get(k) for k in ('counts', 'selected_counts', 'input_references', 'planned_input_references', 'duplicate_copies')}))
        connection.executemany('INSERT OR REPLACE INTO metadata VALUES (?,?)', [(k, json.dumps(v, sort_keys=True)) for k, v in updates.items()])
        connection.commit()
        with phase(f'{domain} ANALYZE main'):
            connection.execute('ANALYZE main')
            connection.commit()
        with phase(f'{domain} INTEGRITY main'):
            if connection.execute('PRAGMA main.integrity_check').fetchall() != [('ok',)]:
                raise ValueError('Invalid output integrity')
        with phase(f'{domain} FOREIGN_KEYS main'):
            if connection.execute('PRAGMA main.foreign_key_check').fetchone():
                raise ValueError('Invalid output foreign keys')
        return dict(domain=domain, counts=counts, selected_counts=selected, verified_rows=verified, integrity='ok')
    finally:
        connection.close()


def split(index, taxonomy, outdir):
    index, taxonomy = Path(index).resolve(strict=True), Path(taxonomy).resolve(strict=True)
    outdir = Path(outdir).absolute()
    if outdir.exists() or outdir.is_symlink():
        raise FileExistsError(outdir)
    if not outdir.parent.is_dir():
        raise ValueError('Output parent must already exist')
    with phase('HASH_INPUTS'):
        original_hash, taxonomy_hash = sha(index), sha(taxonomy)
    classes = read_taxonomy(taxonomy)
    with closing(sqlite3.connect(index.as_uri() + '?mode=ro', uri=True)) as parent:
        metadata = {k: json.loads(v) for k, v in parent.execute('SELECT key,value FROM metadata')}
        if metadata.get('schema') != 'coverage-reference-sqlite-3' or metadata.get('storage_mode') != 'native':
            raise ValueError('Requires native schema-3 index')
        if metadata.get('skipped_samples', 0) != 0 or metadata.get('partial_index', False):
            raise ValueError('Partial indexes require a separate skipped-audit partition design')
        actual = {r[0] for r in parent.execute("SELECT name FROM sqlite_master WHERE type='table' AND name NOT LIKE 'sqlite_%'")}
        if actual - TABLES.keys() or parent.execute("SELECT 1 FROM sqlite_master WHERE type IN ('view','trigger')").fetchone():
            raise ValueError('Unsupported extra SQL objects')
        with phase('PARENT INTEGRITY main'):
            if parent.execute('PRAGMA main.integrity_check').fetchall() != [('ok',)]:
                raise ValueError('Invalid parent index')
        with phase('PARENT FOREIGN_KEYS main'):
            if parent.execute('PRAGMA main.foreign_key_check').fetchone():
                raise ValueError('Invalid parent foreign keys')
        audit = metadata.get('skipped_samples_audit')
        if audit and (audit.get('rows') != 0 or sha(audit['path']) != audit['sha256']):
            raise ValueError('Invalid inherited empty skipped audit')
        rows = [(g, classes.get(g, 'UNMATCHED'), n, selected) for g, n, selected in parent.execute('SELECT g.gcf,COUNT(r.id),COALESCE(SUM(r.selected),0) FROM gcfs g LEFT JOIN references_table r ON r.gcf=g.gcf GROUP BY g.gcf ORDER BY g.gcf')]
    if any(not any(r[1] == d for r in rows) for d in ('bacteria', 'viral')):
        raise ValueError('Both requested domains must have indexed GCFs')
    provenance = dict(parent_index=str(index), parent_sha256=original_hash, taxonomy=str(taxonomy), taxonomy_sha256=taxonomy_hash,
                      matching='exact assembly_accession; bacteria/viral only', payload_validation='inherited; no payload revalidation')
    with tempfile.TemporaryDirectory(prefix='.taxonomy-split-', dir=outdir.parent) as temp:
        staging = Path(temp)
        report = dict(**provenance, outputs={}, groups={}, taxonomy_not_indexed=sorted(classes.keys()-{r[0] for r in rows}))
        for domain in sorted({r[1] for r in rows}):
            subset = [r for r in rows if r[1] == domain]
            report['groups'][domain] = dict(gcfs=len(subset), references=sum(r[2] for r in subset))
        write_mapping(staging/'gcf_classification.tsv', rows)
        write_mapping(staging/'excluded_gcfs.tsv', [r for r in rows if r[1] not in ('bacteria', 'viral')])
        for domain in ('bacteria', 'viral'):
            name = index.stem + '_' + domain + '.sqlite'
            print('PARTITION ' + domain, flush=True)
            result = partition(index, staging/name, domain, rows, metadata, provenance)
            with phase(f'{domain} HASH_OUTPUT'):
                result.update(path=str(outdir/name), sha256=sha(staging/name), bytes=(staging/name).stat().st_size)
            report['outputs'][domain] = result
        with phase('RECHECK_INPUT_HASHES'):
            if sha(index) != original_hash or sha(taxonomy) != taxonomy_hash:
                raise ValueError('Input changed during split')
        report['complete'] = True
        write_json(staging/'split_audit.json', report)
        # Exclusive output-directory reservation; publish only fully validated files.
        # No overwrite is possible, and SUCCESS marks completion of the whole pair.
        with phase('PUBLISH'):
            outdir.mkdir()
            for path in sorted(staging.iterdir()):
                path.chmod(0o444)
                os.link(path, outdir/path.name)
            (outdir/'SUCCESS').write_text('TAXONOMY_SPLIT_COMPLETE\n')
    progress('COMPLETE ' + str(outdir))
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--index', required=True)
    parser.add_argument('--taxonomy', required=True)
    parser.add_argument('--outdir', required=True, help='New output directory; existing paths refused')
    args = parser.parse_args()
    print(json.dumps(split(args.index, args.taxonomy, args.outdir), sort_keys=True, indent=2))


if __name__ == '__main__':
    main()
