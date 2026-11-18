"""Create disjoint 10/90 native reference indexes by exact unique sample ID."""
import argparse
import csv
from contextlib import closing
import itertools
import json
import math
import os
from pathlib import Path
import random
import sqlite3
import sys
import tempfile

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from bin.split_reference_taxonomy import TABLES, sha, write_json


def assign_samples(samples, seed):
    ordered = sorted(set(samples))
    if len(ordered) < 2 or any(not isinstance(s, str) or not s for s in ordered):
        raise ValueError('Need at least two nonempty sample IDs')
    shuffled = ordered.copy()
    random.Random(seed).shuffle(shuffled)
    n = max(1, min(len(ordered)-1, (len(ordered)+5)//10))
    small = set(shuffled[:n])
    return {s: 'sample10' if s in small else 'sample90' for s in ordered}


def tsv(path, columns, rows):
    with Path(path).open('w', newline='') as h:
        writer = csv.writer(h, delimiter='\t', lineterminator='\n')
        writer.writerow(columns)
        writer.writerows(rows)


def summarize(connection):
    connection.execute('DELETE FROM gcf_statistics')
    for (gcf,) in connection.execute('SELECT gcf FROM gcfs ORDER BY gcf').fetchall():
        values = [r[0] for r in connection.execute('SELECT full_mean_depth FROM references_table WHERE gcf=? ORDER BY full_mean_depth,sample', (gcf,))]
        if not values or any(not math.isfinite(x) or x < 0 for x in values):
            raise ValueError('Invalid depth group: ' + gcf)
        n = len(values)
        def quantile(q):
            position = (n-1)*q
            lo, hi = math.floor(position), math.ceil(position)
            return values[lo]+(values[hi]-values[lo])*(position-lo)
        mean = connection.execute('SELECT AVG(full_mean_depth) FROM references_table WHERE gcf=?', (gcf,)).fetchone()[0]
        stats = dict(n=n, min=values[0], q25=quantile(.25), median=quantile(.5), q75=quantile(.75), max=values[-1], mean=mean)
        encoded = json.dumps(stats, sort_keys=True)
        connection.execute('INSERT INTO gcf_statistics VALUES (?,?,?,?,?,?)', (gcf,n,n,1,encoded,encoded))


def write_group(index, destination, samples, group, metadata, provenance):
    connection = sqlite3.connect(destination, uri=True)
    try:
        connection.execute('ATTACH DATABASE ? AS parent', (index.as_uri()+'?mode=ro',))
        connection.execute('CREATE TEMP TABLE keep(sample TEXT PRIMARY KEY)')
        connection.executemany('INSERT INTO keep VALUES (?)', [(s,) for s in samples])
        refs = 'sample IN (SELECT sample FROM keep)'
        ids = f'SELECT id FROM parent.references_table WHERE {refs}'
        gcfs = f'gcf IN (SELECT gcf FROM parent.references_table WHERE {refs})'
        predicates = dict(metadata='1', sources='1', gcfs=gcfs, references_table=refs,
                          depths=f'reference_id IN ({ids})', reference_sources=f'reference_id IN ({ids})',
                          skipped_sample_counts='1')
        for table in TABLES:
            definition = connection.execute("SELECT sql FROM parent.sqlite_master WHERE type='table' AND name=?",(table,)).fetchone()
            if definition is None:
                if table == 'skipped_sample_counts':
                    continue
                raise ValueError('Missing table: '+table)
            connection.execute(definition[0])
            if table != 'gcf_statistics':
                connection.execute(f'INSERT INTO main.{table} SELECT * FROM parent.{table} WHERE {predicates[table]} ORDER BY {TABLES[table][1]}')
        for (definition,) in connection.execute("SELECT sql FROM parent.sqlite_master WHERE type='index' AND sql IS NOT NULL ORDER BY name").fetchall():
            connection.execute(definition)
        verified = {}
        for table, where in predicates.items():
            if table == 'metadata' or not connection.execute("SELECT 1 FROM main.sqlite_master WHERE name=?", (table,)).fetchone():
                continue
            order = TABLES[table][1]
            left = connection.execute(f'SELECT * FROM main.{table} ORDER BY {order}')
            right = connection.execute(f'SELECT * FROM parent.{table} WHERE {where} ORDER BY {order}')
            n = 0
            for a, b in itertools.zip_longest(left, right):
                if a != b:
                    raise ValueError('Copied row mismatch: '+table)
                n += 1
            verified[table] = n
        summarize(connection)
        copies = verified['reference_sources']
        n = verified['references_table']
        counts = dict(gcfs=verified['gcfs'],references=n,depths=verified['depths'],ready=n)
        updates = dict(counts=counts,selected_counts=dict(gcfs=counts['gcfs'],references=n),
                       input_references=copies,planned_input_references=copies,duplicate_copies=copies-n,
                       sample_partition=dict(group=group,samples=len(samples),**provenance,
                           parent_aggregates={k:metadata.get(k) for k in ('counts','selected_counts','input_references','planned_input_references','duplicate_copies')},
                           source_metadata_scope='original full-build provenance; taxonomy_partition describes parent scope'))
        connection.executemany('INSERT OR REPLACE INTO metadata VALUES (?,?)',[(k,json.dumps(v,sort_keys=True)) for k,v in updates.items()])
        connection.commit()
        print('VALIDATE '+group, flush=True)
        connection.execute('ANALYZE main')
        connection.commit()
        if connection.execute('PRAGMA main.quick_check').fetchall() != [('ok',)]:
            raise ValueError('Output SQLite quick_check failed')
        if connection.execute('PRAGMA main.foreign_key_check').fetchone():
            raise ValueError('Output foreign key check failed')
        return dict(counts=counts,samples=len(samples),verified_rows=verified,
                    sqlite_quick_check='ok',foreign_key_check='ok',copied_rows='exact ordered comparison')
    finally:
        connection.close()


def split(index, outdir, seed=20261103, parent_audit=None):
    index, outdir = Path(index).resolve(strict=True), Path(outdir).absolute()
    if outdir.exists() or outdir.is_symlink():
        raise FileExistsError(outdir)
    if not outdir.parent.is_dir():
        raise ValueError('Output parent must already exist')
    print('HASH_PARENT', flush=True)
    source_hash = sha(index)
    validated_parent = None
    if parent_audit is not None:
        audit = json.loads(Path(parent_audit).read_text())
        matches = [r for r in audit.get('outputs',{}).values() if Path(r['path']).resolve()==index and r['sha256']==source_hash and r.get('integrity')=='ok']
        if not audit.get('complete') or len(matches)!=1:
            raise ValueError('Parent validation receipt does not match index')
        validated_parent = dict(path=str(Path(parent_audit).resolve()),sha256=sha(parent_audit))
    with closing(sqlite3.connect(index.as_uri()+'?mode=ro',uri=True)) as c:
        meta={k:json.loads(v) for k,v in c.execute('SELECT * FROM metadata')}
        if meta.get('schema')!='coverage-reference-sqlite-3' or meta.get('storage_mode')!='native':
            raise ValueError('Requires native schema-3 index')
        policy=meta.get('selection',{})
        if meta.get('skipped_samples',0)!=0 or meta.get('partial_index') or policy.get('min_depth')!=0 or policy.get('min_samples')!=1:
            raise ValueError('Requires complete unfiltered index (min-depth0/min-samples1)')
        if c.execute('SELECT 1 FROM references_table WHERE selected!=1 OR eligible_depth!=1 LIMIT 1').fetchone() or c.execute('SELECT 1 FROM gcfs WHERE selected!=1 LIMIT 1').fetchone():
            raise ValueError('Requires all records selected and depth-eligible')
        tables={r[0] for r in c.execute("SELECT name FROM sqlite_master WHERE type='table' AND name NOT LIKE 'sqlite_%'")}
        if tables-TABLES.keys() or c.execute("SELECT 1 FROM sqlite_master WHERE type IN ('view','trigger')").fetchone():
            raise ValueError('Unsupported SQL objects')
        if 'skipped_sample_counts' in tables and c.execute('SELECT 1 FROM skipped_sample_counts WHERE copies!=0 LIMIT 1').fetchone():
            raise ValueError('Nonempty skipped counts')
        audit=meta.get('skipped_samples_audit')
        if audit and (audit.get('rows')!=0 or sha(audit['path'])!=audit['sha256']):
            raise ValueError('Invalid inherited empty skipped audit')
        if validated_parent is None and c.execute('PRAGMA main.quick_check').fetchall()!=[('ok',)]:
            raise ValueError('Parent SQLite quick_check failed')
        assignment=assign_samples([r[0] for r in c.execute('SELECT DISTINCT sample FROM references_table ORDER BY sample')],seed)
        original_n=c.execute('SELECT COUNT(*) FROM references_table').fetchone()[0]
    provenance=dict(parent_index=str(index),parent_sha256=source_hash,seed=seed,ratio='1:9',unit='exact global sample ID',
                    rounding='nearest integer 10% (half-up), minimum1 per group',algorithm='Python random.Random.shuffle over sorted unique sample IDs',
                    validated_parent=validated_parent,payload_validation='inherited; index-only split')
    with tempfile.TemporaryDirectory(prefix='.sample-split-',dir=outdir.parent) as temp:
        staging=Path(temp)
        tsv(staging/'sample_assignment.tsv',['sample','group'],assignment.items())
        report=dict(**provenance,outputs={},total_samples=len(assignment),total_references=original_n)
        for group in ('sample10','sample90'):
            samples=[s for s,g in assignment.items() if g==group]
            (staging/(group+'.samples.txt')).write_text('\n'.join(samples)+'\n')
            name=index.stem+'_'+group+'.sqlite'
            print('COPY '+group,flush=True)
            info=write_group(index,staging/name,samples,group,meta,provenance)
            info.update(path=str(outdir/name),sha256=sha(staging/name),bytes=(staging/name).stat().st_size)
            report['outputs'][group]=info
        small=set(s for s,g in assignment.items() if g=='sample10')
        large=set(s for s,g in assignment.items() if g=='sample90')
        if small&large or small|large!=set(assignment) or sum(r['counts']['references'] for r in report['outputs'].values())!=original_n:
            raise ValueError('Partition conservation failed')
        coverage={}
        for group in ('sample10','sample90'):
            info=report['outputs'][group]
            with closing(sqlite3.connect(staging/Path(info['path']).name)) as c:
                for gcf,n in c.execute('SELECT gcf,COUNT(*) FROM references_table GROUP BY gcf'):
                    coverage.setdefault(gcf,dict(sample10=0,sample90=0))[group]=n
        tsv(staging/'gcf_group_counts.tsv',['gcf','sample10_records','sample90_records'],[(g,counts['sample10'],counts['sample90']) for g,counts in sorted(coverage.items())])
        report['gcf_coverage']=dict(both=sum(v['sample10']>0 and v['sample90']>0 for v in coverage.values()),
                                   sample10_only=sum(v['sample90']==0 for v in coverage.values()),sample90_only=sum(v['sample10']==0 for v in coverage.values()))
        print('RECHECK_PARENT',flush=True)
        if sha(index)!=source_hash or (validated_parent and sha(parent_audit)!=validated_parent['sha256']):
            raise ValueError('Input changed while splitting')
        report.update(complete=True,sample_overlap=0)
        write_json(staging/'sample_split_audit.json',report)
        outdir.mkdir()
        for path in sorted(staging.iterdir()):
            path.chmod(0o444)
            os.link(path,outdir/path.name)
        (outdir/'SUCCESS').write_text('SAMPLE_SPLIT_COMPLETE\n')
    return report


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--index',required=True)
    parser.add_argument('--outdir',required=True)
    parser.add_argument('--seed',type=int,default=20261103)
    parser.add_argument('--parent-audit',help='Prior taxonomy audit; must match actual input SHA and integrity=ok')
    args=parser.parse_args()
    print(json.dumps(split(args.index,args.outdir,args.seed,args.parent_audit),sort_keys=True,indent=2))


if __name__=='__main__':
    main()
