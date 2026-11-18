from contextlib import closing
import importlib.util
import json
from pathlib import Path
import sqlite3
import statistics
import tempfile
import unittest

SPEC=importlib.util.spec_from_file_location('sample_split',Path(__file__).parents[1]/'bin/split_reference_samples.py')
tool=importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(tool)


class SampleSplitTests(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root=Path(self.temp.name)
        self.parent=self.root/'parent.sqlite'
        self.out=self.root/'output'
        self.samples=['s%02d'%n for n in range(10)]
        self.assignment=tool.assign_samples(self.samples,20261103)
        self.rare=next(s for s,g in self.assignment.items() if g=='sample90')
        with closing(sqlite3.connect(self.parent)) as c, c:
            c.executescript('''
                CREATE TABLE metadata(key TEXT PRIMARY KEY,value TEXT NOT NULL);
                CREATE TABLE sources(id INTEGER PRIMARY KEY,root TEXT);
                CREATE TABLE gcfs(gcf TEXT PRIMARY KEY,layout_json TEXT,selected INTEGER);
                CREATE TABLE references_table(id INTEGER PRIMARY KEY,gcf TEXT REFERENCES gcfs(gcf),sample TEXT,selected INTEGER,eligible_depth INTEGER,full_mean_depth REAL,source_id INTEGER REFERENCES sources(id),record_json TEXT);
                CREATE TABLE depths(reference_id INTEGER REFERENCES references_table(id),gcf TEXT,depth_key TEXT,PRIMARY KEY(reference_id,depth_key));
                CREATE TABLE reference_sources(reference_id INTEGER REFERENCES references_table(id),source_id INTEGER REFERENCES sources(id),manifest_path TEXT,PRIMARY KEY(reference_id,source_id));
                CREATE TABLE gcf_statistics(gcf TEXT PRIMARY KEY REFERENCES gcfs(gcf),total_samples INTEGER,eligible_samples INTEGER,selected INTEGER,all_depth_json TEXT,eligible_depth_json TEXT);
                CREATE TABLE skipped_sample_counts(gcf TEXT PRIMARY KEY,copies INTEGER);
                CREATE INDEX selected_depth_lookup ON references_table(gcf,selected,full_mean_depth,sample);
                INSERT INTO sources VALUES (1,'/original/pool');
                INSERT INTO gcfs VALUES ('common','{}',1),('rare','{}',1);
            ''')
            for n,sample in enumerate(self.samples,1):
                c.execute('INSERT INTO references_table VALUES (?,?,?,?,?,?,?,?)',(n,'common',sample,1,1,float(n),1,'{"source":"unchanged"}'))
                c.execute('INSERT INTO reference_sources VALUES (?,1,?)',(n,'original/'+sample))
            c.execute('INSERT INTO references_table VALUES (11,?,?,?,?,?,?,?)',('rare',self.rare,1,1,25.,1,'{"source":"rare"}'))
            c.execute("INSERT INTO reference_sources VALUES (11,1,'original/rare')")
            meta=dict(schema='coverage-reference-sqlite-3',storage_mode='native',selection=dict(min_depth=0,min_samples=1),skipped_samples=0,partial_index=False,counts=dict(gcfs=2,references=11),sources=[dict(id=1,root='/original/pool')])
            c.executemany('INSERT INTO metadata VALUES (?,?)',[(k,json.dumps(v)) for k,v in meta.items()])

    def test_assignment_reproducible_order_independent_and_rounded(self):
        self.assertEqual(tool.assign_samples(self.samples[::-1]+self.samples,20261103),self.assignment)
        self.assertEqual(sum(v=='sample10' for v in tool.assign_samples(['s'+str(n) for n in range(17814)],20261103).values()),1781)
        larger=['s%03d'%n for n in range(100)]
        self.assertNotEqual(tool.assign_samples(larger,42),tool.assign_samples(larger,20261103))
        self.assertEqual(sum(v=='sample10' for v in tool.assign_samples(['s'+str(n) for n in range(15)],1).values()),2)
        with self.assertRaises(ValueError):tool.assign_samples(['one'],1)

    def test_partition_global_disjoint_and_statistics(self):
        before=tool.sha(self.parent)
        result=tool.split(self.parent,self.out)
        self.assertTrue(result['complete'])
        self.assertEqual(result['sample_overlap'],0)
        self.assertEqual(result['gcf_coverage'],dict(both=1,sample10_only=0,sample90_only=1))
        seen={}
        total=[]
        for group,info in result['outputs'].items():
            with closing(sqlite3.connect(info['path'])) as c, c:
                records=c.execute('SELECT * FROM references_table ORDER BY id').fetchall()
                total.extend(records)
                seen[group]={r[2] for r in records}
                for gcf,n,eligible,selected,encoded,_ in c.execute('SELECT * FROM gcf_statistics'):
                    values=[r[5] for r in records if r[1]==gcf]
                    stats=json.loads(encoded)
                    self.assertEqual((n,eligible,selected),(len(values),len(values),1))
                    self.assertEqual(stats['median'],statistics.median(values))
                    self.assertAlmostEqual(stats['mean'],statistics.mean(values))
                self.assertEqual(c.execute('PRAGMA foreign_key_check').fetchall(),[])
                self.assertIsNotNone(c.execute("SELECT sql FROM sqlite_master WHERE name='selected_depth_lookup'").fetchone())
        self.assertFalse(seen['sample10']&seen['sample90'])
        self.assertEqual(seen['sample10']|seen['sample90'],set(self.samples))
        with closing(sqlite3.connect(self.parent)) as c, c:self.assertEqual(sorted(total),c.execute('SELECT * FROM references_table ORDER BY id').fetchall())
        self.assertEqual(tool.sha(self.parent),before)
        self.assertTrue((self.out/'SUCCESS').exists())

    def test_existing_output_and_filtered_input_refused(self):
        self.out.mkdir()
        with self.assertRaises(FileExistsError):tool.split(self.parent,self.out)
        self.out.rmdir()
        with closing(sqlite3.connect(self.parent)) as c, c:c.execute('UPDATE references_table SET selected=0 WHERE id=1')
        with self.assertRaisesRegex(ValueError,'all records selected'):tool.split(self.parent,self.out)
        self.assertFalse(self.out.exists())

    def test_wrong_parent_receipt_refused(self):
        audit=self.root/'audit.json'
        audit.write_text(json.dumps(dict(complete=True,outputs={'bacteria':dict(path=str(self.parent),sha256='wrong',integrity='ok')})))
        with self.assertRaisesRegex(ValueError,'receipt'):tool.split(self.parent,self.out,parent_audit=audit)
        self.assertFalse(self.out.exists())

    def test_all_sqlite_handles_closed_before_cleanup(self):
        original=sqlite3.connect
        connections=[]
        def track(*args,**kwargs):
            connection=original(*args,**kwargs)
            connections.append(connection)
            return connection
        tool.sqlite3.connect=track
        try:
            tool.split(self.parent,self.out)
        finally:
            tool.sqlite3.connect=original
        self.assertTrue(connections)
        for connection in connections:
            with self.assertRaises(sqlite3.ProgrammingError):
                connection.execute('SELECT 1')

    def test_changed_parent_prevents_publication(self):
        original=tool.write_group
        def changing(*args,**kwargs):
            result=original(*args,**kwargs)
            with closing(sqlite3.connect(self.parent)) as c, c:c.execute("INSERT OR REPLACE INTO metadata VALUES ('changed','true')")
            return result
        tool.write_group=changing
        try:
            with self.assertRaisesRegex(ValueError,'Input changed'):tool.split(self.parent,self.out)
        finally:tool.write_group=original
        self.assertFalse(self.out.exists())


if __name__=='__main__':unittest.main()
