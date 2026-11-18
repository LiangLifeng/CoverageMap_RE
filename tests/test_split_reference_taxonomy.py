from contextlib import closing, redirect_stderr
import csv
import importlib.util
import json
import io
from pathlib import Path
import sqlite3
import tempfile
import unittest

SPEC = importlib.util.spec_from_file_location('split_reference_taxonomy', Path(__file__).parents[1]/'bin/split_reference_taxonomy.py')
tool = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(tool)


class TaxonomySplitTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.index = self.root/'all.sqlite'
        self.taxonomy = self.root/'tax.csv'
        self.output = self.root/'split'
        with closing(sqlite3.connect(self.index)) as c, c:
            c.executescript('''
                CREATE TABLE metadata(key TEXT PRIMARY KEY,value TEXT NOT NULL);
                CREATE TABLE sources(id INTEGER PRIMARY KEY,root TEXT);
                CREATE TABLE gcfs(gcf TEXT PRIMARY KEY,layout_json TEXT,selected INTEGER);
                CREATE TABLE references_table(id INTEGER PRIMARY KEY,gcf TEXT REFERENCES gcfs(gcf),sample TEXT,selected INTEGER,full_mean_depth REAL,source_id INTEGER REFERENCES sources(id),record_json TEXT);
                CREATE TABLE depths(reference_id INTEGER REFERENCES references_table(id),gcf TEXT,depth_key TEXT,PRIMARY KEY(reference_id,depth_key));
                CREATE TABLE reference_sources(reference_id INTEGER REFERENCES references_table(id),source_id INTEGER REFERENCES sources(id),manifest_path TEXT,PRIMARY KEY(reference_id,source_id));
                CREATE TABLE gcf_statistics(gcf TEXT PRIMARY KEY REFERENCES gcfs(gcf),total_samples INTEGER);
                CREATE TABLE skipped_sample_counts(gcf TEXT PRIMARY KEY,copies INTEGER);
                CREATE INDEX selected_depth_lookup ON references_table(gcf,selected,full_mean_depth,sample);
                INSERT INTO sources VALUES(1,'/original/source1'),(2,'/original/source2');
            ''')
            for n, gcf in enumerate(['B', 'V', 'A', 'F', 'U'], 1):
                c.execute('INSERT INTO gcfs VALUES (?,?,?)', (gcf, '{"names":["x"]}', int(gcf != 'B')))
                c.execute('INSERT INTO references_table VALUES (?,?,?,?,?,?,?)', (n, gcf, 'sample', int(gcf != 'B'), n*1., 1, '{"untouched":"原始"}'))
                c.execute('INSERT INTO reference_sources VALUES (?,1,?)', (n, 'relative/'+gcf))
                c.execute('INSERT INTO gcf_statistics VALUES (?,1)', (gcf,))
            c.execute("INSERT INTO reference_sources VALUES (1,2,'copy/B')")
            values = dict(schema='coverage-reference-sqlite-3', storage_mode='native', skipped_samples=0,
                          partial_index=False, counts={'gcfs': 5,'references':5}, selected_counts={'gcfs':4,'references':4},
                          selection={'min_depth':0,'min_samples':1}, sources=[{'id':1,'root':'/original/source1'}, {'id':2,'root':'/original/source2'}])
            c.executemany('INSERT INTO metadata VALUES (?,?)', [(k,json.dumps(v)) for k,v in values.items()])
        self.taxonomy.write_text('assembly_accession,domain\nB,bacteria\nV,viral\nA,archaea\nF,fungi\nNOT_IN_DB,bacteria\n')

    def test_partition_and_preservation(self):
        original = tool.sha(self.index)
        report = tool.split(self.index, self.taxonomy, self.output)
        self.assertTrue(report['complete'])
        self.assertEqual(tool.sha(self.index), original)
        self.assertEqual(report['taxonomy_not_indexed'], ['NOT_IN_DB'])
        self.assertEqual(report['outputs']['bacteria']['selected_counts'], {'gcfs':0,'references':0})
        for domain, gcf, identity in [('bacteria','B',1),('viral','V',2)]:
            with closing(sqlite3.connect(report['outputs'][domain]['path'])) as c, c:
                self.assertEqual(c.execute('SELECT id,gcf,record_json FROM references_table').fetchall(), [(identity,gcf,'{"untouched":"原始"}')])
                self.assertEqual(c.execute('PRAGMA foreign_key_check').fetchall(), [])
                self.assertIsNotNone(c.execute("SELECT sql FROM sqlite_master WHERE name='selected_depth_lookup'").fetchone())
                meta={k:json.loads(v) for k,v in c.execute('SELECT * FROM metadata')}
                self.assertEqual(meta['input_references'], 2 if domain=='bacteria' else 1)
                self.assertEqual(meta['duplicate_copies'], 1 if domain=='bacteria' else 0)
                self.assertEqual(meta['selection'], {'min_depth':0,'min_samples':1})
        with (self.output/'excluded_gcfs.tsv').open() as h:
            self.assertEqual([r['gcf'] for r in csv.DictReader(h,delimiter='\t')], ['A','F','U'])
        self.assertTrue((self.output/'SUCCESS').is_file())

    def test_conflict_fails_without_output(self):
        with self.taxonomy.open('a') as h:
            h.write('B,viral\n')
        with self.assertRaisesRegex(ValueError, 'Conflicting'):
            tool.split(self.index, self.taxonomy, self.output)
        self.assertFalse(self.output.exists())

    def test_all_database_handles_close(self):
        original=sqlite3.connect
        connections=[]
        def track(*args,**kwargs):
            connection=original(*args,**kwargs)
            connections.append(connection)
            return connection
        tool.sqlite3.connect=track
        try:
            tool.split(self.index,self.taxonomy,self.output)
        finally:
            tool.sqlite3.connect=original
        for connection in connections:
            with self.assertRaises(sqlite3.ProgrammingError):
                connection.execute('SELECT 1')

    def test_output_integrity_ignores_corrupt_unrelated_attachment(self):
        auxiliary = self.root/'unrelated.sqlite'
        with closing(sqlite3.connect(auxiliary)) as c, c:
            c.executescript("CREATE TABLE unrelated(value TEXT); INSERT INTO unrelated VALUES ('value');")
            page_size = c.execute('PRAGMA page_size').fetchone()[0]
            rootpage = c.execute("SELECT rootpage FROM sqlite_master WHERE name='unrelated'").fetchone()[0]
        data = bytearray(auxiliary.read_bytes())
        data[(rootpage-1)*page_size] = 255  # Invalid b-tree type, schema remains readable.
        auxiliary.write_bytes(data)
        with closing(sqlite3.connect(':memory:')) as c:
            c.execute('ATTACH DATABASE ? AS unrelated', (str(auxiliary),))
            self.assertEqual(c.execute('PRAGMA main.integrity_check').fetchall(), [('ok',)])
            with self.assertRaises(sqlite3.DatabaseError):
                c.execute('PRAGMA integrity_check').fetchall()
        original = sqlite3.connect
        def attach_unrelated(database, *args, **kwargs):
            connection = original(database, *args, **kwargs)
            if not str(database).startswith('file:'):
                connection.execute('ATTACH DATABASE ? AS unrelated', (str(auxiliary),))
            return connection
        tool.sqlite3.connect = attach_unrelated
        log = io.StringIO()
        try:
            with redirect_stderr(log):
                result = tool.split(self.index, self.taxonomy, self.output)
        finally:
            tool.sqlite3.connect = original
        self.assertTrue(result['complete'])
        self.assertEqual(result['outputs']['bacteria']['integrity'], 'ok')
        self.assertIn('DONE bacteria INTEGRITY main', log.getvalue())
        self.assertIn('DONE bacteria CHECK_ROWS references_table', log.getvalue())
        self.assertIn('elapsed_s=', log.getvalue())

    def test_existing_output_refused(self):
        self.output.mkdir()
        marker = self.output/'keep'
        marker.write_text('preserve')
        with self.assertRaises(FileExistsError):
            tool.split(self.index, self.taxonomy, self.output)
        self.assertEqual(marker.read_text(),'preserve')

    def test_partial_and_unknown_schema_fail(self):
        for sql, message in [("UPDATE metadata SET value='1' WHERE key='skipped_samples'", 'Partial'),
                             ('CREATE TABLE unexpected(x)', 'extra SQL')]:
            with closing(sqlite3.connect(self.index)) as c, c:
                c.execute("UPDATE metadata SET value='0' WHERE key='skipped_samples'")
                c.execute(sql)
            with self.assertRaisesRegex(ValueError,message):
                tool.split(self.index,self.taxonomy,self.output)
            self.assertFalse(self.output.exists())

    def test_input_change_prevents_publication(self):
        original = tool.partition
        def changing(*args, **kwargs):
            result = original(*args, **kwargs)
            with self.taxonomy.open('a') as h:
                h.write('EXTRA,bacteria\n')
            return result
        tool.partition = changing
        try:
            with self.assertRaisesRegex(ValueError,'Input changed'):
                tool.split(self.index,self.taxonomy,self.output)
        finally:
            tool.partition = original
        self.assertFalse(self.output.exists())


if __name__ == '__main__':
    unittest.main()
