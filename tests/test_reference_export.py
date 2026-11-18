"""Read-only Sample–GCF export: filter semantics, provenance and safe publication."""
from contextlib import closing
import csv
import hashlib
import json
from pathlib import Path
import sqlite3
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

from bin.reference_export import COLUMNS, export_index


class ExportTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name)
        self.index = self.root / 'source.sqlite'
        with closing(sqlite3.connect(self.index)) as db, db:
            db.executescript('''CREATE TABLE metadata(key TEXT,value TEXT);
                CREATE TABLE sources(id INTEGER,root TEXT);
                CREATE TABLE references_table(sample TEXT,gcf TEXT,full_mean_depth REAL,
                record_json TEXT,source_id INTEGER,selected INTEGER,eligible_depth INTEGER,exclusion_reason TEXT);''')
            db.execute('INSERT INTO metadata VALUES (?,?)', ('schema', json.dumps('coverage-reference-sqlite-3')))
            db.execute('INSERT INTO metadata VALUES (?,?)', ('selection', json.dumps(dict(min_depth=3,min_samples=2))))
            db.executemany('INSERT INTO sources VALUES (?,?)', [(1,'/not-mounted/IH'),(2,'/not-mounted/MAG')])
            rows = [('B','G1',5,2,1,1,''), ('A','G1',1,1,0,0,'BELOW_MIN_DEPTH'),
                    ('A','G2',3,2,0,1,'INSUFFICIENT_SAMPLES_AFTER_DEPTH_FILTER'),
                    ('C','G1',5.001,1,1,1,''), ("x' OR 1=1 --",'G3',0,1,0,0,'BELOW_MIN_DEPTH')]
            for sample,gcf,depth,source,selected,eligible,reason in rows:
                db.execute('INSERT INTO references_table VALUES (?,?,?,?,?,?,?,?)',
                           (sample,gcf,depth,json.dumps({'sources':[{'path':f'/raw/{sample}.cram'}]}),source,selected,eligible,reason))
        self.before = hashlib.sha256(self.index.read_bytes()).hexdigest()

    def tearDown(self):
        self.assertEqual(hashlib.sha256(self.index.read_bytes()).hexdigest(), self.before)
        self.temp.cleanup()

    def export(self, **kwargs):
        out = self.root / 'result.tsv'
        info = export_index(self.index,out,**kwargs)
        with out.open() as handle:
            reader = csv.DictReader(handle,delimiter='\t')
            self.assertEqual(tuple(reader.fieldnames), COLUMNS)
            rows = list(reader)
        return info, rows

    def test_inclusive_depth_order_roots_and_no_original_files(self):
        info, rows = self.export(minimum=1,maximum=5)
        self.assertEqual((info['records'],info['samples'],info['gcfs']), (3,2,2))
        self.assertEqual([(r['Sample'],r['GCFID']) for r in rows], [('A','G1'),('A','G2'),('B','G1')])
        self.assertEqual([r['DatabaseRoot'] for r in rows], ['/not-mounted/IH','/not-mounted/MAG','/not-mounted/MAG'])
        self.assertEqual(rows[0]['CramFilePath'], '/raw/A.cram')

    def test_repeat_ids_or_within_and_across_dimensions(self):
        info, rows = self.export(samples=['B','A','B'],gcfs=['G1','G3'],minimum=1,maximum=5,selected_only=True)
        self.assertEqual(info['records'],1)
        self.assertEqual(rows[0]['Sample'],'B')
        self.assertEqual(info['filters']['samples'],['A','B'])

    def test_exact_ids_parameterized_and_empty_has_header(self):
        info, rows = self.export(samples=["x' OR 1=1 --"],maximum=0)
        self.assertEqual(info['records'],1)
        self.assertEqual(rows[0]['GCFID'],'G3')
        (self.root/'result.tsv').unlink()
        info, rows = self.export(samples=['a'])
        self.assertEqual(info['records'],0)
        self.assertEqual(rows,[])

    def test_selected_is_not_depth_eligibility(self):
        info, rows = self.export(minimum=3,selected_only=True)
        self.assertEqual(info['records'],2)
        self.assertTrue(all(r['Selected']=='1' for r in rows))
        self.assertNotIn(('A','G2'),[(r['Sample'],r['GCFID']) for r in rows])

    def test_invalid_bounds_and_ids_no_output(self):
        for options in [dict(minimum=-1),dict(maximum=float('inf')),dict(minimum=float('nan')),
                        dict(minimum=5,maximum=1),dict(samples=[' '])]:
            with self.subTest(options=options), self.assertRaises(ValueError):
                export_index(self.index,self.root/'out.tsv',**options)
        self.assertFalse((self.root/'out.tsv').exists())

    def test_existing_file_symlink_and_publication_race_preserved(self):
        out = self.root/'old.tsv';out.write_text('keep')
        with self.assertRaises(FileExistsError):export_index(self.index,out)
        link=self.root/'link.tsv';link.symlink_to(self.root/'absent')
        with self.assertRaises(FileExistsError):export_index(self.index,link)
        with patch('bin.reference_export.os.link',side_effect=FileExistsError('race')):
            with self.assertRaises(FileExistsError):export_index(self.index,self.root/'new.tsv')
        self.assertEqual(out.read_text(),'keep')
        self.assertFalse(list(self.root.glob('.reference-export-*')))

    def mutate(self, sql, parameters=()):
        with closing(sqlite3.connect(self.index)) as db, db:db.execute(sql,parameters)
        self.before=hashlib.sha256(self.index.read_bytes()).hexdigest()

    def test_bad_record_or_source_aborts_without_partial(self):
        self.mutate("UPDATE references_table SET record_json='[]' WHERE sample='B'")
        with self.assertRaisesRegex(ValueError,'provenance'):self.export()
        self.assertFalse((self.root/'result.tsv').exists())
        self.assertFalse(list(self.root.glob('.reference-export-*')))
        self.mutate("UPDATE references_table SET source_id=99 WHERE sample='A'")
        with self.assertRaisesRegex(ValueError,'source/depth'):self.export()

    def test_unknown_schema_rejected_and_missing_provenance_counted(self):
        self.mutate("UPDATE references_table SET record_json='{}'")
        info, rows=self.export()
        self.assertEqual(info['missing_alignment_paths'],5)
        self.assertTrue(all(r['CramFilePath']=='' for r in rows))
        (self.root/'result.tsv').unlink()
        self.mutate('UPDATE metadata SET value=? WHERE key=?',(json.dumps('coverage-reference-sqlite-2'),'schema'))
        with self.assertRaisesRegex(ValueError,'sqlite-3'):self.export()

    def test_actual_cli_and_error_exit(self):
        script=Path(__file__).resolve().parents[1]/'bin/reference_export.py'
        args=[sys.executable,'-B',str(script),'--index',str(self.index),'--out',str(self.root/'cli.tsv'),
              '--sample','A','--gcf','G2','--min-depth','3','--max-depth','3']
        result=subprocess.run(args,capture_output=True,text=True,check=True)
        self.assertEqual(json.loads(result.stdout)['records'],1)
        result=subprocess.run(args,capture_output=True,text=True)
        self.assertEqual(result.returncode,2)
        self.assertIn('already exists',result.stderr)
