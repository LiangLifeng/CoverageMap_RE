"""Multi-root resolution, post-depth sample support, reproducible immutable selection."""
import csv
import json
from pathlib import Path
import shutil
import unittest
from unittest.mock import patch

import test_ref_depth_database as fixtures
from bin import concordance_batch as batch, coverage_concordance as cli
from bin.ref_depth_catalog import create_plan, collect, digest
from bin.ref_depth_worker import run_gcf
from bin.reference_index import build_index, ReferenceIndex, _build_single
from bin.reference_registry import filter_registry, export_report


class RegistryTests(unittest.TestCase):
    def setUp(self):
        self.f = fixtures.DepthDatabaseTests('runTest'); self.f.setUp(); self.root = self.f.root
        self.config = dict(self.f.config, storage_mode='native', depths=[], bins=None, window_size=7)
        fixtures.write_csv(self.f.refs, ['GCF_ID','GCF_path'],
            [dict(GCF_ID=g, GCF_path=str(self.f.fai)) for g in (fixtures.GCF1, fixtures.GCF2)])
        self.databases = []
        for number, rows in enumerate(([self.f.rows[0]], [self.f.rows[1], dict(self.f.rows[1], GCFID=fixtures.GCF2)])):
            records = self.root/f'records-{number}.csv'
            fixtures.write_csv(records, ['Sample','CramFilePath','GCFID'], rows)
            database = self.root/f'db-{number}'
            plan = create_plan(records, self.f.refs, database, self.config, fixtures.COLUMNS)
            for gcf in plan['groups']:
                self.assertEqual(run_gcf(database, gcf, self.f.extractor), 0)
            self.assertTrue(collect(database)['complete']); self.databases.append(database)
        self.path = self.root/'all.sqlite'
        self.before = {str(p):digest(p) for db in self.databases for p in db.rglob('*') if p.is_file()}
        build_index(self.databases, self.path, True)

    def tearDown(self):
        if batch.CACHE is not None:
            batch.CACHE.close(); batch.INDEX.close(); batch.CACHE = batch.INDEX = None
        self.f.tearDown()

    def test_depth_then_distinct_sample_filter_and_reports(self):
        out = self.root/'selected.sqlite'
        # a=43/31, b=500/31. Exact boundary is inclusive; only GCF1 has two samples.
        info = filter_registry(self.path, out, self.f.small.total_bases/31, 2)
        self.assertEqual(info['selected_counts'], dict(gcfs=1, references=2))
        index = ReferenceIndex(out)
        try:
            self.assertIsNone(index.group(fixtures.GCF2, True))
            self.assertEqual([r['sample'] for r in index.nearest(fixtures.GCF1, 10, 2)], ['b','a'])
            self.assertEqual(len(list(index.references(fixtures.GCF1, 1))), 2)
        finally: index.close()
        report = self.root/'report'; export_report(out, report)
        with (report/'gcf_depths.tsv').open() as handle: rows = list(csv.DictReader(handle, delimiter='\t'))
        self.assertEqual([r['total_samples'] for r in rows], ['2','1'])
        self.assertEqual(float(rows[0]['all_depth_median']), (self.f.small.total_bases+500)/62)
        with (report/'sample_depths.tsv').open() as handle: samples = list(csv.DictReader(handle, delimiter='\t'))
        self.assertEqual(len(samples), 3)
        self.assertEqual(samples[-1]['exclusion_reason'], 'INSUFFICIENT_SAMPLES_AFTER_DEPTH_FILTER')
        strict = self.root/'strict.sqlite'; filter_registry(out, strict, 2, 2)
        strict_index = ReferenceIndex(strict)
        try: self.assertEqual(strict_index.metadata['selected_counts'], dict(gcfs=0,references=0))
        finally: strict_index.close()
        # Refiltering starts from all audited rows, so relaxing restores previously excluded references.
        relaxed = self.root/'relaxed.sqlite'; filter_registry(strict, relaxed, 0, 1)
        relaxed_index = ReferenceIndex(relaxed)
        try: self.assertEqual(relaxed_index.metadata['selected_counts']['references'], 3)
        finally: relaxed_index.close()
        self.assertEqual(self.before, {str(p):digest(p) for db in self.databases for p in db.rglob('*') if p.is_file()})

    def test_native_query_resolves_each_source_and_filtered_gcf_skips_extraction(self):
        query = self.root/'query.cram'; query.write_bytes(b'independent')
        words = ['--comparison-policy','nearest','--reference-sqlite',str(self.path),'--test-alignment',str(query),'--query-sample','q',
                 '--gcf',fixtures.GCF1,'--depth','.1','--reference-limit','2','--replicates','3','--outdir',str(self.root/'query')]
        with patch('sys.argv',['coverage_concordance.py',*words]): args = cli.parse_args()
        with patch('bin.native_concordance.extract_pool', return_value=(self.f.big,self.f.layout)) as call:
            self.assertEqual(batch.run_indexed(args,cli.depth_values(args)), 0)
        self.assertEqual(call.call_count, 1)
        with (args.outdir/'summary.tsv').open() as handle: rows = list(csv.DictReader(handle,delimiter='\t'))
        self.assertEqual({r['reference_sample'] for r in rows}, {'a','b'})
        self.assertTrue(all(r['status']=='READY' for r in rows))
        self.assertEqual({str(Path(r['sample_manifest']).relative_to(self.root).parts[0]) for r in rows}, {'db-0','db-1'})
        empty=self.root/'empty.sqlite';filter_registry(self.path,empty,0,3)
        args.reference_sqlite=empty;args.outdir=self.root/'empty-query'
        with patch('bin.native_concordance.extract_pool',side_effect=AssertionError('Excluded GCF must skip extraction')):
            self.assertEqual(batch.run_indexed(args,cli.depth_values(args)),0)
        with self.assertRaisesRegex(ValueError, 'cannot relocate'):
            ReferenceIndex(self.path, self.root)

    def test_duplicates_config_and_layout_conflicts_rejected_no_partial_publication(self):
        copy = self.root/'copy'; shutil.copytree(self.databases[0],copy)
        with self.assertRaisesRegex(ValueError,'distinct'):
            build_index([self.databases[0],self.databases[0]],self.root/'alias.sqlite')
        out=self.root/'overlap.sqlite'
        info=build_index([self.databases[0],copy],out,True)
        self.assertEqual(info['counts']['references'],1)
        self.assertEqual(info['duplicate_copies'],1)
        self.assertEqual(info['selected_counts']['references'],1)
        index=ReferenceIndex(out)
        try:self.assertEqual(index.connection.execute('SELECT COUNT(*) FROM reference_sources').fetchone()[0],2)
        finally:index.close()
        # Metadata disagreement is never resolved by picking the deeper/first copy.
        from bin.reference_registry import duplicate_signature
        signature=duplicate_signature(copy,next(copy.glob('*/samples/*/attempt-*/result.json')).relative_to(copy))
        with patch('bin.reference_registry.duplicate_signature',side_effect=[signature,dict(signature,units=-99)]):
            with self.assertRaisesRegex(ValueError,'Conflicting duplicate'):
                build_index([self.databases[0],copy],self.root/'conflicting.sqlite')
        self.assertFalse((self.root/'conflicting.sqlite').exists())
        # Configuration compatibility checked before a conflicting root is merged.
        original = self.config['window_size']
        with patch('bin.reference_registry.scientific_config', side_effect=[{'window':original},{'window':999}]):
            with self.assertRaisesRegex(ValueError,'Incompatible'):
                build_index(self.databases,self.root/'config.sqlite')
        self.assertFalse((self.root/'config.sqlite').exists())

    def test_metadata_filter_does_not_touch_payloads_and_old_schema_still_readable(self):
        with patch('bin.reference_registry._build_single',side_effect=AssertionError('No source reindex')):
            filter_registry(self.path,self.root/'metadata-only.sqlite',100,1)
        old = self.root/'old.sqlite'; _build_single(self.databases[0],old,True)
        index=ReferenceIndex(old)
        try:
            self.assertEqual(index.nearest(fixtures.GCF1,1)[0]['sample'],'a')
            self.assertEqual(index.root_for(index.nearest(fixtures.GCF1,1)[0]),self.databases[0])
        finally:index.close()
        for depth, samples in [(float('nan'),1),(-1,1),(0,0)]:
            with self.assertRaises(ValueError):filter_registry(self.path,self.root/'bad.sqlite',depth,samples)
        with self.assertRaises(FileExistsError):build_index(self.databases,self.path)

    def test_known_m5_conflict_rejected_even_if_all_samples_filtered_out(self):
        from bin.coverage_core import Reference
        records = self.root/'conflict.csv'
        fixtures.write_csv(records,['Sample','CramFilePath','GCFID'],
            [dict(self.f.rows[1],Sample='conflicting')])
        database = self.root/'conflict'
        create_plan(records,self.f.refs,database,self.config,fixtures.COLUMNS)
        wrong = Reference(self.f.layout.names,self.f.layout.lengths,('ccc','bbb'))
        self.assertEqual(run_gcf(database,fixtures.GCF1,lambda *args,**kwargs:(self.f.big,wrong)),0)
        self.assertTrue(collect(database)['complete'])
        with self.assertRaisesRegex(ValueError,'M5'):
            build_index([self.databases[0],database],self.root/'m5.sqlite',min_depth=1000)
        self.assertFalse((self.root/'m5.sqlite').exists())

    def test_build_selection_equals_metadata_filter_and_output_is_exclusive(self):
        built=self.root/'built.sqlite';filtered=self.root/'filtered.sqlite'
        build_index(self.databases,built,True,min_depth=2,min_samples=1)
        filter_registry(self.path,filtered,2,1)
        left=ReferenceIndex(built);right=ReferenceIndex(filtered)
        try:
            for table in ('references_table','gcf_statistics','gcfs','sources'):
                self.assertEqual([tuple(r) for r in left.connection.execute('SELECT * FROM '+table)],
                                 [tuple(r) for r in right.connection.execute('SELECT * FROM '+table)])
            self.assertEqual(left.nearest(fixtures.GCF1,1)[0]['sample'],'b')
        finally:left.close();right.close()
        before=digest(filtered)
        with self.assertRaises(FileExistsError):filter_registry(self.path,filtered,0,1)
        self.assertEqual(digest(filtered),before)
        link=self.root/'link.sqlite';link.symlink_to(self.root/'absent')
        with self.assertRaises(FileExistsError):filter_registry(self.path,link)
        self.assertTrue(link.is_symlink());self.assertFalse(link.exists())


try:
    import pysam
except ImportError:
    pysam=None


@unittest.skipIf(pysam is None,'Real process/CLI multi-root integration requires pysam on NH')
class MultiRootProcessTests(unittest.TestCase):
    def test_native_and_legacy_process_queries_use_all_roots(self):
        import subprocess
        import sys
        import tempfile
        with tempfile.TemporaryDirectory() as temporary:
            root=Path(temporary);fai=root/'one.fai';fai.write_text('one\t31\n')
            bam=root/'source.bam'
            with pysam.AlignmentFile(str(bam),'wb',header={'SQ':[{'SN':'one','LN':31}]}) as handle:
                for i in range(50):
                    read=pysam.AlignedSegment(handle.header);read.query_name=f'q{i}';read.query_sequence='A'*10
                    read.flag=0;read.reference_id=0;read.reference_start=i%20;read.mapping_quality=40;read.cigar=[(0,10)]
                    handle.write(read)
            refs=root/'refs.csv'
            fixtures.write_csv(refs,['GCF_ID','GCF_path'],[dict(GCF_ID=g,GCF_path=str(fai)) for g in (fixtures.GCF1,fixtures.GCF2)])
            query=root/'query.bam';shutil.copyfile(bam,query)
            # Distinct bytes/source hash for normal self-exclusion.
            config=dict(depths=[.1],cpus=1,memory=None,python=sys.executable,replicates=3,bins=7,seed=41,
                batch_size=2,sampler='permutation',geometry_format='packed',matrix_format='npz',insufficient_depth='record',
                policies=dict(mate_policy='both',sampling_unit='read',min_mapq=0,exclude_duplicates=False,exclude_qcfail=False))
            for native in (True,False):
                databases=[]
                for sample in ('a','b'):
                    records=root/f'{native}-{sample}.csv'
                    fixtures.write_csv(records,['Sample','CramFilePath','GCFID'],
                        [dict(Sample=sample,CramFilePath=str(bam),GCFID=g) for g in (fixtures.GCF1,fixtures.GCF2)])
                    database=root/f'{native}-{sample}'
                    cfg=dict(config,storage_mode='native',depths=[],bins=None,window_size=7) if native else config
                    plan=create_plan(records,refs,database,cfg,fixtures.COLUMNS)
                    for gcf in plan['groups']:self.assertEqual(run_gcf(database,gcf),0)
                    self.assertTrue(collect(database)['complete']);databases.append(database)
                index=root/f'{native}.sqlite'
                indexer=Path(__file__).resolve().parents[1]/'bin/reference_index.py'
                subprocess.run([sys.executable,'-B',str(indexer),'build','--database',str(databases[0]),
                    '--database',str(databases[1]),'--out',str(index),'--verify-payloads'],check=True)
                subprocess.run([sys.executable,'-B',str(indexer),'filter','--index',str(index),
                    '--out',str(root/f'{native}-filtered.sqlite'),'--min-depth','0','--min-samples','2'],check=True)
                subprocess.run([sys.executable,'-B',str(indexer),'report','--index',str(root/f'{native}-filtered.sqlite'),
                    '--outdir',str(root/f'{native}-report')],check=True)
                queries=root/f'{native}-queries.csv'
                fixtures.write_csv(queries,['Sample','CramFilePath','GCFID'],
                    [dict(Sample='query',CramFilePath=str(query),GCFID=g) for g in (fixtures.GCF1,fixtures.GCF2)])
                entry=Path(__file__).resolve().parents[1]/'bin/coverage_concordance.py'
                expected=None
                for jobs in (1,2):
                    out=root/f'query-{native}-{jobs}'
                    subprocess.run([sys.executable,'-B',str(entry),'--comparison-policy','nearest','--reference-sqlite',str(index),
                        '--query-manifest',str(queries),'--depth','.1','--include-self','--reference-limit','2',
                        '--replicates','3','--query-jobs',str(jobs),'--outdir',str(out)],check=True)
                    with (out/'summary.tsv').open() as handle:rows=list(csv.DictReader(handle,delimiter='\t'))
                    self.assertEqual(len(rows),4);self.assertEqual({r['status'] for r in rows},{'READY'})
                    self.assertEqual({r['reference_sample'] for r in rows},{'a','b'})
                    observed={(r['gcf'],r['reference_sample']):{k:json.loads(Path(r['result_path']).read_text())[k]
                        for k in ('jaccard','pearson','test_selected_ids_sha256')} for r in rows}
                    if expected is None:expected=observed
                    else:self.assertEqual(observed,expected)
