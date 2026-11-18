"""Original-only persistence, nearest lookup and dynamic-depth numerical checks."""
import csv
import json
from pathlib import Path
import shutil
import subprocess
import sys
import unittest
from unittest.mock import patch
import numpy as np

sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import test_ref_depth_database as fixtures
from bin import ref_depth_database as builder, coverage_concordance as cli, concordance_batch as batch
from bin.ref_depth_catalog import collect, completed_sample, create_plan, digest
from bin.ref_depth_worker import run_gcf
from bin.reference_index import build_index, ReferenceIndex
from bin.native_reference import fixed_windows, covered_sums
from bin.native_concordance import NativeCache, compare
from bin.coverage_core import Reference, BlockPool, dense_coverage, select_units
from bin.concordance_core import task_rng


class NativeTests(unittest.TestCase):
    def setUp(self):
        self.fixture=fixtures.DepthDatabaseTests('runTest');self.fixture.setUp();self.root=self.fixture.root
        self.config=dict(self.fixture.config,storage_mode='native',depths=[],window_size=7,bins=None)
        self.database=self.root/'native';self.plan=create_plan(self.fixture.records,self.fixture.refs,self.database,self.config,fixtures.COLUMNS)
        with patch('bin.ref_depth_worker.simulate_grid',side_effect=AssertionError('No depth simulation allowed during build')):
            self.assertEqual(run_gcf(self.database,fixtures.GCF1,self.fixture.extractor),0)
        self.assertTrue(collect(self.database)['complete'])
        self.index=self.root/'native.sqlite';build_index(self.database,self.index,True)
        self.query=self.root/'query.cram';self.query.write_bytes(b'independent query source')

    def tearDown(self):
        if batch.CACHE is not None:
            batch.CACHE.close();batch.INDEX.close();batch.CACHE=batch.INDEX=None
        self.fixture.tearDown()

    def result(self,sample='a'):
        s=next(s for s in self.plan['groups'][fixtures.GCF1]['samples'] if s['sample']==sample)
        return completed_sample(self.database/fixtures.GCF1/'samples'/s['folder'],digest(self.database/'database.json'))

    def args(self,name='query',depths='.1',**changes):
        words=['coverage_concordance.py','--comparison-policy','nearest','--reference-sqlite',str(self.index),'--test-alignment',str(self.query),
            '--query-sample','new','--gcf',fixtures.GCF1,'--outdir',str(self.root/name),'--replicates','5','--seed','123']
        if depths is not None:words+=['--depths',depths]
        with patch.object(sys,'argv',words):args=cli.parse_args()
        for k,v in changes.items():setattr(args,k,v)
        return args

    def execute(self,args,pool=None):
        with patch('bin.native_concordance.extract_pool',return_value=(pool or self.fixture.big,self.fixture.layout)) as called:
            self.assertEqual(batch.run_indexed(args,cli.depth_values(args)),0)
        with (args.outdir/'summary.tsv').open() as handle:rows=list(csv.DictReader(handle,delimiter='\t'))
        return rows,called.call_count

    def test_original_only_dense_coverage_breadth_and_contig_boundaries(self):
        attempt,record=self.result()
        self.assertNotIn('depths',record);self.assertNotIn('replicates',record)
        self.assertFalse((attempt/'coverage').exists());self.assertFalse(list(attempt.rglob('ids-*.npz')))
        edges=np.load(attempt/'edges.npy')
        np.testing.assert_array_equal(edges,[0,7,14,20,27,31])
        dense=dense_coverage(self.fixture.small,np.arange(len(self.fixture.small)))
        with np.load(attempt/'native_coverage.npz') as data:
            expected=np.array([dense[a:b].sum() for a,b in zip(edges[:-1],edges[1:])])
            breadth=np.array([np.count_nonzero(dense[a:b]) for a,b in zip(edges[:-1],edges[1:])])
            np.testing.assert_array_equal(data['aligned_bases'],expected)
            np.testing.assert_array_equal(data['covered_bases'],breadth)
            np.testing.assert_array_equal(data['mean_depth'],expected/np.diff(edges))
            np.testing.assert_array_equal(data['length'],[7,7,6,7,4])
        self.assertEqual(record['native']['covered_bases'],np.count_nonzero(dense))
        self.assertEqual(record['full_mean_depth'],dense.sum()/31)
        self.assertTrue((self.database/'NATIVE_INDEX.tsv').is_file())
        self.assertFalse((self.database/'READY_DEPTH_INDEX.tsv').exists())

    def test_one_window_and_empty_pool(self):
        layout=Reference(('virus',),(300,),('',))
        edges,_,_,_=fixed_windows(layout,1000);np.testing.assert_array_equal(edges,[0,300])
        empty=BlockPool.from_intervals(np.empty((0,2),dtype=int),300)
        np.testing.assert_array_equal(covered_sums(empty,edges),[0])
        from bin.native_concordance import descriptive
        value=descriptive(np.array([1.]),np.array([2.]))
        self.assertIsNone(value['pearson']['similarity']);self.assertIsNone(value['jaccard']['p_value'])

    def test_index_nearest_exclusions_and_numeric_index(self):
        index=ReferenceIndex(self.index)
        try:
            self.assertEqual(index.metadata['storage_mode'],'native')
            self.assertEqual(index.metadata['counts']['depths'],0)
            self.assertEqual(index.nearest(fixtures.GCF1,1)[0]['sample'],'a')
            self.assertEqual(index.nearest(fixtures.GCF1,15)[0]['sample'],'b')
            self.assertEqual(index.nearest(fixtures.GCF1,1,exclude_sample='a')[0]['sample'],'b')
            self.assertEqual(index.nearest(fixtures.GCF1,1,samples={'b'})[0]['sample'],'b')
            row=index.nearest(fixtures.GCF1,1)[0]
            self.assertEqual(index.nearest(fixtures.GCF1,1,exclude_source=row['source_sha'])[0]['sample'],'b')
            plan=index.connection.execute('EXPLAIN QUERY PLAN SELECT * FROM references_table INDEXED BY native_depth_lookup WHERE gcf=? AND full_mean_depth>=? ORDER BY full_mean_depth,sample LIMIT 1',(fixtures.GCF1,1)).fetchall()
            self.assertIn('SEARCH',' '.join(r[3] for r in plan))
            a=index.nearest(fixtures.GCF1,1)[0]['full_mean_depth'];b=index.nearest(fixtures.GCF1,20)[0]['full_mean_depth']
            self.assertEqual(index.nearest(fixtures.GCF1,(a+b)/2,limit=2)[0]['sample'],'a')
        finally:index.close()

    def test_dynamic_arbitrary_depth_dense_and_query_extract_once(self):
        args=self.args(depths='.1,.2');before={str(p):digest(p) for p in self.database.rglob('*') if p.is_file()}
        rows,calls=self.execute(args);self.assertEqual(calls,1);self.assertEqual(len(rows),2)
        edges=np.load(self.result()[0]/'edges.npy');widths=np.diff(edges)
        for row in rows:
            observed=json.loads(Path(row['result_path']).read_text());depth=float(row['comparison_depth'])
            self.assertEqual(row['comparison_mode'],'dynamic_reference_sampling')
            self.assertIsInstance(observed['jaccard'],dict)
            matrix=[]
            for replicate in range(5):
                ids=select_units(self.fixture.small,task_rng(123,depth,0,replicate),target_bases=depth*31)
                dense=dense_coverage(self.fixture.small,ids)
                matrix.append([dense[a:b].sum() for a,b in zip(edges[:-1],edges[1:])])
            ids=np.load(Path(row['result_path']).parent/'test_ids.npy')
            dense=dense_coverage(self.fixture.big,ids)
            test=np.array([dense[a:b].sum() for a,b in zip(edges[:-1],edges[1:])])
            np.testing.assert_array_equal(test,np.load(Path(row['result_path']).parent/'test.npy'))
            expected=compare(test,np.array(matrix),widths,dict(alpha=.01,similarity_block_size=64))
            self.assertEqual(observed['jaccard'],expected['jaccard']);self.assertEqual(observed['pearson'],expected['pearson'])
        self.assertEqual(before,{str(p):digest(p) for p in self.database.rglob('*') if p.is_file()})

    def test_default_native_query_and_full_depth_descriptive(self):
        args=self.args(depths=None);rows,calls=self.execute(args)
        self.assertEqual(calls,1);self.assertEqual(rows[0]['reference_sample'],'b')
        result=json.loads(Path(rows[0]['result_path']).read_text())
        self.assertEqual(result['simulations'],0);self.assertIsNone(result['jaccard']['p_value'])
        self.assertFalse(result['reference_sampled'])

    def test_exclusion_before_limit_missing_gcf_and_low_query(self):
        args=self.args(query_sample='a');rows,_=self.execute(args)
        self.assertEqual(rows[0]['reference_sample'],'b')
        args=self.args(name='missing',gcf=fixtures.GCF2);rows,calls=self.execute(args)
        self.assertEqual(calls,0);self.assertEqual(rows[0]['status'],'NO_REFERENCE_GCF')
        args=self.args(name='too-deep',depths='1000');rows,_=self.execute(args)
        self.assertEqual(rows[0]['status'],'UNSUPPORTED_QUERY_DEPTH')

    def test_reuse_budget_relocation_and_corruption(self):
        index=ReferenceIndex(self.index);cache=NativeCache(index,10*1024*1024)
        from bin.stage_timing import StageTiming
        try:
            row=index.nearest(fixtures.GCF1,1)[0];cache.load(row,StageTiming());cache.load(row,StageTiming())
            self.assertEqual(cache.hits,1);self.assertLessEqual(cache.bytes,cache.budget)
            attempt,_=self.result();(attempt/'native_coverage.npz').write_bytes(b'bad')
            with self.assertRaises(ValueError):cache.load(row,StageTiming())
            self.assertFalse(collect(self.database)['complete'])
        finally:cache.close();index.close()

    def test_empty_native_build_and_relocated_query_without_original_reference(self):
        empty_db=self.root/'empty-db'
        create_plan(self.fixture.records,self.fixture.refs,empty_db,self.config,fixtures.COLUMNS)
        empty=BlockPool.from_intervals([],31)
        self.assertEqual(run_gcf(empty_db,fixtures.GCF1,lambda *a,**k:(empty,self.fixture.layout)),0)
        self.assertTrue(collect(empty_db)['complete'])
        empty_index=self.root/'empty.sqlite';build_index(empty_db,empty_index,True)
        index=ReferenceIndex(empty_index)
        try:self.assertEqual(index.nearest(fixtures.GCF1,1),[])
        finally:index.close()
        moved=self.root/'moved';shutil.copytree(self.database,moved);shutil.rmtree(self.database)
        self.fixture.a.unlink();self.fixture.b.unlink();self.fixture.fai.unlink()
        args=self.args(name='relocated',database_root=moved)
        rows,calls=self.execute(args)
        self.assertEqual(calls,1);self.assertEqual(rows[0]['status'],'READY')

    def test_native_prepare_offline_defaults_and_no_scheduler(self):
        args=builder.parser().parse_args(['run','--mode','slurm','--dry-run','--records',str(self.fixture.records),
            '--references',str(self.fixture.refs),'--outdir',str(self.root/'dry'),'--jobs','2'])
        with patch('bin.ref_depth_database.subprocess.run') as called,patch('bin.ref_depth_database.subprocess.Popen') as spawned:
            result=builder.run(args)
        called.assert_not_called();spawned.assert_not_called()
        plan=json.loads((args.outdir/'database.json').read_text())
        self.assertEqual(plan['config']['window_size'],1000);self.assertNotIn('depths',plan['config'])
        self.assertIsNone(plan['config']['memory']);self.assertEqual(result['concurrency_limit'],2)


class NativeCliTests(unittest.TestCase):
    def test_synthetic_bam_local_slurm_scripts_and_parallel_queries(self):
        try:import pysam
        except ImportError:self.skipTest('pysam unavailable')
        import tempfile
        with tempfile.TemporaryDirectory() as temporary:
            root=Path(temporary);fai=root/'ref.fai';fai.write_text('one\t2300\ntwo\t700\n')
            bam=root/'pool.bam'
            with pysam.AlignmentFile(bam,'wb',header={'HD':{'VN':'1.6'},'SQ':[{'SN':'one','LN':2300}]}) as out:
                for i in range(160):
                    read=pysam.AlignedSegment();read.query_name=f'r{i}';read.query_sequence='A'*80;read.flag=0
                    read.reference_id=0;read.reference_start=(i*13)%2200;read.cigarstring='80M';read.mapping_quality=60;out.write(read)
            records=root/'records.csv';refs=root/'refs.csv';queries=root/'queries.csv'
            g1,g2=fixtures.GCF1,fixtures.GCF2
            fixtures.write_csv(records,['Sample','CramFilePath','GCFID'],[dict(Sample='ref',CramFilePath=str(bam),GCFID=g) for g in (g1,g2)])
            fixtures.write_csv(refs,['GCF_ID','GCF_path'],[dict(GCF_ID=g,GCF_path=str(fai)) for g in (g1,g2)])
            # Distinct source bytes for independent query/exclusion semantics.
            query=root/'query.cram'
            with pysam.AlignmentFile(bam,'rb') as src,pysam.AlignmentFile(query,'wc',header=src.header,format_options=[b'no_ref=1']) as dst:
                for read in src:
                    read.query_name='q'+read.query_name;dst.write(read)
            fixtures.write_csv(queries,['Sample','CramFilePath','GCFID'],[dict(Sample=s,CramFilePath=str(query),GCFID=g) for g in (g1,g2) for s in ('q1','q2')])
            binroot=Path(__file__).resolve().parents[1]/'bin';summaries=[]
            for mode in ('local','slurm'):
                database=root/mode
                command=[sys.executable,'-B',str(binroot/'ref_depth_database.py'),'run','--mode',mode,'--records',str(records),'--references',str(refs),
                    '--outdir',str(database),'--cpus','2','--jobs','2','--geometry-format','packed']
                if mode=='slurm':command+=['--dry-run']
                subprocess.run(command,check=True,stdout=subprocess.DEVNULL)
                if mode=='slurm':
                    for script in database.glob('*/job.sbatch.sh'):subprocess.run(['bash',str(script)],check=True)
                    self.assertTrue(collect(database)['complete'])
                index=root/(mode+'.sqlite');build_index(database,index,True)
                for jobs,cache in ((1,0),(2,64)):
                    output=root/f'{mode}-q{jobs}'
                    subprocess.run([sys.executable,'-B',str(binroot/'coverage_concordance.py'),'--comparison-policy','nearest','--reference-sqlite',str(index),'--query-manifest',str(queries),
                        '--depths','.1,1','--replicates','5','--query-jobs',str(jobs),'--reference-cache-mb',str(cache),'--outdir',str(output)],check=True)
                    with (output/'summary.tsv').open() as handle:rows=list(csv.DictReader(handle,delimiter='\t'))
                    self.assertEqual(len(rows),8);self.assertTrue(all(r['status']=='READY' for r in rows))
                    summaries.append([{k:r[k] for k in ('query_sample','gcf','depth','jaccard','pearson','test_actual_depth','comparison_depth')} for r in rows])
            for rows in summaries[1:]:self.assertEqual(rows,summaries[0])


if __name__=='__main__':unittest.main()
