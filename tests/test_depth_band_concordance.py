"""Independent statistics, inclusive lookup, layered queries and seeded source draws."""
import csv
import json
from pathlib import Path
from statistics import NormalDist
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import test_native_reference as native_fixture
import test_ref_depth_database as fixtures
from bin import concordance_batch as batch, coverage_concordance as cli
from bin.depth_band_concordance import configure, high_coverage, select_references, selection_size
from bin.native_concordance import compare
from bin.reference_index import ReferenceIndex, build_index
from bin.ref_depth_catalog import create_plan, collect, digest
from bin.ref_depth_worker import run_gcf
from bin.coverage_core import BlockPool


class StatisticTests(unittest.TestCase):
    def test_tools_high_maximum_against_independent_pair_distribution(self):
        refs = np.array([[0,1,4,2,0], [3,0,1,0,5], [0,0,1,3,4], [1,3,0,0,2]], float)
        query = np.array([0,1,3,0,2.])
        observed = high_coverage(query, refs, .01)
        def jaccard(a,b):
            return np.logical_and(a>0,b>0).sum()/np.logical_or(a>0,b>0).sum()
        def pearson(a,b):
            a,b=a-a.mean(),b-b.mean()
            # Exact epsilon convention of unchanged original Tools.py.
            return (a@b)/(np.linalg.norm(a)*np.linalg.norm(b)+1e-12)
        for name, function in (('jaccard',jaccard),('pearson',pearson)):
            distribution = np.array([function(a,b) for i,a in enumerate(refs) for b in refs[i+1:]])
            best = max(function(query,r) for r in refs)
            expected = NormalDist().cdf((best-distribution.mean())/distribution.std())
            self.assertAlmostEqual(observed[name]['similarity'],best,places=11)
            self.assertAlmostEqual(observed[name]['p_value'],expected,places=11)
            self.assertEqual(observed[name]['reference_pairs'],6)

    def test_tools_low_mean_and_off_diagonal_distribution(self):
        matrix=np.array([[0,2,4,3,0],[1,0,3,2,0],[0,2,1,0,2],[3,0,0,1,1]],float)
        query=np.array([0,1,3,1,1.]);widths=np.array([7,7,6,7,4])
        result=compare(query,matrix,widths,dict(alpha=.01,similarity_block_size=2))
        q=query/widths;refs=matrix/widths
        for name in ('jaccard','pearson'):
            def function(a,b):
                if name=='jaccard':return np.count_nonzero((a>0)&(b>0))/np.count_nonzero((a>0)|(b>0))
                return np.corrcoef(a,b)[0,1]
            means=np.array([np.mean([function(a,b) for j,b in enumerate(refs) if i!=j]) for i,a in enumerate(refs)])
            value=np.mean([function(q,b) for b in refs])
            p=NormalDist().cdf((value-means.mean())/means.std())
            self.assertAlmostEqual(result[name]['similarity'],value)
            self.assertAlmostEqual(result[name]['p_value'],p)

    def test_source_sampling_count_percent_all_and_seed_order(self):
        rows=[dict(sample=str(i)) for i in range(11)]
        for spec,count in (('all',11),('3',3),('20%',3),('100%',11),('100',11)):
            first=select_references(rows,spec,401,'GCF')
            self.assertEqual(len(first),count)
            self.assertEqual(first,select_references(rows[::-1],spec,401,'GCF'))
            self.assertEqual(len({r['sample'] for r in first}),count)
        for value in ('0','-2','0%','101%','nan%','bad'):
            with self.assertRaises(ValueError):selection_size(value)

    def test_degenerate_high_distributions_are_not_manufactured(self):
        for refs in (np.array([[1.,2,3]]),np.array([[1.,2,3],[2.,4,6]])):
            result=high_coverage(np.array([1.,3,2]),refs,.01)
            self.assertIsNone(result['pearson']['p_value'])
        result=high_coverage(np.ones(3),np.ones((3,3)),.01)
        self.assertEqual(result['pearson']['status'],'non_estimable')


class QueryTests(unittest.TestCase):
    def setUp(self):
        self.fx=native_fixture.NativeTests('runTest');self.fx.setUp()
        self.root=self.fx.root

    def tearDown(self):
        self.fx.tearDown()

    def args(self,name='band',**changes):
        return self.fx.args(name=name,depths=None,comparison_policy='depth-band',window_size=7,**changes)

    def run_query(self,args,pool=None):
        pool = self.fx.fixture.big if pool is None else pool
        with patch('bin.native_concordance.extract_pool',return_value=(pool,self.fx.fixture.layout)) as called:
            self.assertEqual(batch.run_indexed(args,cli.depth_values(args)),0)
        with (args.outdir/'summary.tsv').open() as handle:
            return list(csv.DictReader(handle,delimiter='\t')),called.call_count

    def test_low_unrestricted_sources_multiple_layers_and_mass_limit(self):
        args=self.args(simulation_depths='.1,.2,2',simulation_references='all')
        rows,calls=self.run_query(args,self.fx.fixture.small)
        self.assertEqual(calls,1)
        self.assertEqual(len(rows),5)  # two references for each supported layer, one unsupported query layer
        self.assertEqual(sum(r['status']=='READY' for r in rows),4)
        self.assertEqual(rows[-1]['status'],'UNSUPPORTED_QUERY_DEPTH')
        for row in rows[:-1]:
            result=json.loads(Path(row['result_path']).read_text())
            self.assertEqual(result['simulations'],5)
            self.assertEqual(len(result['reference_replicates']),5)
        self.assertEqual({r['reference_sample'] for r in rows[:-1]},{'a','b'})
        # Original b depth is >16X, outside the low query's +/-1X band, but usable by default.
        self.assertTrue(any(float(r['reference_full_depth'])>16 for r in rows[:-1]))

    def test_optional_source_band_and_range(self):
        for number,kw in enumerate((dict(simulation_match_depth_band=True),dict(simulation_source_depth_range='0:2'))):
            rows,_=self.run_query(self.args(name=f'limited-{number}',**kw),self.fx.fixture.small)
            self.assertEqual({r['reference_sample'] for r in rows},{'a'})
        rows,_=self.run_query(self.args(name='unsupported',simulation_source_depth_range='0:2'),
                              BlockPool.from_intervals(np.tile([[1,11]],(10,1)),31))
        self.assertEqual(rows[0]['status'],'UNSUPPORTED_REFERENCE_DEPTH')

    def test_high_uses_original_query_and_audits_nearest_outside_band(self):
        rows,_=self.run_query(self.args())
        self.assertEqual(rows[0]['comparison_mode'],'high_depth_reference_pool_maximum')
        result=json.loads(Path(rows[0]['result_path']).read_text())
        self.assertEqual(result['simulations'],0)
        self.assertEqual(result['reference_samples'],['b'])
        # No matched reference at ~8X: compare the nearest, explicitly outside the band.
        pool=BlockPool.from_intervals(np.tile([[1,11]],(25,1)),31)
        rows,_=self.run_query(self.args(name='gap'),pool)
        self.assertEqual(rows[0]['status'],'READY')
        self.assertEqual(rows[0]['selection_mode'],'nearest_depth_fallback')
        self.assertEqual(rows[0]['depth_band_matched'],'False')
        result=json.loads(Path(rows[0]['result_path']).read_text())
        self.assertIsNone(result['pearson']['p_value'])
        self.assertEqual(result['pearson']['reason'],'outside_requested_depth_band')

    def test_zero_missing_and_self_exclusion(self):
        rows,_=self.run_query(self.args(),BlockPool.from_intervals([],31))
        self.assertEqual(rows[0]['status'],'ZERO_QUERY_DEPTH')
        rows,calls=self.run_query(self.args(name='missing',gcf=fixtures.GCF2))
        self.assertEqual(calls,0);self.assertEqual(rows[0]['status'],'NO_REFERENCE_GCF')
        rows,_=self.run_query(self.args(name='exclude',query_sample='a',include_self=False),self.fx.fixture.small)
        self.assertEqual({r['reference_sample'] for r in rows},{'b'})

    def test_same_sample_and_source_are_included_by_default_single_reference_scores(self):
        args=self.args(query_sample='b',test_alignment=self.fx.fixture.b)
        rows,_=self.run_query(args)
        self.assertEqual(rows[0]['status'],'READY')
        self.assertEqual(rows[0]['reference_count'],'1')
        self.assertEqual(rows[0]['comparison_scope'],'self-consistency-included')
        self.assertAlmostEqual(float(rows[0]['pearson']),1.)
        self.assertEqual(rows[0]['pearson_p'],'')
        selection=json.loads(next(args.outdir.glob('queries/*/*/selection.json')).read_text())
        self.assertTrue(selection['include_self'])
        self.assertEqual(selection['selected_samples'],['b'])
        # Explicit exclusion still means exclusion: fall back to a different source, never sneak b back in.
        args=self.args(name='explicit',query_sample='b',test_alignment=self.fx.fixture.b,include_self=False)
        rows,_=self.run_query(args)
        self.assertEqual(rows[0]['reference_count'],'1')
        result=json.loads(Path(rows[0]['result_path']).read_text())
        self.assertEqual(result['reference_samples'],['a'])
        self.assertEqual(result['selection_mode'],'nearest_depth_fallback')

    def test_unestimable_query_is_processed_instead_of_skipped(self):
        args=self.args()
        self.run_query(args,BlockPool.from_intervals([],31))
        report=json.loads((args.outdir/'provenance.json').read_text())
        self.assertEqual(report['skipped_query_count'],0)
        self.assertEqual(report['outcomes'][0]['state'],'COMPLETE_WITH_GAPS')

    def test_filtered_gcf_is_reported_without_skipping_or_restoring_references(self):
        from bin.reference_registry import filter_registry
        target=self.root/'filtered.sqlite'
        filter_registry(self.fx.index,target,10000,1)
        args=self.args(reference_sqlite=target)
        with patch('bin.native_concordance.extract_pool',side_effect=AssertionError('No selected reference')):
            self.assertEqual(batch.run_indexed(args,[]),0)
        report=json.loads((args.outdir/'provenance.json').read_text())
        self.assertEqual(report['skipped_query_count'],0)
        self.assertEqual(report['outcomes'][0]['state'],'COMPLETE_WITH_GAPS')
        self.assertEqual(report['outcomes'][0]['reason'],'NO_SELECTED_REFERENCE_GCF')

    def test_default_1000bp_rebin_and_preoutput_validation(self):
        args=self.args();args.window_size=1000
        rows,_=self.run_query(args,self.fx.fixture.small)
        self.assertTrue(all(r['window_size']=='1000' for r in rows))
        edges=np.load(next(args.outdir.glob('queries/*/*/edges.npy')))
        np.testing.assert_array_equal(edges,[0,20,31])
        for number,kw in enumerate((dict(depth_tolerance=-1),dict(simulation_depths='1,1'),dict(simulation_references='0%'),
                                    dict(simulation_source_depth_range='2:1'),dict(depth_sampler='auto'))):
            args=self.args(name=f'bad-{number}',**kw)
            with self.assertRaises(ValueError):batch.run_indexed(args,[])
            self.assertFalse(args.outdir.exists())

    def test_index_depth_boundaries_exclusions_and_query_plan(self):
        index=ReferenceIndex(self.fx.index)
        try:
            low=self.fx.fixture.small.total_bases/31;high=self.fx.fixture.big.total_bases/31
            self.assertEqual([r['sample'] for r in index.candidates(fixtures.GCF1,low,high)],['a','b'])
            self.assertEqual([r['sample'] for r in index.candidates(fixtures.GCF1,low+1e-8,high)],['b'])
            self.assertEqual(index.candidates(fixtures.GCF1,low,high,exclude_sample='b',samples={'b'}),[])
        finally:index.close()


class NativeProcessTests(unittest.TestCase):
    def test_real_cli_5x_boundary_pool_statistics_and_serial_process(self):
        try:import pysam
        except ImportError:self.skipTest('pysam required for actual BAM CLI integration')
        with tempfile.TemporaryDirectory() as temporary:
            root=Path(temporary);fai=root/'ref.fai';fai.write_text('one\t3000\ntwo\t1000\n')
            # Four full-FAI depths 4,5,6,7X with different nonuniform profiles.
            def bam(name,depth,offset):
                path=root/f'{name}.bam'
                with pysam.AlignmentFile(path,'wb',header={'HD':{'VN':'1.6'},'SQ':[{'SN':'one','LN':3000}]}) as handle:
                    for i in range(depth*40):
                        read=pysam.AlignedSegment();read.query_name=f'{name}-{i}';read.query_sequence='A'*100
                        read.flag=0;read.reference_id=0;read.reference_start=(i*offset)%2900
                        read.mapping_quality=60;read.cigar=[(0,100)];handle.write(read)
                return path
            records=root/'records.csv';refs=root/'refs.csv'
            fixtures.write_csv(refs,['GCF_ID','GCF_path'],[dict(GCF_ID=fixtures.GCF1,GCF_path=str(fai))])
            paths=[bam(f'r{d}',d,17+d) for d in (4,5,6,7)]
            fixtures.write_csv(records,['Sample','CramFilePath','GCFID'],[
                dict(Sample=f'r{d}',CramFilePath=str(p),GCFID=fixtures.GCF1) for d,p in zip((4,5,6,7),paths)])
            database=root/'db';entry=Path(__file__).resolve().parents[1]/'bin'/'coverage_concordance.py'
            config=dict(storage_mode='native',window_size=1000,cpus=1,python=sys.executable,memory=None,geometry_format='packed',
                        policies=dict(mate_policy='both',sampling_unit='read',min_mapq=0,exclude_duplicates=False,exclude_qcfail=False))
            databases=[]
            for number,positions in enumerate(((0,1),(2,3))):
                part=root/f'part-{number}.csv'
                fixtures.write_csv(part,['Sample','CramFilePath','GCFID'],[
                    dict(Sample=f'r{(4,5,6,7)[i]}',CramFilePath=str(paths[i]),GCFID=fixtures.GCF1) for i in positions])
                database=root/f'db-{number}'
                create_plan(part,refs,database,config,fixtures.COLUMNS)
                self.assertEqual(run_gcf(database,fixtures.GCF1),0);self.assertTrue(collect(database)['complete'])
                databases.append(database)
            index=root/'index.sqlite';build_index(databases,index,True)
            queries=root/'queries.csv'
            fixtures.write_csv(queries,['Sample','CramFilePath','GCFID'],[
                dict(Sample='q5',CramFilePath=str(bam('q5',5,13)),GCFID=fixtures.GCF1),
                dict(Sample='q4',CramFilePath=str(bam('q4',4,37)),GCFID=fixtures.GCF1)])
            before={str(p):digest(p) for db in databases for p in db.rglob('*') if p.is_file()}
            science=[]
            for jobs in (1,2):
                out=root/f'out-{jobs}'
                subprocess.run([sys.executable,'-B',str(entry),'--reference-sqlite',str(index),'--query-manifest',str(queries),
                    '--query-jobs',str(jobs),'--query-scheduling','query','--replicates','4','--simulation-depths','1,2',
                    '--simulation-references','50%','--outdir',str(out)],check=True)
                results=[json.loads(p.read_text()) for p in sorted(out.rglob('concordance.json'))]
                high=[r for r in results if r['comparison_mode']=='high_depth_reference_pool_maximum']
                self.assertEqual(high[0]['reference_samples'],['r4','r5','r6'])  # both +/-1 boundaries included, r7 excluded
                low=[r for r in results if r['comparison_mode']=='low_depth_multilayer_simulation']
                self.assertEqual(len(low),4)  # 50% of four sources, two layers, four repeats each
                self.assertEqual({r['requested_depth'] for r in low},{1.,2.})
                science.append([{k:v for k,v in r.items() if k not in ('result_path',)} for r in results])
            self.assertEqual(science[0],science[1])
            self.assertEqual(before,{str(p):digest(p) for db in databases for p in db.rglob('*') if p.is_file()})


if __name__=='__main__':unittest.main()
