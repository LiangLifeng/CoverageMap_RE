import numpy as np
import pytest
from bin.coverage_core import BlockPool, covered_sums, window_sums, canonical_reference


def test_overlap_gaps_units_and_unequal_windows():
    pool = BlockPool(np.array([0, 200, 0, 500]), np.array([150, 250, 150, 650]),
                     np.array([0, 2, 3, 4]), 1050)
    edges = np.array([0, 1000, 1050])
    assert covered_sums(pool, edges).tolist() == [350, 0]
    assert window_sums(pool, np.arange(3), edges).tolist() == [500, 0]
    assert covered_sums(pool, edges, np.array([1, 2])).tolist() == [300, 0]
    assert covered_sums(pool, edges, np.array([], dtype=int)).tolist() == [0, 0]
    with pytest.raises(ValueError):
        covered_sums(pool, np.array([0, 1049]), np.array([], dtype=int))


def test_virtual_contig_boundary_and_short_tail():
    pool = BlockPool.from_intervals(np.array([[60, 70], [70, 90], [135, 150], [65, 70]]), 150)
    assert covered_sums(pool, np.array([0, 50, 100, 145, 150])).tolist() == [0, 30, 10, 5]


def test_union_against_independent_dense_reference():
    rng = np.random.default_rng(88)
    starts = rng.integers(0, 1500, 300)
    ends = np.minimum(starts + rng.integers(1, 200, 300), 1501)
    pool = BlockPool.from_intervals(np.column_stack((starts, ends)), 1501)
    edges = np.array([0, 100, 800, 1501])
    for ids in (np.arange(300), rng.choice(300, 70, replace=False), np.array([], dtype=int)):
        covered = np.zeros(1501, dtype=bool)
        for i in ids:
            covered[starts[i]:ends[i]] = True
        expected = np.array([covered[a:b].sum() for a,b in zip(edges[:-1],edges[1:])])
        assert np.array_equal(covered_sums(pool, edges, ids), expected)

from bin.jaccard_breadth import covered_state, reference_means, thresholds
from bin.query_manifest_stats import similarity, reference_row_means, low_null


def test_threshold_boundaries_and_zero_union():
    lengths=np.array([1000,1000,1000,30])
    covered=np.array([1,10,100,3])
    assert covered_state(covered,lengths,0).tolist()==[True]*4
    assert covered_state(covered,lengths,.01).tolist()==[False,True,True,True]
    assert covered_state(covered,lengths,.1).tolist()==[False,False,True,True]
    assert not covered_state(np.array([0]),np.array([3]),0)[0]
    assert thresholds(0,'default') == [0,.01,.05,.1,.15,.2,.3]
    for bad in (-1,float('nan'),1.1):
        with pytest.raises(ValueError): thresholds(bad)
    score=similarity([0,0],[0,0])
    assert score['jaccard'] is None and score['jaccard_status']=='zero_union'
    with pytest.raises(ValueError): similarity([1,2],[2,3],threshold=.1)


def test_shared_threshold_nulls_and_pearson_exact_regression():
    matrix=np.random.default_rng(41).integers(0,100,(15,20))
    widths=np.full(20,100)
    covered=np.minimum(matrix,30)
    for tau in thresholds(0,'default'):
        expected=low_null(matrix,covered=covered,lengths=widths,threshold=tau)
        actual=reference_row_means(matrix,3,covered=covered,lengths=widths,threshold=tau)
        assert expected==actual
        assert actual['jaccard']==reference_means(covered_state(covered,widths,tau))
        assert actual['pearson']==reference_row_means(matrix,3)['pearson']
    for i in range(1,len(matrix)):
        old=similarity(matrix[0],matrix[i])
        new=similarity(matrix[0],matrix[i],covered_x=covered[0],covered_y=covered[i],lengths=widths)
        assert old==new

@pytest.mark.parametrize('manifest_depth',[1,6])
def test_shared_draw_grid_entrypoint_and_independent_nulls(manifest_depth):
    import csv,gzip,json
    from unittest.mock import patch
    import test_ref_depth_database as fxmod
    from bin.ref_depth_catalog import create_plan,collect
    from bin.ref_depth_worker import run_gcf
    from bin.reference_index import build_index
    from bin.query_manifest import parser,run,table
    fx=fxmod.DepthDatabaseTests('runTest');fx.setUp()
    try:
        config=dict(fx.config,storage_mode='native',window_size=7,depths=[],bins=None)
        db=fx.root/'native'
        create_plan(fx.records,fx.refs,db,config,fxmod.COLUMNS)
        assert run_gcf(db,fxmod.GCF1,fx.extractor)==0
        assert collect(db)['complete']
        index=fx.root/'index.sqlite';build_index(db,index,True)
        query=fx.root/'query.bam';query.write_bytes(b'unique synthetic query')
        manifest=fx.root/'manifest.tsv'
        table(manifest,('Sample','GCFID','Depth','CramPath'),[dict(Sample='external',GCFID=fxmod.GCF1,Depth=manifest_depth,CramPath=str(query))])
        pool=BlockPool.from_intervals(np.array([[0,1],[0,1],[21,23]]),31)
        base=['--manifest',str(manifest),'--index',str(index),'--window-size','7','--min-samples','0','--alpha','.05']
        outputs=[]
        for name,extra in [('single',[]),('grid',['--jaccard-coverage-thresholds','default'])]:
            out=fx.root/name
            args=parser(reference_distribution=True).parse_args(base+['--outdir',str(out)]+extra)
            with patch('bin.query_manifest.extract_pool',return_value=(pool,canonical_reference(fx.layout))):
                status=run(args)
                assert status==0, (out/'queries/q00000001/audit.json').read_text()
            outputs.append(out)
        a,b=[p/'queries/q00000001' for p in outputs]
        for name in ['draws.tsv.gz','null.tsv.gz','query_windows.npy','reference_windows.npy','reference_covered_bases.npy']:
            assert (a/name).read_bytes()==(b/name).read_bytes()
        with (b/'sensitivity.tsv').open() as f: rows=list(csv.DictReader(f,delimiter='\t'))
        assert len(rows)==7
        covered=np.load(b/'reference_covered_bases.npy'); widths=np.load(b/'window_lengths.npy')
        q=np.load(b/'query_covered_bases.npy')
        assert q.sum()==3 and np.load(b/'query_windows.npy').sum()==4
        pearson=[(r['pearson'],r['pearson_p'],r['pearson_status']) for r in rows]
        assert len(set(pearson))==1
        for i,row in enumerate(rows):
            tau=float(row['jaccard_coverage_threshold'])
            states=(covered>0) if tau==0 else covered/widths>=tau
            query_state=(q>0) if tau==0 else q/widths>=tau
            expected=[]
            for j,left in enumerate(states):
                values=[]
                for k,right in enumerate(states):
                    if j!=k:
                        union=(left|right).sum()
                        values.append((left&right).sum()/union if union else None)
                expected.append(np.mean(values) if None not in values else None)
            with gzip.open(b/f'jaccard_thresholds/t{i:02d}/null.tsv.gz','rt') as f:
                actual=[r for r in csv.DictReader(f,delimiter='\t') if r['metric']=='jaccard']
            for record,value in zip(actual,expected):
                if value is None: assert record['value']==''
                else: assert float(record['value'])==pytest.approx(value,abs=1e-15)
            zero=sum(not (query_state|state).any() for state in states)
            assert int(row['jaccard_zero_union_count'])==zero
            if zero: assert row['detection_status']=='NOT_EVALUATED' and row['jaccard_p']==''
    finally:
        fx.tearDown()


def test_collector_retains_grid_and_missing_receipts(tmp_path):
    import json,csv
    from bin.coverage_query_flow import parser,collect
    roster=[dict(query_id=f'q{i:08d}',Sample=f's{i}',GCFID='g',manifest_depth='1') for i in (1,2)]
    rp=tmp_path/'roster.json';rp.write_text(json.dumps(roster))
    receipts=tmp_path/'receipts';receipts.mkdir()
    result=dict(roster[0],status='COMPLETE',detection_status='DETECTED',result_source='TASK_RECEIPT')
    result['threshold_results']=[dict(result,jaccard_coverage_threshold=tau) for tau in thresholds(0,'default')]
    (receipts/'q00000001.json').write_text(json.dumps(result))
    args=parser().parse_args(['collect','--roster',str(rp),'--receipts',str(receipts),'--outdir',str(tmp_path/'out'),'--jaccard-coverage-thresholds','default'])
    collect(args)
    with (tmp_path/'out/sensitivity.tsv').open() as f: rows=list(csv.DictReader(f,delimiter='\t'))
    assert len(rows)==14
    assert all(r['status']=='ERROR' and r['detection_status']=='NOT_EVALUATED' for r in rows[7:])


def test_report_checks_cohort_membership_and_pearson(tmp_path):
    from bin.query_manifest import table
    from bin.jaccard_sensitivity_report import summarize
    fields=('query_id','Sample','GCFID','jaccard_coverage_threshold','pearson','pearson_p','pearson_status','detection_status','status','jaccard_zero_union_count','jaccard_null_zero_union_count')
    rows=[dict(zip(fields,('q1','s1','g',tau,.7,.9,'estimable','DETECTED','COMPLETE',0,0))) for tau in (0,.1)]
    tp=tmp_path/'tp.tsv';fp=tmp_path/'fp.tsv'
    table(tp,fields,rows);table(fp,fields,rows)
    summarize(tp,fp,tmp_path/'report')
    assert (tmp_path/'report/sensitivity.svg').is_file()
    rows[1]['pearson']=.8;table(fp,fields,rows)
    with pytest.raises(ValueError,match='Pearson changed'):summarize(tp,fp,tmp_path/'bad')
