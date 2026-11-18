from pathlib import Path
import sys
import numpy as np
import pytest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import Tools
from bin.coverage_core import BlockPool, window_edges, dense_coverage
from bin.concordance_core import evaluate_windows, simulate_grid, task_rng


def test_blocked_statistics_match_existing_library():
    rng=np.random.default_rng(84)
    matrix=rng.poisson(1.2,(27,81));test=rng.poisson(.9,81)
    for block in (1,7,64):
        actual=evaluate_windows(test,matrix,block_size=block)
        for name,function in [('jaccard',Tools.jaccard_similarity),('pearson',Tools.pearson_corr)]:
            similarities=function(matrix);np.fill_diagonal(similarities,0)
            distribution=similarities.sum(axis=1)/(len(matrix)-1)
            score=float(np.mean(function(test,matrix)))
            assert actual[name]['similarity']==pytest.approx(score,abs=2e-15)
            assert actual[name]['reference_mean']==pytest.approx(distribution.mean(),abs=2e-15)
            assert actual[name]['reference_sd']==pytest.approx(distribution.std(),abs=2e-15)
            assert actual[name]['p_value']==pytest.approx(Tools.z_test(score,distribution),abs=2e-13)


def test_non_estimable_is_distinct_from_zero_similarity():
    result=evaluate_windows([1,1,1],[[2,1,0],[0,1,2]])
    assert result['pearson']['status']=='non_estimable'
    assert result['pearson']['accepted'] is None
    assert result['jaccard']['status']=='estimable'
    with pytest.raises(ValueError,match='Zero coverage'):
        evaluate_windows([0,0,0],[[2,1,0],[0,1,2]])
    with pytest.raises(ValueError,match='N>=2'):
        evaluate_windows([1,0,0],[[2,1,0]])


@pytest.mark.parametrize('sampler',['lazy','permutation'])
def test_worker_count_batch_and_depth_order_do_not_change_draws(tmp_path,sampler):
    rng=np.random.default_rng(920)
    starts=rng.integers(0,470,size=1500)
    pool=BlockPool.from_intervals(np.column_stack([starts,starts+rng.integers(2,30,size=1500)]),500)
    edges=window_edges(500,23);baseline={}
    for name,workers,batch,depths in [('serial',1,7,[.1,1,3]),('thread',3,4,[3,.1,1]),('process',2,5,[1,3,.1])]:
        output=tmp_path/name;output.mkdir();observed={}
        def consume(index,depth,directory,matrix,summaries,calibration,seconds):
            observed[depth]=(np.array(matrix),summaries)
            for path in directory.glob('ids-*.npz'):
                with np.load(path) as ids:
                    for replicate in ids:
                        dense=dense_coverage(pool,ids[replicate])
                        expected=[dense[left:right].sum() for left,right in zip(edges[:-1],edges[1:])]
                        np.testing.assert_array_equal(matrix[int(replicate)],expected)
            assert len(summaries)==13 and all(row['actual_depth']>=depth for row in summaries)
        simulate_grid(pool,edges,depths,13,401,workers,name,batch,sampler,output,consume)
        if name=='serial':baseline=observed
        else:
            for depth,(matrix,summaries) in observed.items():
                np.testing.assert_array_equal(matrix,baseline[depth][0])
                assert summaries==baseline[depth][1]
    ref=task_rng(401,1,0,0).integers(2**32,size=20)
    test=task_rng(401,1,1,0).integers(2**32,size=20)
    assert not np.array_equal(ref,test)
