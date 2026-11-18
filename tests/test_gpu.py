from pathlib import Path
import sys
import numpy as np
import pytest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from bin.coverage_core import BlockPool, window_edges, window_sums, dense_coverage
from bin.concordance_core import evaluate_windows, simulate_grid

@pytest.fixture
def cuda():
    cp=pytest.importorskip('cupy')
    try:
        assert cp.cuda.runtime.getDeviceCount()>0
    except cp.cuda.runtime.CUDARuntimeError:
        pytest.skip('CUDA device unavailable')
    return cp


def test_statistics_backend_rejects_unknown():
    with pytest.raises(ValueError,match='Unknown statistics backend'):
        evaluate_windows([1,2,3],[[3,2,1],[1,3,2]],backend='invalid')

@pytest.mark.parametrize('bins',[1,7,31])
def test_gpu_ragged_blocks_and_dense_exact(cuda,bins):
    from bin.gpu_core import GpuWindowAccumulator
    # Unit 0 spans a wide interval plus a separated block; unequal widths and
    # exact edge ends expose difference-prefix errors and missing CIGAR gaps.
    pool=BlockPool(np.array([0,85,11,0,100]),np.array([73,99,57,101,101]),np.array([0,2,3,4,5]),101)
    edges=window_edges(101,bins);gpu=GpuWindowAccumulator(pool,edges)
    selected=[np.array([0,2]),np.array([],dtype=np.int64),np.array([3,1,0]),np.arange(4)]
    rows=gpu.count_batch(selected)
    np.testing.assert_array_equal(rows,gpu._count_selected_batch(selected))
    for row,ids in zip(rows,selected):
        dense=dense_coverage(pool,ids)
        np.testing.assert_array_equal(row,[dense[s:e].sum() for s,e in zip(edges[:-1],edges[1:])])
        np.testing.assert_array_equal(row,window_sums(pool,ids,edges))
    with pytest.raises(ValueError,match='Invalid'):
        gpu.count_batch([np.array([0,0])])
    assert gpu.metrics()['transfer_bytes']>0

@pytest.mark.parametrize('block',[1,7,64])
def test_gpu_statistics(cuda,block):
    rng=np.random.default_rng(8);x=rng.poisson(.9,(27,81));test=rng.poisson(1.2,81)
    cpu=evaluate_windows(test,x,block_size=block)
    gpu=evaluate_windows(test,x,block_size=block,backend='gpu')
    for name in ['jaccard','pearson']:
        assert cpu[name]['accepted']==gpu[name]['accepted']
        assert cpu[name]['status']==gpu[name]['status']
        for key in ['similarity','reference_mean','reference_sd','p_value']:
            assert gpu[name][key]==pytest.approx(cpu[name][key],abs=2e-11)
    assert evaluate_windows(np.ones(81),x,backend='gpu')['pearson']['status']=='non_estimable'
    with pytest.raises(ValueError,match='Zero coverage'):
        evaluate_windows(np.zeros(81),x,backend='gpu')

@pytest.mark.parametrize('sampler',['lazy','permutation'])
def test_gpu_cpu_sampling_same_ids(cuda,tmp_path,sampler):
    from bin.gpu_core import GpuWindowAccumulator
    starts=np.arange(200)*2
    pool=BlockPool.from_intervals(np.column_stack([starts,starts+7]),407)
    edges=window_edges(407,19);gpu=GpuWindowAccumulator(pool,edges);actual=[]
    for label,accumulator in [('cpu',None),('gpu',gpu)]:
        out=tmp_path/label;out.mkdir();records={}
        def complete(index,depth,directory,matrix,summaries,calibration,seconds):
            records[depth]=(np.array(matrix),summaries)
        simulate_grid(pool,edges,[.1,1,2],9,73,3,'thread',4,sampler,out,complete,gpu=accumulator)
        actual.append(records)
    for depth in actual[0]:
        np.testing.assert_array_equal(actual[0][depth][0],actual[1][depth][0])
        assert actual[0][depth][1]==actual[1][depth][1]
    with pytest.raises(ValueError,match='not process'):
        simulate_grid(pool,edges,[1],2,73,2,'process',2,sampler,tmp_path,lambda *a:None,gpu=gpu)


def test_gpu_edges_reject_unsigned_descending_before_cuda_import():
    from bin.gpu_core import GpuWindowAccumulator
    pool=BlockPool.from_intervals(np.array([[0,101]]),101)
    with pytest.raises(ValueError,match='partition'):
        GpuWindowAccumulator(pool,np.array([0,51,50,101],dtype=np.uint64))
