import hashlib
import importlib.util
from pathlib import Path
import sys
from types import SimpleNamespace

import numpy as np
import pytest

sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import Tools
from bin.coverage_core import (BlockPool, Reference, dense_coverage, extract_pool,
                               read_reference, select_units, window_edges, window_sums)

def make_pool():
    # Unit0 has two separated blocks, unit1 overlaps a window boundary, unit2 spans many windows.
    return BlockPool(np.array([0,8,3,1],dtype=np.uint32),np.array([2,10,7,15],dtype=np.uint32),
                     np.array([0,2,3,4],dtype=np.uint64),17)

@pytest.mark.parametrize("bins",[1,2,3,7,17])
@pytest.mark.parametrize("ids",[[0],[1],[2],[0,2],[0,1,2],[]])
def test_direct_matches_independent_per_base(bins,ids):
    pool=make_pool()
    edges=window_edges(17,bins)
    depth=dense_coverage(pool,ids)
    expected=np.array([depth[left:right].sum() for left,right in zip(edges[:-1],edges[1:])])
    actual=window_sums(pool,ids,edges)
    np.testing.assert_array_equal(actual,expected)
    assert actual.sum()==pool.unit_bases[np.asarray(ids,dtype=int)].sum()

def test_random_geometry():
    rng=np.random.default_rng(181)
    for _ in range(40):
        g=int(rng.integers(1,200))
        starts=rng.integers(0,g,size=40)
        ends=np.array([rng.integers(start+1,g+1) for start in starts])
        pool=BlockPool.from_intervals(np.stack([starts,ends],axis=1),g)
        ids=rng.choice(40,size=int(rng.integers(0,41)),replace=False)
        edges=window_edges(g,int(rng.integers(1,g+1)))
        dense=dense_coverage(pool,ids)
        np.testing.assert_array_equal(window_sums(pool,ids,edges),
            [dense[s:e].sum() for s,e in zip(edges[:-1],edges[1:])])

def test_complete_units_no_recycling():
    pool=BlockPool.from_intervals(np.array([[0,10]]),10)
    selected=select_units(pool,np.random.default_rng(1),target_bases=7)
    assert pool.unit_bases[selected].sum()==10
    with pytest.raises(ValueError):
        select_units(pool,np.random.default_rng(1),target_bases=11)
    assert len(select_units(pool,np.random.default_rng(1),count=0))==0

def test_fixed_count_fraction_seed_independence():
    pool=BlockPool.from_intervals(np.array([[i,i+1] for i in range(20)]),20)
    def draws():
        rng=np.random.default_rng(88)
        return [select_units(pool,rng,count=5).tolist() for _ in range(20)]
    first=draws()
    assert first==draws()
    assert all(len(set(ids))==5 for ids in first)
    assert any(set(a)&set(b) for a,b in zip(first,first[1:]))
    assert len(select_units(pool,np.random.default_rng(2),fraction=.31))==6

@pytest.mark.parametrize("kwargs",[{"count":-1},{"count":4},{"fraction":1.1},{"target_bases":-1},
                                    {"count":1,"fraction":.2},{"count":1.5}])
def test_invalid_targets(kwargs):
    with pytest.raises(ValueError):
        select_units(make_pool(),np.random.default_rng(1),**kwargs)

def test_pool_validation_and_empty(tmp_path):
    with pytest.raises(ValueError):
        BlockPool.from_intervals(np.array([[3,2]]),5)
    with pytest.raises(ValueError):
        window_edges(100,101)
    with pytest.raises(ValueError):
        window_sums(make_pool(),[1,1],window_edges(17,4))
    empty=BlockPool.from_intervals([],5)
    np.testing.assert_array_equal(window_sums(empty,[],window_edges(5,3)),[0,0,0])
    pool=make_pool()
    path=tmp_path/"pool.npz"
    pool.save(path)
    loaded=BlockPool.load(path)
    np.testing.assert_array_equal(loaded.indptr,pool.indptr)
    np.testing.assert_array_equal(window_sums(loaded,[0,2],window_edges(17,4)),window_sums(pool,[0,2],window_edges(17,4)))

def test_legacy_bin_and_statistical_degeneracy():
    assert len(Tools.bin_coverage(np.ones(1001,dtype=int),bin_count=1000))==1000
    assert Tools.bin_coverage(np.array([2_000_000_000]*2,dtype=np.int64),bin_count=1)[0]==4_000_000_000
    assert Tools.bin_coverage(np.ones(10),bin_size=2).tolist()==[2]*5
    assert np.isnan(Tools.pearson_corr([1,1],[2,2]))
    assert Tools.jaccard_similarity([0,0],[0,0])==1
    assert Tools.z_test(.75,[.5,.5])==1
    with pytest.raises(ValueError):
        Tools.concordance_test_low_cov(np.zeros(8),np.empty((0,2),dtype=int),8,n=3)
    np.testing.assert_allclose(Tools.normalize_coverage([.5,1.5],target_sum=4),[1,3])

def test_legacy_state_is_complete_and_independent():
    state=Tools.SubsamplingState(np.array([[0,10]]),np.random.default_rng(4))
    assert Tools.simulate_coverage_subsampling(state,10,7).sum()==10
    with pytest.raises(ValueError):
        Tools.simulate_coverage_subsampling(state,10,11)

def test_fasta_reference(tmp_path):
    path=tmp_path/"ref.fna"
    path.write_text(">a\nACGT\n>b\nTT\n")
    ref=read_reference(path)
    assert ref.offsets=={"a":0,"b":4}
    assert ref.length==6
    assert ref.md5[0]==hashlib.md5(b"ACGT").hexdigest()
    path.write_text(">a\nA\n>a\nC\n")
    with pytest.raises(ValueError):
        read_reference(path)

def test_stream_mock_keeps_mates_gaps_and_read_groups(tmp_path,monkeypatch):
    alignment=tmp_path/"mock.bam"
    alignment.write_bytes(b"BAM1")
    reference=tmp_path/"ref.fna"
    reference.write_text(">a\n"+"A"*20+"\n")
    def record(q,flag,blocks,rg):
        return SimpleNamespace(is_unmapped=False,is_secondary=False,is_supplementary=False,
            is_paired=flag!=0,is_read1=flag==1,is_read2=flag==2,mapping_quality=60,is_duplicate=False,
            is_qcfail=False,reference_name="a",query_name=q,get_blocks=lambda:blocks,
            has_tag=lambda key:True,get_tag=lambda key:rg)
    records=[record("pair",1,[(0,3)],"A"),record("pair",2,[(8,11)],"A"),
             record("pair",1,[(14,16)],"B"),record("gap",0,[(3,5),(7,9)],"A")]
    class FakeBam:
        def __init__(self,*args,**kwargs):
            self.header=SimpleNamespace(to_dict=lambda:{"SQ":[{"SN":"a","LN":20}]})
        def __enter__(self):return self
        def __exit__(self,*args):pass
        def fetch(self,until_eof=False):
            assert until_eof is True
            yield from records
    monkeypatch.setitem(sys.modules,"pysam",SimpleNamespace(AlignmentFile=FakeBam))
    read_pool,_=extract_pool(alignment,reference)
    first_pool,_=extract_pool(alignment,reference,mate_policy="read1")
    fragment_pool,_=extract_pool(alignment,reference,sampling_unit="fragment")
    assert (len(read_pool),len(first_pool),len(fragment_pool))==(4,3,3)
    assert read_pool.total_bases==12
    assert fragment_pool.unit_bases.tolist()==[6,2,4]
    assert read_pool.indptr.tolist()==[0,1,2,3,5]

def test_real_cram_partial_fields_same_as_full_reference(tmp_path,monkeypatch):
    pysam=pytest.importorskip("pysam")
    reference=tmp_path/"ref.fna"
    reference.write_text(">a\n"+"A"*20+"\n>b\n"+"C"*16+"\n")
    pysam.faidx(str(reference))
    header={"HD":{"VN":"1.6"},"SQ":[{"SN":"a","LN":20},{"SN":"b","LN":16}],
            "RG":[{"ID":"A"},{"ID":"B"}]}
    alignment=tmp_path/"synthetic.cram"
    rows=[("pair",99,0,0,"4M","AAAA","A"),("pair",147,0,8,"4M","AAAA","A"),
          ("pair",0,0,14,"2M","AA","B"),("gap",0,0,4,"2M2D2M","AAAA","A"),
          ("skip",0,1,0,"2M2N2M","CCCC","A")]
    with pysam.AlignmentFile(str(alignment),"wc",header=header,reference_filename=str(reference)) as bam:
        for name,flag,tid,pos,cigar,seq,rg in rows:
            read=pysam.AlignedSegment(bam.header)
            read.query_name=name
            read.flag=flag
            read.reference_id=tid
            read.reference_start=pos
            read.mapping_quality=60
            read.cigarstring=cigar
            read.query_sequence=seq
            read.query_qualities=pysam.qualitystring_to_array("I"*len(seq))
            read.set_tag("RG",rg)
            bam.write(read)
    full,_=extract_pool(alignment,reference)
    assert len(full)==5 and full.total_bases==18
    reference.unlink()
    monkeypatch.setenv("REF_PATH","/nonexistent/coverage-map-no-reference")
    monkeypatch.setenv("REF_CACHE","/nonexistent/coverage-map-no-cache")
    partial,_=extract_pool(alignment,reference_index=tmp_path/"ref.fna.fai")
    fragments,_=extract_pool(alignment,reference_index=tmp_path/"ref.fna.fai",sampling_unit="fragment")
    assert len(fragments)==4 and sorted(fragments.unit_bases.tolist())==[2,4,4,8]
    np.testing.assert_array_equal(partial.starts,full.starts)
    np.testing.assert_array_equal(partial.ends,full.ends)
    np.testing.assert_array_equal(partial.indptr,full.indptr)
    expected=np.zeros(36,dtype=np.int64)
    for start,end in [(0,4),(8,12),(14,16),(4,6),(8,10),(20,22),(24,26)]:
        expected[start:end]+=1
    np.testing.assert_array_equal(dense_coverage(partial,np.arange(len(partial))),expected)
