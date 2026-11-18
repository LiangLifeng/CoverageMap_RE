"""Coverage helpers with explicit complete-unit downsampling.

The public names of the original Tools.py are retained. Read-pool extraction now
returns BlockPool, preserving CIGAR blocks and paired sampling units. Legacy (N,2)
integer arrays are still accepted by SubsamplingState. See docs/IO_CONTRACT.md.
"""
from __future__ import annotations
import hashlib
import json
from pathlib import Path
from statistics import NormalDist
import tempfile
import os

import numpy as np
from bin.coverage_core import (BlockPool, read_reference, extract_pool, window_edges,
                               window_sums, select_units)

def build_genome_offset(ref_genome):
    ref=read_reference(ref_genome)
    if "length" in ref.names:
        raise ValueError("The legacy offset API reserves the contig name 'length'")
    return {**ref.offsets,"length":ref.length}

def _iter_primary_reads(bam, mate_policy="both"):
    if mate_policy not in ("both","read1","read2"):
        raise ValueError("Invalid mate policy")
    for read in bam.fetch(until_eof=True):
        if read.is_secondary or read.is_supplementary or read.is_unmapped:
            continue
        if read.is_paired and ((mate_policy=="read1" and not read.is_read1) or (mate_policy=="read2" and not read.is_read2)):
            continue
        yield read

def extract_read_intervals(cram_file, ref_genome, genome_offset, **policies):
    pool, ref=extract_pool(cram_file,ref_genome,**policies)
    if build_genome_offset(ref_genome) != genome_offset:
        raise ValueError("Reference layout and supplied offsets disagree")
    return pool

def _depth(pool, ids):
    starts,ends=pool.selected_blocks(ids)
    diff=np.zeros(pool.genome_len+1,dtype=np.int64)
    np.add.at(diff,starts.astype(np.intp),1)
    np.add.at(diff,ends.astype(np.intp),-1)
    return np.cumsum(diff[:-1],dtype=np.int64)

def cram_to_depth_array(cram_file,ref_genome,genome_offset,dtype=np.int64,**policies):
    pool=extract_read_intervals(cram_file,ref_genome,genome_offset,**policies)
    depth=_depth(pool,np.arange(len(pool)))
    kind=np.dtype(dtype)
    if kind.kind in "iu" and depth.size and depth.max() > np.iinfo(kind).max:
        raise OverflowError("Coverage exceeds requested output dtype")
    return depth.astype(kind)

def _file_hash(path):
    value=hashlib.sha256()
    with open(path,"rb") as handle:
        for chunk in iter(lambda:handle.read(1024*1024),b""):
            value.update(chunk)
    return value.hexdigest()

def _cache_path(sample,cram_file,ref_genome,genome_offset,cache_dir,policies,suffix):
    if not sample or Path(sample).name != sample or sample in (".",".."):
        raise ValueError("sample must be one filename component")
    identity={"version":"coverage-map-1","alignment":_file_hash(cram_file),"reference":_file_hash(ref_genome),
              "layout":genome_offset,"policies":policies}
    digest=hashlib.sha256(json.dumps(identity,sort_keys=True).encode()).hexdigest()
    directory=Path(cache_dir)
    directory.mkdir(parents=True,exist_ok=True)
    return directory/f"{sample}.{digest}.{suffix}"

def get_cached_reads_pool(sample,cram_file,ref_genome,genome_offset,cache_dir,**policies):
    path=_cache_path(sample,cram_file,ref_genome,genome_offset,cache_dir,policies,"reads.npz")
    if path.is_file():
        return BlockPool.load(path)
    pool=extract_read_intervals(cram_file,ref_genome,genome_offset,**policies)
    with tempfile.NamedTemporaryFile(dir=path.parent,suffix=".npz",delete=False) as handle:
        temporary=Path(handle.name)
    try:
        pool.save(temporary)
        os.replace(temporary,path)
    finally:
        temporary.unlink(missing_ok=True)
    return pool

def get_cached_depth_array(sample,cram_file,ref_genome,genome_offset,cache_dir,**policies):
    path=_cache_path(sample,cram_file,ref_genome,genome_offset,cache_dir,policies,"depth.npy")
    if path.is_file():
        return np.load(path,allow_pickle=False)
    depth=cram_to_depth_array(cram_file,ref_genome,genome_offset,**policies)
    with tempfile.NamedTemporaryFile(dir=path.parent,suffix=".npy",delete=False) as handle:
        temporary=Path(handle.name)
        np.save(handle,depth)
    os.replace(temporary,path)
    return depth

def normalize_coverage(values,target_sum=10**6):
    values=np.asarray(values,dtype=np.float64)
    if np.any(~np.isfinite(values)) or np.any(values<0) or not np.isfinite(target_sum) or target_sum<0:
        raise ValueError("Coverage and target sum must be finite and nonnegative")
    total=values.sum()
    if total<=0:
        raise ValueError("Cannot normalize a zero-sum coverage array")
    return values*target_sum/total

def bin_coverage(cov,bin_count=None,bin_size=None,mode="sum",drop_tail=False):
    cov=np.asarray(cov)
    if cov.ndim==0 or mode not in ("sum","mean"):
        raise ValueError("Require an array and sum/mean mode")
    length=cov.shape[-1]
    if length==0:
        return cov.copy()
    if bin_count is not None and bin_size is not None:
        raise ValueError("Specify exactly one window definition")
    if bin_size is not None:
        if not isinstance(bin_size,(int,np.integer)) or bin_size<=0:
            raise ValueError("bin_size must be a positive integer")
        edges=np.arange(0,length+1,bin_size,dtype=np.int64)
        if edges[-1]!=length and not drop_tail:
            edges=np.append(edges,length)
        if len(edges)<2:
            return np.empty((*cov.shape[:-1],0),dtype=np.float64 if mode=="mean" else np.int64)
        cov=cov[...,:edges[-1]]
    else:
        count=min(length,1000) if bin_count is None else bin_count
        edges=window_edges(length,count)
    dtype=np.float64 if cov.dtype.kind=="f" else np.int64
    sums=np.add.reduceat(cov,edges[:-1].astype(np.intp),axis=-1,dtype=dtype)
    return sums/np.diff(edges) if mode=="mean" else sums

def simulate_uniform_coverage(genome_len,depth,read_len=150,region=None,rng=None):
    if genome_len<=0 or read_len<=0 or not np.isfinite(depth) or depth<0:
        raise ValueError("Invalid genome/read length or depth")
    r0,r1=(0,genome_len) if region is None else region
    if not 0<=r0<r1<=genome_len or r1-r0<read_len:
        raise ValueError("Invalid placement region")
    rng=np.random.default_rng() if rng is None else rng
    count=int(round(depth*genome_len/read_len))
    starts=rng.integers(r0,r1-read_len+1,size=count)
    diff=np.zeros(genome_len+1,dtype=np.int64)
    np.add.at(diff,starts,1)
    np.add.at(diff,starts+read_len,-1)
    return np.cumsum(diff[:-1])

def z_test(test_value,ref_value_pool,one_side=True):
    ref=np.asarray(ref_value_pool,dtype=float)
    if not ref.size or np.any(~np.isfinite(ref)) or not np.isfinite(test_value):
        raise ValueError("Non-estimable or empty similarity distribution")
    mu,sigma=float(ref.mean()),float(ref.std())
    if sigma==0:
        return float(test_value>=mu) if one_side else float(test_value==mu)
    z=(test_value-mu)/sigma
    normal=NormalDist()
    return normal.cdf(z) if one_side else 2*normal.cdf(-abs(z))

def pearson_corr(x,y=None,eps=1e-12):
    x=np.asarray(x,dtype=float)
    if x.ndim==1 and y is None:
        raise ValueError("y required for vector input")
    x=np.atleast_2d(x)
    if x.ndim!=2 or np.any(~np.isfinite(x)):
        raise ValueError("Expected finite vectors/matrices")
    y=x if y is None else np.atleast_2d(np.asarray(y,dtype=float))
    if y.ndim!=2 or y.shape[1]!=x.shape[1] or np.any(~np.isfinite(y)):
        raise ValueError("Incompatible vectors")
    xc=x-x.mean(axis=1,keepdims=True)
    yc=y-y.mean(axis=1,keepdims=True)
    denominator=np.linalg.norm(xc,axis=1)[:,None]*np.linalg.norm(yc,axis=1)[None,:]
    corr=np.full(denominator.shape,np.nan)
    np.divide(xc@yc.T,denominator,out=corr,where=denominator>0)
    corr=np.clip(corr,-1,1)
    return _squeeze_similarity(corr)

def _squeeze_similarity(matrix):
    if matrix.shape==(1,1):
        return float(matrix[0,0])
    if matrix.shape[0]==1:
        return matrix[0]
    if matrix.shape[1]==1:
        return matrix[:,0]
    return matrix

def jaccard_similarity(x,y=None,threshold=0):
    x=np.asarray(x)
    if x.ndim==1 and y is None:
        raise ValueError("y required for vector input")
    x=np.atleast_2d(x)
    y=x if y is None else np.atleast_2d(np.asarray(y))
    if x.ndim!=2 or y.ndim!=2 or x.shape[1]!=y.shape[1] or np.any(~np.isfinite(x)) or np.any(~np.isfinite(y)):
        raise ValueError("Expected compatible finite vectors")
    xb,yb=x>threshold,y>threshold
    sim=np.empty((len(xb),len(yb)),dtype=float)
    for i,row in enumerate(xb):
        union=np.logical_or(row,yb).sum(axis=1)
        intersection=np.logical_and(row,yb).sum(axis=1)
        sim[i]=np.divide(intersection,union,out=np.ones(len(yb)),where=union>0)
    return _squeeze_similarity(sim)

class SubsamplingState:
    """RNG holder: each call independently samples the full original pool."""
    def __init__(self,intervals,rng):
        self.intervals,self.rng=intervals,rng
        self.last_selected_ids=None
        self._cached_pool=None
    def reset(self):
        self.last_selected_ids=None
    def pool(self,genome_len):
        if isinstance(self.intervals,BlockPool):
            if self.intervals.genome_len!=genome_len:
                raise ValueError("Pool reference length mismatch")
            return self.intervals
        if self._cached_pool is None:
            self._cached_pool=BlockPool.from_intervals(self.intervals,genome_len)
        elif self._cached_pool.genome_len!=genome_len:
            raise ValueError("Pool reference length changed")
        return self._cached_pool

def simulate_coverage_subsampling(state,genome_len,target_bases):
    pool=state.pool(genome_len)
    ids=select_units(pool,state.rng,target_bases=target_bases)
    state.last_selected_ids=ids
    return _depth(pool,ids)

def concordance_test_low_cov(test_cov,ref_reads_pool,genome_len,n=100,alpha=.01,*,seed=None):
    if n<2 or not 0<alpha<1:
        raise ValueError("Require n>=2 and 0<alpha<1")
    raw=np.asarray(test_cov)
    if raw.shape!=(genome_len,) or np.any(raw<0) or np.any(~np.isfinite(raw)):
        raise ValueError("test_cov must be finite nonnegative per-base coverage")
    if raw.sum()==0:
        raise ValueError("Zero coverage is descriptive only; Pearson concordance is non-estimable")
    pool=ref_reads_pool if isinstance(ref_reads_pool,BlockPool) else BlockPool.from_intervals(ref_reads_pool,genome_len)
    if pool.genome_len!=genome_len or len(pool)==0:
        raise ValueError("Missing/mismatched reference evidence")
    edges=window_edges(genome_len,min(1000,genome_len))
    test=bin_coverage(raw,bin_count=len(edges)-1)
    rng=np.random.default_rng(seed)
    simulated=np.empty((n,len(edges)-1),dtype=np.int64)
    for i in range(n):
        ids=select_units(pool,rng,target_bases=float(raw.sum()))
        simulated[i]=window_sums(pool,ids,edges)
    m_j=jaccard_similarity(simulated)
    m_r=pearson_corr(simulated)
    np.fill_diagonal(m_j,0)
    np.fill_diagonal(m_r,0)
    p_j=z_test(float(np.mean(jaccard_similarity(test,simulated))),m_j.sum(axis=1)/(n-1))
    p_r=z_test(float(np.mean(pearson_corr(test,simulated))),m_r.sum(axis=1)/(n-1))
    return p_j>=alpha,p_j,p_r>=alpha,p_r

def concordance_test_high_cov(test_cov,ref_cov_pool,ref_simi_jcd,ref_simi_r,alpha=.01,n_samples=1000,*,seed=None):
    if not 0<alpha<1 or n_samples<1:
        raise ValueError("Invalid alpha/n_samples")
    rng=np.random.default_rng(seed)
    pools=[]
    for values in (ref_simi_jcd,ref_simi_r):
        values=np.asarray(values,dtype=float)
        if values.size>n_samples:
            values=values[rng.choice(values.size,n_samples,replace=False)]
        pools.append(values)
    test=bin_coverage(test_cov)
    ref=bin_coverage(ref_cov_pool)
    p_j=z_test(float(np.max(jaccard_similarity(test,ref))),pools[0])
    p_r=z_test(float(np.max(pearson_corr(test,ref))),pools[1])
    return p_j>=alpha,p_j,p_r>=alpha,p_r

_WORKER_REF_READS_POOL=None
_WORKER_GENOME_LEN=None

def init_worker(ref_reads_pool,genome_len):
    global _WORKER_REF_READS_POOL,_WORKER_GENOME_LEN
    _WORKER_REF_READS_POOL,_WORKER_GENOME_LEN=ref_reads_pool,genome_len

def parallel_subsampling_test(args):
    if _WORKER_REF_READS_POOL is None or _WORKER_GENOME_LEN is None:
        raise RuntimeError("init_worker must be called first")
    sample,coverage_test=args
    return (sample,)+concordance_test_low_cov(coverage_test,_WORKER_REF_READS_POOL,_WORKER_GENOME_LEN)
