#!/usr/bin/env python3
"""Validate direct window coverage against independent dense geometry on NH."""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import resource
import time
import numpy as np
from coverage_core import BlockPool, dense_coverage, extract_pool, select_units, window_edges, window_sums

def compare(pool,ids,edges):
    t=time.perf_counter()
    fast=window_sums(pool,ids,edges)
    fast_seconds=time.perf_counter()-t
    t=time.perf_counter()
    dense=dense_coverage(pool,ids)
    expected=np.array([dense[s:e].sum() for s,e in zip(edges[:-1],edges[1:])])
    dense_seconds=time.perf_counter()-t
    np.testing.assert_array_equal(fast,expected)
    assert fast.sum()==pool.unit_bases[ids].sum()
    return {"units":len(ids),"aligned_bases":int(fast.sum()),"direct_seconds":fast_seconds,
            "dense_seconds":dense_seconds,"equal":True}

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--alignment",type=Path,required=True)
    parser.add_argument("--reference-index",type=Path,required=True)
    parser.add_argument("--output",type=Path,required=True)
    parser.add_argument("--seed",type=int,default=401)
    args=parser.parse_args()
    if args.output.exists():raise FileExistsError(args.output)
    rng=np.random.default_rng(args.seed)
    records={}
    for unit,mate in (("read","both"),("read","read1"),("fragment","both")):
        t=time.perf_counter()
        pool,ref=extract_pool(args.alignment,reference_index=args.reference_index,sampling_unit=unit,mate_policy=mate)
        parse_seconds=time.perf_counter()-t
        edges=window_edges(pool.genome_len,1000)
        comparisons=[compare(pool,np.arange(len(pool)),edges)]
        if pool.total_bases>=pool.genome_len:
            comparisons += [compare(pool,select_units(pool,rng,target_bases=pool.genome_len),edges) for _ in range(5)]
        records[f"{unit}_{mate}"]={"units_total":len(pool),"blocks_total":len(pool.starts),
            "total_aligned_bases":pool.total_bases,"reference_length":pool.genome_len,
            "pool_numeric_bytes":pool.nbytes,"parse_seconds":parse_seconds,"comparisons":comparisons}
    # Declared synthetic scale experiment, separate from the real genome evidence.
    g=5_000_000
    starts=rng.integers(0,g-150,size=100_000,dtype=np.int64)
    stress=BlockPool.from_intervals(np.stack([starts,starts+150],axis=1),g)
    selected=select_units(stress,rng,count=10_000)
    records["synthetic_5Mb"]=compare(stress,selected,window_edges(g,1000))
    records["peak_rss_kib"]=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    records["scope"]="geometry/mass/performance validation; no concordance p-value calibration"
    args.output.write_text(json.dumps(records,indent=2)+"\n")
    print(json.dumps(records,indent=2))
if __name__=="__main__":main()
