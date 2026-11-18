import csv
from collections import Counter
import importlib
import itertools
import json
from pathlib import Path
import sys

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "bin"))
from coverage_core import BlockPool, Reference, select_units, choose_depth_sampler
cli = importlib.import_module("coverage_map")


def pool_fixture():
    return BlockPool.from_intervals(np.array([[0,1],[3,5],[7,10]]), 10)


def test_lazy_uniform_permutation_and_weighted_stopping():
    pool = pool_fixture()
    class RankedRNG:
        def __init__(self, ranks): self.ranks = ranks
        def integers(self, low, high, dtype):
            assert low == 0
            values = np.array(self.ranks, dtype=dtype)
            assert np.all(values < high)
            return values
    # Enumerate every equally likely Fisher-Yates rank path, not a noisy frequency test.
    for target in (1, 2, 3, 5, 6):
        actual = Counter()
        expected = Counter()
        for ranks in itertools.product(range(3), range(2), range(1)):
            ids = select_units(pool, RankedRNG(ranks), target_bases=target, depth_sampler="lazy")
            actual[tuple(ids)] += 1
        for order in itertools.permutations(range(3)):
            prefix=[]; mass=0
            for unit in order:
                prefix.append(unit); mass += int(pool.unit_bases[unit])
                if mass >= target: break
            expected[tuple(prefix)] += 1
        assert actual == expected


def test_lazy_reproducibility_batches_and_complete_units():
    pool = BlockPool.from_intervals(np.stack([np.arange(300), np.arange(300)+2], axis=1), 302)
    for target in (0, 1, 257, 600):
        a=select_units(pool,np.random.default_rng(81),target_bases=target,depth_sampler="lazy")
        b=select_units(pool,np.random.default_rng(81),target_bases=target,depth_sampler="lazy")
        np.testing.assert_array_equal(a,b)
        assert len(a) == len(set(a))
        assert int(pool.unit_bases[a].sum()) >= target
        if len(a): assert int(pool.unit_bases[a[:-1]].sum()) < target
    rng=np.random.default_rng(55)
    draws=[tuple(select_units(pool,rng,target_bases=20,depth_sampler="lazy")) for _ in range(5)]
    assert len(set(draws)) > 1
    with pytest.raises(ValueError):select_units(pool,rng,target_bases=601,depth_sampler="lazy")
    with pytest.raises(ValueError):select_units(pool,rng,target_bases=2,depth_sampler="missing")


def test_historical_permutation_seed_ids_unchanged():
    pool=pool_fixture()
    rng=np.random.default_rng(401)
    oracle=np.random.default_rng(401)
    for _ in range(20):
        order=oracle.permutation(len(pool))
        stop=int(np.searchsorted(np.cumsum(pool.unit_bases[order], dtype=np.int64),4))+1
        np.testing.assert_array_equal(select_units(pool,rng,target_bases=4),order[:stop])


@pytest.mark.parametrize("sampler",["permutation","lazy"])
def test_cli_text_binary_same_ids_counts_edges_and_stages(tmp_path,monkeypatch,sampler):
    pool=pool_fixture();ref=Reference(("synthetic",),(10,),("",))
    monkeypatch.setattr(cli,"extract_pool",lambda *a,**kw:(pool,ref))
    alignment=tmp_path/"fixture";alignment.write_bytes(b"synthetic fixture only")
    fai=tmp_path/"fixture.fai";fai.write_text("synthetic\t10\t0\t0\t0\n")
    records=[]
    for fmt in ("tsv","npy"):
        out=tmp_path/fmt
        monkeypatch.setattr(sys,"argv",["coverage_map.py","--alignment",str(alignment),"--reference-index",str(fai),
            "--outdir",str(out),"--depth","0.4","--bins","3","--replicates","7","--seed","401",
            "--depth-sampler",sampler,"--output-format",fmt,"--save-pool"])
        cli.main()
        p=json.loads((out/"provenance.json").read_text())
        assert p["depth_sampler"] == sampler and p["output_format"] == fmt
        assert all(v>=0 for v in p["stages_seconds"].values())
        records.append(p["replicates"])
    assert records[0] == records[1]
    edges=np.load(tmp_path/"npy/edges.npy")
    matrix=np.load(tmp_path/"npy/simulations.npy",mmap_mode="r")
    assert matrix.shape==(7,3) and matrix.dtype==np.dtype("int64")
    for filename,values in (("simulations.tsv",matrix),("full.tsv",np.load(tmp_path/"npy/full.npy")[None,:])):
        with (tmp_path/"tsv"/filename).open() as f: rows=list(csv.DictReader(f,delimiter="\t"))
        for i,r in enumerate(rows):
            replicate=i//3;window=i%3
            assert int(r["start"])==edges[window] and int(r["end"])==edges[window+1]
            assert int(r["aligned_bases"])==values[replicate,window]
            assert float(r["mean_depth"])==pytest.approx(values[replicate,window]/(edges[window+1]-edges[window]))


def test_invalid_target_does_not_create_partial_output(tmp_path,monkeypatch):
    pool=pool_fixture();ref=Reference(("synthetic",),(10,),("",))
    monkeypatch.setattr(cli,"extract_pool",lambda *a,**kw:(pool,ref))
    out=tmp_path/"bad"
    monkeypatch.setattr(sys,"argv",["coverage_map.py","--alignment",str(tmp_path/"unused"),
        "--reference-index",str(tmp_path/"unused.fai"),"--outdir",str(out),"--depth","100","--depth-sampler","lazy"])
    with pytest.raises(ValueError):cli.main()
    assert not out.exists()


@pytest.mark.parametrize("lazy_cost,perm_cost,chosen",[(1.,3.,"lazy"),(3.,1.,"permutation"),(0.95,1.,"permutation")])
def test_auto_paired_calibration_and_rng_isolation(monkeypatch,lazy_cost,perm_cost,chosen):
    import coverage_core as core
    original = core.select_units
    clock = [0.]
    calls = []
    def timed(pool,rng,**kw):
        calls.append(kw["depth_sampler"])
        ids = original(pool,rng,**kw)
        clock[0] += lazy_cost if kw["depth_sampler"] == "lazy" else perm_cost
        return ids
    monkeypatch.setattr(core,"select_units",timed)
    monkeypatch.setattr(core.time,"perf_counter",lambda:clock[0])
    actual,info=choose_depth_sampler(pool_fixture(),4)
    assert actual == chosen and info["selected"] == chosen
    assert calls == ["lazy","permutation","permutation","lazy","lazy","permutation"]
    assert info["median_seconds"] == pytest.approx({"lazy":lazy_cost,"permutation":perm_cost})
    rng=np.random.default_rng(401);oracle=np.random.default_rng(401)
    np.testing.assert_array_equal(original(pool_fixture(),rng,target_bases=4,depth_sampler=actual),
                                  original(pool_fixture(),oracle,target_bases=4,depth_sampler=actual))


def test_auto_zero_and_invalid_no_calibration():
    chosen,info=choose_depth_sampler(pool_fixture(),0)
    assert chosen=="permutation" and not info["trials_seconds"]["lazy"]
    with pytest.raises(ValueError):choose_depth_sampler(pool_fixture(),7)


@pytest.mark.parametrize("selected",["lazy","permutation"])
def test_cli_auto_resolves_once_and_matches_explicit(tmp_path,monkeypatch,selected):
    pool=pool_fixture();ref=Reference(("synthetic",),(10,),("",))
    monkeypatch.setattr(cli,"extract_pool",lambda *a,**kw:(pool,ref))
    calls=[]
    def choose(pool,target):calls.append(target);return selected,{"selected":selected,"policy":"fixture"}
    monkeypatch.setattr(cli,"choose_depth_sampler",choose)
    alignment=tmp_path/"fixture";alignment.write_bytes(b"synthetic")
    fai=tmp_path/"fixture.fai";fai.write_text("synthetic\t10\t0\t0\t0\n")
    records=[]
    for sampler in ("auto",selected):
        out=tmp_path/sampler
        monkeypatch.setattr(sys,"argv",["coverage_map.py","--alignment",str(alignment),"--reference-index",str(fai),
            "--outdir",str(out),"--depth","0.4","--bins","3","--replicates","7","--seed","401",
            "--depth-sampler",sampler,"--output-format","npy"])
        cli.main();p=json.loads((out/"provenance.json").read_text());records.append(p["replicates"])
        assert p["depth_sampler"]==selected and p["depth_sampler_requested"]==sampler
    assert calls==[4.] and records[0]==records[1]
    np.testing.assert_array_equal(np.load(tmp_path/"auto/simulations.npy"),np.load(tmp_path/selected/"simulations.npy"))


@pytest.mark.parametrize("option,value",[("--reads","1"),("--fraction","0.5")])
def test_cli_auto_count_modes_skip_probes(tmp_path,monkeypatch,option,value):
    pool=pool_fixture();ref=Reference(("synthetic",),(10,),("",))
    monkeypatch.setattr(cli,"extract_pool",lambda *a,**kw:(pool,ref))
    def unexpected(*a):raise AssertionError("Depth probes must not run in count mode")
    monkeypatch.setattr(cli,"choose_depth_sampler",unexpected)
    alignment=tmp_path/"a";alignment.write_bytes(b"synthetic")
    fai=tmp_path/"r";fai.write_bytes(b"synthetic")
    out=tmp_path/"out"
    monkeypatch.setattr(sys,"argv",["coverage_map.py","--alignment",str(alignment),"--reference-index",str(fai),
        "--outdir",str(out),option,value,"--bins","3","--replicates","2","--depth-sampler","auto"])
    cli.main();p=json.loads((out/"provenance.json").read_text())
    assert p["depth_sampler_requested"]=="auto" and "unused" in p["sampler_calibration"]["reason"]
    assert all(r["units"]==1 for r in p["replicates"])
