"""Reference CSV-depth floor: source quality is distinct from Query routing."""
import base64
import csv
import json
import shutil
import sqlite3
from pathlib import Path
from types import SimpleNamespace

import pytest

from bin.raw_cram_query import candidates, calculate, plan
from test_raw_cram_query import fixture
from test_raw_query_vnext import _fixture


def _sources(db, query, low, floor):
    with sqlite3.connect(db) as conn:
        conn.row_factory = sqlite3.Row
        audit = {}
        count, refs = candidates(conn, query, low, floor, audit)
    return count, refs, audit


def test_reference_floor_applies_to_low_and_high_and_is_inclusive(tmp_path):
    fai, source, refs, db, manifest = fixture(tmp_path)
    with sqlite3.connect(db) as con:
        for i, depth in enumerate(("4.99", "5", "8"), 1):
            con.execute("UPDATE reference_records SET Depth=? WHERE id=?", (depth, i))
    query = dict(GCFID="GCF_1", Sample="test", CramPath=str(source), Depth="2")
    for low in (True, False):
        available, selected, audit = _sources(db, query, low, 5)
        assert available == 3
        assert {float(r["Depth"]) for r in selected} == {5, 8}
        assert audit["excluded_below_min_reference_depth"] == 1
        assert audit["excluded_below_query_depth"] == 0
        assert audit["eligible_distinct_samples"] == 2
    query["Depth"] = "6"
    _, selected, audit = _sources(db, query, True, 5)
    assert [float(r["Depth"]) for r in selected] == [8]
    assert audit["excluded_below_min_reference_depth"] == 1
    assert audit["excluded_below_query_depth"] == 1
    query["Depth"] = "10"
    _, selected, _ = _sources(db, query, False, 0)
    assert {float(r["Depth"]) for r in selected} == {4.99, 5, 8}


def _opts(root, db, manifest, *, floor=5, policy="unscored", catalog=False):
    return SimpleNamespace(
        manifest=manifest, index=db, output=root / "plan",
        fai_map=None, fai_dir=None if catalog else root,
        reference_draws=10, alpha=.05, seed=401, depth_threshold=5,
        min_reference_depth=floor, bin_count=None if catalog else 7,
        window_size=1 if catalog else None,
        sampling_unit="read", min_mapq=0, min_samples=0, threads=1,
        jaccard_coverage_threshold=0, sampling_tolerance_bp=150,
        insufficient_reference_policy=policy,
        geometry_cache_dir=root / "geometry-cache" if catalog else None,
        coverage_unit=None, coverage_alignment_policy="primary_only",
        jaccard_methods="normal", bootstrap_replicates=100,
        bootstrap_confidence=.95, jaccard_coverage_thresholds=None)


def _planned(path):
    with (path / "tasks.tsv").open() as handle:
        return [json.loads(base64.b64decode(r["payload"]))
                for r in csv.DictReader(handle, delimiter="\t")]


def test_planner_excludes_shallow_sources_from_both_arms_and_records_qc(tmp_path):
    _, _, _, db, manifest = fixture(tmp_path)
    with sqlite3.connect(db) as con:
        for i, depth in enumerate(("4.99", "5", "8"), 1):
            con.execute("UPDATE reference_records SET Depth=? WHERE id=?", (depth, i))
    opts = _opts(tmp_path, db, manifest)
    plan(opts)
    specs = _planned(opts.output)
    assert [spec["arm"] for spec in specs] == ["low", "high"]
    for spec in specs:
        assert spec["config"]["min_reference_depth"] == 5
        assert spec["available_count"] == 3
        assert spec["eligible_count"] == 2
        assert spec["reference_depth_audit"]["excluded_below_min_reference_depth"] == 1
        assert {float(r["Depth"]) for r in spec["references"]} <= {5, 8}
        assert set(spec["choices"]) <= {r["id"] for r in spec["references"]}
    target = tmp_path / "scored"
    target.mkdir()
    result = calculate(specs[1], target)
    assert result["min_reference_depth"] == 5
    assert result["excluded_below_min_reference_depth"] == 1
    assert json.loads((target / "reference_depth_audit.json").read_text())["eligible_after_depth_filter"] == 2
    opts.output = tmp_path / "unrestricted"
    opts.min_reference_depth = 0
    plan(opts)
    for spec in _planned(opts.output):
        assert spec["eligible_count"] == 3
        assert spec["reference_depth_audit"]["excluded_below_min_reference_depth"] == 0


def test_replace_pool_never_contains_below_floor_reference(tmp_path):
    first, _ = _fixture(tmp_path)
    db = tmp_path / "reference.sqlite"
    shallow = tmp_path / "shallow.bam"
    qualified = tmp_path / "qualified.bam"
    shutil.copyfile(first["references"][0]["CramFilePath"], shallow)
    shutil.copyfile(first["references"][0]["CramFilePath"], qualified)
    with sqlite3.connect(db) as con:
        con.execute("INSERT INTO reference_records "
                    "(Sample,CramFilePath,GCFID,Coverage,Depth) VALUES (?,?,?,?,?)",
                    ("shallow", str(shallow), "GCF_1.1", "50", "4.5"))
        con.execute("INSERT INTO reference_records "
                    "(Sample,CramFilePath,GCFID,Coverage,Depth) VALUES (?,?,?,?,?)",
                    ("qualified", str(qualified), "GCF_1.1", "50", "5"))
    opts = _opts(tmp_path, db, tmp_path / "queries.tsv",
                 policy="replace", catalog=True)
    opts.output = tmp_path / "filtered-plan"
    plan(opts)
    spec, = _planned(opts.output)
    assert spec["eligible_count"] == 2
    assert {float(r["Depth"]) for r in spec["references"]} == {5, 6}
    assert shallow.as_posix() not in spec["stamps"]
    assert spec["reference_depth_audit"]["excluded_below_min_reference_depth"] == 1
    assert set(spec["choices"]) <= {r["id"] for r in spec["references"]}


@pytest.mark.parametrize("floor", ["-1", "nan", "inf"])
def test_invalid_reference_floor_is_rejected_before_planning(tmp_path, floor):
    _, _, _, db, manifest = fixture(tmp_path)
    opts = _opts(tmp_path, db, manifest, floor=floor)
    with pytest.raises(ValueError, match="min_reference_depth"):
        plan(opts)
