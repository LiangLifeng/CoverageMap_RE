import csv
import gzip
import json
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pysam
import pytest

from bin.coverage_core import canonical_reference, read_reference_index, extract_pool
from bin.query_decision_methods import evaluate, COLUMNS
from bin.raw_cram_pool import extract_once
from bin.raw_cram_query import calculate, plan
from bin.raw_query_outputs import read_rows
from test_raw_cram_query import specs, fixture


def test_all_methods_and_thresholds_share_draws(tmp_path):
    planned, _ = specs(tmp_path)
    spec = planned[0]
    spec["config"].update(jaccard_methods="normal,empirical,minimum,bootstrap",
                          bootstrap_replicates=100, bootstrap_confidence=.9,
                          jaccard_coverage_thresholds="default")
    output = tmp_path / "parity"
    output.mkdir()
    calls = []
    def decode(path, layout, **policies):
        calls.append(path)
        return extract_once(path, layout, **policies)
    summary = calculate(spec, output, extractor=decode)
    assert len(calls) == 1 + len(set(spec["choices"]))
    methods = read_rows(output / "decision_methods.tsv")
    assert [r["method"] for r in methods] == ["normal", "empirical", "minimum", "bootstrap"]
    with gzip.open(output / "draws.tsv.gz", "rt") as handle:
        draws = list(csv.DictReader(handle, delimiter="\t"))
    with gzip.open(output / "null.tsv.gz", "rt") as handle:
        nulls = list(csv.DictReader(handle, delimiter="\t"))
    expected = evaluate(summary, draws, nulls, methods=spec["config"]["jaccard_methods"],
                        alpha=.05, bootstrap_replicates=100, bootstrap_confidence=.9, seed=401)
    for actual, oracle in zip(methods, expected):
        for key in COLUMNS:
            assert actual[key] == ("" if oracle.get(key) is None else str(oracle.get(key)))
    levels = read_rows(output / "sensitivity.tsv")
    assert [float(r["jaccard_coverage_threshold"]) for r in levels] == [0, .01, .05, .1, .15, .2, .3]
    assert all(r["pearson_p"] == str(summary["pearson_p"]) for r in levels)
    assert (output / "breadth.tsv.gz").stat().st_size > 0
    assert np.array_equal(np.load(output / "window_lengths.npy"), np.diff(np.load(output / "edges.npy")))


@pytest.mark.parametrize("unit", ["read", "fragment"])
@pytest.mark.parametrize("reverse", [False, True])
@pytest.mark.parametrize("encoding", ["bam", "cram"])
def test_supplementary_matches_original_projection(tmp_path, unit, reverse, encoding):
    fai = tmp_path / "ref.fai"
    fai.write_text("chr1\t100\n")
    bam = tmp_path / "split.bam"
    header = {"SQ": [{"SN": "chr1", "LN": 100}]}
    with pysam.AlignmentFile(str(bam), "wb", header=header) as output:
        for role in (64, 128):
            for supplementary in (False, True):
                read = pysam.AlignedSegment()
                read.query_name = "pair1"
                read.query_sequence = "A" * 10
                read.flag = 1 | role | (16 if reverse else 0) | (2048 if supplementary else 0)
                read.reference_id = 0
                read.reference_start = (10 if role == 64 else 40) + (10 if supplementary else 0)
                read.mapping_quality = 60
                read.cigarstring = "5S5M" if supplementary else "5M5S"
                output.write(read)
    source = bam
    if encoding == "cram":
        source = tmp_path / "split.cram"
        with pysam.AlignmentFile(str(bam), "rb") as original:
            with pysam.AlignmentFile(str(source), "wc", header=original.header,
                                     format_options=[b"no_ref=1"]) as target:
                for read in original.fetch(until_eof=True):
                    target.write(read)
    layout = canonical_reference(read_reference_index(fai))
    pool, _, _ = extract_once(source, layout, sampling_unit=unit,
                              coverage_alignment_policy="primary_plus_supplementary")
    oracle, _ = extract_pool(bam, reference_index=fai, sort_contigs=True,
                             sampling_unit=unit, coverage_alignment_policy="primary_plus_supplementary",
                             capture_geometry=tmp_path / "oracle")
    for name in ("starts", "ends", "indptr"):
        np.testing.assert_array_equal(getattr(pool, name), getattr(oracle, name))
    assert pool.total_bases == 20
    assert int(pool.unit_read_counts.sum()) == 2
    assert int(pool.unit_segment_counts.sum()) == 4


@pytest.mark.parametrize("supp_cigar,primary", [("4S6M", True), ("5S5M", False)])
def test_bad_split_rejected(tmp_path, supp_cigar, primary):
    fai = tmp_path / "ref.fai"
    fai.write_text("chr1\t100\n")
    bam = tmp_path / "bad.bam"
    with pysam.AlignmentFile(str(bam), "wb", header={"SQ": [{"SN": "chr1", "LN": 100}]}) as output:
        for kind in ([0, 2048] if primary else [2048]):
            read = pysam.AlignedSegment()
            read.query_name = "r"
            read.query_sequence = "A" * 10
            read.flag = kind
            read.reference_id = 0
            read.reference_start = 10
            read.mapping_quality = 60
            read.cigarstring = supp_cigar if kind else "5M5S"
            output.write(read)
    with pytest.raises(ValueError, match="split|Overlapping"):
        extract_once(bam, canonical_reference(read_reference_index(fai)),
                     coverage_alignment_policy="primary_plus_supplementary")


def test_coverage_unit_validation(tmp_path):
    planned, _ = specs(tmp_path)
    manifest = tmp_path / "manifest.tsv"
    rows = read_rows(manifest)
    with manifest.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, [*rows[0], "Coverage"], delimiter="\t")
        writer.writeheader()
        writer.writerows(dict(r, Coverage="44.03") for r in rows)
    args = SimpleNamespace(fai_map=None, fai_dir=tmp_path, reference_draws=17, alpha=.05,
        seed=401, depth_threshold=5, bin_count=7, window_size=None, sampling_unit="read",
        min_mapq=0, jaccard_coverage_threshold=0, threads=2, min_samples=0,
        manifest=manifest, index=tmp_path/"index.sqlite", output=tmp_path/"with-unit",
        coverage_unit="percent")
    plan(args)
    import base64
    spec = json.loads(base64.b64decode(read_rows(tmp_path/"with-unit/tasks.tsv")[0]["payload"]))
    output = tmp_path/"coverage"
    output.mkdir()
    assert calculate(spec, output)["manifest_coverage"] == "44.03"
    args.output, args.coverage_unit = tmp_path/"missing-unit", None
    plan(args)
    spec = json.loads(base64.b64decode(read_rows(tmp_path/"missing-unit/tasks.tsv")[0]["payload"]))
    assert "coverage-unit" in spec["error"]


def test_legacy_parameter_inventory():
    root = Path(__file__).resolve().parents[1]
    import re
    original = (root/"nextflow.config").read_text().split("params {", 1)[1].split("}", 1)[0]
    current = (root/"raw_query/nextflow.config").read_text().split("params {", 1)[1].split("}", 1)[0]
    names = lambda text: set(re.findall(r"^\s*([a-z_]+)\s*=", text, re.M))
    assert names(original) - {"cache_mib", "work_mib"} <= names(current)


def test_missing_receipt_uses_legacy_error(tmp_path):
    from bin.raw_cram_query import collect
    roster = tmp_path / "roster.json"
    roster.write_text(json.dumps([dict(query_id="q00000001", Sample="s", GCFID="GCF_1", manifest_depth="1")]))
    collect(SimpleNamespace(roster=roster, results_root=tmp_path/"missing", output=tmp_path/"collected"))
    row, = read_rows(tmp_path/"collected/summary.tsv")
    assert row["status"] == "ERROR"
    assert row["result_source"] == "NO_RECEIPT"
    assert row["manifest_depth"] == "1"
    assert len(read_rows(tmp_path/"collected/errors.tsv")) == 1
    assert json.loads((tmp_path/"collected/completion.json").read_text())["ERROR"] == 1


@pytest.mark.parametrize("state", ["samples", "depth", "zero"])
def test_unscored_is_successful_and_keeps_outputs(tmp_path, state):
    import base64
    import sqlite3
    from bin.raw_cram_query import score
    fai, query, refs, index, manifest = fixture(tmp_path)
    if state == "depth":
        with sqlite3.connect(index) as connection:
            connection.execute("UPDATE reference_records SET Depth='0.25'")
    if state == "zero":
        with pysam.AlignmentFile(str(query), "wb", header={"SQ": [{"SN":"chr2","LN":100},{"SN":"chr10","LN":300}]}):
            pass
    args = SimpleNamespace(fai_map=None, fai_dir=tmp_path, reference_draws=17, alpha=.05,
        seed=401, depth_threshold=5, bin_count=7, window_size=None, sampling_unit="read",
        min_mapq=0, jaccard_coverage_threshold=0, threads=2, min_samples=3 if state=="samples" else 0,
        manifest=manifest, index=index, output=tmp_path/"plan")
    plan(args)
    spec = json.loads(base64.b64decode(read_rows(tmp_path/"plan/tasks.tsv")[0]["payload"]))
    saved = tmp_path/"spec.json"
    saved.write_text(json.dumps(spec))
    assert score(SimpleNamespace(spec=saved, output=tmp_path/"output", stub=False)) == 0
    summary = json.loads((tmp_path/"output/summary.json").read_text())
    assert summary["status"] == {"samples":"INSUFFICIENT_REFERENCE_SAMPLES",
        "depth":"INSUFFICIENT_REFERENCE_DEPTH","zero":"ZERO_QUERY_DEPTH"}[state]
    assert summary["detection_status"] == "NOT_EVALUATED"
    assert summary["requested_draws"] == 17
    assert (tmp_path/"output/decision_methods.tsv").is_file()
    assert (tmp_path/"output/sensitivity.tsv").is_file()

