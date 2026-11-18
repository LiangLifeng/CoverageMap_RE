"""Focused vNext regressions: FASTA catalog, mate splits, tolerance and cache."""
import csv
import hashlib
import json
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pysam
import pytest

from bin.build_cram_index import build_index
from bin.raw_cram_query import plan, calculate
from bin.raw_query_sampling import select_nearest_units


def _bam(path, m5, *, short_last=False):
    sq = [{"SN": "chr1", "LN": 400, "M5": m5}]
    with pysam.AlignmentFile(str(path), "wb", header={"HD": {"VN": "1.6"}, "SQ": sq}) as writer:
        for i in range(10):
            for mate in (1, 2):
                read = pysam.AlignedSegment()
                read.query_name = f"pair{i}"
                read.flag = 1 | (64 if mate == 1 else 128)
                read.reference_id = 0
                read.reference_start = 10 + 4 * i + 5 * (mate - 1)
                read.mapping_quality = 60
                size = 8 if short_last and i == 9 and mate == 2 else 10
                read.cigarstring = f"{size}M"
                read.query_sequence = "A" * size
                writer.write(read)


def _fixture(tmp_path):
    gcf = "GCF_1.1"
    fna = tmp_path / (gcf + ".fna")
    fai = tmp_path / (gcf + ".fai")
    fna.write_text(">chr1\n" + ("A" * 400) + "\n")
    fai.write_text("chr1\t400\t6\t400\t401\n")
    m5 = hashlib.md5(("A" * 400).encode()).hexdigest()
    qpath, rpath = tmp_path / "query.bam", tmp_path / "reference.bam"
    _bam(qpath, m5)
    _bam(rpath, m5, short_last=True)
    build = tmp_path / "build.csv"
    with build.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["Sample", "CramFilePath", "GCFID", "Coverage", "Depth"])
        writer.writerow(["reference", str(rpath), gcf, "50", "6"])
    mapping = []
    for name, column, path in (("fna.csv", "GCF_path", fna), ("fai.csv", "FAI_path", fai)):
        filename = tmp_path / name
        with filename.open("w", newline="") as handle:
            writer = csv.writer(handle)
            writer.writerow(["GCF_ID", column])
            writer.writerow([gcf, str(path)])
        mapping.append(filename)
    sql = tmp_path / "reference.sqlite"
    build_index(build, sql, *mapping)
    manifest = tmp_path / "queries.tsv"
    with manifest.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["Sample", "GCFID", "Depth", "CramPath"])
        writer.writerow(["query", gcf, "0.5", str(qpath)])
    opts = SimpleNamespace(manifest=manifest, index=sql, output=tmp_path / "plan",
        fai_map=None, fai_dir=None, bin_count=None, window_size=1, alpha=.05,
        reference_draws=10, seed=401, depth_threshold=5, min_samples=0,
        sampling_unit="read", sampling_tolerance_bp=150,
        geometry_cache_dir=tmp_path / "shared_geometry", min_mapq=0, threads=1,
        jaccard_coverage_threshold=0, coverage_unit=None,
        coverage_alignment_policy="primary_only", jaccard_methods="normal,minimum",
        bootstrap_replicates=100, bootstrap_confidence=.95,
        jaccard_coverage_thresholds=None)
    plan(opts)
    import base64
    with (opts.output / "tasks.tsv").open() as handle:
        row, = list(csv.DictReader(handle, delimiter="\t"))
    return json.loads(base64.b64decode(row["payload"])), m5


def test_nearest_prefix_allows_small_undershoot():
    class Pool:
        unit_bases = np.array([100, 98])
        def __len__(self):
            return 2
    result = select_nearest_units(Pool(), [0, 1], np.random.default_rng(7), 200, 150)
    assert result.actual_bases == 198
    assert result.delta_bases == -2
    assert result.within_tolerance and result.sampling_status == "UNDERSHOOT"


def test_sql_fasta_mate_split_and_cross_query_cache(tmp_path):
    spec, _ = _fixture(tmp_path)
    assert spec["fna"].endswith("GCF_1.1.fna")
    assert spec["fai"].endswith("GCF_1.1.fai")
    first = tmp_path / "score_1"
    first.mkdir()
    result = calculate(spec, first)
    assert result["status"] in ("COMPLETE", "COMPLETE_WITH_GAPS")
    assert result["valid_draws"] == 10
    assert result["query_read1_bases"] == 100
    assert result["query_read2_bases"] == 100
    assert result["tolerated_undershoot_draws"] == 5
    with (first / "draws.tsv").open() as handle:
        draws = list(csv.DictReader(handle, delimiter="\t"))
    assert [d["mate_side"] for d in draws] == ["read1"] * 5 + ["read2"] * 5
    assert all(d["status"] == "READY" for d in draws)
    assert all(d["sampling_delta_bases"] == "-2" for d in draws[5:])
    second = tmp_path / "score_2"
    second.mkdir()
    cached = calculate(spec, second)
    assert cached["valid_draws"] == 10
    assert cached["decoded_reference_files"] == 0
    assert cached["geometry_cache_hits"] == 2


def test_fasta_m5_disagreement_is_strict(tmp_path):
    spec, _ = _fixture(tmp_path)
    source = Path(spec["query"]["CramPath"])
    _bam(source, "f" * 32)
    # Planner's file stamp intentionally detects the changed source.
    output = tmp_path / "invalid"
    output.mkdir()
    with pytest.raises(ValueError, match="Input changed"):
        calculate(spec, output)
