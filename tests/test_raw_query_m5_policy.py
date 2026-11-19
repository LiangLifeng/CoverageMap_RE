"""Strict/geometry-only M5 policies and persisted Query/Reference header audit."""
import json
from pathlib import Path

import pysam
import pytest

from bin.coverage_core import Reference
from bin.raw_cram_pool import M5ValidationError, extract_once
from bin.raw_cram_query import calculate, stamp
from test_raw_query_vnext import _bam, _fixture


@pytest.mark.parametrize("unit", ["read", "fragment"])
def test_geometry_only_keeps_fasta_layout_and_records_mismatch(tmp_path, unit):
    spec, digest = _fixture(tmp_path)
    source = Path(spec["query"]["CramPath"])
    _bam(source, "f" * 32)
    spec["stamps"][str(source)] = stamp(source)
    spec["config"]["sampling_unit"] = unit
    strict_dir = tmp_path / "strict"
    strict_dir.mkdir()
    with pytest.raises(M5ValidationError, match="M5 mismatch"):
        calculate(spec, strict_dir)
    strict_audit = json.loads((strict_dir / "m5_audit.json").read_text())
    assert strict_audit[0]["source_role"] == "query"
    assert strict_audit[0]["status"] == "ERROR"
    assert strict_audit[0]["contigs"][0]["status"] == "MISMATCH"
    spec["config"]["m5_policy"] = "geometry_only"
    loose_dir = tmp_path / "geometry"
    loose_dir.mkdir()
    summary = calculate(spec, loose_dir)
    assert summary["valid_draws"] == 10
    assert summary["status"] in ("COMPLETE", "COMPLETE_WITH_GAPS")
    audit = json.loads((loose_dir / "m5_audit.json").read_text())
    assert audit[0]["source_role"] == "query"
    assert audit[0]["status"] == "MISMATCH_ALLOWED"
    assert audit[0]["mismatch_contigs"] == 1
    assert audit[0]["contigs"][0]["fna_m5"] == digest
    assert audit[0]["contigs"][0]["cram_m5"] == "f" * 32
    assert audit[1]["source_role"] == "reference"
    assert audit[1]["status"] == "PASS"
    warm_dir = tmp_path / "warm"
    warm_dir.mkdir()
    warm = calculate(spec, warm_dir)
    assert warm["decoded_reference_files"] == 0
    warm_audit = json.loads((warm_dir / "m5_audit.json").read_text())
    assert [x["contigs"] for x in warm_audit] == [x["contigs"] for x in audit]


def test_reference_mismatch_rejected_in_strict_accepted_in_geometry_only(tmp_path):
    spec, _ = _fixture(tmp_path)
    reference = Path(spec["references"][0]["CramFilePath"])
    _bam(reference, "0" * 32, short_last=True)
    spec["stamps"][str(reference)] = stamp(reference)
    spec["config"]["insufficient_reference_policy"] = "replace"
    strict_dir = tmp_path / "strictref"
    strict_dir.mkdir()
    with pytest.raises(M5ValidationError, match="M5 mismatch"):
        calculate(spec, strict_dir)
    strict_audit = json.loads((strict_dir / "m5_audit.json").read_text())
    assert strict_audit[-1]["status"] == "ERROR"
    spec["config"]["m5_policy"] = "geometry_only"
    loose_dir = tmp_path / "loose_ref"
    loose_dir.mkdir()
    result = calculate(spec, loose_dir)
    assert result["valid_draws"] == 10
    records = json.loads((loose_dir / "m5_audit.json").read_text())
    assert records[-1]["source_role"] == "reference"
    assert records[-1]["mismatch_contigs"] == 1
    assert records[-1]["status"] == "MISMATCH_ALLOWED"


def test_missing_m5_audited_but_missing_contig_not_mismatch(tmp_path):
    spec, digest = _fixture(tmp_path)
    source = Path(spec["query"]["CramPath"])
    with pysam.AlignmentFile(str(source), "rb") as bam:
        reads = list(bam.fetch(until_eof=True))
        header = bam.header.to_dict()
    header["SQ"][0].pop("M5")
    with pysam.AlignmentFile(str(source), "wb", header=header) as out:
        for read in reads:
            out.write(read)
    layout = Reference(("chr1", "chr2"), (400, 50), (digest, "a" * 32))
    with pytest.raises(M5ValidationError, match="M5 mismatch"):
        extract_once(source, layout, strict_reference=True)
    pool, observed, missing = extract_once(
        source, layout, strict_reference=True, m5_policy="geometry_only")
    assert missing == ["chr2"]
    assert observed.md5 == layout.md5
    assert [r["status"] for r in pool.m5_audit] == ["MISSING_M5", "MISSING_CONTIG"]


def test_geometry_only_never_disables_sn_ln_validation(tmp_path):
    spec, digest = _fixture(tmp_path)
    source = spec["query"]["CramPath"]
    wrong_length = Reference(("chr1",), (401,), (digest,))
    with pytest.raises(ValueError, match="length differs"):
        extract_once(source, wrong_length, strict_reference=True, m5_policy="geometry_only")
    wrong_name = Reference(("other",), (400,), (digest,))
    with pytest.raises(ValueError, match="duplicate/extra contigs"):
        extract_once(source, wrong_name, strict_reference=True, m5_policy="geometry_only")
    with pytest.raises(ValueError, match="requires authoritative SQL FASTA"):
        extract_once(source, wrong_length, strict_reference=False, m5_policy="geometry_only")
    with pytest.raises(ValueError, match="strict or geometry_only"):
        extract_once(source, wrong_length, strict_reference=True, m5_policy="skip")


def test_policy_defaults_strict_and_native_interface_wired(tmp_path):
    spec, _ = _fixture(tmp_path)
    assert spec["config"]["m5_policy"] == "strict"
    root = Path(__file__).resolve().parents[1]
    assert "m5_policy = 'strict'" in (root / "raw_query/nextflow.config").read_text()
    plan = (root / "modules/local/raw_query_plan/main.nf").read_text()
    assert plan.count("--m5-policy") == 2
