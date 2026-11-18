import base64
from collections import Counter
import csv
import json
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pysam
import pytest

from bin.build_cram_index import build_index
from bin.coverage_core import dense_coverage, canonical_reference, read_reference_index, select_units
from bin.raw_cram_pool import extract_once
from bin.raw_cram_query import plan, calculate, score, rng
from bin.query_manifest_stats import low_null, aggregate


def alignment(path, seed=1, pairs=40, ambiguous=False):
    generator = np.random.default_rng(seed)
    header = {"HD": {"VN": "1.6"}, "SQ": [{"SN": "chr2", "LN": 100}, {"SN": "chr10", "LN": 300}]}
    with pysam.AlignmentFile(str(path), "wb", header=header) as output:
        for i in range(pairs):
            start = int(generator.integers(0, 270))
            for mate in (1, 2):
                read = pysam.AlignedSegment()
                read.query_name = f"pair{i}"
                read.flag = 1 | (64 if mate == 1 else 128)
                read.reference_id = 1
                read.reference_start = start + (mate - 1) * 2
                read.mapping_quality = 60
                read.cigarstring = "5M3N5M"
                read.query_sequence = "A" * 10
                output.write(read)
                if ambiguous and i == 0 and mate == 1:
                    output.write(read)


def fixture(root):
    root = Path(root)
    root.mkdir(parents=True, exist_ok=True)
    fai = root / "GCF_1.fai"
    fai.write_text("chr2\t100\nchr10\t300\n")
    query = root / "query.bam"
    alignment(query, 93, 12)
    refs = [root / f"ref{i}.bam" for i in range(3)]
    for i, ref in enumerate(refs):
        alignment(ref, 10 + i * 7)
    records = root / "records.csv"
    with records.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["Sample", "CramFilePath", "GCFID", "Coverage", "Depth"])
        for i, ref in enumerate(refs):
            writer.writerow([f"ref{i}", str(ref), "GCF_1", "90", "20"])
    index = root / "index.sqlite"
    build_index(records, index)
    manifest = root / "manifest.tsv"
    with manifest.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["Sample", "GCFID", "Depth", "CramPath"])
        writer.writerow(["test", "GCF_1", "1", str(query)])
        writer.writerow(["test2", "GCF_1", "6", str(query)])
    return fai, query, refs, index, manifest


def specs(root, *, unit="read", missing=False):
    fai, query, refs, index, manifest = fixture(root)
    if missing:
        refs[0].unlink()
    args = SimpleNamespace(fai_map=None, fai_dir=root, reference_draws=17, alpha=.05, seed=401,
                           depth_threshold=5, bin_count=7, window_size=None, sampling_unit=unit,
                           min_mapq=0, jaccard_coverage_threshold=0., threads=2, min_samples=0,
                           manifest=manifest, index=index, output=root / "plan")
    plan(args)
    with (root / "plan/tasks.tsv").open() as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    return [json.loads(base64.b64decode(row["payload"])) for row in rows], refs


@pytest.mark.parametrize("unit", ["read", "fragment"])
def test_dense_sampling_and_single_decode(tmp_path, unit):
    planned, refs = specs(tmp_path, unit=unit)
    spec = planned[0]
    calls = Counter()
    def decode(path, layout, **policies):
        calls[str(path)] += 1
        return extract_once(path, layout, **policies)
    output = tmp_path / "result"
    output.mkdir()
    summary = calculate(spec, output, extractor=decode)
    assert summary["status"] == "COMPLETE"
    assert all(value == 1 for value in calls.values())
    assert summary["decoded_reference_files"] == len(set(spec["choices"]))
    assert summary["duplicate_draws"] > 0
    matrices = np.load(output / "reference_windows.npy")
    cov = np.load(output / "reference_covered_bases.npy")
    edges = np.load(output / "edges.npy")
    layout = canonical_reference(read_reference_index(tmp_path / "GCF_1.fai"))
    seen = {}
    for ref in spec["references"]:
        pool, _, _ = extract_once(ref["CramFilePath"], layout, sampling_unit=unit)
        if unit == "fragment":
            assert (pool.unit_read_counts == 2).all()
        for i, refid in enumerate(spec["choices"]):
            if refid != ref["id"]:
                continue
            generator, _ = rng(401, spec["query_id"], "test", "GCF_1", "read-units", i)
            ids = select_units(pool, generator, target_bases=240)
            dense = dense_coverage(pool, ids)
            expected = np.array([dense[a:b].sum() for a, b in zip(edges[:-1], edges[1:])])
            coverage = np.array([(dense[a:b] > 0).sum() for a, b in zip(edges[:-1], edges[1:])])
            np.testing.assert_array_equal(matrices[i], expected)
            np.testing.assert_array_equal(cov[i], coverage)
            seen.setdefault(refid, []).append(tuple(ids))
    assert any(len(set(values)) > 1 for values in seen.values())
    with (output / "draws.tsv").open() as handle:
        draws = list(csv.DictReader(handle, delimiter="\t"))
    # Independent scalar pair-loop oracle agrees with the reused v2 null.
    nulls = low_null(matrices, covered=cov, lengths=np.diff(edges))
    saved = json.loads((output / "statistics.json").read_text())
    numeric = [dict(d, jaccard=float(d["jaccard"]), pearson=float(d["pearson"])) for d in draws]
    oracle = aggregate(numeric, nulls, "low", 17, stable_tail=True)
    for metric in ("jaccard", "pearson"):
        assert saved[metric]["p"] == pytest.approx(oracle[metric]["p"], abs=1e-12)


def test_high_full_coverage_once(tmp_path):
    planned, _ = specs(tmp_path)
    output = tmp_path / "high"
    output.mkdir()
    summary = calculate(planned[1], output)
    assert summary["status"] in ("COMPLETE", "COMPLETE_WITH_GAPS")
    assert summary["valid_draws"] == 17
    audit = json.loads((output / "decode_audit.json").read_text())
    assert all(r["decode_attempts"] == r["full_coverage_calculations"] == r["query_similarity_calculations"] == 1 for r in audit)
    matrix = np.load(output / "reference_windows.npy")
    for row in audit:
        for slot in row["draw_slots"]:
            np.testing.assert_array_equal(matrix[slot], matrix[row["draw_slots"][0]])


def test_missing_reference_retains_all_slots(tmp_path):
    planned, _ = specs(tmp_path, missing=True)
    output = tmp_path / "failed"
    output.mkdir()
    summary = calculate(planned[0], output)
    assert summary["status"] == "ERROR"
    assert 0 < summary["valid_draws"] < 17
    assert summary["detection_status"] == "NOT_EVALUATED"
    with (output / "draws.tsv").open() as handle:
        draws = list(csv.DictReader(handle, delimiter="\t"))
    assert len(draws) == 17
    assert {d["status"] for d in draws} == {"ERROR", "READY"}
    # A failed score exits nonzero, so Nextflow cannot reuse it as success.
    spec = tmp_path / "spec.json"
    spec.write_text(json.dumps(planned[0]))
    assert score(SimpleNamespace(spec=spec, output=tmp_path / "failed-score", stub=False)) == 1


@pytest.mark.parametrize("unit", ["read", "fragment"])
def test_identity_collision_rejected(tmp_path, unit):
    fai, *_ = fixture(tmp_path)
    path = tmp_path / "ambiguous.bam"
    alignment(path, ambiguous=True)
    with pytest.raises(ValueError, match="Ambiguous"):
        extract_once(path, canonical_reference(read_reference_index(fai)), sampling_unit=unit)


def test_fai_full_denominator_and_gap(tmp_path):
    fai, query, *_ = fixture(tmp_path)
    pool, layout, _ = extract_once(query, canonical_reference(read_reference_index(fai)), sampling_unit="fragment")
    assert layout.length == 400
    assert len(pool) == 12
    assert pool.total_bases == 240
    dense = dense_coverage(pool, np.arange(len(pool)))
    assert dense[:100].sum() == 0
    assert dense.sum() == 240


if __name__ == "__main__":
    import sys
    fixture(Path(sys.argv[1]).resolve())

