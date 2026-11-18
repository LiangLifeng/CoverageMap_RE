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


@pytest.mark.parametrize('encoding', ['bam', 'cram'])
@pytest.mark.parametrize('unit', ['read', 'fragment'])
@pytest.mark.parametrize('branch', ['low', 'high'])
def test_vnext_frozen_outputs_dense_and_cache(tmp_path, unit, branch, encoding):
    import gzip
    from bin.coverage_core import canonical_reference, read_reference, dense_coverage
    from bin.raw_cram_pool import extract_once
    from bin.raw_cram_query import rng
    from bin.raw_query_sampling import query_arms, compatible_reference_ids
    from bin.query_manifest_stats import low_null, aggregate, similarity
    spec, _ = _fixture(tmp_path)
    if encoding == 'cram':
        for source in [spec['query']['CramPath'], spec['references'][0]['CramFilePath']]:
            with pysam.AlignmentFile(source, 'rb') as bam:
                reads = list(bam.fetch(until_eof=True))
                header = bam.header
            with pysam.AlignmentFile(source, 'wc', header=header, format_options=[b'no_ref=1']) as cram:
                for read in reads:
                    cram.write(read)
            from bin.raw_cram_query import stamp
            spec['stamps'][source] = stamp(source)
    spec['arm'] = branch
    spec['config'].update(sampling_unit=unit, reference_draws=11,
                          jaccard_methods='normal,empirical,minimum,bootstrap',
                          jaccard_coverage_thresholds='default')
    spec['choices'].append(spec['choices'][0])
    root = tmp_path / 'cold'
    root.mkdir()
    summary = calculate(spec, root)
    layout = canonical_reference(read_reference(spec['fna']))
    query, _, _ = extract_once(spec['query']['CramPath'], layout, sampling_unit=unit,
                               strict_reference=True)
    ref, _, _ = extract_once(spec['references'][0]['CramFilePath'], layout,
                             sampling_unit=unit, strict_reference=True)
    arms = query_arms(query, unit, 11)
    assigned = [a for a in arms for _ in range(a.slots)]
    matrix = np.load(root / 'reference_windows.npy')
    covered = np.load(root / 'reference_covered_bases.npy')
    with gzip.open(root / 'draws.tsv.gz', 'rt') as handle:
        draws = list(csv.DictReader(handle, delimiter='\t'))
    numeric = []
    for i, draw in enumerate(draws):
        if branch == 'low':
            arm = assigned[i]
            generator, _ = rng(401, spec['query_id'], 'query', 'GCF_1.1', 'read-units', i)
            selected = select_nearest_units(ref, compatible_reference_ids(ref, arm.mate_side),
                                             generator, arm.target_bases, 150)
            ids, qids = selected.indices, arm.query_ids
        else:
            ids, qids = np.arange(len(ref)), np.arange(len(query))
        dense = dense_coverage(ref, ids)
        qdense = dense_coverage(query, qids)
        np.testing.assert_array_equal(matrix[i], dense)
        np.testing.assert_array_equal(covered[i], dense > 0)
        values = similarity(qdense, dense)
        assert float(draw['jaccard']) == pytest.approx(values['jaccard'])
        assert float(draw['pearson']) == pytest.approx(values['pearson'])
        numeric.append(dict(draw, jaccard=float(draw['jaccard']), pearson=float(draw['pearson'])))
    expected = aggregate(numeric, low_null(matrix, covered=covered, lengths=np.ones(400, dtype=np.int64)),
                         'low', 11, stable_tail=True)
    actual = json.loads((root / 'statistics.json').read_text())
    for metric in expected:
        if expected[metric]['p'] is None:
            assert actual[metric]['p'] is None
        else:
            assert actual[metric]['p'] == pytest.approx(expected[metric]['p'])
    warm = tmp_path / 'warm'
    warm.mkdir()
    cached = calculate(spec, warm)
    assert cached['decoded_reference_files'] == 0 and cached['geometry_cache_hits'] == 2
    for name in ['reference_windows.npy', 'reference_covered_bases.npy', 'statistics.json',
                 'sensitivity.tsv', 'decision_methods.tsv', 'query_windows.npy']:
        assert (root / name).read_bytes() == (warm / name).read_bytes()
    for row in draws:
        assert row['geometry_cache_hit'] == 'False'
    assert summary['valid_draws'] == 11
    for directory in (root / 'jaccard_thresholds').iterdir():
        with gzip.open(directory / 'draws.tsv.gz', 'rt') as handle:
            threshold_draws = list(csv.DictReader(handle, delimiter='\t'))
        assert [r['mate_side'] for r in threshold_draws] == [r['mate_side'] for r in draws]


def test_out_of_tolerance_never_produces_a_decision(tmp_path):
    spec, _ = _fixture(tmp_path)
    spec['config'].update(sampling_tolerance_bp=0,
                          jaccard_methods='normal,empirical,minimum,bootstrap')
    output = tmp_path / 'result'
    output.mkdir()
    summary = calculate(spec, output)
    assert summary['status'] == 'COMPLETE_WITH_GAPS'
    assert summary['valid_draws'] == 5 and summary['out_of_tolerance_draws'] == 5
    assert summary['detection_status'] == 'NOT_EVALUATED'
    with (output / 'decision_methods.tsv').open() as handle:
        assert all(r['detection_status'] == 'NOT_EVALUATED' for r in csv.DictReader(handle, delimiter='\t'))


def test_fasta_m5_validation_after_replanning(tmp_path):
    from bin.raw_cram_query import stamp
    spec, _ = _fixture(tmp_path)
    source = Path(spec['query']['CramPath'])
    _bam(source, 'f' * 32)
    spec['stamps'][str(source)] = stamp(source)
    output = tmp_path / 'invalid-m5'
    output.mkdir()
    with pytest.raises(ValueError, match='M5 mismatch'):
        calculate(spec, output)


def test_cache_lock_corruption_and_source_change(tmp_path):
    from concurrent.futures import ThreadPoolExecutor
    from bin.coverage_core import canonical_reference, read_reference
    from bin.raw_cram_pool import extract_once
    from bin.raw_geometry_cache import cached_extract
    spec, _ = _fixture(tmp_path)
    source = spec['query']['CramPath']
    layout = canonical_reference(read_reference(spec['fna']))
    cache = tmp_path / 'cache'
    calls = []
    def extract(path, ref, **policies):
        calls.append(path)
        return extract_once(path, ref, **policies)
    def fetch(_):
        return cached_extract(source, layout, extract, cache_dir=cache,
                               sampling_unit='read', strict_reference=True)
    with ThreadPoolExecutor(max_workers=4) as executor:
        results = list(executor.map(fetch, range(4)))
    assert len(calls) == 1
    assert sum(r[3] for r in results) == 3
    assert all(r[0].sampling_unit == 'read' for r in results)
    archive, = cache.glob('*.npz')
    archive.write_bytes(b'broken disposable cache')
    assert not fetch(0)[3] and len(calls) == 2
    assert fetch(0)[3]
    _bam(Path(source), hashlib.md5(b'A' * 400).hexdigest(), short_last=True)
    assert not fetch(0)[3] and len(calls) == 3
    changed_cache = tmp_path / 'changed-cache'
    def changing(path, ref, **policies):
        result = extract_once(path, ref, **policies)
        Path(path).touch()
        return result
    with pytest.raises(ValueError, match='Input changed during'):
        cached_extract(source, layout, changing, cache_dir=changed_cache,
                       sampling_unit='read', strict_reference=True)
    assert not list(changed_cache.glob('*.npz'))


def test_vnext_diagnostics_preserves_mate_coverage(tmp_path):
    from bin.raw_cram_query import table, SUMMARY
    from bin.diagnostics import query_similarity_plot as plot
    spec, _ = _fixture(tmp_path)
    spec['config']['jaccard_methods'] = 'normal,empirical,minimum,bootstrap'
    root = tmp_path / 'queries' / spec['query_id']
    root.mkdir(parents=True)
    summary = calculate(spec, root)
    table(root / 'summary.tsv', SUMMARY, [summary])
    row, = plot.read_table(root / 'summary.tsv', plot.IDENTITY)
    coverage, reasons = plot.coverage_data(row, (root / 'null.tsv.gz', None))
    assert not reasons
    assert coverage['query_mate_mean'] == pytest.approx(46/400)
    assert 'mate mean' in plot.coverage_label(coverage)
    out = tmp_path / 'diagnostics'
    assert plot.run(SimpleNamespace(summary=root/'summary.tsv', decision_methods=root/'decision_methods.tsv',
        methods=spec['config']['jaccard_methods'], jobs=1, query_id=None, outdir=out,
        null_root=root.parent, work_dir=None, staged_nulls=None, trace=None, stub=False)) == 0
    assert len(list(out.rglob('*.png'))) == 4


@pytest.mark.parametrize('unit', ['read', 'fragment'])
def test_strict_subset_header_keeps_full_reference(tmp_path, unit):
    from bin.raw_cram_pool import extract_once
    from bin.coverage_core import Reference, dense_coverage
    digest = hashlib.md5(b'A' * 400).hexdigest()
    source = tmp_path / 'subset.bam'
    _bam(source, digest)
    layout = Reference(('chr1', 'chr2'), (400, 100),
                       (digest, hashlib.md5(b'C' * 100).hexdigest()))
    pool, observed, missing = extract_once(source, layout, sampling_unit=unit, strict_reference=True)
    assert missing == ['chr2']
    assert observed.names == layout.names and observed.length == 500
    dense = dense_coverage(pool, np.arange(len(pool)))
    assert dense[:400].sum() == 200
    assert not dense[400:].any()
    _bam(source, 'f' * 32)
    with pytest.raises(ValueError, match='M5 mismatch'):
        extract_once(source, layout, sampling_unit=unit, strict_reference=True)


def test_failed_score_preserves_small_audit_only(tmp_path):
    from bin.raw_cram_query import score
    spec, _ = _fixture(tmp_path)
    spec['reference_errors'] = {str(spec['references'][0]['id']): 'Synthetic reference unavailable'}
    saved = tmp_path/'spec.json'
    saved.write_text(json.dumps(spec))
    output, audit = tmp_path/'failed', tmp_path/'work'/'failure_audit'
    assert score(SimpleNamespace(spec=saved, output=output, stub=False, failure_audit_dir=audit)) == 1
    for name in ['summary.json', 'draws.tsv', 'decode_audit.json', 'statistics.json']:
        assert (audit/name).read_bytes() == (output/name).read_bytes()
    assert not list(audit.glob('*.npy'))
    assert json.loads((audit/'summary.json').read_text())['status'] == 'ERROR'
