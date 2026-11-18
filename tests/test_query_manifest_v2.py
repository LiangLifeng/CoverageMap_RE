"""Exercise both v2 arms through the CLI and reconstruct their distributions."""

import csv
import gzip
import json
import math
from pathlib import Path
import subprocess
import sys
from threading import Barrier

import numpy as np
import pytest

from bin.coverage_core import Reference
from bin.gcf_windows import GCF_CONCAT
from bin.query_manifest import fill_distribution_draws, table
from bin.query_manifest_stats import low_null, normal_test, reference_row_means
from bin.ref_depth_catalog import collect, create_plan, digest
from bin.ref_depth_worker import run_gcf
from bin.reference_index import build_index
from bin.stage_timing import StageTiming


def bam(pysam, root, sample, counts):
    path = root / (sample + '.bam')
    header = {'HD': {'VN': '1.6'}, 'SQ': [{'SN': 'one', 'LN': 80}]}
    with pysam.AlignmentFile(path, 'wb', header=header) as handle:
        for bin_id, count in enumerate(counts):
            for number in range(count):
                read = pysam.AlignedSegment()
                read.query_name = f'{sample}-{bin_id}-{number}'
                read.query_sequence = 'AAAA'
                read.flag = 0
                read.reference_id = 0
                read.reference_start = 10 * bin_id
                read.mapping_quality = 60
                read.cigar = [(0, 4)]
                handle.write(read)
    return path


def expected_distributions(query, references):
    all_maps = np.vstack((query, references))
    support = all_maps > 0
    intersections = np.logical_and(support[:, None], support[None, :]).sum(axis=2)
    unions = np.logical_or(support[:, None], support[None, :]).sum(axis=2)
    for metric, similarities in (('jaccard', intersections / unions), ('pearson', np.corrcoef(all_maps))):
        reference_pairs = similarities[1:, 1:]
        n = len(references)
        off_diagonal = ~np.eye(n, dtype=bool)
        yield metric, similarities[0, 1:].mean(), reference_pairs[off_diagonal].reshape(n, n-1).mean(axis=1)


def test_normal_left_tail_preserves_tiny_p_value_and_reports_zero_variance():
    # Reference mean .8 and SD .01 gives a query of .6 a z-score of -20.
    result = normal_test(.6, [.79, .81] * 50, stable_tail=True)
    assert result['p'] == pytest.approx(2.7536241186062337e-89, rel=1e-12, abs=0)
    result = normal_test(.6, [.8] * 100, stable_tail=True)
    assert result['p'] is None and result['status'] == 'zero_variance_null'


@pytest.mark.parametrize('windows', [1, 60, 131_073])
def test_threaded_reference_means_match_scalar_pairs_exactly(windows):
    refs = np.random.default_rng(71).integers(0, 3, (12, windows))
    refs[1] = refs[0]
    for constant_reference in (False, True):
        if constant_reference:
            refs[2] = 1
        expected = low_null(refs)
        assert reference_row_means(refs, 1) == reference_row_means(refs, 4) == expected


def test_preparation_and_simulations_both_use_multiple_threads(monkeypatch):
    prepared = Barrier(2, timeout=5)
    simulated = Barrier(2, timeout=5)
    layout = Reference(('one',), (3,), ('',))
    records = dict(policies={}, units=1)

    class PoolToken(int):
        pass

    class Cache:
        def stamps(self, row):
            return [row['id']]

        def load(self, row, timing):
            prepared.wait()
            return PoolToken(row['id']), layout, records, None

    def sample(ref, generator, target, edges, audit=None):
        simulated.wait()
        return np.array([ref, 0, generator]), 1, str(generator)

    monkeypatch.setattr('bin.query_manifest.sample_reference_windows', sample)
    groups = {}
    for ref_id in (1, 2):
        groups[ref_id] = [
            (dict(id=ref_id), dict(draw_id=slot), slot+1)
            for slot in range(ref_id-1, 4, 2)
        ]
    matrix = np.zeros((4, 3), dtype=np.int64)
    verified = fill_distribution_draws(
        groups, matrix, np.array([1, 0, 2]), layout, np.arange(4), Cache(), StageTiming(),
        arm='low', target_bases=1, policies={}, threads=2)
    np.testing.assert_array_equal(matrix, [[1, 0, 1], [2, 0, 2], [1, 0, 3], [2, 0, 4]])
    assert sorted(verified) == [1, 2]
    assert all(entry['status'] == 'READY' for group in groups.values() for _, entry, _ in group)


@pytest.mark.parametrize('draw_count', [None, 2, 17, 137])
def test_v2_cli_uses_the_same_draws_for_query_and_null(tmp_path, draw_count):
    requested = 100 if draw_count is None else draw_count
    pysam = pytest.importorskip('pysam')
    gcfs = ('GCF_000000001.1', 'GCF_000000002.1')
    fai = tmp_path / 'ref.fai'
    fai.write_text('one\t80\n')
    refs = tmp_path / 'refs.csv'
    refs.write_text('GCF_ID,GCF_path\n' + ''.join(f'{gcf},{fai}\n' for gcf in gcfs))
    patterns = {
        'a': [12, 6, 0, 0, 2],
        'b': [10, 8, 2, 0, 0],
        'c': [9, 6, 1, 0, 2],
        'd': [14, 8, 2, 0, 0, 0, 1],
        'shallow': [1, 1],
    }
    inputs = {sample: bam(pysam, tmp_path, sample, counts) for sample, counts in patterns.items()}
    records = tmp_path / 'records.csv'
    records.write_text('Sample,CramFilePath,GCFID\n' + ''.join(
        f'{sample},{path},{gcf}\n' for gcf in gcfs for sample, path in inputs.items()
        if gcf == gcfs[0] or sample != 'shallow'))
    config = dict(storage_mode='native', window_size=10, window_policy=GCF_CONCAT,
                  cpus=1, python=sys.executable, memory=None,
                  geometry_format='packed', policies=dict(mate_policy='both', sampling_unit='read',
                                                         min_mapq=0, exclude_duplicates=False, exclude_qcfail=False))
    columns = dict(sample='Sample', cram='CramFilePath', gcf='GCFID', ref_gcf='GCF_ID', fai='GCF_path')
    database = tmp_path / 'db'
    create_plan(records, refs, database, config, columns)
    for gcf in gcfs:
        assert run_gcf(database, gcf) == 0
    assert collect(database)['complete']
    index = tmp_path / 'index.sqlite'
    build_index(database, index, True)
    queries = (
        ('low', gcfs[0], 5, bam(pysam, tmp_path, 'low', [5, 3, 1, 0, 1])),
        ('alias', gcfs[1], 6, inputs['a']),
        ('wrong_shape', gcfs[1], 6, bam(pysam, tmp_path, 'wrong_shape', [1, 1, 18, 0, 20])),
        ('constant', gcfs[1], 6, bam(pysam, tmp_path, 'constant', [6] * 8)),
        ('zero', gcfs[0], 0, bam(pysam, tmp_path, 'zero', [])),
    )
    manifest = tmp_path / 'query.tsv'
    fields = ('Sample', 'GCFID', 'Depth', 'CramPath')
    table(manifest, fields, [dict(zip(fields, query)) for query in queries])
    before = {path: digest(path) for path in tmp_path.rglob('*') if path.is_file()}
    entry = Path(__file__).resolve().parents[1] / 'bin' / 'query_manifest_v2.py'
    out = tmp_path / 'out'
    command = [sys.executable, '-B', str(entry), '--manifest', str(manifest), '--index', str(index),
               '--outdir', str(out), '--window-size', '10', '--min-samples', '2', '--jobs', '2']
    missing_alpha = subprocess.run(command, capture_output=True, text=True)
    assert missing_alpha.returncode == 2 and '--alpha' in missing_alpha.stderr and not out.exists()
    command += ['--alpha', '0.05']
    bad_threads = subprocess.run(command + ['--threads', '1'], capture_output=True, text=True)
    assert bad_threads.returncode == 2 and '--threads must be at least 2' in bad_threads.stderr
    if draw_count is not None:
        command += ['--reference-draws', str(draw_count)]
    command += ['--threads', '4']
    result = subprocess.run(command, capture_output=True, text=True)
    assert result.returncode == 0, result.stderr
    with (out / 'summary.tsv').open() as handle:
        rows = list(csv.DictReader(handle, delimiter='\t'))
    assert [row['Sample'] for row in rows] == [query[0] for query in queries]
    assert [row['arm'] for row in rows[:3]] == ['low', 'high', 'high']
    assert rows[0]['eligible_count'] == '4'  # The shallow reference cannot match actual query mass.
    assert rows[1]['excluded_self_count'] == '1'  # Identical file under a different Sample.
    if requested == 100:
        assert rows[0]['detection_status'] == 'DETECTED'
        assert rows[1]['detection_status'] == 'NOT_DETECTED'
        assert (rows[2]['jaccard_decision'], rows[2]['pearson_decision']) == ('PASS', 'FAIL')
        assert rows[2]['detection_status'] == 'NOT_DETECTED'
        assert rows[3]['pearson_decision'] == 'NOT_EVALUATED'
        assert rows[3]['detection_status'] == rows[4]['detection_status'] == 'NOT_EVALUATED'
        assert rows[4]['status'] == 'ZERO_QUERY_DEPTH'

    for row in rows[:3]:
        path = out / 'queries' / row['query_id']
        with gzip.open(path / 'draws.tsv.gz', 'rt') as handle:
            draws = list(csv.DictReader(handle, delimiter='\t'))
        with gzip.open(path / 'null.tsv.gz', 'rt') as handle:
            null = list(csv.DictReader(handle, delimiter='\t'))
        query_map = np.load(path / 'query_windows.npy')
        ref_maps = np.load(path / 'reference_windows.npy')
        assert row['requested_draws'] == row['valid_draws'] == str(requested)
        assert ref_maps.shape == (requested, 8) and len(draws) == requested and len(null) == 2*requested
        assert all(d['same_sample'] == d['same_source'] == 'False' for d in draws)
        if row['arm'] == 'low':
            assert (ref_maps.sum(axis=1) == query_map.sum()).all()
        for metric, red, blue in expected_distributions(query_map, ref_maps):
            assert float(row[metric]) == pytest.approx(red, abs=1e-12)
            metric_rows = [item for item in null if item['metric'] == metric]
            assert [item['left_reference_id'] for item in metric_rows] == [d['reference_id'] for d in draws]
            np.testing.assert_allclose([float(item['value']) for item in metric_rows], blue, rtol=0, atol=1e-12)
            # Two draws have the same single unordered reference pair exactly;
            # corrcoef may differ by one rounding bit across the mirrored cells.
            if requested == 2 or min(blue) == max(blue):
                assert row[metric + '_p'] == ''
            else:
                p = 0.5 * math.erfc((blue.mean() - red) / (math.sqrt(2) * blue.std()))
                assert float(row[metric + '_p']) == pytest.approx(p, rel=1e-8)

    assert json.loads((out / 'run.json').read_text())['reference_pool_size'] == requested
    assert subprocess.run(command + ['--resume', '--reference-draws', str(requested+1)], capture_output=True).returncode == 2
    completion = json.loads((out / 'completion.json').read_text())
    assert completion['threads'] == 4 and completion['requested_cpu_slots'] == 8
    assert subprocess.run(command + ['--resume'], capture_output=True).returncode == 0
    assert subprocess.run(command + ['--resume', '--alpha', '0.01'], capture_output=True).returncode == 2
    two_threads = tmp_path / 'two_threads'
    result = subprocess.run(command + ['--outdir', str(two_threads), '--threads', '2', '--jobs', '1',
                                     '--reference-draws', str(requested)],
                            capture_output=True, text=True)
    assert result.returncode == 0, result.stderr
    assert (out / 'summary.tsv').read_bytes() == (two_threads / 'summary.tsv').read_bytes()
    for row in rows[:4]:
        for name in ('draws.tsv.gz', 'null.tsv.gz', 'query_windows.npy', 'reference_windows.npy'):
            path = Path('queries') / row['query_id'] / name
            assert (out / path).read_bytes() == (two_threads / path).read_bytes()

    single = tmp_path / 'single'
    single_command = command + ['--outdir', str(single), '--jobs', '1']
    missing = subprocess.run(single_command + ['--query-id', 'q99999999'], capture_output=True, text=True)
    assert missing.returncode == 2 and 'Query ID not found' in missing.stderr and not single.exists()
    result = subprocess.run(single_command + ['--query-id', rows[2]['query_id']], capture_output=True, text=True)
    assert result.returncode == 0, result.stderr
    with (single / 'summary.tsv').open() as handle:
        assert list(csv.DictReader(handle, delimiter='\t')) == [rows[2]]
    for name in ('draws.tsv.gz', 'null.tsv.gz', 'query_windows.npy', 'reference_windows.npy'):
        path = Path('queries') / rows[2]['query_id'] / name
        assert (out / path).read_bytes() == (single / path).read_bytes()
    assert before == {path: digest(path) for path in before}


@pytest.mark.parametrize('value', ['0', '1', '-1', '2.5', '10%'])
def test_v2_rejects_invalid_reference_draws_before_reading_inputs(tmp_path, value):
    from bin.query_manifest import parser, run
    command = ['--manifest', str(tmp_path/'missing.tsv'), '--index', str(tmp_path/'missing.sqlite'),
               '--outdir', str(tmp_path/'out'), '--alpha', '.05', '--bin-count', '60',
               '--reference-draws', value]
    with pytest.raises((ValueError, SystemExit)):
        run(parser(reference_distribution=True).parse_args(command))
    assert not (tmp_path/'out').exists()
