"""Prove legacy geometry produces exactly the same GCF-concat query results."""
import csv
import gzip
from pathlib import Path
import subprocess
import sys

import numpy as np
import pytest

from bin.coverage_core import BlockPool, Reference, canonical_reference
from bin.gcf_windows import GCF_CONCAT, PER_CONTIG
from bin.query_manifest import align_reference_geometry, table
from bin.ref_depth_catalog import collect, create_plan, digest
from bin.ref_depth_worker import run_gcf
from bin.reference_index import build_index


def test_rebase_preserves_units_and_rejects_cross_contig_blocks():
    original = Reference(('z', 'a'), (40, 40), ('', ''))
    target = canonical_reference(original)
    pool = BlockPool(np.array([0, 38, 40, 72]), np.array([3, 40, 44, 80]),
                     np.array([0, 2, 4]), 80)
    remapped, layout = align_reference_geometry(pool, original, target)
    np.testing.assert_array_equal(remapped.starts, [40, 78, 0, 32])
    np.testing.assert_array_equal(remapped.ends, [43, 80, 4, 40])
    np.testing.assert_array_equal(remapped.indptr, pool.indptr)
    np.testing.assert_array_equal(remapped.unit_bases, pool.unit_bases)
    assert layout == target and remapped.total_bases == pool.total_bases
    invalid = BlockPool.from_intervals(np.array([[39, 42]]), 80)
    with pytest.raises(ValueError, match='crosses a contig'):
        align_reference_geometry(invalid, original, target)


def bam(root, name, multiplier):
    pysam = pytest.importorskip('pysam')
    path = root / (name + '.bam')
    with pysam.AlignmentFile(path, 'wb', header={
            'SQ': [{'SN': 'z', 'LN': 19}, {'SN': 'a', 'LN': 61}]}) as handle:
        for contig, positions in ((0, (0, 4, 11)), (1, (0, 10, 25, 49))):
            for pos in positions:
                for repeat in range(multiplier + (pos % 3)):
                    read = pysam.AlignedSegment()
                    read.query_name = f'{name}-{contig}-{pos}-{repeat}'
                    read.flag = 0
                    read.reference_id = contig
                    read.reference_start = pos
                    read.mapping_quality = 60
                    read.cigar = [(0, 4)]
                    handle.write(read)
    return path


def test_legacy_and_concat_indexes_have_identical_low_high_scores(tmp_path):
    gcf = 'GCF_000000001.1'
    fai = tmp_path / 'ref.fai'
    fai.write_text('z\t19\na\t61\n')
    refs = tmp_path / 'refs.csv'
    refs.write_text(f'GCF_ID,GCF_path\n{gcf},{fai}\n')
    paths = [bam(tmp_path, f'ref{i}', i+3) for i in range(4)]
    records = tmp_path / 'records.csv'
    records.write_text('Sample,CramFilePath,GCFID\n' + ''.join(
        f'ref{i},{path},{gcf}\n' for i, path in enumerate(paths)))
    config = dict(storage_mode='native', window_size=10, cpus=1, python=sys.executable,
                  memory=None, geometry_format='packed',
                  policies=dict(mate_policy='both', sampling_unit='read', min_mapq=0,
                                exclude_duplicates=False, exclude_qcfail=False))
    columns = dict(sample='Sample', cram='CramFilePath', gcf='GCFID', ref_gcf='GCF_ID', fai='GCF_path')
    indexes = []
    for policy in (PER_CONTIG, GCF_CONCAT):
        database = tmp_path / policy
        create_plan(records, refs, database, dict(config, window_policy=policy), columns)
        assert run_gcf(database, gcf) == 0
        assert collect(database)['complete']
        index = tmp_path / (policy + '.sqlite')
        build_index(database, index, True)
        indexes.append(index)
    query = bam(tmp_path, 'query', 1)
    manifest = tmp_path / 'query.tsv'
    table(manifest, ('Sample', 'GCFID', 'Depth', 'CramPath'), [
        dict(Sample=f'query-{depth}', GCFID=gcf, Depth=depth, CramPath=query) for depth in (1, 6)])
    unchanged = {path: digest(path) for path in tmp_path.rglob('*') if path.is_file()}
    outputs = []
    for i, index in enumerate(indexes):
        output = tmp_path / f'output{i}'
        result = subprocess.run([sys.executable, '-B', str(Path(__file__).resolve().parents[1] / 'bin/query_manifest_v2.py'),
                                 '--manifest', str(manifest), '--index', str(index), '--outdir', str(output),
                                 '--bin-count', '60', '--min-samples', '2', '--alpha', '.05', '--threads', '4'],
                                capture_output=True, text=True)
        assert result.returncode == 0, result.stderr
        outputs.append(output)
    assert (outputs[0] / 'summary.tsv').read_bytes() == (outputs[1] / 'summary.tsv').read_bytes()
    with (outputs[0] / 'summary.tsv').open() as handle:
        rows = list(csv.DictReader(handle, delimiter='\t'))
    assert [r['arm'] for r in rows] == ['low', 'high']
    assert all(r['windows'] == '60' and r['window_policy'] == GCF_CONCAT and r['valid_draws'] == '100' for r in rows)
    for query_id in ('q00000001', 'q00000002'):
        for name in ('query_windows.npy', 'reference_windows.npy', 'windows.tsv', 'null.tsv.gz'):
            relative = Path('queries') / query_id / name
            assert (outputs[0] / relative).read_bytes() == (outputs[1] / relative).read_bytes()
        draws = []
        for output in outputs:
            with gzip.open(output / 'queries' / query_id / 'draws.tsv.gz', 'rt') as handle:
                draws.append(list(csv.DictReader(handle, delimiter='\t')))
        for left, right in zip(*draws):
            for key in ('reference_sample', 'selected_ids_sha256', 'rng_sha256', 'units', 'actual_depth', 'jaccard', 'pearson'):
                assert left[key] == right[key]
    assert unchanged == {path: digest(path) for path in unchanged}
