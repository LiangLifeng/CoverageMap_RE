"""Natural-order virtual coordinates, immutable migration, and v2 query smoke."""

import csv
import json
from pathlib import Path
import subprocess
import sys

import numpy as np
import pytest

from bin.coverage_core import (BlockPool, Reference, canonical_reference, dense_coverage, extract_pool,
                               layout_metadata, merge_references, reference_from_layout,
                               reproject_pool, window_sums)
from bin.gcf_windows import GCF_CONCAT, PER_CONTIG, concat_edges, coordinate_rows
from bin.native_concordance import NativeCache
from bin.native_reference import covered_sums, validate_native
from bin.ref_depth_catalog import collect, create_plan, digest, sample_folder
from bin.ref_depth_database import parser, prepare
from bin.ref_depth_worker import run_gcf
from bin.reference_index import ReferenceIndex, build_index

GCF = 'GCF_000000001.1'
LENGTHS = {'contig_1': 1300, 'contig_2': 800, 'contig_10': 700}
OLD_ORDER = ('contig_10', 'contig_2', 'contig_1')
COLUMNS = dict(sample='Sample', cram='CramFilePath', gcf='GCFID', ref_gcf='GCF_ID', fai='GCF_path')
POLICIES = dict(mate_policy='both', sampling_unit='read', min_mapq=0,
                exclude_duplicates=False, exclude_qcfail=False)


def test_natural_order_and_cross_contig_windows():
    old = Reference(OLD_ORDER, tuple(LENGTHS[n] for n in OLD_ORDER), ('',)*3)
    layout = canonical_reference(old)
    assert layout.names == tuple(LENGTHS)
    assert list(layout.offsets.values()) == [0, 1300, 2100]
    ids = ('ctg_10a', 'ctg', 'ctg_2a', 'ctg_02a', 'ctg_alpha')
    mixed = Reference(ids, (1,)*len(ids), ('',)*len(ids))
    assert canonical_reference(mixed).names == ('ctg', 'ctg_02a', 'ctg_2a', 'ctg_10a', 'ctg_alpha')
    saved = layout_metadata(layout)
    assert reference_from_layout(saved) == layout
    for altered in (dict(saved, offsets=[0, 1300, 2200]),
                    dict(saved, names=list(OLD_ORDER)),
                    dict(saved, lengths=[1300, 700, 800])):
        with pytest.raises(ValueError, match='order or offsets'):
            reference_from_layout(altered)
    with pytest.raises(ValueError, match='layouts differ'):
        merge_references(layout, old)

    size = concat_edges(layout.length, size=1000)
    count = concat_edges(layout.length, count=3)
    np.testing.assert_array_equal(size, [0, 1000, 2000, 2800])
    np.testing.assert_array_equal(count, [0, 933, 1866, 2800])
    pieces = list(coordinate_rows(layout, size, 1000))
    assert [(r['bin_index'], r['contig'], r['local_start'], r['local_end']) for r in pieces] == [
        (0, 'contig_1', 0, 1000), (1, 'contig_1', 1000, 1300), (1, 'contig_2', 0, 700),
        (2, 'contig_2', 700, 800), (2, 'contig_10', 0, 700)]
    assert len(list(coordinate_rows(layout, concat_edges(2800, count=1), None))) == 3


@pytest.mark.parametrize('length, kwargs, expected', [
    (1, {'size': 1}, [0, 1]), (3, {'size': 10}, [0, 3]),
    (4, {'count': 4}, [0, 1, 2, 3, 4]), (5, {'size': 2}, [0, 2, 4, 5]),
    (2**63-1, {'count': 2}, [0, (2**63-1)//2, 2**63-1]),
    (2**63-1, {'size': 2**63}, [0, 2**63-1]),
])
def test_edge_limits(length, kwargs, expected):
    np.testing.assert_array_equal(concat_edges(length, **kwargs), expected)


@pytest.mark.parametrize('length, kwargs', [
    (3, {'count': 4}), (3, {'count': 0}), (3, {'size': 0}),
    (3, {'size': -1}), (3, {'count': True}), (3, {'size': 2, 'count': 2}),
    (3, {}), (0, {'size': 2}), (2**63, {'count': 1}),
])
def test_invalid_edges(length, kwargs):
    with pytest.raises(ValueError):
        concat_edges(length, **kwargs)


def test_mass_and_covered_union_after_reprojection():
    old = Reference(('contig_2', 'contig_1', 'contig_3'), (4, 5, 1), ('',)*3)
    layout = canonical_reference(old)
    old_pool = BlockPool(np.array([0, 1, 7, 9]), np.array([2, 3, 9, 10]),
                         np.array([0, 2, 3, 4]), 10)
    pool = reproject_pool(old_pool, old, layout)
    np.testing.assert_array_equal(pool.starts, [5, 6, 3, 9])
    np.testing.assert_array_equal(pool.ends, [7, 8, 5, 10])
    np.testing.assert_array_equal(pool.unit_bases, old_pool.unit_bases)
    np.testing.assert_array_equal(pool.indptr, old_pool.indptr)
    edges = concat_edges(layout.length, size=6)
    dense = dense_coverage(pool, np.arange(len(pool)))
    np.testing.assert_array_equal(dense, [0, 0, 0, 1, 1, 1, 2, 1, 0, 1])
    sums = window_sums(pool, np.arange(len(pool)), edges)
    np.testing.assert_array_equal(sums, [3, 4])
    assert sums.sum() == old_pool.total_bases == 7
    np.testing.assert_array_equal(covered_sums(pool, edges), [3, 3])
    np.testing.assert_array_equal(covered_sums(BlockPool.from_intervals([], 10), edges), [0, 0])
    with pytest.raises(ValueError, match='crosses a source contig boundary'):
        reproject_pool(BlockPool.from_intervals([[3, 5]], 10), old, layout)


def write_bam(pysam, path, order, blocks):
    header = {'SQ': [dict(SN=name, LN=LENGTHS[name]) for name in order]}
    with pysam.AlignmentFile(path, 'wb', header=header) as handle:
        for i, (name, start, length) in enumerate(blocks):
            read = pysam.AlignedSegment()
            read.query_name = f'{path.stem}-{i}'
            read.query_sequence = 'A'*length
            read.reference_id = order.index(name)
            read.reference_start = start
            read.mapping_quality = 60
            read.cigar = [(0, length)]
            handle.write(read)


def result(database, sample):
    attempt = database/GCF/'samples'/sample_folder(sample)/'attempt-000001'
    return attempt, json.loads((attempt/'result.json').read_text())


def test_native_build_migration_and_v2_query(tmp_path):
    pysam = pytest.importorskip('pysam')
    fai, refs, records = (tmp_path/name for name in ('ref.fai', 'ref.csv', 'records.csv'))
    fai.write_text(''.join(f'{n}\t{LENGTHS[n]}\n' for n in OLD_ORDER))
    refs.write_text(f'GCF_ID,GCF_path\n{GCF},{fai}\n')
    inputs = {name: tmp_path/f'{name}.bam' for name in ('a', 'b', 'query')}
    common = [('contig_1', 1250, 50), ('contig_2', 0, 40), ('contig_2', 650, 100),
              ('contig_2', 700, 70), ('contig_10', 0, 30)]
    write_bam(pysam, inputs['a'], ('contig_2', 'contig_10', 'contig_1'),
              [*common, ('contig_1', 880, 40)])
    write_bam(pysam, inputs['b'], OLD_ORDER, [*common, ('contig_2', 680, 35)])
    write_bam(pysam, inputs['query'], tuple(LENGTHS),
              [('contig_1', 1260, 40), ('contig_2', 695, 40)])
    fasta = tmp_path/'reference.fa'
    fasta.write_text(''.join(f'>{n}\n'+'A'*LENGTHS[n]+'\n' for n in reversed(tuple(LENGTHS))))
    query_pool, query_layout = extract_pool(inputs['query'], reference_path=fasta, sort_contigs=True)
    assert query_layout.names == tuple(LENGTHS)
    np.testing.assert_array_equal(query_pool.starts, [1260, 1995])
    records.write_text('Sample,CramFilePath,GCFID\n'+''.join(
        f'{n},{inputs[n]},{GCF}\n' for n in ('a','b')))

    old, fresh, migrated = (tmp_path/name for name in ('old', 'fresh', 'migrated'))
    old_config = dict(storage_mode='native', window_size=1000, cpus=1, python=sys.executable,
                      geometry_format='packed', policies=POLICIES)
    create_plan(records, refs, old, old_config, COLUMNS)
    assert run_gcf(old, GCF) == 0 and collect(old)['complete']
    old_index = tmp_path/'old.sqlite'
    build_index(old, old_index, True)
    old_snapshot = {path: digest(path) for path in [old_index, *old.rglob('*')] if path.is_file()}
    old_attempt, _ = result(old, 'a')
    np.testing.assert_array_equal(np.load(old_attempt/'edges.npy'), [0, 700, 1500, 2500, 2800])

    def args(database, *geometry):
        return parser().parse_args(['prepare', '--records', str(records), '--references', str(refs),
                                    '--outdir', str(database), '--cpus', '1', *geometry])

    with pytest.raises(ValueError, match='Too many bins'):
        prepare(args(tmp_path/'bad', '--bin-count', '2801'))
    prepare(args(fresh, '--bin-size', '1000'))
    assert run_gcf(fresh, GCF) == 0 and collect(fresh)['complete']
    prepare(args(migrated, '--bin-count', '3', '--reuse-geometry-from', str(old), '--alignment-capture', 'legacy'))

    def reject_extraction(*_args, **_kwargs):
        raise AssertionError('Verified old geometry must be reprojected')

    assert run_gcf(migrated, GCF, reject_extraction) == 0 and collect(migrated)['complete']
    assert {path: digest(path) for path in old_snapshot} == old_snapshot
    attempt, record = result(fresh, 'a')
    np.testing.assert_array_equal(np.load(attempt/'edges.npy'), [0, 1000, 2000, 2800])
    assert record['layout']['names'] == list(LENGTHS)
    assert record['layout']['offsets'] == [0, 1300, 2100]
    with np.load(attempt/'native_coverage.npz') as coverage:
        np.testing.assert_array_equal(coverage['aligned_bases'], [40, 140, 150])
        np.testing.assert_array_equal(coverage['covered_bases'], [40, 140, 100])
        np.testing.assert_array_equal(coverage['segment_ptr'], [0, 1, 3, 5])
        np.testing.assert_allclose(coverage['mean_depth'], [40/1000, 140/1000, 150/800])
        np.testing.assert_allclose(coverage['covered_fraction'], [40/1000, 140/1000, 100/800])
    migrated_attempt, migrated_record = result(migrated, 'a')
    assert migrated_record['reused_geometry']['source_layout']['names'] == list(OLD_ORDER)
    assert migrated_record['aligned_bases'] == record['aligned_bases'] == 330
    np.testing.assert_array_equal(np.load(migrated_attempt/'edges.npy'), [0, 933, 1866, 2800])
    with np.load(migrated_attempt/'native_coverage.npz') as coverage:
        np.testing.assert_array_equal(coverage['aligned_bases'], [40, 90, 200])
        assert coverage['covered_bases'].sum() == 280
    original_edges = (attempt/'edges.npy').read_bytes()
    np.save(attempt/'edges.npy', [0, 999, 2000, 2800])
    with pytest.raises(ValueError, match='edges mismatch'):
        validate_native(attempt, record)
    (attempt/'edges.npy').write_bytes(original_edges)
    index = tmp_path/'native.sqlite'
    build_index(fresh, index, True)
    indexed = ReferenceIndex(index)
    try:
        cache = NativeCache(indexed, None)
        with pytest.raises(ValueError, match='policy or parameters'):
            cache.check_window_policy(dict(record, native=dict(record['native'], window_policy=PER_CONTIG)))
        with pytest.raises(ValueError, match='policy or parameters'):
            cache.check_window_policy(dict(record, window_size=999))
    finally:
        indexed.close()

    manifest = tmp_path/'queries.tsv'
    manifest.write_text('Sample\tGCFID\tDepth\tCramPath\n'+''.join(
        f'{name}\t{GCF}\t{depth}\t{inputs["query"]}\n' for name, depth in [('low', 1), ('high', 6)]))
    entry = Path(__file__).resolve().parents[1]/'bin/query_manifest_v2.py'
    base = [sys.executable, '-B', str(entry), '--manifest', str(manifest), '--alpha', '.05',
            '--min-samples', '1', '--threads', '2']
    out = tmp_path/'scores'
    legacy_out = tmp_path/'legacy-rebinned'
    legacy = subprocess.run([*base, '--index', str(old_index), '--outdir', str(legacy_out),
                             '--bin-size', '1000'], capture_output=True, text=True)
    assert legacy.returncode == 0, legacy.stderr
    completed = subprocess.run([*base, '--index', str(index), '--outdir', str(out), '--bin-size', '1000'],
                               capture_output=True, text=True)
    assert completed.returncode == 0, completed.stderr
    assert (legacy_out/'summary.tsv').read_bytes() == (out/'summary.tsv').read_bytes()
    assert {path: digest(path) for path in old_snapshot} == old_snapshot
    with (out/'summary.tsv').open() as handle:
        rows = list(csv.DictReader(handle, delimiter='\t'))
    assert [r['arm'] for r in rows] == ['low', 'high']
    assert all(r['window_policy'] == GCF_CONCAT and r['valid_draws'] == '100' for r in rows)
    for row in rows:
        directory = out/'queries'/row['query_id']
        np.testing.assert_array_equal(np.load(directory/'edges.npy'), [0, 1000, 2000, 2800])
        np.testing.assert_array_equal(np.load(directory/'query_windows.npy'), [0, 45, 35])
        assert np.load(directory/'reference_windows.npy').shape == (100, 3)
        for filename in ('query_windows.npy', 'reference_windows.npy', 'edges.npy', 'null.tsv.gz'):
            assert (legacy_out/'queries'/row['query_id']/filename).read_bytes() == (directory/filename).read_bytes()
    count_out = tmp_path/'count-scores'
    count_index = tmp_path/'count.sqlite'
    build_index(migrated, count_index, True)
    completed = subprocess.run([*base, '--index', str(count_index), '--outdir', str(count_out),
                                '--bin-count', '3'], capture_output=True, text=True)
    assert completed.returncode == 0, completed.stderr
    np.testing.assert_array_equal(np.load(count_out/'queries/q00000001/edges.npy'), [0, 933, 1866, 2800])
