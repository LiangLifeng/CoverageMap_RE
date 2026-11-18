"""Detection rules and query-to-center parity with the Bin Count grid."""

import csv
import json
from pathlib import Path
import sqlite3
import subprocess
import sys

import numpy as np
import pytest

from bin.bin_count_eval_math import METRICS, distances, score
from bin.query_detection import read_rules, select_rule
from bin.query_manifest import table
from bin.ref_depth_catalog import collect, create_plan
from bin.ref_depth_worker import run_gcf
from bin.reference_index import build_index

ROOT = Path(__file__).resolve().parents[1]
pysam = pytest.importorskip('pysam')


def rule(lower, upper, distance='euclidean', threshold=None):
    return dict(min_genome_bp=lower, max_genome_bp_exclusive=upper, bin_count=3,
                distance=distance, max_distance=threshold)


def test_rules_validate_both_domains_and_exact_boundaries(tmp_path):
    path = tmp_path/'rules.json'
    path.write_text(json.dumps(dict(bacteria=[rule(40, None), rule(0, 20)], viral=[rule(0, None)])))
    rules = read_rules(path, 'bacteria')
    assert [select_rule(rules, bp) for bp in (19, 20, 39, 40)] == [rules[0], None, None, rules[1]]
    assert len(read_rules(path, 'viral')) == 1
    for invalid in (rule(10, 50), rule(21, None, 'pearson'), rule(21, None, threshold=True),
                    rule(21, None, threshold=float('nan')), rule(21, None, threshold=-1)):
        path.write_text(json.dumps(dict(bacteria=[rule(0, 20), invalid], viral=[rule(0, None)])))
        with pytest.raises(ValueError):
            read_rules(path, 'viral')
    path.write_text('{"viral":[],"viral":[]}')
    with pytest.raises(ValueError, match='Duplicate'):
        read_rules(path, 'viral')


def bam(root, sample, positions):
    path = root/(sample+'.bam')
    with pysam.AlignmentFile(path, 'wb', header={'HD': {'VN': '1.6'}, 'SQ': [{'SN': 'one', 'LN': 15}]}) as handle:
        for number, position in enumerate(positions):
            read = pysam.AlignedSegment()
            read.query_name = f'{sample}-{number}'
            read.query_sequence = 'A'*4
            read.flag = 0
            read.reference_id = 0
            read.reference_start = position
            read.mapping_quality = 60
            read.cigar = [(0, 4)]
            handle.write(read)
    return path


def dense_profile(positions):
    depth = np.zeros(31)
    for position in positions:
        depth[position:position+4] += 1
    depths = np.array([depth[:10].mean(), depth[10:20].mean(), depth[20:].mean()])
    return depths * (1e6 / depths.sum())


def test_detection_cli_grid_parity_thresholds_and_self_exclusion(tmp_path):
    gcf = 'GCF_000000001.1'
    columns = dict(sample='Sample', cram='CramFilePath', gcf='GCFID', ref_gcf='GCF_ID', fai='GCF_path')
    ref_positions = {'a': [1, 8, 10], 'b': [2, 3, 4, 8], 'self': [8, 9, 11]}
    inputs = {sample: bam(tmp_path, sample, starts) for sample, starts in ref_positions.items()}
    fai = tmp_path/'ref.fai'
    fai.write_text('one\t15\ntwo\t16\n')
    records, refs = tmp_path/'records.csv', tmp_path/'refs.csv'
    with records.open('w') as handle:
        writer = csv.writer(handle)
        writer.writerow([columns[key] for key in ('sample', 'cram', 'gcf')])
        writer.writerows((sample, path, gcf) for sample, path in inputs.items())
    refs.write_text(f'GCF_ID,GCF_path\n{gcf},{fai}\n')
    config = dict(storage_mode='native', window_size=7, cpus=1, python=sys.executable, memory=None,
                  geometry_format='packed', policies=dict(mate_policy='both', sampling_unit='read', min_mapq=0,
                                                         exclude_duplicates=False, exclude_qcfail=False))
    database = tmp_path/'db'
    create_plan(records, refs, database, config, columns)
    assert run_gcf(database, gcf) == 0
    assert collect(database)['complete']
    index = tmp_path/'index.sqlite'
    build_index(database, index, True)
    index.chmod(0o644)
    with sqlite3.connect(index) as db:
        db.execute('INSERT INTO metadata VALUES (?, ?)', ('taxonomy_partition', '{"domain": "bacteria"}'))
    index.chmod(0o444)
    manifest = tmp_path/'query.tsv'
    table(manifest, ('Sample', 'GCFID', 'Depth', 'CramPath'), [
        dict(Sample='renamed', GCFID=gcf, Depth=1, CramPath=str(inputs['self'])),
        dict(Sample='missing', GCFID='GCF_999999999.1', Depth=1, CramPath=str(inputs['self'])),
    ])
    profile = dense_profile(ref_positions['self'])
    matrix = np.stack([dense_profile(ref_positions[s]) for s in ('a', 'b')])
    rules = tmp_path/'rules.json'
    command = [sys.executable, '-B', str(ROOT/'bin/query_manifest.py'), '--manifest', str(manifest),
               '--index', str(index), '--min-samples', '1', '--jobs', '2', '--detection-config', str(rules)]

    def run(out, domain='bacteria', resume=False):
        flags = ['--outdir', str(out), '--domain', domain] + (['--resume'] if resume else [])
        return subprocess.run(command+flags, capture_output=True, text=True)

    for metric in METRICS:
        _, _, centers, *_ = score(matrix, [gcf, gcf], metric)
        expected = float(distances(profile, centers, metric)[0, 0])
        rules.write_text(json.dumps(dict(bacteria=[rule(0, 31), rule(31, None, metric, expected)],
                                        viral=[rule(0, None)])))
        out = tmp_path/metric
        result = run(out)
        assert result.returncode == 0, result.stderr
        with (out/'summary.tsv').open() as handle:
            actual, missing = csv.DictReader(handle, delimiter='\t')
        assert actual['distance'] == metric
        assert float(actual['distance_value']) == expected
        assert actual['detection_status'] == 'DETECTED'
        assert actual['center_reference_count'] == '2' and actual['excluded_self_count'] == '1'
        assert actual['min_genome_bp'] == '31'
        assert missing['status'] == 'GCF_NOT_INDEXED'
        query = out/'queries/q00000001'
        np.testing.assert_array_equal(np.load(query/'bin_edges.npy'), [0, 10, 20, 31])
        np.testing.assert_array_equal(np.load(query/'query_profile.npy'), profile)
        np.testing.assert_array_equal(np.load(query/'reference_center.npy'), centers[0])
        assert run(out, resume=True).returncode == 0
    assert run(tmp_path/'wrong-domain', 'viral').returncode == 2
    for label, threshold, detection in (('score-only', None, 'THRESHOLD_UNSET'), ('negative', 0, 'NOT_DETECTED')):
        rules.write_text(json.dumps(dict(bacteria=[rule(31, None, threshold=threshold)])))
        out = tmp_path/label
        assert run(out).returncode == 0
        assert json.loads((out/'queries/q00000001/summary.json').read_text())['detection_status'] == detection
    assert run(tmp_path/'euclidean', resume=True).returncode == 2
