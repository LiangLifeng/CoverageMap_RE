"""Exercise geometry-only BAM sampling, checkpoint replay, and failure recovery."""
import argparse
import csv
import json
from pathlib import Path
import subprocess
import sys

import pytest

from bin import sample_cram_modes as sampler
from bin.ref_depth_catalog import collect, create_plan
from bin.ref_depth_worker import run_gcf
from bin.reference_index import build_index


def fixture(root):
    import pysam
    gcf = 'GCF_000000001.1'
    fai = root / 'ref.fai'
    fai.write_text('one\t80\n')
    bam = root / 'input.bam'
    with pysam.AlignmentFile(bam, 'wb', header={'SQ': [{'SN': 'one', 'LN': 80}]}) as handle:
        for pair in range(12):
            for mate in range(2):
                read = pysam.AlignedSegment()
                read.query_name = f'pair{pair}'
                read.flag = 65 if mate == 0 else 129
                read.reference_id = 0
                read.reference_start = (pair * 4 + mate * 2) % 70
                read.mapping_quality = 60
                read.cigar = [(0, 4)]
                handle.write(read)  # No SEQ or QUAL, matching the benchmark output.
    records, refs = root / 'records.csv', root / 'refs.csv'
    records.write_text(f'Sample,CramFilePath,GCFID\nref,{bam},{gcf}\n')
    refs.write_text(f'GCF_ID,GCF_path\n{gcf},{fai}\n')
    config = dict(storage_mode='native', window_size=10, window_policy='gcf-concat-v1',
                  cpus=1, python=sys.executable, memory=None, geometry_format='packed',
                  policies=dict(mate_policy='both', sampling_unit='read', min_mapq=0,
                                exclude_duplicates=False, exclude_qcfail=False))
    columns = dict(sample='Sample', cram='CramFilePath', gcf='GCFID', ref_gcf='GCF_ID', fai='GCF_path')
    database = root / 'database'
    create_plan(records, refs, database, config, columns)
    assert run_gcf(database, gcf) == 0
    assert collect(database)['complete']
    index = root / 'index.sqlite'
    build_index(database, index, True)
    manifest = root / 'input.tsv'
    manifest.write_text(f'Sample\tGCFID\tDepth\tCramPath\nquery\t{gcf}\t1.2\t{bam}\n')
    output = root / 'output'
    command = [sys.executable, '-B', sampler.__file__, '--manifest', str(manifest),
               '--index', str(index), '--outdir', str(output), '--target-depth', '0.5x',
               '--reference-fai', str(fai), '--threads', '2', '--seed', '401']
    return output, command


@pytest.mark.parametrize('value', ['0.5', '0.5x', '0.5×'])
def test_depth_units(value):
    assert sampler.parse_depth(value) == .5


@pytest.mark.parametrize('value', ['0', '-1x', 'NaN', 'inf', '0.5xx'])
def test_invalid_depth(value):
    with pytest.raises(argparse.ArgumentTypeError):
        sampler.parse_depth(value)


def test_sample_resume_repair_and_corruption(tmp_path):
    import pysam
    output, command = fixture(tmp_path)
    result = subprocess.run(command, capture_output=True, text=True)
    assert result.returncode == 0, result.stderr
    summary = json.loads((output / 'summary.json').read_text())
    assert summary['read_manifest_rows'] == summary['fragment_manifest_rows'] == 1
    bams = sorted(output.glob('*/*.bam'))
    original = {path: (path.read_bytes(), path.stat().st_mtime_ns) for path in bams}
    fragment = next(output.glob('fragment/*.bam'))
    with pysam.AlignmentFile(fragment, 'rb') as alignment:
        names = [read.query_name for read in alignment.fetch(until_eof=True)]
    assert len(names) == 10 and len(set(names)) == 5
    result = subprocess.run(command + ['--resume'], capture_output=True, text=True)
    assert result.returncode == 0 and result.stdout.count('resumed verified BAM') == 2
    assert original == {path: (path.read_bytes(), path.stat().st_mtime_ns) for path in bams}
    assert subprocess.run(command, capture_output=True).returncode == 2
    assert subprocess.run(command + ['--resume', '--seed', '402'], capture_output=True).returncode == 2
    # Simulate an interruption before the second mode completed.
    fragment.unlink()
    (output / '.checkpoints/q00000001_fragment.json').unlink()
    result = subprocess.run(command + ['--resume'], capture_output=True, text=True)
    assert result.returncode == 0 and result.stdout.count('resumed verified BAM') == 1
    read_bam = next(output.glob('read/*.bam'))
    assert (read_bam.read_bytes(), read_bam.stat().st_mtime_ns) == original[read_bam]
    assert fragment.read_bytes() == original[fragment][0]
    with read_bam.open('ab') as handle:
        handle.write(b'corruption')
    result = subprocess.run(command + ['--resume'], capture_output=True, text=True)
    assert result.returncode == 1 and 'checksum mismatch' in result.stderr
    for mode in ('read', 'fragment'):
        with (output / f'manifest_{mode}.tsv').open() as handle:
            assert list(csv.DictReader(handle, delimiter='\t')) == []


def test_missing_fai_fails_before_output(tmp_path):
    output, command = fixture(tmp_path)
    result = subprocess.run(command + ['--reference-fai', str(tmp_path / 'missing.fai')],
                            capture_output=True, text=True)
    assert result.returncode == 2 and not output.exists()
