"""Native cache compatibility and per-file isolation for strict mate indexing."""
import csv
import os
from pathlib import Path
import sqlite3
import subprocess
import sys

import pysam
import pytest

from test_computed_cram_index import fixture

NEXTFLOW = os.environ.get('COVERAGEMAP_TEST_NEXTFLOW')
ROOT = Path(__file__).resolve().parents[1]
LEGACY = Path(os.environ.get('COVERAGEMAP_TEST_LEGACY_INDEX_CODE', ROOT))


@pytest.mark.skipif(not NEXTFLOW, reason='Requires native Nextflow')
def test_legacy_resume_and_strict_mates(tmp_path):
    records, fna, fai, mixed = fixture(tmp_path)
    good, bad = tmp_path/'good.bam', tmp_path/'bad.bam'
    with pysam.AlignmentFile(mixed, 'rb') as src:
        with pysam.AlignmentFile(good, 'wb', header=src.header) as dst:
            for read in src:
                if read.is_read1 or read.is_read2:
                    dst.write(read)
    with pysam.AlignmentFile(good, 'rb') as src:
        with pysam.AlignmentFile(bad, 'wb', header=src.header) as dst:
            for read in src:
                dst.write(read)
                dst.write(read)
    with records.open('w', newline='') as stream:
        writer = csv.writer(stream)
        writer.writerow(('Sample','CramFilePath','GCFID','Coverage','Depth'))
        writer.writerows([('a',good,'GCF_1.1',0,0),('b',mixed,'GCF_1.1',0,0)])
    config = tmp_path/'local.config'
    config.write_text(f'conda.enabled=false\nprocess.executor="local"\n'
                      f'env.PATH="{Path(sys.executable).parent}:"+System.getenv("PATH")\n')
    recovery = tmp_path/'recovery.config'
    recovery.write_text("workflow.failOnIgnore=false\n"
                        "process { withName: 'CRAM_INDEX_BUILD:CRAM_METRICS_CHUNK' { errorStrategy = 'ignore' } }\n")
    env = dict(os.environ, NXF_SYNTAX_PARSER='v2', NXF_OFFLINE='true')

    def launch(index, entry, extra=()):
        configs = [ROOT/'nextflow.config', config]
        if index > 1:
            configs.append(recovery)
        command = [NEXTFLOW,'-C',','.join(map(str,configs)),'run',str((LEGACY if index == 1 else ROOT)/entry),
                   '-work-dir',str(tmp_path/'work'),'-with-trace',str(tmp_path/f'trace-{index}.tsv'),
                   '--records',str(records),'--fna',str(fna),'--fai',str(fai),
                   '--outdir',str(tmp_path/f'result-{index}'),'--index_chunk_buckets','1',
                   '--index_chunk_size','2','--index_chunk_cpus','2',*extra]
        completed = subprocess.run(command, cwd=tmp_path, env=env,
                                   capture_output=True, text=True, timeout=150)
        (tmp_path/f'console-{index}.log').write_text(completed.stdout+completed.stderr)
        assert completed.returncode == 0, completed.stdout+completed.stderr

    launch(1, 'build_cram_index.nf')
    with records.open('a', newline='') as stream:
        csv.writer(stream).writerows([('c',bad,'GCF_1.1',0,0),('d',good,'GCF_1.1',0,0)])
    launch(2, 'build_cram_index_mates.nf', ('-resume','--cached_index_code',str(LEGACY)))
    with (tmp_path/'trace-2.tsv').open() as stream:
        rows = list(csv.DictReader(stream, delimiter='\t'))
    assert sorted(r['status'] for r in rows if 'CRAM_METRICS_CHUNK' in r['name']) == ['CACHED','FAILED']
    assert [r['status'] for r in rows if 'CRAM_METRICS_RECOVER' in r['name']] == ['COMPLETED']
    out = tmp_path/'result-2/cram_index'
    with sqlite3.connect(out/'reference_cram.sqlite') as db:
        for table in ('reference_records','source_records','reference_metrics','reference_mate_audit'):
            assert db.execute(f'SELECT id FROM {table} ORDER BY id').fetchall() == [(1,),(4,)]
    with (out/'excluded_crams.tsv').open() as stream:
        exclusions = list(csv.DictReader(stream, delimiter='\t'))
    assert [(r['OriginalCSVRow'],r['Reason']) for r in exclusions] == [
        ('2','UNLABELLED_READS'), ('3','AMBIGUOUS_READ_IDENTITY')]
