"""Real BAM, bounded process workers, checkpoint resume and multi-root CLI smoke."""
import csv
import json
import os
from pathlib import Path
import subprocess
import sys

import pytest


def make_fixture(root):
    import pysam
    root=Path(root)
    root.mkdir(parents=True,exist_ok=False)
    fai=root/'reference.fai'
    fai.write_text('c1\t120\nc2\t80\n')
    header=dict(HD=dict(VN='1.6',SO='coordinate'),SQ=[dict(SN='c1',LN=120),dict(SN='c2',LN=80)],
                RG=[dict(ID='g',SM='fixture')])
    rows=[]
    for group,count in enumerate((3,2,1),1):
        for sample in range(count):
            bam=root/f'g{group}-s{sample}.bam'
            with pysam.AlignmentFile(bam,'wb',header=header) as handle:
                for contig in (0,1):
                    for number in range(4+sample):
                        read=pysam.AlignedSegment()
                        read.query_name=f'g{group}s{sample}-c{contig}-r{number}'
                        read.query_sequence='A'*20
                        read.flag=0
                        read.reference_id=contig
                        read.reference_start=number*9
                        read.mapping_quality=60
                        read.cigartuples=[(0,20)]
                        read.query_qualities=pysam.qualitystring_to_array('I'*20)
                        read.set_tag('RG','g')
                        handle.write(read)
            rows.append(dict(Sample=f'g{group}s{sample}',CramFilePath=str(bam),GCFID=f'GCF_{group:09d}.1'))
    records=root/'records.csv'
    with records.open('w',newline='') as handle:
        writer=csv.DictWriter(handle,list(rows[0]));writer.writeheader();writer.writerows(reversed(rows))
    references=root/'references.csv'
    references.write_text('GCF_ID,GCF_path\n'+''.join(f'GCF_{g:09d}.1,{fai}\n' for g in (1,2,3)))
    return records,references


def test_real_bam_cli_batches_resume_and_merge(tmp_path):
    pytest.importorskip('pysam')
    records,refs=make_fixture(tmp_path/'inputs')
    script=Path(__file__).resolve().parents[1]/'bin/ref_build_pipeline.py'
    env=dict(os.environ,PYTHONDONTWRITEBYTECODE='1')
    def run(*args):
        completed=subprocess.run([sys.executable,'-B',str(script),*map(str,args)],cwd=tmp_path,
            env=env,stdout=subprocess.PIPE,stderr=subprocess.PIPE,text=True)
        assert completed.returncode==0,completed.stdout+'\n'+completed.stderr
    plan=tmp_path/'plan'
    out=tmp_path/'out'
    run('plan','--records',records,'--references',refs,'--outdir',plan,'--chunk-size',2,'--cpus',2)
    receipts=tmp_path/'receipts';receipts.mkdir()
    tasks=[]
    for chunk in sorted((plan/'chunks').glob('*.json')):
        root=out/'batches'/chunk.stem
        args=('chunk','--chunk',chunk,'--checkpoints',root/'checkpoints','--database',root/'database',
              '--receipt',receipts/chunk.name,'--cpus',2)
        run(*args)
        row=json.loads((receipts/chunk.name).read_text())
        assert row['workers']==min(2,row['completed'])
        assert not any(record['reused'] for record in row['records'])
        tasks.append(args)
    run(*tasks[0])
    replay=json.loads((receipts/'chunk-000000.json').read_text())
    assert all(row['reused'] for row in replay['records'])
    run('collect','--plan',plan,'--receipts',receipts,'--database',out/'database',
        '--receipt',tmp_path/'collected.json','--cpus',2)
    run('index','--database',out/'database','--receipt',tmp_path/'indexed.json','--cpus',2)
    summary=json.loads((out/'database/SUMMARY.json').read_text())
    assert summary['complete'] and summary['records_completed']==6 and summary['gcfs']==3
    assert (out/'database/SUCCESS').is_file()
    assert (out/'reference.sqlite').is_file()
