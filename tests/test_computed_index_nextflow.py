"""Opt-in native strict-parser chunk/index/resume acceptance, synthetic files only."""
import csv
import json
import os
from pathlib import Path
import sqlite3
import subprocess
import sys

import pytest
from test_computed_cram_index import fixture

NEXTFLOW=os.environ.get('COVERAGEMAP_TEST_NEXTFLOW')
pytestmark=pytest.mark.skipif(not NEXTFLOW,reason='Set COVERAGEMAP_TEST_NEXTFLOW')
ROOT=Path(__file__).resolve().parents[1]


def test_native_chunks_resume_and_source_dependency(tmp_path):
    records,fna,fai,_=fixture(tmp_path)
    config=tmp_path/'local.config'
    config.write_text('conda.enabled=false\nprocess.executor="local"\n' +
        f'env.PATH="{Path(sys.executable).parent}:"+System.getenv("PATH")\n')
    base=[NEXTFLOW,'-C',f'{ROOT}/nextflow.config,{config}','run',str(ROOT/'build_cram_index.nf'),
          '-work-dir',str(tmp_path/'work'),'--records',str(records),'--fna',str(fna),'--fai',str(fai),
          '--index_chunk_size','2','--index_chunk_cpus','2']
    env=dict(os.environ,NXF_SYNTAX_PARSER='v2',NXF_OFFLINE='true')
    for n,extra in [(1,[]),(2,['-resume']),(3,['-resume'])]:
        if n==3:
            src=tmp_path/'unlabelled.bam';st=src.stat()
            os.utime(src,ns=(st.st_atime_ns,st.st_mtime_ns+2_000_000_000))
        trace=tmp_path/f'trace-{n}.tsv'
        out=tmp_path/('output' if n<3 else 'output-source-change')
        cmd=base+['--outdir',str(out),'-with-trace',str(trace)]+extra
        result=subprocess.run(cmd,cwd=tmp_path,env=env,text=True,capture_output=True,timeout=180)
        (tmp_path/f'command-{n}.json').write_text(json.dumps(cmd,indent=2))
        (tmp_path/f'console-{n}.log').write_text(result.stdout+result.stderr)
        assert result.returncode==0,result.stdout+result.stderr
        rows=list(csv.DictReader(trace.open(),delimiter='\t'))
        chunks=[r for r in rows if 'CRAM_METRICS_CHUNK' in r['name']]
        assert len(chunks)==2
        expected=['COMPLETED','COMPLETED'] if n==1 else ['CACHED','CACHED'] if n==2 else ['CACHED','COMPLETED']
        assert sorted(r['status'] for r in chunks)==sorted(expected)
        with sqlite3.connect(out/'cram_index/reference_cram.sqlite') as db:
            assert db.execute('SELECT COUNT(*) FROM reference_metrics').fetchone()[0]==3
            assert db.execute('SELECT read1_bases,read2_bases FROM reference_metrics WHERE id=1').fetchone()==(30,20)
            assert float(db.execute('SELECT Depth FROM reference_records WHERE id=1').fetchone()[0])==57/150
    stub=subprocess.run(base+['--outdir',str(tmp_path/'stub'),'-stub-run'],cwd=tmp_path,
                        env=env,text=True,capture_output=True,timeout=180)
    (tmp_path/'stub.log').write_text(stub.stdout+stub.stderr)
    assert stub.returncode==0,stub.stdout+stub.stderr
