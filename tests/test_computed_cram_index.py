import csv
import hashlib
import json
from pathlib import Path
import sqlite3

import pytest

from bin.cram_index_metrics import compute, merge
from bin.raw_cram_query import candidates


def fixture(tmp_path, kind='bam'):
    pysam = pytest.importorskip('pysam')
    fasta = tmp_path / 'GCF_1.1.fna'
    fasta.write_text('>chr1\n' + 'A'*100 + '\n>chr2\n' + 'C'*50 + '\n')
    pysam.faidx(str(fasta))
    header = {'HD': {'VN': '1.6'}, 'SQ': [{'SN': 'chr1', 'LN': 100,
              'M5': hashlib.md5(b'A'*100).hexdigest()}]}
    bam = tmp_path/'known.bam'
    # Independent expected geometry: R1=30, R2=20, unlabelled=7 bases;
    # overall union48, full denominator150 (missing header chr2 contributes zero).
    entries = [('a',65,0,[(0,10),(2,5),(0,5)]), ('a',129,5,[(0,10)]),
               ('b',65,30,[(0,4),(3,3),(0,6)]), ('b',129,38,[(0,4)]),
               ('single',0,60,[(0,7)]), ('dup',1089,80,[(0,5)]),
               ('qc',641,85,[(0,6)]), ('secondary',256,10,[(0,12)]),
               ('supp',2048,10,[(0,12)])]
    with pysam.AlignmentFile(bam,'wb',header=header) as out:
        for name,flag,pos,cigar in entries:
            a=pysam.AlignedSegment();a.query_name=name;a.flag=flag;a.reference_id=0
            a.reference_start=pos;a.mapping_quality=60;a.cigartuples=cigar
            a.query_sequence='A'*sum(n for op,n in cigar if op in (0,1,4,7,8))
            out.write(a)
    if kind=='cram':
        path=tmp_path/'known.cram'
        with pysam.AlignmentFile(bam,'rb') as src, pysam.AlignmentFile(
                path,'wc',header=src.header,reference_filename=str(fasta)) as out:
            for a in src:out.write(a)
    else:path=bam
    unlabelled=tmp_path/'unlabelled.bam'
    with pysam.AlignmentFile(unlabelled,'wb',header=header) as out:
        a=pysam.AlignedSegment();a.query_name='single';a.flag=0;a.reference_id=0
        a.reference_start=10;a.mapping_quality=60;a.cigartuples=[(0,10)];a.query_sequence='A'*10
        out.write(a)
    records=tmp_path/'records.csv'
    with records.open('w',newline='') as handle:
        w=csv.writer(handle);w.writerow(['Sample','CramFilePath','GCFID','Coverage','Depth'])
        w.writerows([['s',str(path),'GCF_1.1','bogus','nan'],
                     ['u',str(unlabelled),'GCF_1.1','999','-1'],
                     ['s',str(path),'GCF_1.1','','not-used']])
    fna,fai=tmp_path/'fna.csv',tmp_path/'fai.csv'
    fna.write_text(f'GCF_ID,GCF_path\nGCF_1.1,{fasta}\n')
    fai.write_text(f'GCF_ID,FAI_path\nGCF_1.1,{fasta}.fai\n')
    return records,fna,fai,path


def build(tmp_path, jobs=1, kind='bam'):
    records,fna,fai,path=fixture(tmp_path,kind)
    bundle=tmp_path/'rows'
    bundle.mkdir()
    rows=list(csv.DictReader(records.open()))
    unique={tuple(row[k] for k in ('Sample','CramFilePath','GCFID')) for row in rows}
    chunks=[]
    for number,identity in enumerate(sorted(unique),1):
        row_file=bundle/f'chunk-{number:06d}.json'
        row_file.write_text(json.dumps([[*identity,str(tmp_path/'GCF_1.1.fna'),str(tmp_path/'GCF_1.1.fna.fai')]]))
        output=tmp_path/(row_file.stem+'.sqlite')
        compute(row_file,output,jobs,0,'geometry_only');chunks.append(output)
    output=tmp_path/'final.sqlite'
    merge(records,chunks,fna,fai,output,0,'geometry_only')
    return output,chunks,bundle,fna,fai


@pytest.mark.parametrize('kind',['bam','cram'])
def test_dense_expected_metrics_schema_and_old_query(tmp_path,kind):
    output,_,_,_,_=build(tmp_path,kind=kind)
    with sqlite3.connect(output) as db:
        db.row_factory=sqlite3.Row
        assert [tuple(r) for r in db.execute('PRAGMA table_info(reference_records)')]==[
            (0,'id','INTEGER',0,None,1),(1,'Sample','TEXT',1,None,0),
            (2,'CramFilePath','TEXT',1,None,0),(3,'GCFID','TEXT',1,None,0),
            (4,'Coverage','TEXT',1,None,0),(5,'Depth','TEXT',1,None,0)]
        assert [r[0] for r in db.execute('SELECT id FROM reference_records ORDER BY id')]==[1,3]
        r=dict(db.execute('SELECT * FROM reference_metrics WHERE id=1').fetchone())
        assert (r['genome_length'],r['aligned_bases'],r['covered_bases'],r['read_count'],r['fragment_count'])==(150,57,48,7,5)
        assert (r['read1_bases'],r['read2_bases'],r['unlabelled_bases'])==(30,20,7)
        assert (r['read1_count'],r['read2_count'],r['unlabelled_count'])==(3,3,1)
        assert r['Depth']==pytest.approx(57/150)
        assert r['Coverage']==pytest.approx(32)
        assert r['read1_depth']==pytest.approx(.2)
        assert r['read2_depth']==pytest.approx(20/150)
        row=db.execute('SELECT Depth,Coverage FROM reference_records WHERE id=1').fetchone()
        assert tuple(map(float,row))==pytest.approx((57/150,32))
        assert db.execute('SELECT Depth FROM source_records WHERE id=1').fetchone()[0]=='nan'
        assert db.execute('SELECT COUNT(*) FROM reference_metrics').fetchone()[0]==3
        assert db.execute('SELECT scan_complete FROM reference_mate_audit WHERE id=2').fetchone()[0]==1
        assert json.loads(dict(db.execute('SELECT * FROM metadata'))['schema'])=='coverage-raw-cram-index-1'
        available,refs=candidates(db,{'Sample':'query','GCFID':'GCF_1.1','CramPath':'/not-self','Depth':'.1'},True,0)
        assert available==1 and len(refs)==1  # duplicate source path collapsed as before
        assert float(refs[0]['Depth'])==pytest.approx(57/150)
        assert db.execute('SELECT COUNT(*) FROM reference_cram_catalog').fetchone()[0]==2
        assert db.execute('PRAGMA integrity_check').fetchone()[0]=='ok'


def test_parallel_and_shards_exact(tmp_path):
    serial=tmp_path/'serial';parallel=tmp_path/'parallel';serial.mkdir();parallel.mkdir()
    a,*_=build(serial,1);b,*_=build(parallel,2)
    for table in ('reference_metrics','reference_mate_audit'):
        with sqlite3.connect(a) as c,sqlite3.connect(b) as d:
            fields=[r[1] for r in c.execute(f'PRAGMA table_info({table})') if not r[1].endswith('_mtime_ns')]
            selected=','.join(fields)
            ar=c.execute(f'SELECT {selected} FROM {table} ORDER BY id').fetchall()
            br=d.execute(f'SELECT {selected} FROM {table} ORDER BY id').fetchall()
            if table=='reference_mate_audit':
                ar=[r[:3]+r[4:] for r in ar];br=[r[:3]+r[4:] for r in br]
            assert ar==br


def test_failure_and_missing_shards_protect_output(tmp_path):
    output,chunks,bundle,fna,fai=build(tmp_path)
    before=output.read_bytes()
    with pytest.raises(FileExistsError):merge(tmp_path/'records.csv',chunks,fna,fai,output,0,'geometry_only')
    assert output.read_bytes()==before
    for shards in [chunks[:1],[chunks[0],chunks[0]]]:
        with pytest.raises(ValueError):merge(tmp_path/'records.csv',shards,fna,fai,tmp_path/'bad.sqlite',0,'geometry_only')
        assert not (tmp_path/'bad.sqlite').exists()
    with pytest.raises(ValueError):merge(tmp_path/'records.csv',chunks,fna,fai,tmp_path/'bad.sqlite',1,'geometry_only')
    assert not (tmp_path/'bad.sqlite').exists()
    chunk=next(bundle.glob('*.json'))
    rows=json.loads(chunk.read_text());rows[0][1]='/missing.cram'
    broken=tmp_path/'broken.json';broken.write_text(json.dumps(rows))
    with pytest.raises(OSError):compute(broken,tmp_path/'bad.sqlite',1,0,'geometry_only')
    assert not (tmp_path/'bad.sqlite').exists()
