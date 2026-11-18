"""Full-input preprocessing, independent complete batches and exact final merge."""
import argparse
import csv
import json
import os
from pathlib import Path
import shutil
import sys

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from bin import ref_build_pipeline as pipeline
from bin.ref_build_batches import batches, publish_batch
from bin.ref_depth_catalog import atomic_json, digest
from bin.reference_index import ReferenceIndex
from test_alignment_geometry import alignments
from test_plan_parallel_validation import fixture as roster_fixture


def prepare(tmp_path, alignments, *, target=2, eager=False, out='plan'):
    bam, fai, _ = alignments
    records = tmp_path/'records.csv'
    rows = []
    for group, count in enumerate((3, 2, 1), 1):
        for number in range(count):
            path = tmp_path/f'g{group}s{number}.bam'
            path.write_bytes(bam.read_bytes()+bytes([group, number]))
            rows.append(dict(Sample=f's{number}', CramFilePath=str(path), GCFID=f'GCF_{group:09d}.1'))
    with records.open('w', newline='') as handle:
        writer=csv.DictWriter(handle, list(rows[0]));writer.writeheader();writer.writerows(reversed(rows))
    refs=tmp_path/'refs.csv'
    refs.write_text('GCF_ID,GCF_path\n'+''.join(f'GCF_{group:09d}.1,{fai}\n' for group in (1,2,3)))
    args=argparse.Namespace(records=records,references=refs,outdir=tmp_path/out,
        chunk_size=target,cpus=2,smoke_records=0,window_size=20,bin_count=None)
    pipeline.plan_stage(args,eager_sources=eager)
    return args


def build_manual(tmp_path,args, *, output='out'):
    root=tmp_path/output
    receipts=root/'task-receipts';receipts.mkdir(parents=True)
    for path in sorted((args.outdir/'chunks').glob('*.json')):
        plan=json.loads(path.read_text())
        batch_root=root/'batches'/path.stem
        rows=[pipeline.build_record((plan,g,s,str(batch_root/'checkpoints')))
              for g,group in plan['groups'].items() for s in group['samples']]
        assert all(row['state']=='SUCCESS' for row in rows), rows
        published=publish_batch(plan,rows,batch_root/'database',2)
        atomic_json(receipts/path.name,dict(chunk_id=path.stem,chunk_sha256=digest(path),records=rows,
            completed=len(rows),failed=0,skipped=0,state='SUCCESS',**published))
    return argparse.Namespace(plan=args.outdir,receipts=receipts,database=root/'database',
        cpus=2,receipt=root/'collected.json')


def index(collect):
    pipeline.collect_stage(collect)
    args=argparse.Namespace(database=collect.database,cpus=2,receipt=collect.database.parent/'indexed.json')
    pipeline.index_stage(args)
    return collect.database.parent/'reference.sqlite'


def contents(path):
    db=ReferenceIndex(path)
    try:
        result={}
        for row in db.connection.execute('SELECT * FROM references_table ORDER BY gcf,sample'):
            row=dict(row)
            manifest=db.root_for(row)/row['manifest_path']
            record=json.loads(manifest.read_text())
            arrays={}
            with np.load(manifest.parent/'native_coverage.npz',allow_pickle=False) as archive:
                for key in archive.files:
                    arrays[key]=archive[key].copy()
            result[row['gcf'],row['sample']]=(row['source_sha'],record['aligned_bases'],record['full_mean_depth'],arrays)
        return result,db.metadata
    finally:
        db.close()


def test_preprocessing_never_reads_full_source_and_keeps_gcfs(tmp_path,monkeypatch):
    records,refs,_=roster_fixture(tmp_path,12)
    original=pipeline.digest
    def no_source_hash(path):
        assert Path(path).suffix not in ('.cram','.bam','.dat')
        return original(path)
    monkeypatch.setattr(pipeline,'digest',no_source_hash)
    args=argparse.Namespace(records=records,references=refs,outdir=tmp_path/'plan',
        cpus=4,chunk_size=3,smoke_records=0,window_size=20,bin_count=None)
    pipeline.plan_stage(args)
    plan=json.loads((args.outdir/'database.json').read_text())
    chunks=list(batches(plan,3))
    assert len(chunks)==2  # two complete GCFs, even though target suggests four batches
    assert sum(c['counts']['sample_gcf_pairs'] for c in chunks)==12
    assert all(len(c['groups'])==1 and c['counts']['sample_gcf_pairs']==6 for c in chunks)
    assert all('source_sha256' not in s for group in plan['groups'].values() for s in group['samples'])
    report=json.loads((args.outdir/'planning.json').read_text())
    assert report['source_bytes_hashed']==0
    assert report['source_bytes_total']==sum(c['source_bytes'] for c in chunks)
    assert (args.outdir/'batches.tsv').is_file()
    args.outdir=tmp_path/'again';args.cpus=1
    pipeline.plan_stage(args)
    assert {p.name:p.read_bytes() for p in (args.outdir/'chunks').glob('*.json')} == {
        p.name:p.read_bytes() for p in (tmp_path/'plan/chunks').glob('*.json')}


def test_byte_balancing_not_equal_record_slicing():
    groups={f'GCF_{i:09d}.1':dict(samples=[dict(sample='s',source_stat=dict(bytes=size))])
            for i,size in enumerate((1000,900,100,90),1)}
    plan=dict(schema='test',config={},groups=groups,counts=dict(sample_gcf_pairs=4),scope=dict(mode='full'))
    values=list(batches(plan,2))
    assert sorted(c['source_bytes'] for c in values)==[1000,1090]
    assert len({g for c in values for g in c['groups']})==4
    reverse=dict(plan,groups=dict(reversed(list(groups.items()))))
    assert list(batches(reverse,2))==values


def test_batches_equal_eager_geometry_and_global_index(tmp_path,alignments):
    args=prepare(tmp_path,alignments)
    collect=build_manual(tmp_path,args)
    actual,meta=contents(index(collect))
    assert len(actual)==6 and meta['counts']['gcfs']==3
    assert len(meta['sources'])==3
    assert meta['geometry_capabilities']['sampling_units']==['read','fragment']
    assert meta['selected_counts']['references']==6
    # Aggregation retains references to batch payloads rather than another copy.
    assert not list(collect.database.glob('GCF_*'))
    assert len(json.loads((collect.database/'BATCHES.json').read_text()))==3
    eager=prepare(tmp_path,alignments,eager=True,out='eager-plan')
    from test_reference_build_pipeline import build_chunks
    receipts=build_chunks(tmp_path,eager)
    old=argparse.Namespace(plan=eager.outdir,receipts=receipts,database=tmp_path/'eager-out/database',
        cpus=2,receipt=tmp_path/'eager-out/collected.json')
    expected,_=contents(index(old))
    assert actual.keys()==expected.keys()
    for key in actual:
        assert actual[key][:3]==expected[key][:3]
        assert actual[key][3].keys()==expected[key][3].keys()
        for name,value in actual[key][3].items():
            np.testing.assert_array_equal(value,expected[key][3][name])
    before=digest(collect.database.parent/'reference.sqlite')
    index(collect)
    assert digest(collect.database.parent/'reference.sqlite')==before


def test_source_hash_bound_before_checkpoint_reuse(tmp_path,alignments):
    args=prepare(tmp_path,alignments)
    chunk=json.loads(next((args.outdir/'chunks').glob('*.json')).read_text())
    gcf=next(iter(chunk['groups']));sample=chunk['groups'][gcf]['samples'][0]
    task=(chunk,gcf,sample,str(tmp_path/'checkpoints'))
    first=pipeline.build_record(task)
    assert first['state']=='SUCCESS' and not first['reused']
    assert pipeline.build_record(task)['reused']
    source=Path(sample['cram']);before=source.stat()
    value=source.read_bytes();source.write_bytes(value[:-1]+bytes([value[-1]^1]))
    os.utime(source,ns=(before.st_atime_ns,before.st_mtime_ns))
    changed=pipeline.build_record(task)
    assert changed['state']=='SUCCESS' and not changed['reused']
    assert first['source_sha256']!=changed['source_sha256'] and first['key']!=changed['key']
    assert pipeline.build_record(task)['reused']


def test_changed_source_fails_inside_batch_and_preserves_good_checkpoint(tmp_path,alignments):
    args=prepare(tmp_path,alignments)
    chunk=next(value for path in sorted((args.outdir/'chunks').glob('*.json'))
        if (value:=json.loads(path.read_text()))['counts']['sample_gcf_pairs']>=2)
    gcf=next(iter(chunk['groups']));first,second=chunk['groups'][gcf]['samples'][:2]
    checkpoints=str(tmp_path/'checkpoints')
    good=pipeline.build_record((chunk,gcf,first,checkpoints))
    source=Path(second['cram']);source.write_bytes(source.read_bytes()+b'changed')
    bad=pipeline.build_record((chunk,gcf,second,checkpoints))
    assert good['state']=='SUCCESS' and bad['state']=='FAILED' and 'changed' in bad['error'].lower()
    assert pipeline.build_record((chunk,gcf,first,checkpoints))['reused']


@pytest.mark.parametrize('damage',['missing','duplicate','failed','foreign-path','source-sha','key','plan','summary','extra-planned'])
def test_global_merge_rejects_incomplete_or_foreign_batches(tmp_path,alignments,damage):
    args=prepare(tmp_path,alignments)
    collect=build_manual(tmp_path,args)
    path=next(collect.receipts.glob('*.json'))
    row=json.loads(path.read_text())
    if damage=='missing':path.unlink()
    elif damage=='duplicate':shutil.copy2(path,collect.receipts/'copy.json')
    elif damage=='extra-planned':shutil.copy2(next((args.outdir/'chunks').glob('*.json')),args.outdir/'chunks/extra.json')
    elif damage=='plan':
        target=Path(row['database'])/'database.json'
        value=json.loads(target.read_text());value['counts']['sample_gcf_pairs']+=1;atomic_json(target,value)
    elif damage=='summary':
        target=Path(row['database'])/'SUMMARY.json'
        value=json.loads(target.read_text());value['complete']=False;atomic_json(target,value)
    else:
        if damage=='failed':row.update(state='FAILED',failed=1)
        elif damage=='foreign-path':row['database']=str(tmp_path/'foreign')
        elif damage=='source-sha':row['records'][0]['source_sha256']='0'*64
        elif damage=='key':row['records'][0]['key']='0'*64
        atomic_json(path,row)
    with pytest.raises((ValueError,FileNotFoundError)):
        pipeline.collect_stage(collect)
    assert not collect.database.exists()
    assert not (collect.database.parent/'reference.sqlite').exists()


def test_final_index_rejects_payload_corruption(tmp_path,alignments):
    args=prepare(tmp_path,alignments)
    collect=build_manual(tmp_path,args)
    pipeline.collect_stage(collect)
    batch=next((collect.database.parent/'batches').iterdir())
    payload=next((batch/'database').glob('GCF_*/samples/*/attempt-*/native_coverage.npz'))
    payload.write_bytes(payload.read_bytes()+b'corrupted')
    with pytest.raises(ValueError):
        pipeline.index_stage(argparse.Namespace(database=collect.database,cpus=2,receipt=tmp_path/'indexed.json'))
    assert not (collect.database/'SUCCESS').exists()
    assert not (collect.database.parent/'reference.sqlite').exists()


def test_content_change_during_chunk_validation_is_rejected(tmp_path,alignments,monkeypatch):
    args=prepare(tmp_path,alignments)
    chunk=json.loads(next((args.outdir/'chunks').glob('*.json')).read_text())
    gcf=next(iter(chunk['groups']));sample=chunk['groups'][gcf]['samples'][0]
    source=Path(sample['cram'])
    original=pipeline.digest
    def mutate(path):
        value=original(path)
        if Path(path)==source:
            stat=source.stat();data=source.read_bytes()
            source.write_bytes(data[:-1]+bytes([data[-1]^1]))
            os.utime(source,ns=(stat.st_atime_ns,stat.st_mtime_ns))
        return value
    monkeypatch.setattr(pipeline,'digest',mutate)
    row=pipeline.build_record((chunk,gcf,sample,str(tmp_path/'checkpoints')))
    assert row['state']=='FAILED' and 'cram content changed' in row['error'].lower()


def test_reference_modules_pin_staged_code_and_disable_chunk_wrapper_cache():
    root=Path(__file__).resolve().parents[1]
    for name in ('plan','chunk','collect','index'):
        source=(root/f'modules/local/reference_{name}/main.nf').read_text()
        assert "path code, stageAs: 'reference_code/bin/*'" in source
        assert source.count('export PATH=')==2 and 'reference_code/bin:' in source
    config=(root/'conf/reference_build.config').read_text()
    block=config.split("withName: 'REFERENCE_CHUNK'",1)[1].split("withName:",1)[0]
    assert 'cache = false' in block
