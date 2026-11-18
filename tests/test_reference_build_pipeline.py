import argparse
import json
from pathlib import Path
import shutil
import sys

import pytest

sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from bin import ref_build_pipeline as pipeline
from bin.ref_depth_catalog import atomic_json, digest
from bin.reference_index import ReferenceIndex
from test_alignment_geometry import alignments


def prepare(tmp_path,alignments):
    bam,fai,_=alignments
    records=tmp_path/'records.csv'
    rows=[]
    for i in range(4):
        source=tmp_path/f'source{i}.bam'
        source.write_bytes(bam.read_bytes()+str(i).encode())
        rows.append(f's{i},{source},GCF_000000001.1\n')
    records.write_text('Sample,CramFilePath,GCFID\n'+''.join(rows))
    refs=tmp_path/'refs.csv'
    refs.write_text(f'GCF_ID,GCF_path\nGCF_000000001.1,{fai}\n')
    args=argparse.Namespace(records=records,references=refs,outdir=tmp_path/'plan',chunk_size=2,
        cpus=1,smoke_records=0,window_size=20,bin_count=None)
    pipeline.plan_stage(args, eager_sources=True)
    return args


def build_chunks(tmp_path,args):
    receipt_dir=tmp_path/'receipts';receipt_dir.mkdir()
    for path in sorted((args.outdir/'chunks').glob('*.json')):
        plan=json.loads(path.read_text())
        results=[pipeline.build_record((plan,g,s,str(tmp_path/'checkpoints')))
            for g,group in plan['groups'].items() for s in group['samples']]
        assert all(r['state']=='SUCCESS' for r in results),results
        atomic_json(receipt_dir/path.name,dict(chunk_id=path.stem,chunk_sha256=digest(path),
            records=results,completed=len(results),failed=0,skipped=0,state='SUCCESS'))
    return receipt_dir


def test_global_plan_conservation_and_scoped_input_change(tmp_path,alignments):
    args=prepare(tmp_path,alignments)
    original={p.name:json.loads(p.read_text()) for p in (args.outdir/'chunks').glob('*.json')}
    assert len(original)==2
    keys=[(g,s['sample']) for chunk in original.values() for g,group in chunk['groups'].items() for s in group['samples']]
    assert len(keys)==len(set(keys))==4
    # Replan on unchanged inputs is byte-stable even under different CPU allocation.
    args.outdir=tmp_path/'plan2';args.cpus=2
    pipeline.plan_stage(args, eager_sources=True)
    assert {p.name:json.loads(p.read_text()) for p in (args.outdir/'chunks').glob('*.json')}==original
    source=tmp_path/'source0.bam';source.write_bytes(source.read_bytes()+b'changed')
    args.outdir=tmp_path/'plan3'
    pipeline.plan_stage(args, eager_sources=True)
    updated={p.name:json.loads(p.read_text()) for p in (args.outdir/'chunks').glob('*.json')}
    assert updated['chunk-000000.json']!=original['chunk-000000.json']
    assert updated['chunk-000001.json']==original['chunk-000001.json']


def test_record_reuse_collect_single_writer_and_index(tmp_path,alignments):
    args=prepare(tmp_path,alignments)
    receipts=build_chunks(tmp_path,args)
    chunk=json.loads((args.outdir/'chunks/chunk-000000.json').read_text())
    gcf=next(iter(chunk['groups']));sample=chunk['groups'][gcf]['samples'][0]
    replay=pipeline.build_record((chunk,gcf,sample,str(tmp_path/'checkpoints')))
    assert replay['reused'] and replay['state']=='SUCCESS'
    collect=argparse.Namespace(plan=args.outdir,receipts=receipts,database=tmp_path/'database',
        cpus=1,receipt=tmp_path/'collected.json')
    pipeline.collect_stage(collect)
    index_args=argparse.Namespace(database=collect.database,cpus=1,receipt=tmp_path/'indexed.json')
    pipeline.index_stage(index_args)
    with_index=ReferenceIndex(collect.database.parent/'reference.sqlite')
    try:
        assert len(with_index.candidates(gcf))==4
        assert with_index.metadata['geometry_capabilities']['sampling_units']==['read','fragment']
    finally:
        with_index.close()
    sha=digest(collect.database.parent/'reference.sqlite')
    pipeline.index_stage(index_args)
    assert digest(collect.database.parent/'reference.sqlite')==sha
    pipeline.collect_stage(collect)
    assert digest(collect.database.parent/'reference.sqlite')==sha


@pytest.mark.parametrize('damage',['missing','duplicate','failed'])
def test_collector_rejects_missing_duplicate_failed(tmp_path,alignments,damage):
    args=prepare(tmp_path,alignments)
    receipts=build_chunks(tmp_path,args)
    path=receipts/'chunk-000000.json'
    if damage=='missing':
        path.unlink()
    elif damage=='duplicate':
        shutil.copy2(path,receipts/'duplicate.json')
    else:
        record=json.loads(path.read_text());record.update(state='FAILED',failed=1)
        atomic_json(path,record)
    with pytest.raises(ValueError):
        pipeline.collect_stage(argparse.Namespace(plan=args.outdir,receipts=receipts,
            database=tmp_path/'database',cpus=1,receipt=tmp_path/'collect.json'))
    assert not (tmp_path/'database').exists()


def test_failed_record_retains_verified_successes(tmp_path,alignments):
    args=prepare(tmp_path,alignments)
    chunk=json.loads((args.outdir/'chunks/chunk-000000.json').read_text())
    gcf=next(iter(chunk['groups']));first,second=chunk['groups'][gcf]['samples']
    good=pipeline.build_record((chunk,gcf,first,str(tmp_path/'checkpoints')))
    assert good['state']=='SUCCESS'
    root=tmp_path/'checkpoints'/pipeline.fingerprint(pipeline.record_plan(chunk,gcf,second))
    root.mkdir(parents=True);(root/'.record.lock').mkdir()
    failed=pipeline.build_record((chunk,gcf,second,str(tmp_path/'checkpoints')))
    assert failed['state']=='FAILED' and 'Locked' in failed['error']
    (root/'.record.lock').rmdir()
    assert pipeline.build_record((chunk,gcf,first,str(tmp_path/'checkpoints')))['reused']
    assert pipeline.build_record((chunk,gcf,second,str(tmp_path/'checkpoints')))['state']=='SUCCESS'
