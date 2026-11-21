import json
import os
from pathlib import Path
import sqlite3
from types import SimpleNamespace

from bin.indexed_reference_metrics import bind_metrics, capacities
from test_computed_cram_index import build


def test_sql_capacity_and_policy_or_stale_fallback(tmp_path):
    index,*_=build(tmp_path)
    config=dict(min_mapq=0,m5_policy='geometry_only',coverage_alignment_policy='primary_only',sampling_unit='read')
    with sqlite3.connect(index) as db:
        db.row_factory=sqlite3.Row
        refs=[dict(r) for r in db.execute('SELECT * FROM reference_records')]
        bind_metrics(db,refs,config)
        assert len(refs)==2 and all(r.get('sql_capacity') for r in refs)
        spec=dict(config=config,fna=str(tmp_path/'GCF_1.1.fna'),fai=str(tmp_path/'GCF_1.1.fna.fai'))
        arms=[SimpleNamespace(name='read1',mate_side='read1'),SimpleNamespace(name='read2',mate_side='read2')]
        assert capacities(refs[0],arms,spec)=={'read1':(30,3),'read2':(20,3)}
        pair=[SimpleNamespace(name='paired',mate_side='both')]
        spec['config']=dict(config,sampling_unit='fragment')
        assert capacities(refs[0],pair,spec)=={'paired':(57,5)}
        changed=[dict(r) for r in db.execute('SELECT * FROM reference_records')]
        bind_metrics(db,changed,dict(config,min_mapq=1))
        assert all('sql_capacity' not in r for r in changed)
        bam=Path(refs[0]['CramFilePath']);st=bam.stat();os.utime(bam,ns=(st.st_atime_ns,st.st_mtime_ns+1_000_000))
        assert capacities(refs[0],arms,spec) is None


def test_old_index_uses_existing_decode_fallback(tmp_path):
    with sqlite3.connect(':memory:') as db:
        db.execute('CREATE TABLE metadata(key TEXT,value TEXT)')
        db.execute('INSERT INTO metadata VALUES (?,?)',('schema',json.dumps('coverage-raw-cram-index-1')))
        refs=[dict(id=1,GCFID='g')]
        bind_metrics(db,refs,{})
        assert 'sql_capacity' not in refs[0]
        assert capacities(refs[0],[],{}) is None


def test_sql_preflight_reduces_decodes_without_changing_science(tmp_path):
    from collections import Counter
    from copy import deepcopy
    import shutil
    import pysam
    from bin.cram_index_metrics import measure
    from bin.indexed_reference_metrics import FIELDS
    from bin.raw_cram_pool import extract_once
    from bin.raw_cram_query import calculate,stamp
    from test_raw_query_vnext import _fixture
    spec,_=_fixture(tmp_path)
    insufficient=Path(spec['references'][0]['CramFilePath'])
    modified=tmp_path/'modified.bam'
    with pysam.AlignmentFile(insufficient,'rb') as src,pysam.AlignmentFile(modified,'wb',header=src.header) as out:
        for read in src:
            read.query_sequence='A'*8;read.cigarstring='8M';out.write(read)
    modified.replace(insufficient)
    spec['stamps'][str(insufficient)]=stamp(insufficient)
    sufficient=tmp_path/'sufficient.bam';shutil.copyfile(spec['query']['CramPath'],sufficient)
    spec['references'].append(dict(id=999,Sample='sufficient',GCFID='GCF_1.1',CramFilePath=str(sufficient),Depth='6'))
    spec['stamps'][str(sufficient)]=stamp(sufficient)
    spec['config'].update(insufficient_reference_policy='replace',sampling_tolerance_bp=0)
    original=deepcopy(spec);optimized=deepcopy(spec)
    for ref in optimized['references']:
        _,metrics=measure((dict(id=ref['id'],row=ref,fna=spec['fna'],fai=spec['fai']),0,'strict'))
        ref['sql_capacity']={key:metrics[key] for key in FIELDS}
    calls=[]
    for label,data in [('original',original),('optimized',optimized)]:
        data['config']['geometry_cache_dir']=str(tmp_path/(label+'-cache'))
        count=Counter()
        def decode(path,layout,**options):
            count[str(path)]+=1
            return extract_once(path,layout,**options)
        output=tmp_path/label;output.mkdir()
        result=calculate(data,output,extractor=decode)
        assert result['valid_draws']==10 and result['replaced_draws']==10
        calls.append(count)
    assert calls[0][str(insufficient)]==1 and calls[1][str(insufficient)]==0
    assert calls[0][str(sufficient)]==calls[1][str(sufficient)]==1
    for name in ['reference_windows.npy','reference_covered_bases.npy','statistics.json','replacement_audit.json']:
        assert (tmp_path/'original'/name).read_bytes()==(tmp_path/'optimized'/name).read_bytes()
    audit=json.loads((tmp_path/'optimized/replacement_capacity_audit.json').read_text())
    assert audit and all(r['capacity_source']=='sqlite' for r in audit)


def test_sql_capacity_does_not_bypass_self_source_guard(tmp_path):
    import pytest
    from bin.cram_index_metrics import measure
    from bin.indexed_reference_metrics import FIELDS
    from bin.raw_cram_query import calculate,stamp
    from test_raw_query_vnext import _fixture
    spec,_=_fixture(tmp_path)
    alias=tmp_path/'hardlink-query.bam';os.link(spec['query']['CramPath'],alias)
    ref=spec['references'][0];ref['CramFilePath']=str(alias)
    spec['stamps'][str(alias)]=stamp(alias)
    _,metric=measure((dict(id=ref['id'],row=ref,fna=spec['fna'],fai=spec['fai']),0,'strict'))
    ref['sql_capacity']={key:metric[key] for key in FIELDS}
    spec['config']['insufficient_reference_policy']='replace'
    output=tmp_path/'guard';output.mkdir()
    with pytest.raises(ValueError,match='aliases query source'):
        calculate(spec,output)
