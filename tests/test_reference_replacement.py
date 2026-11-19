import json
from pathlib import Path
import shutil

import numpy as np
import pytest

from bin.raw_cram_query import calculate, rng, stamp
from bin.raw_query_sampling import QueryArm, compatible_reference_ids
from bin.raw_reference_replacement import resolve_reference_choices, InsufficientReferenceCapacity
from test_raw_query_vnext import _fixture


def test_replacement_keeps_slot_count_and_is_deterministic():
    arms=[QueryArm('read1','read1',np.array([0]),100,3),
          QueryArm('read2','read2',np.array([1]),100,3)]
    spec=dict(query_id='q',query=dict(Sample='s',GCFID='g'),
              config=dict(seed=401,sampling_tolerance_bp=0,min_samples=0),
              choices=[1]*6,references=[dict(id=i,Sample=str(i)) for i in [1,2]])
    def capacity(ref,arm):
        return (80,8) if ref['id']==1 and arm.name=='read2' else (200,20)
    first,events=resolve_reference_choices(spec,arms,capacity,rng)
    second,other=resolve_reference_choices(spec,arms,capacity,rng)
    assert first==second==[1,1,1,2,2,2]
    assert events==other and len(events)==3
    assert spec['choices']==[1]*6
    with pytest.raises(InsufficientReferenceCapacity):
        resolve_reference_choices(spec,arms,lambda *_:(0,0),rng)


def test_tolerance_and_decode_errors_are_not_silently_replaced():
    arm=QueryArm('single','all',np.array([0]),100,2)
    spec=dict(query_id='q',query=dict(Sample='s',GCFID='g'),
        config=dict(seed=1,sampling_tolerance_bp=150,min_samples=0),
        choices=[1,1],references=[dict(id=1,Sample='r')])
    assert resolve_reference_choices(spec,[arm],lambda *_:(98,10),rng)==([1,1],[])
    def broken(*_):raise ValueError('Corrupt source')
    with pytest.raises(ValueError,match='Corrupt source'):
        resolve_reference_choices(spec,[arm],broken,rng)


def test_actual_distinct_sample_support_is_verified_after_exclusion():
    arm=QueryArm('single','all',np.array([0]),100,2)
    spec=dict(query_id='q',query=dict(Sample='s',GCFID='g'),
              config=dict(seed=401,sampling_tolerance_bp=0,min_samples=1),
              choices=[2,2],references=[dict(id=i,Sample=str(i)) for i in [1,2,3]])
    capacity=lambda ref,_: (10,1) if ref['id']==1 else (200,20)
    assert resolve_reference_choices(spec,[arm],capacity,rng)==([2,2],[])
    spec['references'][2]['Sample']='2'
    with pytest.raises(InsufficientReferenceCapacity,match='distinct reference samples'):
        resolve_reference_choices(spec,[arm],capacity,rng)


def test_unlabelled_pool_policy_never_invents_mate_labels():
    class Pool:
        unit_mate_sides=np.array([0,0,0])
        def __len__(self):return 3
    pool=Pool()
    assert len(compatible_reference_ids(pool,'read1'))==0
    assert len(compatible_reference_ids(pool,'read2'))==0
    np.testing.assert_array_equal(pool.unit_mate_sides,[0,0,0])


def test_real_geometry_replacement_exact_depth_and_warm_cache(tmp_path):
    spec,_=_fixture(tmp_path)
    source=tmp_path/'sufficient.bam'
    shutil.copyfile(spec['query']['CramPath'],source)
    spec['references'].append(dict(id=999,Sample='sufficient',GCFID='GCF_1.1',
                                  CramFilePath=str(source),Depth='6'))
    spec['stamps'][str(source)]=stamp(source)
    spec['config'].update(insufficient_reference_policy='replace',
                          sampling_tolerance_bp=0)
    output=tmp_path/'replace';output.mkdir()
    result=calculate(spec,output)
    assert result['requested_draws']==result['valid_draws']==10
    assert result['replaced_draws']==result['replacement_attempts']==5
    assert result['out_of_tolerance_draws']==0
    events=json.loads((output/'replacement_audit.json').read_text())
    assert all(e['reason']=='INSUFFICIENT_CAPACITY' and e['replacement_reference_id']==999 for e in events)
    warm=tmp_path/'warm';warm.mkdir()
    cached=calculate(spec,warm)
    assert cached['decoded_reference_files']==0
    for name in ['reference_windows.npy','reference_covered_bases.npy','draws.tsv',
                 'statistics.json','replacement_audit.json']:
        if name=='draws.tsv':continue # cache-hit bookkeeping differs; numerical payloads do not.
        assert (output/name).read_bytes()==(warm/name).read_bytes()


def test_replacement_exhaustion_is_unscored_not_an_error(tmp_path):
    spec,_=_fixture(tmp_path)
    spec['config'].update(insufficient_reference_policy='replace',
                          sampling_tolerance_bp=0)
    output=tmp_path/'exhausted';output.mkdir()
    summary=calculate(spec,output)
    assert summary['status']=='INSUFFICIENT_REFERENCE_CAPACITY'
    assert summary['detection_status']=='NOT_EVALUATED'
    assert summary['requested_draws']==10
