import json
from pathlib import Path
import sys
from types import SimpleNamespace

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from bin.alignment_geometry import (capture, Geometry, query_interval, SAME_PAIR,
    CROSS_PAIR, SPLIT_CROSS, BOTH_UNMAPPED, SINGLE_MAPPED, MISSING_RG, INCOMPLETE)
from bin.coverage_core import Reference, extract_pool, window_edges, window_sums, runtime_pool, select_units
from bin.native_reference import covered_sums
from bin.reference_store import build_bundle, load_reference_bundle


class Read:
    def __init__(self, name, flag, ci=-1, start=-1, cigar=None, mate=-1, mate_start=-1, rg='g'):
        self.query_name, self.flag, self.reference_id = name,flag,ci
        self.reference_name = ('c1','c2')[ci] if ci >= 0 else None
        self.reference_start, self.cigartuples = start,cigar
        self.reference_end = start+sum(n for op,n in cigar or [] if op in (0,2,3,7,8)) if ci>=0 else None
        self.next_reference_id, self.next_reference_start = mate,mate_start
        self.next_reference_name = ('c1','c2')[mate] if mate>=0 else None
        self.cigarstring = ''.join(str(n)+'MIDNSHP=X'[op] for op,n in cigar or [])
        self.mapping_quality = 42
        self.tags = {'RG':rg,'SA':'c2,50,+,10S10M,42,0;','AS':20,'NM':0} if rg is not None else {}
        for key, bit in (('is_paired',1),('is_unmapped',4),('mate_is_unmapped',8),('is_reverse',16),
                         ('mate_is_reverse',32),('is_read1',64),('is_read2',128),('is_secondary',256),
                         ('is_qcfail',512),('is_duplicate',1024),('is_supplementary',2048)):
            setattr(self,key,bool(flag & bit))
    def has_tag(self,key):
        return key in self.tags
    def get_tag(self,key):
        return self.tags[key]
    def get_blocks(self):
        position=self.reference_start
        blocks=[]
        for op,n in self.cigartuples or []:
            if op in (0,7,8):
                blocks.append((position,position+n))
            if op in (0,2,3,7,8):
                position+=n
        return blocks


@pytest.fixture
def alignments(tmp_path,monkeypatch):
    reads = [
        Read('same',65,0,2,[(0,20)],0,30), Read('same',129,0,30,[(0,20)],0,2),
        Read('cross',65,0,6,[(0,10),(4,10)],1,40),
        Read('cross',65|2048,1,20,[(4,10),(0,10)],1,40),
        Read('cross',129,1,40,[(0,20)],0,6),
        Read('cross',65|256,1,60,[(0,20)],1,40),
        Read('single',65|8,0,80,[(0,10)],-1,-1), Read('single',129|4,-1,-1,None,0,80),
        Read('unmapped',65|4|8),Read('unmapped',129|4|8),
        Read('orphan',65,1,85,[(0,10)],0,5,rg=None),
        Read('unpaired',0,0,55,[(0,5),(2,3),(0,5)]),
    ]
    class Alignment:
        header=SimpleNamespace(to_dict=lambda:{'SQ':[{'SN':'c1','LN':100},{'SN':'c2','LN':100}]})
        def __enter__(self):return self
        def __exit__(self,*args):pass
        def fetch(self,until_eof):return iter(reads)
    monkeypatch.setitem(sys.modules,'pysam',SimpleNamespace(AlignmentFile=lambda *a,**k:Alignment()))
    bam=tmp_path/'reads.bam';bam.write_bytes(b'BAM synthetic fixture')
    fai=tmp_path/'ref.fai';fai.write_text('c1\t100\nc2\t100\n')
    return bam,fai,reads


def test_cram_capture_preserves_used_fields_without_tlen(alignments,tmp_path,monkeypatch):
    bam,_,_=alignments
    bam.write_bytes(b'CRAM synthetic fixture')
    factory=sys.modules['pysam'].AlignmentFile
    observed=[]
    def open_alignment(*args,**kwargs):
        observed.append(kwargs)
        return factory(*args,**kwargs)
    monkeypatch.setattr(sys.modules['pysam'],'AlignmentFile',open_alignment)
    geometry,_=capture(bam,Reference(('c1','c2'),(100,100),('','')),tmp_path/'capture')
    mask=int(observed[0]['format_options'][0].split(b'=')[1])
    assert mask & 0x100 == 0
    assert mask & 0x600 == 0
    assert mask & 0x18ff == 0x18ff
    assert b'decode_md=0' in observed[0]['format_options']
    assert geometry.project('fragment').total_bases==geometry.project('read').total_bases


def test_scratch_location_does_not_change_geometry(alignments,tmp_path,monkeypatch):
    bam,_,_=alignments
    reference=Reference(('c1','c2'),(100,100),('',''))
    scratch=tmp_path/'scratch';scratch.mkdir()
    monkeypatch.setenv('COVERAGEMAP_SPOOL_DIR',str(scratch))
    first,_=capture(bam,reference,tmp_path/'first')
    assert not list(scratch.iterdir())
    other=tmp_path/'other-scratch';other.mkdir()
    monkeypatch.setenv('COVERAGEMAP_SPOOL_DIR',str(other))
    second,_=capture(bam,reference,tmp_path/'second')
    assert not list(other.iterdir())
    assert {p.name:p.read_bytes() for p in first.directory.iterdir()} == {
        p.name:p.read_bytes() for p in second.directory.iterdir()}


def test_capture_rejects_missing_scratch_before_output(alignments,tmp_path,monkeypatch):
    bam,_,_=alignments
    monkeypatch.setenv('COVERAGEMAP_SPOOL_DIR',str(tmp_path/'missing'))
    with pytest.raises(ValueError,match='existing absolute directory'):
        capture(bam,Reference(('c1','c2'),(100,100),('','')),tmp_path/'capture')
    assert not (tmp_path/'capture').exists()


def test_capture_scratch_reserve_and_error_cleanup(alignments,tmp_path,monkeypatch):
    from bin import alignment_geometry as module
    bam,_,_=alignments
    scratch=tmp_path/'scratch';scratch.mkdir()
    monkeypatch.setenv('COVERAGEMAP_SPOOL_DIR',str(scratch))
    monkeypatch.setenv('COVERAGEMAP_SPOOL_MIN_FREE_MIB','1')
    with monkeypatch.context() as scope:
        scope.setattr(module.shutil,'disk_usage',lambda p:SimpleNamespace(free=0))
        with pytest.raises(OSError,match='Capture scratch'):
            capture(bam,Reference(('c1','c2'),(100,100),('','')),tmp_path/'no-space')
    assert not (tmp_path/'no-space').exists()
    def fail_export(*args):
        raise OSError('injected export failure')
    monkeypatch.setattr(module,'_export',fail_export)
    with pytest.raises(OSError,match='injected export failure'):
        capture(bam,Reference(('c1','c2'),(100,100),('','')),tmp_path/'failed')
    assert not list(scratch.iterdir())


def test_capture_primary_mass_and_dynamic_units(alignments,tmp_path):
    bam,fai,_=alignments
    old,layout=extract_pool(bam,reference_index=fai)
    geometry,_=capture(bam,layout,tmp_path/'capture')
    reads=geometry.project('read')
    fragments=geometry.project('fragment')
    assert len(geometry.segments)==12 and len(geometry.reads)==10
    assert reads.total_bases==fragments.total_bases==old.total_bases
    edges=window_edges(200,13)
    for pool in (reads,fragments):
        np.testing.assert_array_equal(window_sums(pool,np.arange(len(pool)),edges),window_sums(old,np.arange(len(old)),edges))
        np.testing.assert_array_equal(covered_sums(pool,edges),covered_sums(old,edges))
    assert len(reads)==7 and len(fragments)==5
    assert reads.unit_bases[:2].tolist()==[20,20]
    assert fragments.unit_bases[0]==40
    assert geometry.fragments['flags'][0] & SAME_PAIR
    assert geometry.fragments['flags'][1] & CROSS_PAIR
    assert geometry.fragments['flags'][1] & SPLIT_CROSS
    assert geometry.fragments['flags'][2] & SINGLE_MAPPED
    assert geometry.fragments['flags'][3] & BOTH_UNMAPPED
    assert geometry.fragments['flags'][4] & MISSING_RG
    assert geometry.fragments['flags'][4] & INCOMPLETE
    full=geometry.project('read','primary_plus_supplementary')
    assert len(full)==len(reads) and full.total_bases==reads.total_bases+10
    selected=select_units(fragments,np.random.default_rng(1),target_bases=1)
    assert len(selected)==1 and fragments.unit_bases[selected].sum()>=1
    assert len(geometry.string(int(geometry.segments['cigar'][0])))>0


def test_streamed_audit_preserves_roles_across_update_batches(alignments,tmp_path):
    bam,_,reads=alignments
    reads[:]=[r for i in range(6001) for r in (
        Read(f'pair-{i}',65,0,2,[(0,20)],0,30),
        Read(f'pair-{i}',129,0,30,[(0,20)],0,2))]
    # Exercise mate mismatch on both sides of a 10,000-read update boundary.
    reads[9999].next_reference_start=7
    reads[10000].next_reference_start=8
    geometry,_=capture(bam,Reference(('c1','c2'),(100,100),('','')),tmp_path/'capture')
    from bin.alignment_geometry import MATE_DISAGREEMENT
    assert len(geometry.fragments)==6001 and len(geometry.reads)==12002
    assert np.all(geometry.fragments['flags'] & SAME_PAIR)
    assert np.flatnonzero(geometry.fragments['flags'] & MATE_DISAGREEMENT).tolist()==[4999,5000]
    assert np.flatnonzero(geometry.reads['audit_flags'] & MATE_DISAGREEMENT).tolist()==[9999,10000]
    assert geometry.project('read').total_bases==geometry.project('fragment').total_bases==240040


def test_string_ids_survive_bounded_cache_eviction():
    import sqlite3
    from bin.alignment_geometry import _string_interner
    with sqlite3.connect(':memory:') as db:
        db.execute('CREATE TABLE strings(id INTEGER PRIMARY KEY,value TEXT UNIQUE NOT NULL)')
        intern=_string_interner(db,budget_bytes=1024,capacity=2)
        assert [intern(f'tag-{i}') for i in range(20)]==list(range(20))
        assert intern('tag-0')==0 and intern('tag-19')==19
        long_tag='S'*2048
        assert intern(long_tag)==20 and intern(long_tag)==20
        assert intern('after')==21
        assert db.execute('SELECT count(*) FROM strings').fetchone()[0]==22


@pytest.mark.parametrize('reverse',[False,True])
def test_query_interval_clips_and_gaps(reverse):
    cigar=[(5,3),(4,5),(0,10),(1,2),(2,4),(3,7),(0,8),(4,4),(5,6)]
    assert query_interval(cigar,reverse)==((10,30,38) if reverse else (8,28,38))


def test_overlap_is_saved_but_explicit_coverage_fails(alignments,tmp_path):
    bam,_,reads=alignments
    reads[3].cigartuples=[(4,5),(0,10),(4,5)]
    geometry,_=capture(bam,Reference(('c1','c2'),(100,100),('','')),tmp_path/'capture')
    assert geometry.project().total_bases>0
    with pytest.raises(ValueError,match='Overlapping'):
        geometry.project('read','primary_plus_supplementary')


def test_bundle_roundtrip_capability_and_checksum(alignments,tmp_path):
    bam,_,_=alignments
    geometry,layout=capture(bam,Reference(('c1','c2'),(100,100),('','')),tmp_path/'capture')
    pool=geometry.project()
    policies=dict(mate_policy='both',min_mapq=0,exclude_duplicates=False,exclude_qcfail=False,
                  coverage_alignment_policy='primary_only')
    target=tmp_path/'bundle'
    build_bundle(pool,layout,target,'s',policies,[],workers=1)
    replay,_,manifest=load_reference_bundle(target,workers=1)
    fragment=runtime_pool(replay,policies,dict(policies,sampling_unit='fragment'))
    assert fragment.total_bases==pool.total_bases and len(fragment)<len(pool)
    assert manifest['capabilities']['supplementary_available']
    path=target/'alignment/segments.npy'
    with path.open('r+b') as handle:
        handle.seek(-1,2);handle.write(b'!')
    with pytest.raises(ValueError,match='Checksum'):
        load_reference_bundle(target,workers=1)


def test_legacy_capabilities_are_not_guessed(alignments):
    bam,fai,_=alignments
    pool,_=extract_pool(bam,reference_index=fai)
    policies=dict(mate_policy='both',sampling_unit='read',min_mapq=0,exclude_duplicates=False,exclude_qcfail=False)
    assert runtime_pool(pool,policies,policies) is pool
    with pytest.raises(ValueError,match='Legacy'):
        runtime_pool(pool,policies,dict(policies,sampling_unit='fragment'))


def test_complete_policies_work_with_legacy_concordance_adapter(alignments,tmp_path):
    from bin.concordance_inputs import resolve_policies
    bam,fai,_=alignments
    pool,_=extract_pool(bam,reference_index=fai,capture_geometry=tmp_path/'geometry')
    saved=dict(mate_policy='both',min_mapq=0,exclude_duplicates=False,exclude_qcfail=False,
               coverage_alignment_policy='primary_only')
    resolved=resolve_policies(SimpleNamespace(sampling_unit='fragment'),saved)
    assert resolved['sampling_unit']=='fragment'
    assert runtime_pool(pool,saved,resolved).total_bases==pool.total_bases
    legacy={k:v for k,v in dict(saved,sampling_unit='read').items() if k!='coverage_alignment_policy'}
    with pytest.raises(ValueError,match='conflicts'):
        resolve_policies(SimpleNamespace(sampling_unit='fragment'),legacy)


def test_empty_capture_and_fasta_header_strictness(alignments,tmp_path):
    bam,fai,reads=alignments
    with pytest.raises(ValueError,match='FASTA mode'):
        capture(bam,Reference(('c1','c2','extra'),(100,100,10),('','','')),tmp_path/'strict',reference_path='reference.fa')
    reads.clear()
    pool,layout=extract_pool(bam,reference_index=fai,capture_geometry=tmp_path/'empty')
    build_bundle(pool,layout,tmp_path/'empty-bundle','s',dict(mate_policy='both'),[],workers=1)
    loaded,_,_=load_reference_bundle(tmp_path/'empty-bundle',workers=1)
    assert len(loaded)==0 and loaded.total_bases==0


def test_query_fragment_projection_and_reference_draw_keep_observed_mates(alignments, tmp_path):
    from bin.query_manifest import extract_query, sample_reference_windows
    from bin.alignment_geometry import selection_audit
    from bin.stage_timing import StageTiming

    bam, fai, _ = alignments
    layout = Reference(('c1', 'c2'), (100, 100), ('', ''))
    policies = dict(mate_policy='both', sampling_unit='fragment',
                    coverage_alignment_policy='primary_only',
                    min_mapq=0, exclude_duplicates=False, exclude_qcfail=False)
    output = tmp_path/'query';output.mkdir()
    pool, got, inputs = extract_query(
        dict(CramPath=str(bam)), layout, output, policies, StageTiming())
    assert got.length == 200 and inputs['CramPath']['sha256']
    assert pool.sampling_unit == 'fragment'
    assert len(pool) == 5  # same, cross, single, orphan, unpaired
    first = selection_audit(pool, np.array([0]), target_bases=20)
    assert first['selected_units'] == first['selected_fragments'] == 1
    assert first['selected_reads'] == first['selected_segments'] == 2
    assert first['actual_bases'] == 40 and first['overshoot_bases'] == 20

    ref_pool, _ = extract_pool(bam, reference_index=fai,
                               capture_geometry=tmp_path/'reference-capture')
    fragment = runtime_pool(ref_pool, dict(policies, sampling_unit='read'), policies)
    edges = np.arange(0, 201, dtype=np.int64)  # 1 bp windows
    np.testing.assert_array_equal(
        window_sums(fragment, np.arange(len(fragment)), edges),
        window_sums(ref_pool, np.arange(len(ref_pool)), edges))
    audit={}
    sampled, count, selected_sha = sample_reference_windows(
        fragment, np.random.default_rng(17), 1, edges, audit=audit)
    assert count == audit['selected_units'] == audit['selected_fragments'] == 1
    assert audit['selected_reads'] in (1, 2)
    assert audit['actual_bases'] == sampled.sum() and selected_sha
    assert audit['coverage_alignment_policy'] == 'primary_only'


def test_contig_rebase_preserves_fragment_identity_and_audit_handles():
    from bin.coverage_core import BlockPool, canonical_reference
    from bin.query_manifest import align_reference_geometry

    source = Reference(('c2','c1'), (100,100), ('',''))
    target = canonical_reference(source)
    pool = BlockPool(np.array([2,110]),np.array([8,119]),np.array([0,2]),200)
    handle = SimpleNamespace(nbytes=1)
    pool.alignment_geometry = handle
    pool.geometry_unit_ids = np.array([9])
    pool.sampling_unit = 'fragment'
    pool.coverage_alignment_policy = 'primary_only'
    rebased, layout = align_reference_geometry(pool, source, target)
    assert layout.names == target.names
    np.testing.assert_array_equal(rebased.starts,[102,10])
    np.testing.assert_array_equal(rebased.ends,[108,19])
    assert rebased.alignment_geometry is handle
    np.testing.assert_array_equal(rebased.geometry_unit_ids,[9])
    assert rebased.sampling_unit == 'fragment'
    assert rebased.coverage_alignment_policy == 'primary_only'
