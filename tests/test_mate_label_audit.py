import csv
import sqlite3
from pathlib import Path

import pysam
import pytest

from bin.audit_cram_mate_labels import audit, inspect_record
from bin.build_cram_index import build_index


def write_bam(path, flags):
    with pysam.AlignmentFile(str(path), 'wb', header={'SQ': [{'SN': 'c', 'LN': 100}]}) as bam:
        for i, flag in enumerate(flags):
            read = pysam.AlignedSegment()
            read.query_name = f'r{i}'
            read.flag = flag
            read.reference_id = 0
            read.reference_start = i
            read.mapping_quality = 60
            read.cigarstring = '10M'
            read.query_sequence = 'A'*10
            bam.write(read)


@pytest.mark.parametrize('flags,status,complete', [
    ([65,129], 'BOTH_LABELS_PRESENT', False),
    ([0,16], 'NO_MATE_LABELS', True),
    ([65], 'ONLY_READ1', True), ([129], 'ONLY_READ2', True),
    ([1], 'INVALID_MATE_FLAGS', True),
    ([256], 'NO_ELIGIBLE_PRIMARY', True),
    ([2048,65,129], 'BOTH_LABELS_PRESENT', False)])
def test_mate_presence(tmp_path, flags, status, complete):
    path = tmp_path/'input.bam'
    write_bam(path, flags)
    result = inspect_record(dict(id=1, Sample='s', GCFID='g', CramFilePath=str(path)))
    assert result['status'] == status and result['scan_complete'] == complete


def test_range_determinism_missing_and_overwrite(tmp_path):
    source = tmp_path/'input.bam'
    write_bam(source, [0])
    csvfile = tmp_path/'records.csv'
    with csvfile.open('w',newline='') as h:
        writer=csv.writer(h)
        writer.writerow(['Sample','CramFilePath','GCFID','Coverage','Depth'])
        writer.writerow(['s',str(source),'g','1','1'])
        writer.writerow(['t',str(tmp_path/'missing.bam'),'g','1','1'])
    index=tmp_path/'index.sqlite'
    build_index(csvfile,index)
    first,second=tmp_path/'a.tsv',tmp_path/'b.tsv'
    assert audit(index,first,1,2,1,0)['errors']==1
    assert audit(index,second,1,2,2,0)['errors']==1
    assert first.read_bytes()==second.read_bytes()
    with pytest.raises(FileExistsError): audit(index,first,1,2,1,0)


def test_filtered_index_preserves_original_ids_and_all_audit_rows(tmp_path):
    source = tmp_path/'records.csv'
    paths = [tmp_path/f'{i}.bam' for i in range(3)]
    for path,flags in zip(paths, [[0],[65,129],[65]]):write_bam(path,flags)
    with source.open('w',newline='') as h:
        writer=csv.writer(h);writer.writerow(['Sample','CramFilePath','GCFID','Coverage','Depth'])
        writer.writerows([[f's{i}',str(p),'g','1','1'] for i,p in enumerate(paths)])
    original=tmp_path/'original.sqlite';build_index(source,original)
    inventory=tmp_path/'audit.tsv';audit(original,inventory,1,3,2,0)
    filtered=tmp_path/'filtered.sqlite'
    metadata=build_index(source,filtered,mate_audit=inventory)
    assert metadata['records']==2 and metadata['source_records']==3
    assert metadata['excluded_mate_records']==1
    with sqlite3.connect(filtered) as db:
        assert db.execute('SELECT id FROM reference_records ORDER BY id').fetchall()==[(2,),(3,)]
        assert db.execute('SELECT id,included FROM reference_mate_audit ORDER BY id').fetchall()==[(1,0),(2,1),(3,1)]
    changed=inventory.read_text().replace('s0','different')
    inventory.write_text(changed)
    invalid=tmp_path/'bad.sqlite'
    with pytest.raises(ValueError,match='identity differs'):
        build_index(source,invalid,mate_audit=inventory)
    assert not invalid.exists()


def test_incomplete_or_error_audit_cannot_silently_filter(tmp_path):
    source=tmp_path/'records.csv'
    with source.open('w',newline='') as h:
        writer=csv.writer(h);writer.writerow(['Sample','CramFilePath','GCFID','Coverage','Depth'])
        writer.writerow(['s',str(tmp_path/'missing.bam'),'g','1','1'])
    original=tmp_path/'original.sqlite';build_index(source,original)
    inventory=tmp_path/'audit.tsv';audit(original,inventory,1,1,1,0)
    output=tmp_path/'filtered.sqlite'
    with pytest.raises(ValueError,match='input errors'):
        build_index(source,output,mate_audit=inventory)
    assert not output.exists()
