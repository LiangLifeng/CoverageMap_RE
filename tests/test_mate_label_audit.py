import csv
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
