#!/usr/bin/env python3
"""Generate portable synthetic alignment relationships for real executor tests."""
import argparse
import csv
from pathlib import Path
import sys

sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from test_alignment_geometry import Read


def write(path,repeats,prefix):
    import pysam
    header={'SQ':[{'SN':'c1','LN':100},{'SN':'c2','LN':100}], 'RG':[{'ID':'g','SM':'synthetic'}]}
    rows=[Read('same',65,0,2,[(0,20)],0,30),Read('same',129,0,30,[(0,20)],0,2),
          Read('cross',65,0,6,[(0,10),(4,10)],1,40),Read('cross',65|2048,1,20,[(4,10),(0,10)],1,40),
          Read('cross',129,1,40,[(0,20)],0,6),Read('cross',65|256,1,60,[(0,20)],1,40),
          Read('single',65|8,0,80,[(0,10)],-1,-1),Read('single',129|4,-1,-1,None,0,80),
          Read('unmapped',65|4|8),Read('unmapped',129|4|8),Read('orphan',65,1,85,[(0,10)],0,5,rg=None),
          Read('unpaired',0,0,55,[(0,5),(2,3),(0,5)])]
    with pysam.AlignmentFile(path,'wb',header=header) as handle:
        for repeat in range(repeats):
            for row in rows:
                read=pysam.AlignedSegment(handle.header)
                read.query_name=f'{prefix}-{repeat}-{row.query_name}'
                read.flag=row.flag;read.reference_id=row.reference_id;read.reference_start=row.reference_start
                read.cigartuples=row.cigartuples;read.mapping_quality=row.mapping_quality
                read.next_reference_id=row.next_reference_id;read.next_reference_start=row.next_reference_start
                read.query_sequence='A'*sum(n for op,n in row.cigartuples or [] if op in (0,1,4,7,8)) or None
                for key,value in row.tags.items():
                    read.set_tag(key,value)
                handle.write(read)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--outdir',type=Path,required=True)
    parser.add_argument('--repeats',type=int,default=1000)
    args=parser.parse_args()
    if args.repeats<1:parser.error('repeats must be positive')
    args.outdir.mkdir(parents=True,exist_ok=False)
    fai=args.outdir/'ref.fai';fai.write_text('c1\t100\nc2\t100\n')
    with (args.outdir/'records.csv').open('w',newline='') as stream:
        writer=csv.writer(stream);writer.writerow(['Sample','CramFilePath','GCFID'])
        for i in range(6):
            path=args.outdir/f's{i}.bam';write(path,args.repeats,f's{i}')
            writer.writerow([f's{i}',str(path),f'GCF_00000000{1+i//3}.1'])
    (args.outdir/'references.csv').write_text('GCF_ID,GCF_path\n'+''.join(f'GCF_00000000{i}.1,{fai}\n' for i in (1,2)))
    write(args.outdir/'query-low.bam',1,'query-low')
    write(args.outdir/'query-high.bam',12,'query-high')
    (args.outdir/'queries.tsv').write_text('Sample\tGCFID\tDepth\tCramPath\n'+''.join(
        f'query-{arm}\tGCF_000000001.1\t{depth}\t{args.outdir}/query-{arm}.bam\n' for arm,depth in [('low',.5),('high',6.0)]))


if __name__=='__main__':main()
