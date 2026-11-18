#!/usr/bin/env python3
"""Compare experimental Jaccard decision methods on frozen TP and FP query sets.

Counts include failed and NOT_EVALUATED queries in denominators. The units are
Sample/GCF query records, not unique patients. No method is auto-selected.
"""
import argparse
from collections import Counter
import csv
from pathlib import Path

OUT_COLUMNS = ('method', 'tp_total', 'tp_detected', 'tp_not_detected', 'tp_not_evaluated',
               'tpr_total', 'tpr_evaluable', 'fp_total', 'fp_detected', 'fp_not_detected',
               'fp_not_evaluated', 'fpr_total', 'fpr_evaluable', 'delta_tpr_vs_normal',
               'delta_fpr_vs_normal')
STATES = ('DETECTED', 'NOT_DETECTED', 'NOT_EVALUATED')


def load(path):
    """Require identical per-query rosters for every method in this cohort."""
    with Path(path).open(encoding='utf-8', newline='') as handle:
        reader = csv.DictReader(handle, delimiter='\t')
        required = {'query_id','Sample','GCFID','method','detection_status'}
        if not required <= set(reader.fieldnames or []):
            raise ValueError('Missing decision method columns in ' + str(path))
        rows = list(reader)
    if not rows:
        raise ValueError('Empty decision table: ' + str(path))
    grouped = {}
    for row in rows:
        method = row['method']
        key = (row['query_id'], row['Sample'], row['GCFID'])
        if not all(key) or not method or row['detection_status'] not in STATES:
            raise ValueError('Incomplete identity/method/outcome in ' + str(path))
        bucket = grouped.setdefault(method, {})
        if key in bucket:
            raise ValueError('Duplicate method/query row: ' + str((method, key)))
        bucket[key] = row['detection_status']
    rosters = list(grouped.values())
    if any(set(bucket) != set(rosters[0]) for bucket in rosters[1:]):
        raise ValueError('Different query rosters across methods in ' + str(path))
    return grouped


def summarize(tp, fp):
    if set(tp) != set(fp):
        raise ValueError('TP/FP method sets differ')
    result=[]
    normal_tpr=normal_fpr=None
    for method in sorted(tp, key=lambda k:(k!='normal', k)):
        true=Counter(tp[method].values())
        false=Counter(fp[method].values())
        tp_n,fp_n=len(tp[method]),len(fp[method])
        tp_eval=true['DETECTED']+true['NOT_DETECTED']
        fp_eval=false['DETECTED']+false['NOT_DETECTED']
        tpr=true['DETECTED']/tp_n
        fpr=false['DETECTED']/fp_n
        row=dict(method=method, tp_total=tp_n, tp_detected=true['DETECTED'],
                 tp_not_detected=true['NOT_DETECTED'], tp_not_evaluated=true['NOT_EVALUATED'],
                 tpr_total=tpr, tpr_evaluable=true['DETECTED']/tp_eval if tp_eval else None,
                 fp_total=fp_n, fp_detected=false['DETECTED'],
                 fp_not_detected=false['NOT_DETECTED'], fp_not_evaluated=false['NOT_EVALUATED'],
                 fpr_total=fpr, fpr_evaluable=false['DETECTED']/fp_eval if fp_eval else None,
                 delta_tpr_vs_normal=None, delta_fpr_vs_normal=None)
        if method == 'normal':
            normal_tpr,normal_fpr=tpr,fpr
        result.append(row)
    if normal_tpr is not None:
        for row in result:
            row['delta_tpr_vs_normal']=row['tpr_total']-normal_tpr
            row['delta_fpr_vs_normal']=row['fpr_total']-normal_fpr
    return result


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--tp',type=Path,required=True,help='TP decision_methods.tsv')
    parser.add_argument('--fp',type=Path,required=True,help='FP decision_methods.tsv')
    parser.add_argument('--out',type=Path,required=True,help='New summary TSV (must not exist)')
    args=parser.parse_args()
    if args.out.exists():
        parser.error('Output already exists: '+str(args.out))
    result=summarize(load(args.tp), load(args.fp))
    args.out.parent.mkdir(parents=True,exist_ok=True)
    with args.out.open('w',encoding='utf-8',newline='') as handle:
        writer=csv.DictWriter(handle,OUT_COLUMNS,delimiter='\t',lineterminator='\n')
        writer.writeheader()
        writer.writerows(result)


if __name__ == '__main__':
    main()
