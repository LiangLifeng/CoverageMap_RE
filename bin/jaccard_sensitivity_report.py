#!/usr/bin/env python3
"""Summarize paired TP/FP threshold results; retain all unscored/error rows."""
import argparse
from collections import Counter,defaultdict
import csv
import html
import json
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from bin.query_manifest import table


def summarize(tp,fp,outdir):
    groups=defaultdict(list)
    for label,path in [('TP',tp),('FP',fp)]:
        with path.open() as handle:
            for row in csv.DictReader(handle,delimiter='\t'):
                tau=float(row['jaccard_coverage_threshold'])
                groups[label,tau].append(row)
    universes={}
    pearson={}
    for (label,tau),rows in groups.items():
        identities=[(r['query_id'],r['Sample'],r['GCFID']) for r in rows]
        if len(identities)!=len(set(identities)): raise ValueError('Duplicate query within threshold')
        membership=set(identities)
        if label in universes and membership!=universes[label]: raise ValueError('Threshold cohorts differ')
        universes[label]=membership
        for row in rows:
            key=label,row['query_id']
            signature=tuple(row.get(k,'') for k in ('pearson','pearson_p','pearson_status','pearson_null_mean','pearson_null_sd'))
            if key in pearson and pearson[key]!=signature: raise ValueError('Pearson changed across thresholds')
            pearson[key]=signature
    if not groups or {t for l,t in groups if l=='TP'}!={t for l,t in groups if l=='FP'}:
        raise ValueError('TP and FP require matching nonempty threshold grids')
    outdir.mkdir(parents=True,exist_ok=False)
    rows=[]
    for (label,tau),items in sorted(groups.items()):
        counts=Counter(r['detection_status'] for r in items)
        if set(counts)-{'DETECTED','NOT_DETECTED','NOT_EVALUATED'}: raise ValueError('Unknown decision')
        rows.append(dict(cohort=label,threshold=tau,total=len(items),detected=counts['DETECTED'],
            not_detected=counts['NOT_DETECTED'],not_evaluated=counts['NOT_EVALUATED'],
            errors=sum(r['status']=='ERROR' for r in items),detected_fraction=counts['DETECTED']/len(items),
            query_zero_union=sum(int(r['jaccard_zero_union_count'] or 0) for r in items),
            reference_zero_union=sum(int(r['jaccard_null_zero_union_count'] or 0) for r in items)))
    table(outdir/'comparison.tsv',tuple(rows[0]),rows)
    svg=['<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 720 400" role="img" aria-label="Detected fraction over all query rows by threshold">',
         '<rect width="720" height="400" fill="white"/><path d="M60 30V340H680" fill="none" stroke="#222"/>']
    for value in (0,.25,.5,.75,1):
        y=340-value*280
        svg.append(f'<text x="15" y="{y+4}" font-size="12">{value:.0%}</text><path d="M60 {y}H680" stroke="#ddd"/>')
    grid=sorted({r['threshold'] for r in rows}); maximum=max(grid) or 1
    for label,color in [('TP','#2878a3'),('FP','#bf5b32')]:
        items=sorted((r for r in rows if r['cohort']==label),key=lambda r:r['threshold'])
        points=' '.join(f'{60+600*r["threshold"]/maximum},{340-280*r["detected_fraction"]}' for r in items)
        svg.append(f'<polyline points="{points}" fill="none" stroke="{color}" stroke-width="2"/>')
        for r in items:
            svg.append(f'<circle cx="{60+600*r["threshold"]/maximum}" cy="{340-280*r["detected_fraction"]}" r="4" fill="{color}"/>')
        svg.append(f'<text x="{80 if label=="TP" else 340}" y="25" fill="{color}">{label} detected / all query rows</text>')
    for tau in grid: svg.append(f'<text x="{60+600*tau/maximum}" y="365" text-anchor="middle" font-size="12">{tau:g}</text>')
    svg.append('<text x="300" y="390">Coverage breadth threshold</text></svg>')
    (outdir/'sensitivity.svg').write_text(''.join(svg))
    headers=tuple(rows[0]); cells=['<table><tr>'+''.join(f'<th>{html.escape(k)}</th>' for k in headers)+'</tr>']
    for row in rows: cells.append('<tr>'+''.join(f'<td>{html.escape(str(row[k]))}</td>' for k in headers)+'</tr>')
    cells.append('</table>')
    (outdir/'index.html').write_text('<!doctype html><meta charset="utf-8"><title>Jaccard breadth sensitivity</title><style>body{font:16px system-ui;max-width:1200px;margin:30px auto}td,th{border:1px solid #ddd;padding:7px}table{border-collapse:collapse}svg{max-width:850px}</style><h1>Paired threshold sensitivity</h1><p>All query rows remain in denominators. Counts describe sample–GCF queries, not unique patients. Unscored results remain NOT_EVALUATED. No optimal threshold is selected.</p>'+''.join(svg)+''.join(cells))
    (outdir/'validation.json').write_text(json.dumps(dict(cohorts={k:len(v) for k,v in universes.items()},thresholds=grid,paired_membership=True,pearson_invariant=True),indent=2)+'\n')


def main():
    cli=argparse.ArgumentParser(description=__doc__)
    cli.add_argument('--version',action='version',version='1.0.0')
    for key in ('tp','fp','outdir'): cli.add_argument('--'+key,type=Path,required=True)
    args=cli.parse_args()
    try: summarize(args.tp,args.fp,args.outdir)
    except (ValueError,OSError,KeyError) as error: cli.exit(2,f'ERROR: {error}\n')

if __name__=='__main__': main()
