"""Persist and score breadth thresholds without drawing or decoding again."""
from pathlib import Path
import numpy as np
from bin.jaccard_breadth import covered_state, packed, pair, reference_means, thresholds
from bin.query_manifest_stats import aggregate

EXTRA = ('jaccard_coverage_threshold', 'query_covered_bins', 'jaccard_zero_union_count',
         'jaccard_null_zero_union_count', 'jaccard_zero_union_fraction', 'jaccard_null_zero_union_fraction')


def write_sensitivity(output, summary, draws, test, matrix, query_covered, covered,
                      edges, config, query, primary_nulls, *, draw_columns=None):
    from bin.query_manifest import table, DRAW_COLUMNS, NULL_COLUMNS, DISTRIBUTION_COLUMNS, SUMMARY_COLUMNS, high_null
    widths = np.diff(edges)
    grid = thresholds(config.get('jaccard_coverage_threshold',0.), config.get('jaccard_coverage_thresholds'))
    distribution = config.get('alpha') is not None
    summaries=[]
    root=output/'jaccard_thresholds'
    root.mkdir()
    for number,tau in enumerate(grid):
        directory=root/f't{number:02d}'
        directory.mkdir()
        qstate=covered_state(query_covered,widths,tau)
        rstates=covered_state(covered,widths,tau)
        split=qstate.ndim==2
        if split and qstate.shape[0]!=2:
            raise ValueError('Mate-stratified query breadth must contain Read1 and Read2')
        qpacked=([packed(row) for row in qstate] if split else packed(qstate))
        profiles=[packed(row) for row in rstates]
        def qprofile(draw):
            if not split:
                return qpacked, int(qstate.sum())
            side=draw.get('mate_side')
            if side not in ('read1','read2'):
                raise ValueError('Mate-stratified draw lacks read1/read2 identity')
            index=0 if side=='read1' else 1
            return qpacked[index], int(qstate[index].sum())
        qbin_summary=None if split else int(qstate.sum())
        scored=[]
        for draw in draws:
            query_bits, query_bins=qprofile(draw)
            entry=dict(draw, jaccard_coverage_threshold=tau,
                       query_covered_bins=query_bins, reference_covered_bins=int(rstates[draw['draw_id']].sum()))
            if draw['status']=='READY':
                entry.update(pair(query_bits,profiles[draw['draw_id']]))
            scored.append(entry)
        # Pearson values are copied from the primary run, without threshold-dependent arithmetic.
        if all(d['status']=='READY' for d in draws):
            if distribution or summary['arm']=='low':
                jnull=primary_nulls.get('jaccard',[]) if number==0 else reference_means(rstates)
                nulls=dict(jaccard=jnull,pearson=primary_nulls.get('pearson',[]))
                records=[dict(metric=k,left_draw=i,right_draw='',value=v,rng_sha256='',
                              left_reference_id=draws[i]['reference_id'],
                              left_reference_sample=draws[i]['reference_sample'])
                         for k,values in nulls.items() for i,v in enumerate(values)]
            else:
                nulls,records=high_null(matrix,draws,config,query,covered=covered,lengths=widths,threshold=tau)
        else:
            nulls,records={},[]
        row=dict(summary,jaccard_coverage_threshold=tau,query_covered_bins=qbin_summary)
        stats=aggregate(scored,nulls,'low' if distribution else summary['arm'],len(draws),stable_tail=distribution)['jaccard']
        if distribution and stats['valid_draws']!=len(draws):
            stats.update(value=None,mean=None)
        for source,suffix in (('value',''),('mean','_mean'),('max','_max'),('p','_p'),('status','_status'),
                              ('null_count','_null_count'),('null_mean','_null_mean'),('null_sd','_null_sd')):
            row['jaccard'+suffix]=stats[source]
        if distribution:
            p=stats['p']
            row['jaccard_decision']='NOT_EVALUATED' if p is None else 'PASS' if p>=config['alpha'] else 'FAIL'
            decisions={row['jaccard_decision'],row['pearson_decision']}
            row['detection_status']='NOT_EVALUATED' if 'NOT_EVALUATED' in decisions else 'DETECTED' if decisions=={'PASS'} else 'NOT_DETECTED'
            row['status']=summary['status'] if row['detection_status']!='NOT_EVALUATED' else ('ERROR' if summary['status']=='ERROR' else 'COMPLETE_WITH_GAPS')
        ready=[i for i,d in enumerate(draws) if d['status']=='READY']
        qzero=[i for i in ready if not (qprofile(draws[i])[0]|profiles[i])]
        # Full row-mean null uses every off-diagonal pair. Legacy high null retains its fixed pair design.
        if distribution or summary['arm']=='low':
            pairs=((i,j) for i in ready for j in ready if i<j)
        else:
            pairs=((r['left_draw'],r['right_draw']) for r in records if r['metric']=='jaccard')
        zero_count=[0]
        pair_count=[0]
        def zero_rows():
            for i in qzero:
                yield dict(kind='query_reference',left_draw='query',right_draw=i,
                           left_id=query['query_id'],right_id=draws[i]['reference_id'])
            for i,j in pairs:
                pair_count[0]+=1
                if not (profiles[i]|profiles[j]):
                    zero_count[0]+=1
                    yield dict(kind='reference_reference',left_draw=i,right_draw=j,
                               left_id=draws[i]['reference_id'],right_id=draws[j]['reference_id'])
        table(directory/'zero_unions.tsv.gz',('kind','left_draw','right_draw','left_id','right_id'),zero_rows())
        row.update(jaccard_zero_union_count=len(qzero),jaccard_null_zero_union_count=zero_count[0],
                   jaccard_zero_union_fraction=len(qzero)/len(ready) if ready else None,
                   jaccard_null_zero_union_fraction=zero_count[0]/pair_count[0] if pair_count[0] else None)
        table(directory/'draws.tsv.gz',draw_columns or DRAW_COLUMNS,scored)
        table(directory/'null.tsv.gz',NULL_COLUMNS,records)
        columns=DISTRIBUTION_COLUMNS if distribution else SUMMARY_COLUMNS
        table(directory/'summary.tsv',columns,[row])
        summaries.append(row)
        if number==0:
            for old,new in zip(draws,scored): old.update(new)
    table(output/'sensitivity.tsv',columns,summaries)
    # Store breadth explicitly once, shared by every threshold, streaming by bin.
    def bins():
        query_rows=([('query_read1',query_covered[0]),('query_read2',query_covered[1])]
                    if np.ndim(query_covered)==2 else [('query',query_covered)])
        for draw_id,vector in [*query_rows,*enumerate(covered)]:
            for i,(count,width) in enumerate(zip(vector,widths)):
                yield dict(draw_id=draw_id,bin=i,covered_bases=int(count),window_length=int(width),coverage_breadth=float(count/width))
    table(output/'breadth.tsv.gz',('draw_id','bin','covered_bases','window_length','coverage_breadth'),bins())
    return [{key:row.get(key) for key in EXTRA} for row in summaries]
