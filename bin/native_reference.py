"""Original window coverage without retained depth simulations."""
import json
import time
from pathlib import Path
import numpy as np
from bin.coverage_core import BlockPool, reference_from_layout, window_sums, covered_sums
from bin.gcf_windows import GCF_CONCAT, PER_CONTIG, edges_for, segments, window_policy
from bin.ref_depth_catalog import (atomic_json, atomic_table, completed_sample, digest,
    exclusive, file_stat, validate_cram_index)


def fixed_windows(layout, size=1000):
    """Decode the frozen per-contig payload coordinates."""
    edges = edges_for(layout, PER_CONTIG, size=size)
    mapping = segments(layout, edges)
    return (edges, mapping['segment_contig_index'], mapping['segment_local_start'],
            mapping['segment_local_end'])


def native_coordinates(layout, size, count, policy):
    if policy == PER_CONTIG:
        if count is not None:
            raise ValueError('Legacy native coverage requires a fixed per-contig window size')
        edges, contigs, starts, ends = fixed_windows(layout, size)
        return edges, dict(contig_index=contigs, start=starts, end=ends)
    edges = edges_for(layout, policy, size=size, count=count)
    return edges, dict(global_start=edges[:-1], global_end=edges[1:], **segments(layout, edges))


def native_edges(layout, record):
    return edges_for(layout, window_policy(record['native']), size=record['window_size'],
                     count=record.get('bin_count'))


def save_native(pool, layout, directory, size, *, count=None, policy=PER_CONTIG):
    if getattr(pool, 'coverage_alignment_policy', 'primary_only') != 'primary_only':
        raise ValueError('Native reference windows require primary_only; supplementary coverage is a runtime choice')
    edges, coordinates = native_coordinates(layout, size, count, policy)
    sums = window_sums(pool, np.arange(len(pool)), edges)
    covered = covered_sums(pool, edges)
    widths = np.diff(edges)
    means = sums / widths
    if int(sums.sum()) != pool.total_bases:
        raise AssertionError('Original coverage mass mismatch')
    np.save(directory/'edges.npy', edges, allow_pickle=False)
    np.savez_compressed(directory/'native_coverage.npz', **coordinates,
        length=widths, aligned_bases=sums, mean_depth=means,
        covered_bases=covered, covered_fraction=covered/widths)
    return edges, dict(coverage_path='native_coverage.npz', window_size=size, windows=len(sums),
        covered_bases=int(covered.sum()), covered_fraction=float(covered.sum()/layout.length),
        mean_window_depth=float(means.mean()), median_window_depth=float(np.median(means)),
        window_depth_quantiles=dict(zip(('min','p25','p50','p75','max'), map(float,np.quantile(means,[0,.25,.5,.75,1])))),
        window_policy=policy, **({'bin_count': count} if policy == GCF_CONCAT else {}),
        full_mean_depth=pool.total_bases/layout.length)


def validate_native(directory, record):
    layout = reference_from_layout(record['layout'])
    policy = window_policy(record['native'])
    count = record.get('bin_count')
    expected, coordinates = native_coordinates(layout, record['window_size'], count, policy)
    widths = np.diff(expected)
    edges = np.load(directory/'edges.npy', allow_pickle=False)
    if not np.array_equal(edges, expected):
        raise ValueError('Native window edges mismatch')
    with np.load(directory/'native_coverage.npz', allow_pickle=False) as data:
        for name, values in dict(coordinates, length=widths).items():
            if data[name].dtype.kind not in 'iu' or not np.array_equal(data[name], values):
                raise ValueError('Native window coordinates mismatch')
        sums, covered = data['aligned_bases'], data['covered_bases']
        if (sums.shape != widths.shape or covered.shape != widths.shape or sums.dtype.kind not in 'iu'
            or covered.dtype.kind not in 'iu' or np.any(sums < 0) or np.any(covered < 0)
            or np.any(covered > widths) or np.any(covered > sums)
            or int(sums.sum()) != record['aligned_bases']
            or int(covered.sum()) != record['native']['covered_bases']
            or not np.array_equal(data['mean_depth'], sums/widths)
            or not np.array_equal(data['covered_fraction'],covered/widths)):
            raise ValueError('Native coverage content mismatch')
    if (record['genome_len'] != layout.length or record['bins'] != len(widths)
        or record['native']['windows'] != len(widths)
        or record['native']['window_size'] != record['window_size']
        or (policy == GCF_CONCAT and record['native'].get('bin_count') != count)
        or record['full_mean_depth'] != record['aligned_bases']/layout.length
        or record['native']['covered_fraction'] != record['native']['covered_bases']/layout.length
        or record.get('depths') or record.get('replicates',0) != 0):
        raise ValueError('Native metadata mismatch')
    return edges


def verify_native_sample(item):
    database, fingerprint, gcf, sample, reference_error, submitted = item
    parent=database/gcf/'samples'/sample['folder']
    row=dict(gcf=gcf,sample=sample['sample'],cram=sample['cram'],state='PENDING')
    resource=None;layout=None
    try:
        if reference_error:
            raise ValueError(reference_error)
        if file_stat(sample['cram']) != sample['source_stat']:
            raise ValueError('Source changed since plan')
        validate_cram_index(sample)
        value=completed_sample(parent,fingerprint,verify=submitted['verify'])
        if value:
            attempt,result=value
            if result['gcf'] != gcf or result['sample'] != sample['sample']:
                raise ValueError('Sample result identity differs from plan')
            if (any(result['layout'].get(key) != submitted['layout'].get(key)
                    for key in ('names','lengths','sort_policy','offsets'))
                    or result['policies'] != submitted['policies']):
                raise ValueError('Sample layout/policies differ from plan')
            if json.loads((attempt/'geometry/manifest.json').read_text())['layout'] != result['layout']:
                raise ValueError('Geometry/sample layout identities differ')
            row.update(state='COMPLETE',sample_manifest=str(attempt/'result.json'),geometry_path=str(attempt/'geometry'),
                coverage_path=str(attempt/'native_coverage.npz'),**{k:result[k] for k in ('full_mean_depth','genome_len','aligned_bases','units','blocks','window_size')},
                **{k:result['native'][k] for k in ('covered_bases','covered_fraction','windows')})
            sizes=[(p.relative_to(attempt),p.stat().st_size) for p in attempt.rglob('*') if p.is_file()]
            resource=dict(gcf=gcf,sample=sample['sample'],**{k:result[k] for k in ('elapsed_seconds','user_cpu_seconds','system_cpu_seconds','process_peak_rss_bytes')},
                geometry_bytes=sum(size for path,size in sizes if path.parts[0]=='geometry'),
                sample_output_bytes=sum(size for _,size in sizes))
            layout=result['layout']
            if file_stat(sample['cram']) != sample['source_stat']:
                raise ValueError('Source changed during collection')
        else:
            state=json.loads((parent/'progress.json').read_text()).get('state') if (parent/'progress.json').is_file() else ''
            if list(parent.glob('attempt-*/status/FAILED')) or state in ('FAILED','INTERRUPTED') or submitted['state']=='SUBMIT_FAILED':
                raise ValueError('Build or submission failed; inspect sample/GCF worker.log')
    except (ValueError,OSError,KeyError) as error:
        row.update(state='FAILED',error=str(error));resource=None;layout=None
    return row,resource,layout


def collect_native(database, plan, verify=True, jobs=4, progress=None, publish_root=None):
    """Bounded verification workers; deterministic single-writer catalogs."""
    from bin.reference_store import bounded_map
    from bin.ref_depth_worker import merge_layout
    started=time.monotonic()
    fingerprint = (database/'PLAN.sha256').read_text().strip()
    samples, resources, groups = [], [], []
    submitted = json.loads((database/'submission.json').read_text()) if (database/'submission.json').exists() else {}
    columns = ['gcf','sample','state','sample_manifest','geometry_path','coverage_path','cram',
        'full_mean_depth','genome_len','aligned_bases','units','blocks','covered_bases','covered_fraction','windows','window_size','error']
    resource_columns = ['gcf','sample','elapsed_seconds','user_cpu_seconds','system_cpu_seconds','process_peak_rss_bytes','geometry_bytes','sample_output_bytes']
    total=sum(len(group['samples']) for group in plan['groups'].values());processed=failures=0
    if progress:progress(0,total,0)
    with exclusive(database/'.catalog.lock'):
        for gcf, group in plan['groups'].items():
            rows=[];group_resources=[];layout=None;reference_error='';reference_stat=None
            try:
                reference_stat=file_stat(group['fai'])
                if digest(group['fai']) != group['fai_sha256']:
                    raise ValueError('Reference FAI changed since plan')
            except (ValueError,OSError) as error:reference_error=str(error)
            work=((database,fingerprint,gcf,sample,reference_error,
                   dict(state=submitted.get(gcf,{}).get('state'),verify=verify,layout=group['layout'],
                        policies=plan['config']['policies'])) for sample in group['samples'])
            for row,resource,observed in bounded_map(verify_native_sample,work,jobs):
                if publish_root is not None:
                    for key in ('sample_manifest','geometry_path','coverage_path'):
                        if row.get(key):
                            row[key] = str(Path(publish_root)/Path(row[key]).relative_to(database))
                if observed is not None:
                    try:layout=merge_layout(layout,observed)
                    except ValueError as error:
                        row.update(state='FAILED',error=str(error));resource=None
                rows.append(row);samples.append(row)
                if resource:group_resources.append(resource)
                processed+=1;failures+=row['state']=='FAILED'
                if progress:progress(processed,total,failures)
            try:
                if reference_error or file_stat(group['fai']) != reference_stat or digest(group['fai']) != group['fai_sha256']:
                    raise ValueError(reference_error or 'Reference FAI changed during collection')
            except (ValueError,OSError) as error:
                for row in rows:
                    row.update(state='FAILED',error=str(error))
                group_resources=[]
            resources.extend(group_resources)
            atomic_table(database/gcf/'NATIVE_INDEX.tsv',columns,rows)
            groups.append(dict(gcf=gcf,state='FAILED' if any(r['state']=='FAILED' for r in rows) else 'COMPLETE' if all(r['state']=='COMPLETE' for r in rows) else 'PENDING',
                samples=len(rows),ready_samples=sum(r['state']=='COMPLETE' for r in rows)))
        atomic_table(database/'GCF_INDEX.tsv',list(groups[0]),groups)
        atomic_table(database/'SAMPLE_INDEX.tsv',columns,samples)
        atomic_table(database/'NATIVE_INDEX.tsv',columns,samples)
        atomic_table(database/'RESOURCE_INDEX.tsv',resource_columns,resources)
        complete=all(r['state']=='COMPLETE' for r in samples)
        summary=dict(schema=plan['schema'],storage_mode='native',counts=plan['counts'],complete=complete,
            state='COMPLETE' if complete else 'INCOMPLETE',ready_samples=sum(r['state']=='COMPLETE' for r in samples),
            failed_samples=sum(r['state']=='FAILED' for r in samples),pending_samples=sum(r['state']=='PENDING' for r in samples),payloads_verified=verify)
        summary.update(collection_jobs=jobs,collection_elapsed_seconds=time.monotonic()-started)
        atomic_json(database/'SUMMARY.json',summary)
        if progress:progress(processed,total,summary['failed_samples'])
    return summary
