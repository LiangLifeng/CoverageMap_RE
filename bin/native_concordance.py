"""Nearest original-depth lookup and on-demand matched-depth comparisons."""
from collections import OrderedDict
import hashlib
import json
import os
from pathlib import Path
import resource
import sys
from threading import Lock
import time
import traceback
import csv

import numpy as np
from bin.coverage_core import Reference, extract_pool, merge_references, compatible_policies, runtime_pool, select_units, window_sums
from bin.concordance_core import task_rng, evaluate_windows
from bin.gcf_windows import GCF_CONCAT, edges_for, window_policy
from bin.reference_store import load_reference_bundle
from bin.ref_depth_catalog import atomic_json, digest, file_stat, safe_file, sample_folder, validate_sample_payloads
from bin.stage_timing import StageTiming

NATIVE_COLUMNS=['query_full_depth','reference_full_depth','requested_depth','comparison_depth',
    'reference_depth_difference','reference_depth_ratio','comparison_mode','reference_actual_depth',
    'jaccard_status','pearson_status','reference_count','candidate_count','depth_metric','window_size',
    'selection_mode','depth_band_matched']


class NativeCache:
    """One budget for pools/simulations; None retains entries without eviction."""
    def __init__(self,index,budget):
        self.index,self.budget=index,budget
        self.entries=OrderedDict();self.bytes=self.hits=self.misses=self.evictions=0
        self.lock = Lock()

    def stats(self):
        return dict(hits=self.hits,misses=self.misses,evictions=self.evictions,accounted_bytes=self.bytes,budget_bytes=self.budget)

    def close(self):
        self.entries.clear();self.bytes=0

    def put(self,key,data,size,stamps):
        if self.budget is not None and size > self.budget:return
        with self.lock:
            old = self.entries.pop(key, None)
            if old is not None:
                self.bytes -= old[1]
            while self.budget is not None and self.entries and self.bytes+size > self.budget:
                _,entry=self.entries.popitem(last=False);self.bytes-=entry[1];self.evictions+=1
            self.entries[key]=(data,size,stamps);self.bytes+=size

    def stamps(self,row):
        attempt=safe_file(self.index.root_for(row),row['manifest_path']).parent
        if (attempt/'status/FAILED').exists() or (attempt/'status/RUNNING').exists():
            raise ValueError('Reference no longer successful')
        if not (attempt/'status/SUCCESS').is_file() or (attempt/'status/exit_code.txt').read_text().strip()!='0':
            raise ValueError('Reference status invalid')
        # Check registered file stamps even on cache hits; initial load verifies hashes.
        checks=json.loads((attempt/'checksums.json').read_text())
        paths=[attempt/'checksums.json',attempt/'status/SUCCESS',attempt/'status/exit_code.txt']
        paths += [safe_file(attempt,n) for n in sorted(checks)]
        return [(str(p),file_stat(p)) for p in paths]

    def get(self,key,stamps):
        with self.lock:
            if key in self.entries:
                data,size,old=self.entries[key]
                if old==stamps:
                    self.entries.move_to_end(key);self.hits+=1;return data
                self.bytes-=size;del self.entries[key]
            self.misses+=1
            return None

    def check_window_policy(self, record):
        config = self.index.metadata['config']
        if (window_policy(record['native']) != self.index.window_policy
                or record['window_size'] != config['window_size']
                or record.get('bin_count') != config.get('bin_count')):
            raise ValueError('Native reference window policy or parameters differ from SQLite')

    def load(self,row,timing):
        with timing.measure('reference_metadata'):
            stamps=self.stamps(row)
            manifest=safe_file(self.index.root_for(row),row['manifest_path']);attempt=manifest.parent
            if digest(manifest)!=row['manifest_sha'] or digest(attempt/'checksums.json')!=row['checksums_sha']:
                raise ValueError('Reference checksum differs from SQLite index')
            cached=self.get(('pool',row['id']),stamps)
        if cached is not None:return cached
        record=json.loads(manifest.read_text())
        self.check_window_policy(record)
        with timing.measure('reference_payload_sha256'):
            checks=json.loads((attempt/'checksums.json').read_text())
            validate_sample_payloads(attempt,record,checks,verify=True)
        with timing.measure('reference_load_decompress'):
            pool,layout,metadata=load_reference_bundle(attempt/'geometry',workers=1)
            if (record['gcf']!=row['gcf'] or record['sample']!=row['sample'] or record['full_mean_depth']!=row['full_mean_depth']
                or pool.total_bases!=record['aligned_bases'] or metadata['policies']!=record['policies']
                or metadata['layout']!=record['layout']):
                raise ValueError('Native reference metadata mismatch')
            edges=np.load(attempt/'edges.npy',allow_pickle=False)
            if stamps!=self.stamps(row):raise ValueError('Reference changed during loading')
            data=(pool,layout,record,edges)
        self.put(('pool',row['id']),data,pool.nbytes+edges.nbytes+len(json.dumps(record).encode()),stamps)
        return data

    def simulate(self,row,pool,edges,depth,config,timing):
        stamps=self.stamps(row)
        key=('simulation',row['id'],hashlib.sha256(edges.astype('<i8').tobytes()).hexdigest(),
             float(depth).hex(),config['seed'],config['replicates'],config['sampler'],config.get('policies',{}).get('sampling_unit','read'),
             config.get('policies',{}).get('coverage_alignment_policy','primary_only'))
        value=self.get(key,stamps)
        if value is not None:return value
        with timing.measure('reference_sample_count'):
            matrix=np.empty((config['replicates'],len(edges)-1),dtype=np.int64);audits=[]
            for replicate in range(config['replicates']):
                ids=select_units(pool,task_rng(config['seed'],depth,0,replicate),target_bases=min(depth*pool.genome_len,pool.total_bases),depth_sampler=config['sampler'])
                matrix[replicate]=window_sums(pool,ids,edges)
                audits.append(dict(replicate=replicate,units=len(ids),aligned_bases=int(matrix[replicate].sum()),actual_depth=int(matrix[replicate].sum())/pool.genome_len,
                    selected_ids_sha256=hashlib.sha256(ids.astype('<i8').tobytes()).hexdigest()))
                if int(matrix[replicate].sum())!=int(pool.unit_bases[ids].sum()):raise AssertionError('Reference sampled mass mismatch')
            value=(matrix,audits)
        self.put(key,value,matrix.nbytes+len(json.dumps(audits).encode()),stamps)
        return value


def descriptive(test,reference):
    positive=test>0;other=reference>0
    union=np.count_nonzero(positive|other)
    j=float(np.count_nonzero(positive&other)/union) if union else None
    estimable=len(test)>1 and np.ptp(test)>0 and np.ptp(reference)>0
    pearson=float(np.clip(np.corrcoef(test,reference)[0,1],-1,1)) if estimable else None
    def entry(value):
        return dict(similarity=value,status='descriptive_only' if value is not None else 'non_estimable',p_value=None,accepted=None)
    return dict(jaccard=entry(j),pearson=entry(pearson),model='original-reference-descriptive-only',simulations=0)


def compare(test,matrix,widths,config):
    # Terminal short windows contribute mean depth, not artificially low total mass.
    test=test/widths;matrix=matrix/widths
    with np.errstate(divide='ignore',invalid='ignore'):
        result=evaluate_windows(test,matrix,config['alpha'],config['similarity_block_size'])
    for metric in ('jaccard','pearson'):
        value=result[metric]
        if value.get('reference_sd')==0 or (value.get('p_value') is not None and not np.isfinite(value['p_value'])):
            value.update(status='non_estimable',reason='zero_variance_reference_distribution',p_value=None,accepted=None)
    return result


def perform_native(query,config,index,cache):
    if config.get('comparison_policy') == 'depth-band':
        from bin.depth_band_concordance import perform_depth_band
        return perform_depth_band(query,config,index,cache)
    from bin.concordance_batch import COLUMNS
    directory=Path(config['outdir'])/'queries'/query['gcf']/sample_folder(query['sample'])
    (directory/'status').mkdir(parents=True);status=directory/'status';(status/'RUNNING').touch()
    timing=StageTiming();started=time.perf_counter();before=resource.getrusage(resource.RUSAGE_SELF)
    hits=cache.stats();rows=[];ready=failed=gaps=0;calls=0;full_depth=None;geometry_loads=0
    try:
        group=index.group(query['gcf'],missing_ok=True)
        if group is None:
            rows.append(dict(query_sample=query['sample'],gcf=query['gcf'],status='NO_REFERENCE_GCF'))
            gaps=1
        else:
            layout_file=directory/'query_layout.fai'
            layout_file.write_text(''.join(f'{n}\t{l}\n' for n,l in zip(group['names'],group['lengths'])))
            with timing.measure('query_sha256'):
                source_stat=file_stat(query['alignment']);source_hash=digest(query['alignment'])
            with timing.measure('query_extract'):
                order = {'sort_contigs': True} if index.window_policy == GCF_CONCAT else {}
                pool,layout=extract_pool(query['alignment'],reference_index=layout_file,**order,**config['policies']);calls=1
            full_depth=pool.total_bases/pool.genome_len
            edges=edges_for(layout,index.window_policy,size=config['window_size'],count=config.get('bin_count'))
            widths=np.diff(edges)
            np.save(directory/'edges.npy',edges,allow_pickle=False)
            depths=config['depths'] or [full_depth]
            for number,target in enumerate(depths):
                base=dict(query_sample=query['sample'],gcf=query['gcf'],depth=target,requested_depth=target,query_full_depth=full_depth)
                if target<=0 or target>full_depth:
                    rows.append(dict(base,status='UNSUPPORTED_QUERY_DEPTH'));gaps+=1;continue
                with timing.measure('reference_sql'):
                    selected=index.nearest(query['gcf'],target,config['reference_limit'],
                        exclude_sample=None if config['include_self'] else query['sample'],
                        exclude_source=None if config['include_self'] else source_hash,samples=config['reference_samples'])
                if not selected:
                    rows.append(dict(base,status='NO_ELIGIBLE_REFERENCE'));gaps+=1;continue
                for row in selected:
                    entry=dict(base,reference_sample=row['sample'],reference_full_depth=row['full_mean_depth'],
                        reference_depth_difference=abs(row['full_mean_depth']-target),reference_depth_ratio=row['full_mean_depth']/target,
                        sample_manifest=str(index.root_for(row)/row['manifest_path']))
                    try:
                        old_misses=cache.misses
                        ref,ref_layout,record,ref_edges=cache.load(row,timing)
                        ref=runtime_pool(ref,record['policies'],config['policies'])
                        geometry_loads+=cache.misses-old_misses
                        merge_references(ref_layout,layout)
                        if not compatible_policies(record['policies'],config['policies']) or not np.array_equal(ref_edges,edges):raise ValueError('Reference/query layout/policy mismatch')
                        common=min(target,row['full_mean_depth'])
                        with timing.measure('query_sample_count'):
                            ids=select_units(pool,task_rng(config['seed'],common,1,0),target_bases=min(common*pool.genome_len,pool.total_bases),depth_sampler=config['sampler'])
                            test=window_sums(pool,ids,edges)
                            if int(test.sum())!=int(pool.unit_bases[ids].sum()):raise AssertionError('Query sampled mass mismatch')
                        output=directory/'comparisons'/sample_folder(row['sample'])/f'target-{number:05d}'
                        output.mkdir(parents=True)
                        if common==row['full_mean_depth']:
                            with timing.measure('reference_full_coverage'):
                                with np.load(index.root_for(row)/row['coverage_path'],allow_pickle=False) as saved:original=saved['aligned_bases']
                            with timing.measure('statistics'):result=descriptive(test/widths,original/widths)
                            audits=[];matrix=None;mode='original_reference_at_common_depth'
                            ref_actual=row['full_mean_depth']
                        else:
                            matrix,audits=cache.simulate(row,ref,edges,common,config,timing)
                            with timing.measure('statistics'):result=compare(test,matrix,widths,config)
                            mode='dynamic_reference_sampling';ref_actual=float(np.mean([a['actual_depth'] for a in audits]))
                        entry.update(comparison_depth=common,comparison_mode=mode,test_actual_depth=int(test.sum())/pool.genome_len,
                            reference_actual_depth=ref_actual,status='READY',jaccard=result['jaccard']['similarity'],pearson=result['pearson']['similarity'],
                            jaccard_status=result['jaccard']['status'],pearson_status=result['pearson']['status'],
                            jaccard_p=result['jaccard']['p_value'],pearson_p=result['pearson']['p_value'],
                            jaccard_accepted=result['jaccard']['accepted'],pearson_accepted=result['pearson']['accepted'],result_path=str(output/'concordance.json'))
                        scope='same-source-self-consistency' if source_hash==row['source_sha'] else 'same-sample-self-consistency' if query['sample']==row['sample'] else 'separate-input-pools'
                        entry['comparison_scope']=scope
                        result.update(**{k:v for k,v in entry.items() if k not in ('jaccard','pearson')},selection_policy='nearest-original-mean-depth-absolute-difference',window_size=config['window_size'],
                            window_policy=index.window_policy,bin_count=config.get('bin_count'),
                            statistics_calibrated=False,reference_replicates=audits,reference_geometry_loaded=True,
                            reference_sampled=matrix is not None,seed=config['seed'],sampler=config['sampler'],
                            test_selected_ids_sha256=hashlib.sha256(ids.astype('<i8').tobytes()).hexdigest(),
                            reference_missing_header_contigs=record.get('missing_header_contigs',[]),test_missing_header_contigs=list(layout.missing_header_contigs),
                            identity_unverified_contigs=[n for n,a,b in zip(layout.names,layout.md5,ref_layout.md5) if not a or not b])
                        with timing.measure('query_output'):
                            np.save(output/'test.npy',test,allow_pickle=False);np.save(output/'test_ids.npy',ids,allow_pickle=False)
                            atomic_json(output/'concordance.json',result)
                        ready+=1
                    except Exception as error:
                        failed+=1;entry.update(status='ERROR',error=str(error))
                    rows.append(entry)
                atomic_json(directory/'progress.json',dict(state='RUNNING',ready_comparisons=ready,depth_gaps=gaps,failed_comparisons=failed))
            with timing.measure('query_sha256'):
                if digest(query['alignment'])!=source_hash or file_stat(query['alignment'])!=source_stat:raise ValueError('Query changed during calculation')
        state='FAILED' if failed else 'SKIPPED' if not ready else 'COMPLETE_WITH_GAPS' if gaps else 'COMPLETE'
        reason='' if ready or failed else 'NO_REFERENCE_GCF' if group is None else 'NO_COMPARABLE_DEPTH'
    except BaseException as error:
        if isinstance(error,(KeyboardInterrupt,SystemExit)):raise
        failed+=1;state='FAILED';reason=str(error)
        rows.append(dict(query_sample=query['sample'],gcf=query['gcf'],status='ERROR',error=str(error)))
        (directory/'error.txt').write_text(traceback.format_exc())
    with (directory/'summary.tsv').open('w',newline='') as handle:
        writer=csv.DictWriter(handle,COLUMNS,delimiter='\t');writer.writeheader();writer.writerows(rows)
    after=resource.getrusage(resource.RUSAGE_SELF);now=cache.stats()
    report=dict(query=query,state=state,reason=reason,query_full_depth=full_depth,extraction_calls=calls,reference_geometry_loads=geometry_loads,
        ready_comparisons=ready,failed_comparisons=failed,depth_gaps=gaps,cache={k:now[k]-hits[k] for k in ('hits','misses','evictions')},cache_after=now,
        worker_pid=os.getpid(),stage_timing=timing.snapshot(),elapsed_seconds=time.perf_counter()-started,
        user_cpu_seconds=after.ru_utime-before.ru_utime,system_cpu_seconds=after.ru_stime-before.ru_stime,
        process_peak_rss_bytes=int(after.ru_maxrss if sys.platform=='darwin' else after.ru_maxrss*1024),statistics_calibrated=False)
    atomic_json(directory/'provenance.json',report);atomic_json(directory/'progress.json',dict(state=state))
    (status/'RUNNING').unlink();(status/'exit_code.txt').write_text('1\n' if failed else '0\n')
    (status/('FAILED' if failed else 'SKIPPED' if state=='SKIPPED' else 'SUCCESS')).touch()
    return dict(query=query,state=state,reason=reason,report=str(directory/'provenance.json'),summary=str(directory/'summary.tsv'),exit_code=int(bool(failed)))


def run_native(args,depths):
    from bin.concordance_batch import run_indexed
    if args.reference_limit<1:raise ValueError('reference-limit must be positive')
    if args.depth_sampler=='auto':
        raise ValueError('Native queries require an explicit permutation/lazy sampler for reproducible cache keys')
    return run_indexed(args,depths)
