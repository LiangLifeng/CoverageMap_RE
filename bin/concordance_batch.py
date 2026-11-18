"""Read-only indexed reference comparisons with one extraction per query CRAM."""
from collections import OrderedDict
from concurrent.futures import ProcessPoolExecutor, wait, FIRST_COMPLETED
import csv
import hashlib
import json
import multiprocessing
import os
from pathlib import Path
import resource
import sys
import time
import traceback

import numpy as np

from bin.coverage_core import (Reference, extract_pool, merge_references, select_units,
                               window_edges, window_sums)
from bin.concordance_core import evaluate_windows, task_rng, RNG_POLICY
from bin.concordance_inputs import file_hash, resolve_policies
from bin.ref_depth_catalog import atomic_json, digest, file_stat, input_file, read_rows, safe_file, sample_folder
from bin.ref_depth_worker import load_depth
from bin.reference_index import ReferenceIndex, depth_key
from bin.stage_timing import StageTiming
from bin.native_concordance import NATIVE_COLUMNS

INDEX = None
CACHE = None


class MatrixCache:
    """Per-process verified LRU; budget is accounted payload bytes, not RSS."""
    def __init__(self, index, budget):
        self.index, self.budget = index, budget
        self.entries = OrderedDict()
        self.bytes = self.hits = self.misses = self.evictions = 0

    def stats(self):
        return dict(hits=self.hits, misses=self.misses, evictions=self.evictions,
                    accounted_bytes=self.bytes, budget_bytes=self.budget)

    def discard(self, key):
        value = self.entries.pop(key)
        self.bytes -= value['bytes']
        matrix = value['data'][0]
        if isinstance(matrix, np.memmap):
            matrix._mmap.close()

    def close(self):
        while self.entries:
            self.discard(next(iter(self.entries)))

    def load(self, row, depth, timing=None):
        timing = timing or StageTiming()
        with timing.measure('reference_metadata'):
            root = self.index.root_for(row)
            manifest = safe_file(root, row['manifest_path'])
            attempt = manifest.parent
            paths = [manifest, attempt / 'checksums.json', attempt / 'status/SUCCESS',
                     attempt / 'status/exit_code.txt', safe_file(root, row['matrix_path']),
                     attempt / 'edges.npy', safe_file(root, row['metadata_path'])]
            if (attempt / 'status/FAILED').exists() or (attempt / 'status/RUNNING').exists():
                raise ValueError('Reference status is no longer successful')
            stamps = [file_stat(path) for path in paths]
            key = (row['id'], depth_key(depth))
            if key in self.entries:
                value = self.entries[key]
                if value['stamps'] == stamps:
                    self.entries.move_to_end(key); self.hits += 1
                    return value['data']
                self.discard(key)
            self.misses += 1
            if digest(manifest) != row['manifest_sha'] or digest(attempt / 'checksums.json') != row['checksums_sha']:
                raise ValueError('Reference manifest/checksum registry differs from SQLite index')
            checks = json.loads((attempt / 'checksums.json').read_text())
            for path, expected in ((paths[4], row['matrix_sha']), (paths[5], row['edges_sha']), (paths[6], row['metadata_sha'])):
                if checks.get(path.relative_to(attempt).as_posix()) != expected:
                    raise ValueError('Indexed payload checksum differs from registry')
            data = load_depth(manifest, depth, timing=timing)
            if stamps != [file_stat(path) for path in paths]:
                raise ValueError('Reference changed during load')
        data[0].setflags(write=False); data[1].setflags(write=False)
        size = data[0].nbytes + data[1].nbytes + len(json.dumps(data[2]).encode())
        if size <= self.budget:
            while self.entries and self.bytes + size > self.budget:
                self.discard(next(iter(self.entries))); self.evictions += 1
            self.entries[key] = dict(data=data, bytes=size, stamps=stamps)
            self.bytes += size
        return data


def initialize(index_path, database_root, budget, native=False):
    global INDEX, CACHE
    if CACHE is not None:
        CACHE.close(); INDEX.close()
    INDEX = ReferenceIndex(index_path, database_root)
    if native:
        from bin.native_concordance import NativeCache
        CACHE=NativeCache(INDEX,budget)
    else:
        CACHE = MatrixCache(INDEX, budget)


def read_queries(args):
    if args.query_manifest is not None:
        if args.gcf or args.query_sample:
            raise ValueError('GCF/sample must come from query manifest in manifest mode')
        source = args.query_manifest.resolve()
        rows = read_rows(source, ['Sample', 'CramFilePath', 'GCFID'])
        queries = [dict(sample=r['Sample'], gcf=r['GCFID'],
                        alignment=str(input_file(r['CramFilePath'], source))) for r in rows]
    else:
        if not args.test_alignment or not args.gcf or not args.query_sample:
            raise ValueError('SQLite single query requires --test-alignment, --gcf and --query-sample')
        queries = [dict(sample=args.query_sample, gcf=args.gcf, alignment=str(args.test_alignment.resolve()))]
    keys = [(q['sample'], q['gcf']) for q in queries]
    if len(set(keys)) != len(keys):
        raise ValueError('Duplicate query sample/GCF')
    for query in queries:
        if (not query['sample'].strip() or not query['gcf'].strip()
                or any(c in query['sample'] + query['gcf'] for c in '\r\n\t\x00')
                or not Path(query['alignment']).is_file()):
            raise ValueError('Invalid query identity/path')
        query['folder'] = sample_folder(query['sample']) + '-' + sample_folder(query['gcf'])
    return sorted(queries, key=lambda q: (q['gcf'], q['sample']))


COLUMNS = ['query_sample', 'gcf', 'reference_sample', 'depth', 'status', 'comparison_scope',
           'jaccard', 'pearson', 'jaccard_p', 'pearson_p', 'jaccard_accepted', 'pearson_accepted',
           'test_actual_depth', 'sample_manifest', 'result_path', 'error'] + NATIVE_COLUMNS


def missing_gcf(query, config, directory, status, timing, started):
    reason = 'NO_REFERENCE_GCF'
    message = 'Query GCF is absent from reference SQLite index'
    with (directory / 'summary.tsv').open('w', newline='') as handle:
        writer = csv.DictWriter(handle, COLUMNS, delimiter='\t'); writer.writeheader()
        for depth in config['depths']:
            writer.writerow(dict(query_sample=query['sample'], gcf=query['gcf'], depth=depth,
                                 status=reason, error=message))
    report = dict(query=query, state='SKIPPED', reason=reason, message=message,
        ready_comparisons=0, failed_comparisons=0, depth_gaps=len(config['depths']),
        extraction_calls=0, reference_geometry_loads=0, worker_pid=os.getpid(),
        cache=dict(hits=0, misses=0, evictions=0), cache_after=CACHE.stats(),
        stage_timing=timing.snapshot(), elapsed_seconds=time.perf_counter()-started)
    atomic_json(directory / 'provenance.json', report)
    atomic_json(directory / 'progress.json', dict(state='SKIPPED', reason=reason))
    (status / 'RUNNING').unlink(); (status / 'SKIPPED').touch(); (status / 'exit_code.txt').write_text('0\n')
    return dict(query=query, state='SKIPPED', reason=reason, report=str(directory / 'provenance.json'),
                summary=str(directory / 'summary.tsv'), exit_code=0)


def perform_query(query, config):
    if config.get('native'):
        from bin.native_concordance import perform_native
        return perform_native(query,config,INDEX,CACHE)
    directory = Path(config['outdir']) / 'queries' / query['folder']
    directory.mkdir(parents=True)
    status = directory / 'status'; status.mkdir(); (status / 'RUNNING').touch()
    started = time.perf_counter(); before = resource.getrusage(resource.RUSAGE_SELF)
    timing = StageTiming()
    cache_before = CACHE.stats()
    compared = gaps = failed = excluded = eligible = 0
    extraction_calls = 0
    try:
        group = INDEX.group(query['gcf'], missing_ok=True)
        if group is None:
            return missing_gcf(query, config, directory, status, timing, started)
        policies = config['policies']
        layout_file = directory / 'query_layout.fai'
        layout_file.write_text(''.join(f'{n}\t{length}\n' for n, length in zip(group['names'], group['lengths'])))
        source_stat = file_stat(query['alignment'])
        with timing.measure('query_sha256'):
            source_hash = file_hash(query['alignment'])
        with timing.measure('query_extract'):
            tick = time.perf_counter()
            extraction_calls += 1
            pool, query_layout = extract_pool(query['alignment'], reference_index=layout_file, **policies)
            extraction_seconds = time.perf_counter() - tick
        with timing.measure('query_output'):
            edges = window_edges(query_layout.length, config['bins'])
            np.save(directory / 'edges.npy', edges, allow_pickle=False)
        draws = {}
        with (directory / 'summary.tsv').open('w', newline='') as table:
            writer = csv.DictWriter(table, COLUMNS, delimiter='\t'); writer.writeheader()
            for depth in config['depths']:
                depth_eligible = 0
                for row in timed_references(query['gcf'], depth, config['reference_batch_size'], timing):
                    if config['reference_samples'] and row['sample'] not in config['reference_samples']:
                        continue
                    self_source = row['source_sha'] == source_hash
                    self_sample = row['sample'] == query['sample']
                    if not config['include_self'] and (self_source or self_sample):
                        excluded += 1; continue
                    eligible += 1; depth_eligible += 1
                    scope = 'same-source-self-consistency' if self_source else 'same-sample-self-consistency' if self_sample else 'separate-input-pools'
                    entry = dict(query_sample=query['sample'], gcf=query['gcf'], reference_sample=row['sample'],
                        depth=depth, comparison_scope=scope, sample_manifest=str(INDEX.root_for(row) / row['manifest_path']))
                    try:
                        if depth * pool.genome_len > pool.total_bases:
                            entry.update(status='UNSUPPORTED_QUERY_DEPTH', error='Query aligned mass below target')
                        elif row['depth_status'] != 'READY':
                            entry.update(status='UNSUPPORTED_REFERENCE_DEPTH' if row['depth_status'] == 'UNSUPPORTED_DEPTH'
                                         else 'REFERENCE_DEPTH_NOT_BUILT', error='No exact READY reference depth')
                        else:
                            matrix, saved_edges, metadata = CACHE.load(row, depth, timing=timing)
                            with timing.measure('reference_compatibility'):
                                record, ref_meta = metadata['sample'], metadata['reference']
                                if record['gcf'] != query['gcf'] or record['sample'] != row['sample']:
                                    raise ValueError('Indexed reference identity mismatch')
                                reference_layout = Reference(*(tuple(record['layout'][k]) for k in ('names', 'lengths', 'md5')),
                                                             tuple(record.get('missing_header_contigs', [])))
                                merge_references(reference_layout, query_layout)
                                if record['policies'] != policies or not np.array_equal(saved_edges, edges):
                                    raise ValueError('Reference/query policy or window layout mismatch')
                                sampler = ref_meta['sampler']['sampler']
                                draw_key = (depth_key(depth), sampler)
                            if draw_key not in draws:
                                with timing.measure('query_sample_count'):
                                    ids = select_units(pool, task_rng(config['seed'], depth, 1, 0),
                                        target_bases=depth * pool.genome_len, depth_sampler=sampler)
                                    test = window_sums(pool, ids, edges)
                                    if int(test.sum()) != int(pool.unit_bases[ids].sum()):
                                        raise AssertionError('Test mass mismatch')
                                with timing.measure('query_output'):
                                    draw_dir = directory / 'draws' / f'draw-{len(draws):05d}-{sampler}'
                                    draw_dir.mkdir(parents=True)
                                    np.save(draw_dir / 'test.npy', test, allow_pickle=False)
                                    np.save(draw_dir / 'test_ids.npy', ids, allow_pickle=False)
                                    draws[draw_key] = dict(test=test, units=len(ids), actual_depth=int(test.sum()) / pool.genome_len,
                                        ids_sha256=hashlib.sha256(ids.astype('<i8').tobytes()).hexdigest(), directory=str(draw_dir))
                            draw = draws[draw_key]
                            with timing.measure('statistics'):
                                tick = time.perf_counter()
                                result = evaluate_windows(draw['test'], matrix, config['alpha'], config['similarity_block_size'])
                                result.update(depth=depth, reference_replicates=ref_meta['reference_replicates'],
                                    depth_sampler=ref_meta['sampler'], test_units=draw['units'], test_actual_depth=draw['actual_depth'],
                                    test_selected_ids_sha256=draw['ids_sha256'], test_draw_directory=draw['directory'],
                                    reference_sample=row['sample'], sample_manifest=entry['sample_manifest'],
                                    reference_seed=record['seed'], comparison_scope=scope,
                                    reference_missing_header_contigs=record.get('missing_header_contigs', []),
                                    test_missing_header_contigs=list(query_layout.missing_header_contigs),
                                    identity_unverified_contigs=[n for n, a, b in zip(reference_layout.names, reference_layout.md5, query_layout.md5) if not a or not b],
                                    reference_simulations_reused=True, reference_geometry_loaded=False,
                                    statistics_seconds=time.perf_counter() - tick)
                            with timing.measure('query_output'):
                                target = directory / 'comparisons' / sample_folder(row['sample']) / ('depth-' + depth_key(depth))
                                target.mkdir(parents=True)
                                atomic_json(target / 'concordance.json', result)
                            entry.update(status='READY', jaccard=result['jaccard']['similarity'],
                                pearson=result['pearson'].get('similarity'), jaccard_p=result['jaccard']['p_value'],
                                pearson_p=result['pearson'].get('p_value'), jaccard_accepted=result['jaccard']['accepted'],
                                pearson_accepted=result['pearson'].get('accepted'), test_actual_depth=draw['actual_depth'],
                                result_path=str(target / 'concordance.json'))
                            compared += 1
                            del matrix
                        if entry['status'] != 'READY':
                            gaps += 1
                    except Exception as error:
                        failed += 1
                        entry.update(status='ERROR', error=str(error))
                    with timing.measure('query_output'):
                        writer.writerow(entry)
                        if (compared + gaps + failed) % config['reference_batch_size'] == 0:
                            table.flush()
                            atomic_json(directory / 'progress.json', dict(state='RUNNING',
                                ready_comparisons=compared, depth_gaps=gaps, failed_comparisons=failed,
                                reference_sample=row['sample'], depth=depth, updated=time.time()))
                if not depth_eligible:
                    gaps += 1
                    writer.writerow(dict(query_sample=query['sample'], gcf=query['gcf'], depth=depth,
                        status='NO_ELIGIBLE_REFERENCE', error='No reference remains after sample filters and self/source exclusion'))
        with timing.measure('query_sha256'):
            if file_hash(query['alignment']) != source_hash or file_stat(query['alignment']) != source_stat:
                raise ValueError('Query changed during comparison')
        after = resource.getrusage(resource.RUSAGE_SELF)
        cache_after = CACHE.stats()
        report = dict(query=query, state='FAILED' if failed else 'SKIPPED' if not compared else 'COMPLETE_WITH_GAPS' if gaps else 'COMPLETE',
            reason=('NO_ELIGIBLE_REFERENCE' if not eligible else 'NO_COMPARABLE_DEPTH') if not failed and not compared else '',
            failed_comparisons=failed, ready_comparisons=compared, depth_gaps=gaps, excluded_reference_depths=excluded,
            query_sha256=source_hash, query_source_stat=source_stat, extraction_calls=extraction_calls,
            extraction_seconds=extraction_seconds, query_pool_numeric_bytes=pool.nbytes,
            query_full_depth=pool.total_bases / pool.genome_len, reference_geometry_loads=0,
            test_draws=len(draws), seed=config['seed'], rng_policy=RNG_POLICY, policies=policies,
            query_missing_header_contigs=list(query_layout.missing_header_contigs),
            cache={k: cache_after[k] - cache_before[k] for k in ('hits', 'misses', 'evictions')},
            worker_pid=os.getpid(), stage_timing=timing.snapshot(),
            cache_after=cache_after, elapsed_seconds=time.perf_counter() - started,
            user_cpu_seconds=after.ru_utime - before.ru_utime, system_cpu_seconds=after.ru_stime - before.ru_stime,
            process_peak_rss_bytes=int(after.ru_maxrss if sys.platform == 'darwin' else after.ru_maxrss * 1024),
            rss_scope='worker cumulative peak, not exclusive query allocation', statistics_calibrated=False)
        atomic_json(directory / 'provenance.json', report)
        atomic_json(directory / 'progress.json', dict(state=report['state'], ready_comparisons=compared,
                    depth_gaps=gaps, failed_comparisons=failed, updated=time.time()))
        success = report['state'].startswith('COMPLETE') or report['state']=='SKIPPED'
        (status / 'RUNNING').unlink(); (status / 'exit_code.txt').write_text('0\n' if success else '1\n')
        (status / ('SKIPPED' if report['state']=='SKIPPED' else 'SUCCESS' if success else 'FAILED')).touch()
        return dict(query=query, state=report['state'], reason=report['reason'], report=str(directory / 'provenance.json'),
                    summary=str(directory / 'summary.tsv'), exit_code=0 if success else 1)
    except BaseException as error:
        (directory / 'error.txt').write_text(traceback.format_exc())
        (status / 'RUNNING').unlink(missing_ok=True); (status / 'FAILED').touch(); (status / 'exit_code.txt').write_text('1\n')
        atomic_json(directory / 'provenance.json', dict(query=query, state='FAILED', error=str(error), extraction_calls=extraction_calls,
            worker_pid=os.getpid(), stage_timing=timing.snapshot()))
        if isinstance(error, (KeyboardInterrupt, SystemExit)):
            raise
        return dict(query=query, state='FAILED', report=str(directory / 'provenance.json'),
                    summary=str(directory / 'summary.tsv'), exit_code=1, error=str(error))


def timed_references(gcf, depth, batch_size, timing):
    iterator = iter(INDEX.references(gcf, depth, batch_size))
    while True:
        with timing.measure('reference_sql'):
            row = next(iterator, None)
        if row is None:
            return
        yield row


def query_tasks(queries, scheduling):
    if scheduling == 'query':
        return [[(number, query)] for number, query in enumerate(queries)]
    groups = {}
    for number, query in enumerate(queries):
        groups.setdefault(query['gcf'], []).append((number, query))
    return sorted(groups.values(), key=lambda group: (-len(group), group[0][1]['gcf']))


def perform_group(group, config):
    return [(number, perform_query(query, config)) for number, query in group]


def run_indexed(args, depths):
    if (args.query_jobs < 1 or args.reference_batch_size < 1 or not np.isfinite(args.reference_cache_mb)
            or args.reference_cache_mb < 0 or args.workers != 1):
        raise ValueError('Invalid indexed budgets; --query-jobs controls concurrency and --workers must remain 1')
    if args.coverage_backend != 'cpu' or args.statistics_backend != 'cpu' or args.extraction_python:
        raise ValueError('SQLite batch mode currently uses the CPU pysam runtime without an extraction bridge')
    if args.reference or args.reference_index:
        raise ValueError('SQLite mode derives full FAI geometry from indexed metadata')
    if len(set(depths)) != len(depths):
        raise ValueError('Duplicate depths')
    from bin.reference_index import depth_key as canonical
    for depth in depths:
        canonical(depth)
    index_path = args.reference_sqlite.resolve()
    index_stat, index_hash = file_stat(index_path), digest(index_path)
    index = ReferenceIndex(index_path, args.database_root)
    try:
        queries = read_queries(args)
        policies = resolve_policies(args, index.metadata['config']['policies'])
        config = dict(index=str(index_path), database_root=str(index.root) if index.root else None, outdir=str(args.outdir.resolve()),
            depths=depths, seed=args.seed, bins=index.metadata['config'].get('bins'), policies=policies,
            reference_batch_size=args.reference_batch_size,
            reference_samples=set(args.reference_samples.split(',')) if args.reference_samples else set(),
            include_self=bool(args.include_self), alpha=args.alpha, similarity_block_size=args.similarity_block_size)
        native=index.metadata.get('storage_mode')=='native'
        config.update(native=native,window_size=index.metadata['config'].get('window_size'),
            bin_count=index.metadata['config'].get('bin_count'),
            replicates=args.replicates,sampler=args.depth_sampler,reference_limit=args.reference_limit)
        if native:
            from bin.depth_band_concordance import configure
            config.update(configure(args, depths))
        index_metadata = index.metadata
    finally:
        index.close()
    args.outdir.mkdir(parents=True, exist_ok=False)
    status = args.outdir / 'status'; status.mkdir(); (status / 'RUNNING').touch()
    atomic_json(args.outdir / 'queries.json', queries)
    started = time.perf_counter(); outcomes = {}
    budget = int(args.reference_cache_mb * 1024 * 1024)
    scheduling = getattr(args, 'query_scheduling', 'gcf')
    tasks = query_tasks(queries, scheduling)
    active_workers = min(args.query_jobs, len(tasks))
    try:
        if active_workers == 1:
            initialize(index_path, config['database_root'], budget,native)
            try:
                for group in tasks:
                    outcomes.update(perform_group(group, config))
            finally:
                CACHE.close(); INDEX.close()
                # Do not reuse the closed singleton in another in-process invocation.
                globals()['CACHE'] = globals()['INDEX'] = None
        else:
            with ProcessPoolExecutor(max_workers=active_workers, mp_context=multiprocessing.get_context('spawn'),
                    initializer=initialize, initargs=(index_path, config['database_root'], budget,native)) as executor:
                pending = {}; iterator = iter(tasks); exhausted = False
                while pending or not exhausted:
                    while len(pending) < 2 * active_workers and not exhausted:
                        try:
                            group = next(iterator)
                        except StopIteration:
                            exhausted = True; break
                        pending[executor.submit(perform_group, group, config)] = True
                    if pending:
                        done, _ = wait(pending, return_when=FIRST_COMPLETED)
                        for future in done:
                            pending.pop(future)
                            outcomes.update(future.result())
        if file_stat(index_path) != index_stat or digest(index_path) != index_hash:
            raise ValueError('SQLite index changed during query run')
        with (args.outdir / 'summary.tsv').open('w', newline='') as handle:
            writer = csv.DictWriter(handle, COLUMNS, delimiter='\t'); writer.writeheader()
            for number in range(len(queries)):
                path = Path(outcomes[number]['summary'])
                if path.exists():
                    with path.open() as source:
                        writer.writerows(csv.DictReader(source, delimiter='\t'))
        with (args.outdir / 'stage_timing.tsv').open('w', newline='') as handle:
            columns = ['query_sample', 'gcf', 'worker_pid', 'stage', 'wall_seconds',
                       'user_cpu_seconds', 'system_cpu_seconds', 'calls']
            writer = csv.DictWriter(handle, columns, delimiter='\t'); writer.writeheader()
            for number in range(len(queries)):
                value = json.loads(Path(outcomes[number]['report']).read_text())
                for stage, metrics in value.get('stage_timing', {}).get('stages', {}).items():
                    writer.writerow(dict(query_sample=queries[number]['sample'], gcf=queries[number]['gcf'],
                        worker_pid=value.get('worker_pid'), stage=stage, **metrics))
        failed = any(r['exit_code'] for r in outcomes.values())
        report = dict(schema='indexed-coverage-concordance-1', state='FAILED' if failed else
            'COMPLETE_WITH_GAPS' if any(r['state'] in ('COMPLETE_WITH_GAPS','SKIPPED') for r in outcomes.values()) else 'COMPLETE',
            ready_query_count=sum(r['state'].startswith('COMPLETE') for r in outcomes.values()),
            skipped_query_count=sum(r['state']=='SKIPPED' for r in outcomes.values()),
            failed_query_count=sum(bool(r['exit_code']) for r in outcomes.values()),
            query_count=len(queries), outcomes=[outcomes[i] for i in range(len(queries))],
            reference_index=str(index_path), index_sha256=index_hash, plan_sha256=index_metadata['plan_sha256'],
            database_root=config['database_root'], query_jobs=args.query_jobs, native_threads_per_worker=1,
            query_scheduling=scheduling, active_workers=active_workers, gcf_task_count=len(tasks),
            cache_budget_per_worker_bytes=budget, reference_batch_size=args.reference_batch_size,
            elapsed_seconds=time.perf_counter() - started, statistics_calibrated=False,
            reference_storage_unchanged=True, query_configuration={k:v for k,v in config.items() if k!='reference_samples'},
            reference_samples=sorted(config['reference_samples']), python=sys.version, numpy=np.__version__)
        atomic_json(args.outdir / 'provenance.json', report)
        (status / 'RUNNING').unlink(); (status / ('FAILED' if failed else 'SUCCESS')).touch()
        (status / 'exit_code.txt').write_text('1\n' if failed else '0\n')
        return 1 if failed else 0
    except BaseException:
        (args.outdir / 'error.txt').write_text(traceback.format_exc())
        (status / 'RUNNING').unlink(missing_ok=True); (status / 'FAILED').touch(); (status / 'exit_code.txt').write_text('1\n')
        raise
