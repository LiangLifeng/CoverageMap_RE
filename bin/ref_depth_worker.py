"""Per-GCF, per-sample depth simulation and verified matrix reading."""
from pathlib import Path
import json
import os
import resource
import sys
import time
import traceback
import shutil
from contextlib import nullcontext
from bin.stage_timing import StageTiming

import numpy as np

from bin.ref_depth_catalog import (atomic_json, completed_sample, digest, exclusive,
                                   file_stat, load_plan, safe_file, validate_cram_index, validate_sample_payloads,
                                   required_sample_payloads)
from bin.coverage_core import (canonical_reference, extract_pool, layout_metadata, merge_references,
                               reference_from_layout, reproject_pool, window_edges)
from bin.gcf_windows import GCF_CONCAT, window_policy
from bin.concordance_core import RNG_POLICY, simulate_grid
from bin.reference_store import build_bundle, load_reference_bundle, payload_sha256, retry_payload_io


def resources():
    usage = resource.getrusage(resource.RUSAGE_SELF)
    return dict(user=usage.ru_utime, system=usage.ru_stime,
                peak_rss=int(usage.ru_maxrss if sys.platform == 'darwin' else usage.ru_maxrss * 1024))


def merge_layout(expected, observed):
    if expected is None:
        return observed
    refs = [reference_from_layout(item) for item in (expected, observed)]
    combined = merge_references(*refs)
    return layout_metadata(combined)


def reuse_geometry(source, gcf, sample, source_hash, policies):
    root = Path(source['database'])
    if digest(root/'database.json') != source['plan_sha256']:
        raise ValueError('Source geometry database plan changed')
    parent = root/gcf/'samples'/sample['folder']
    saved = completed_sample(parent, source['plan_sha256'], verify=True)
    if saved is None:
        raise ValueError('No verified reusable geometry for this sample; re-extract into a new database')
    attempt, record = saved
    if ((record['gcf'], record['sample'], record['sources'][0]['sha256'], record['policies'])
            != (gcf, sample['sample'], source_hash, policies)):
        raise ValueError('Reusable geometry source identity or read policies differ')
    pool, original, bundle = load_reference_bundle(attempt/'geometry', workers=1)
    if (record['layout'] != bundle['layout'] or record['policies'] != bundle['policies']
            or pool.total_bases != record['aligned_bases']):
        raise ValueError('Reusable geometry and original result differ')
    canonical = canonical_reference(original)
    if hasattr(pool,'alignment_geometry'):
        if original != canonical or original.sort_policy != canonical.sort_policy:
            raise ValueError('Complete geometry layout needs re-extraction, not a legacy coordinate-only migration')
        rebased = pool
    else:
        rebased = reproject_pool(pool, original, canonical)
    return rebased, canonical, dict(
        manifest=str(attempt/'result.json'), sha256=digest(attempt/'result.json'),
        source_layout=layout_metadata(original))


def build_sample(database, plan, gcf, sample, expected_layout=None, extractor=extract_pool):
    group, config = plan['groups'][gcf], plan['config']
    parent = database / gcf / 'samples' / sample['folder']
    previous = sorted(parent.glob('attempt-*'))
    number = max([int(p.name.split('-')[-1]) for p in previous] + [0]) + 1
    attempt = parent / f'attempt-{number:06d}'
    (attempt / 'status').mkdir(parents=True)
    (attempt / 'status/RUNNING').touch()
    started, before = time.perf_counter(), resources()
    timing = StageTiming()
    try:
        if file_stat(sample['cram']) != sample['source_stat']:
            raise ValueError('CRAM stat changed since plan; create a new plan')
        if digest(group['fai']) != group['fai_sha256']:
            raise ValueError('FAI changed since plan')
        with timing.measure('source_sha256'):
            validate_cram_index(sample)
            source_hash = digest(sample['cram'])
            if sample.get('source_sha256') is not None and source_hash != sample['source_sha256']:
                raise ValueError('CRAM content changed since immutable plan')
        policy = window_policy(config)
        ordering = {'sort_contigs': True} if policy == GCF_CONCAT else {}
        migrated = {}
        if config.get('reuse_geometry'):
            with timing.measure('geometry_reprojection'):
                pool, ref, reused = reuse_geometry(config['reuse_geometry'], gcf, sample,
                                                  source_hash, config['policies'])
                migrated['reused_geometry'] = reused
        else:
            with timing.measure('alignment_extract'):
                capture = {'capture_geometry': attempt/'alignment-capture', 'geometry_namespace': source_hash} if config.get('capture_alignment') else {}
                pool, ref = extractor(sample['cram'], reference_index=group['fai'],
                                      **ordering, **config['policies'], **capture)
        layout = layout_metadata(ref)
        merge_references(reference_from_layout(group['layout']), ref)
        merge_layout(expected_layout, layout)
        available = [d for d in config.get('depths', []) if d * ref.length <= pool.total_bases]
        if len(available) != len(config.get('depths', [])) and config['insufficient_depth'] == 'error':
            raise ValueError('Requested depth exceeds sample aligned mass')
        sources = [dict(path=sample['cram'], sha256=source_hash),
                   dict(path=group['fai'], sha256=group['fai_sha256'])]
        if sample.get('cram_index'):
            with timing.measure('source_index_archive'):
                index = sample['cram_index']
                shutil.copyfile(index['path'], attempt / 'source.cram.crai')
                if digest(attempt / 'source.cram.crai') != index['sha256']:
                    raise ValueError('Archived CRAM index differs from plan')
                sources.append(dict(path=index['path'], sha256=index['sha256']))
        with timing.measure('geometry_storage'):
            build_bundle(pool, ref, attempt / 'geometry', sample['sample'], config['policies'], sources,
                         workers=config['cpus'], storage_format=config['geometry_format'])
            if config.get('capture_alignment'):
                from bin.alignment_geometry import Geometry
                pool.alignment_geometry = Geometry(attempt/'geometry/alignment')
                if (attempt/'alignment-capture').exists():
                    shutil.rmtree(attempt/'alignment-capture')
        native = config.get('storage_mode') == 'native'
        native_info = {}
        rows = {}
        if native:
            from bin.native_reference import save_native
            with timing.measure('native_window_coverage'):
                edges, native_info = save_native(pool, ref, attempt, config['window_size'],
                                                count=config.get('bin_count'), policy=policy)
        else:
            edges = window_edges(ref.length, config['bins'])
            np.save(attempt / 'edges.npy', edges, allow_pickle=False)
            coverage = attempt / 'coverage'
            coverage.mkdir()
            rows = {d: dict(depth=d, status='UNSUPPORTED_DEPTH', matrix='') for d in config.get('depths', [])}
            converted = []
            def complete(index, depth, directory, matrix, summaries, calibration, seconds):
                with timing.measure('matrix_write_compress_verify'):
                    path = directory / 'simulations.npy'
                    if config['matrix_format'] == 'npz':
                        compressed = directory / 'simulations.npz'
                        with timing.measure('matrix_compress_write'):
                            np.savez_compressed(compressed, counts=matrix)
                        with timing.measure('matrix_readback_verify'):
                            with np.load(compressed, allow_pickle=False) as saved:
                                if not np.array_equal(matrix, saved['counts']):
                                    raise AssertionError('Compressed matrix readback differs')
                        converted.append(path)
                        path = compressed
                    else:
                        with path.open('rb') as handle:
                            if not np.array_equal(matrix, np.load(handle, allow_pickle=False)):
                                raise AssertionError('NPY matrix readback differs')
                metadata = dict(depth=depth, reference_replicates=summaries,
                                sampler=calibration, worker_simulation_seconds_sum=seconds)
                atomic_json(directory / 'reference.json', metadata)
                rows[depth] = dict(depth=depth, status='READY', matrix=str(path.relative_to(attempt)),
                                   reference_metadata=str((directory / 'reference.json').relative_to(attempt)))
            if available:
                with timing.measure('sampling_counting_and_dispatch'):
                    simulate_grid(pool, edges, available, config['replicates'], config['seed'], config['cpus'],
                                  'thread', config['batch_size'], config['sampler'], coverage, complete, timing=timing)
            # simulate_grid has flushed and closed all maps before removing our own converted scratch NPYs.
            for path in converted:
                path.unlink()
        with timing.measure('source_sha256'):
            validate_cram_index(sample)
            if digest(sample['cram']) != source_hash or file_stat(sample['cram']) != sample['source_stat']:
                raise ValueError('CRAM changed while building sample')
            if digest(group['fai']) != group['fai_sha256']:
                raise ValueError('FAI changed while building sample')
        after = resources()
        missing_header = set(ref.missing_header_contigs)
        result = dict(schema='gcf-native-sample-2' if native else 'gcf-depth-sample-1',
            storage_mode='native' if native else 'legacy-depths', native=native_info, plan_sha256=(database / 'PLAN.sha256').read_text().strip(),
            gcf=gcf, sample=sample['sample'], sources=sources, source_stat=sample['source_stat'],
            cram_index=sample.get('cram_index'), cram_index_archive='source.cram.crai' if sample.get('cram_index') else None,
            layout=layout, genome_len=ref.length, units=len(pool), blocks=len(pool.starts),
            aligned_bases=pool.total_bases, full_mean_depth=pool.total_bases / ref.length,
            empty_eligible_pool=len(pool) == 0, missing_m5=[n for n, m in zip(ref.names, ref.md5) if not m],
            missing_header_contigs=list(ref.missing_header_contigs),
            missing_header_contig_count=len(ref.missing_header_contigs),
            missing_header_bases=sum(length for name, length in zip(ref.names, ref.lengths) if name in missing_header),
            reference_layout_policy='full-fai-subset-header-zero-coverage-v1',
            policies=config['policies'], bins=len(edges)-1,
            window_size=config.get('window_size') if native else None, **migrated,
            **({'bin_count': config.get('bin_count')} if native and policy == GCF_CONCAT else {}),
            elapsed_seconds=time.perf_counter() - started,
            user_cpu_seconds=after['user'] - before['user'], system_cpu_seconds=after['system'] - before['system'],
            process_peak_rss_bytes=after['peak_rss'],
            rss_scope='worker process cumulative high-water across tasks, not exclusive sample memory')
        if not native:
            result.update(rng_policy=RNG_POLICY, seed=config['seed'],replicates=config['replicates'],
                          matrix_format=config['matrix_format'],depths=[rows[d] for d in config['depths']])
        result['capabilities'] = json.loads((attempt/'geometry/manifest.json').read_text()).get('capabilities', {})
        atomic_json(attempt / 'result.json', result)
        with timing.measure('final_payload_sha256'):
            required = retry_payload_io(lambda: required_sample_payloads(attempt, result), attempt)
            # Discovery retains optional artifacts; required files never depend on rglob/is_file.
            discovered = retry_payload_io(lambda: {
                p.relative_to(attempt).as_posix() for p in attempt.rglob('*')
                if p.is_file() and 'status' not in p.relative_to(attempt).parts}, attempt)
            names = (required | discovered) - {'timing.json', 'checksums.json'}
            checks = {name: payload_sha256(attempt, name) for name in sorted(names)}
        atomic_json(attempt / 'timing.json', timing.snapshot())
        checks['timing.json'] = payload_sha256(attempt, 'timing.json')
        atomic_json(attempt / 'checksums.json', checks)
        retry_payload_io(lambda: validate_sample_payloads(attempt, result, checks, verify=False), attempt)
        (attempt / 'status/RUNNING').unlink()
        (attempt / 'status/exit_code.txt').write_text('0\n')
        (attempt / 'status/SUCCESS').touch()
        return attempt, result
    except BaseException:
        (attempt / 'error.txt').write_text(traceback.format_exc())
        (attempt / 'status/RUNNING').unlink(missing_ok=True)
        (attempt / 'status/exit_code.txt').write_text('1\n')
        (attempt / 'status/FAILED').touch()
        raise


def execute_sample(database, plan, gcf, sample, expected_layout=None, extractor=extract_pool, locked=False):
    """One exclusive attempt allocator shared by GCF and sample workers."""
    parent = database/gcf/'samples'/sample['folder']
    with nullcontext() if locked else exclusive(parent/'.sample.lock'):
        if file_stat(sample['cram']) != sample['source_stat']:
            raise ValueError('Source stat changed; use a new plan')
        if digest(plan['groups'][gcf]['fai']) != plan['groups'][gcf]['fai_sha256']:
            raise ValueError('FAI changed or unavailable since plan')
        validate_cram_index(sample)
        try:
            previous = completed_sample(parent, (database/'PLAN.sha256').read_text().strip(), verify=True)
        except (ValueError, KeyError, OSError):
            previous = None
        if previous:
            attempt, result = previous
        else:
            attempt, result = build_sample(database, plan, gcf, sample, expected_layout, extractor)
        return attempt, result, previous is not None


def run_sample(database, plan, gcf, sample, extractor=extract_pool):
    """Sample task; its persistent process may subsequently execute any GCF."""
    with exclusive(database/gcf/'samples'/sample['folder']/'.sample.lock'):
        return run_locked_sample(database,plan,gcf,sample,extractor)


def run_locked_sample(database, plan, gcf, sample, extractor):
    parent=database/gcf/'samples'/sample['folder']
    started=time.perf_counter()
    value=dict(gcf=gcf,sample=sample['sample'],state='RUNNING',worker_pid=os.getpid(),updated=time.time())
    atomic_json(parent/'progress.json',value)
    print(time.strftime('[%Y-%m-%d %H:%M:%S] ')+f"SAMPLE START {gcf}/{sample['sample']}",flush=True)
    try:
        attempt,result,reused=execute_sample(database,plan,gcf,sample,extractor=extractor,locked=True)
        value.update(state='COMPLETE',reused=reused,manifest=str(attempt/'result.json'),exit_code=0)
    except (KeyboardInterrupt,InterruptedError) as error:
        value.update(state='INTERRUPTED',error=str(error),exit_code=1)
    except Exception as error:
        value.update(state='FAILED',error=f'{type(error).__name__}: {error}',exit_code=1)
        traceback.print_exc()
    value.update(elapsed_seconds=time.perf_counter()-started,updated=time.time())
    atomic_json(parent/'progress.json',value)
    print(time.strftime('[%Y-%m-%d %H:%M:%S] ')+f"SAMPLE END {gcf}/{sample['sample']} state={value['state']} elapsed={value['elapsed_seconds']:.2f}s",flush=True)
    return value


def run_gcf(database, gcf, extractor=extract_pool):
    database, plan = load_plan(database, check_code=True)
    if gcf not in plan['groups']:
        raise ValueError('GCF is not in the plan')
    directory = database / gcf
    fingerprint = (database / 'PLAN.sha256').read_text().strip()
    with exclusive(directory / '.worker.lock'):
        failures, finished, reused = [], 0, 0
        expected_layout = None
        started = time.perf_counter()
        def log(message):
            print(time.strftime('[%Y-%m-%d %H:%M:%S] ')+gcf+' '+message, flush=True)
        log(f"START samples={len(plan['groups'][gcf]['samples'])}")
        try:
            reference_valid = digest(plan['groups'][gcf]['fai']) == plan['groups'][gcf]['fai_sha256']
        except OSError:
            reference_valid = False
        if not reference_valid:
            log('FAILED: FAI changed or unavailable since plan')
            atomic_json(directory / 'progress.json', dict(state='FAILED', samples_finished=0,
                failed_samples=len(plan['groups'][gcf]['samples']),
                failures=[dict(sample=r['sample'], error='FAI changed or unavailable since plan') for r in plan['groups'][gcf]['samples']]))
            return 1
        for sample in plan['groups'][gcf]['samples']:
            log(f"SAMPLE START {sample['sample']} processed={finished+len(failures)}/{len(plan['groups'][gcf]['samples'])}")
            atomic_json(directory / 'progress.json', dict(state='RUNNING', current_sample=sample['sample'],
                samples_finished=finished, samples_reused=reused, failed_samples=len(failures),
                samples_total=len(plan['groups'][gcf]['samples']), updated=time.time()))
            try:
                attempt,result,was_reused=execute_sample(database,plan,gcf,sample,expected_layout,extractor)
                reused += int(was_reused)
                expected_layout = merge_layout(expected_layout, result['layout'])
                finished += 1
                log(f"SAMPLE OK {sample['sample']} ok={finished} failed={len(failures)} reused={reused}")
            except (KeyboardInterrupt, InterruptedError):
                log(f"INTERRUPTED sample={sample['sample']}")
                atomic_json(directory / 'progress.json', dict(state='INTERRUPTED', samples_finished=finished,
                    failed_samples=len(failures), updated=time.time()))
                raise
            except Exception as error:
                failures.append(dict(sample=sample['sample'], error=str(error)))
                log(f"SAMPLE FAILED {sample['sample']}: {type(error).__name__}: {error}")
            atomic_json(directory / 'progress.json', dict(state='RUNNING',
                samples_finished=finished, samples_reused=reused, failed_samples=len(failures),
                samples_total=len(plan['groups'][gcf]['samples']), failures=failures, updated=time.time()))
        state = 'FAILED' if failures else 'COMPLETE'
        atomic_json(directory / 'progress.json', dict(state=state, samples_finished=finished,
            samples_reused=reused, samples_total=len(plan['groups'][gcf]['samples']),
            failed_samples=len(failures), failures=failures, elapsed_seconds=time.perf_counter() - started,
            process_peak_rss_bytes=resources()['peak_rss'], updated=time.time()))
        log(f"END state={state} ok={finished} failed={len(failures)} elapsed={time.perf_counter()-started:.1f}s")
        return 1 if failures else 0


def load_depth(sample_manifest, depth, timing=None):
    """Return (N-by-B count matrix, shared edges, metadata) for an exact stored depth."""
    manifest = Path(sample_manifest).resolve()
    root = manifest.parent
    if (not (root / 'status/SUCCESS').is_file() or (root / 'status/FAILED').exists()
            or (root / 'status/RUNNING').exists()
            or (root / 'status/exit_code.txt').read_text().strip() != '0'):
        raise ValueError('Sample build is not successful')
    checks = json.loads((root / 'checksums.json').read_text())
    if checks.get('result.json') != digest(manifest):
        raise ValueError('Sample manifest checksum mismatch')
    result = json.loads(manifest.read_text())
    rows = [row for row in result['depths'] if row['depth'] == float(depth)]
    if len(rows) != 1 or rows[0]['status'] != 'READY':
        raise ValueError('Exact depth is unavailable; no interpolation is performed')
    row = rows[0]
    with timing.measure('reference_payload_sha256') if timing else nullcontext():
        for name in (row['matrix'], 'edges.npy', row['reference_metadata']):
            if digest(safe_file(root, name)) != checks.get(name):
                raise ValueError(f'Payload checksum mismatch: {name}')
    with timing.measure('reference_load_decompress') if timing else nullcontext():
        path = safe_file(root, row['matrix'])
        if path.suffix == '.npz':
            with np.load(path, allow_pickle=False) as archive:
                matrix = archive['counts']
        else:
            matrix = np.load(path, mmap_mode='r', allow_pickle=False)
        edges = np.load(root / 'edges.npy', allow_pickle=False)
    if matrix.shape != (result['replicates'], result['bins']) or matrix.dtype.kind not in 'iu':
        raise ValueError('Invalid matrix dimensions/dtype')
    if len(edges) != result['bins'] + 1 or edges[0] != 0 or edges[-1] != result['genome_len'] or np.any(np.diff(edges) <= 0):
        raise ValueError('Invalid window geometry')
    if np.any(matrix < 0):
        raise ValueError('Negative coverage')
    metadata = json.loads((root / row['reference_metadata']).read_text())
    if not np.array_equal(matrix.sum(axis=1), [r['aligned_bases'] for r in metadata['reference_replicates']]):
        raise ValueError('Matrix mass differs from recorded replicates')
    return matrix, edges, dict(sample=result, reference=metadata)
