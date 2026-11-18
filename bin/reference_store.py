#!/usr/bin/env python3
"""Build and verify separate reusable reference geometry bundles."""
from __future__ import annotations

import argparse
from collections import deque
from concurrent.futures import ThreadPoolExecutor
import errno
import hashlib
import json
from pathlib import Path
import resource
import shutil
import stat
import subprocess
import sys
import tempfile
import time
import zlib

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from bin.coverage_core import (BlockPool, Reference, extract_pool, layout_metadata,
                               reference_from_layout, select_units, window_edges, window_sums)

SCHEMA = 'coverage-reference-1'
ARRAYS = ('starts', 'ends', 'indptr')
PAYLOAD_RETRY_DELAYS = (.2, .5, 1., 2.)
TRANSIENT_PAYLOAD_ERRNOS = {errno.ENOENT, errno.ESTALE, errno.EAGAIN, errno.ETIMEDOUT, errno.EIO}


def retry_payload_io(operation, path):
    """Retry transient filesystem access only; never retry validation failures."""
    for attempt in range(len(PAYLOAD_RETRY_DELAYS) + 1):
        try:
            return operation()
        except OSError as error:
            if error.errno not in TRANSIENT_PAYLOAD_ERRNOS or attempt == len(PAYLOAD_RETRY_DELAYS):
                raise
            delay = PAYLOAD_RETRY_DELAYS[attempt]
            print(f'PAYLOAD RETRY path={path} errno={error.errno} '
                  f'attempt={attempt + 1}/{len(PAYLOAD_RETRY_DELAYS) + 1} delay={delay}s',
                  file=sys.stderr, flush=True)
            time.sleep(delay)


def payload_file(directory, name):
    """Preserve path safety while allowing missing/stale stat errors to propagate."""
    directory, relative = Path(directory), Path(name)
    if relative.is_absolute() or '..' in relative.parts or not relative.parts:
        raise ValueError('Invalid bundle path')
    path = directory / relative
    if not path.resolve().is_relative_to(directory.resolve()):
        raise ValueError('Escaping bundle file')
    info = path.stat()
    if not stat.S_ISREG(info.st_mode) or info.st_size == 0:
        raise ValueError(f'Payload is not a nonempty regular file: {path}')
    return path


def payload_sha256(directory, name):
    return retry_payload_io(lambda: sha256(payload_file(directory, name)), Path(directory) / name)


def storage_payload_names(storage):
    extra = set(storage.get('alignment_geometry_files', []))
    if storage['codec'] == 'packed':
        names = [c['file'] for c in storage['chunks']]
        if len(names) != len(set(names)):
            raise ValueError('Duplicate storage payload')
        return set(names) | extra
    if storage['codec'] == 'npz':
        return {'pool.npz'} | extra
    if storage['codec'] == 'npy':
        return {f'{name}.npy' for name in ARRAYS} | extra
    raise ValueError('Unknown storage codec')


def validate_layout(layout):
    if (not layout.names or len(set(layout.names)) != len(layout.names)
            or any(not isinstance(n, str) or not n for n in layout.names)
            or len(layout.names) != len(layout.lengths) or len(layout.md5) != len(layout.names)
            or any(not isinstance(n, int) or n <= 0 for n in layout.lengths)):
        raise ValueError('Invalid reference layout')
    if (len(set(layout.missing_header_contigs)) != len(layout.missing_header_contigs)
            or not set(layout.missing_header_contigs) <= set(layout.names)):
        raise ValueError('Invalid missing header contigs')


def write_json(path, value):
    path.write_text(json.dumps(value, sort_keys=True, indent=2, allow_nan=False) + '\n')


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def bounded_map(function, jobs, workers):
    """Ordered results, no more than 2*workers pending jobs/buffers."""
    if workers < 1:
        raise ValueError('workers must be positive')
    with ThreadPoolExecutor(max_workers=workers) as executor:
        pending = deque()
        for job in jobs:
            pending.append(executor.submit(function, job))
            if len(pending) >= 2 * workers:
                yield pending.popleft().result()
        while pending:
            yield pending.popleft().result()


def safe_file(directory, name):
    relative = Path(name)
    if relative.is_absolute() or '..' in relative.parts or not relative.parts:
        raise ValueError('Invalid bundle path')
    path = directory / relative
    if not path.resolve().is_relative_to(directory.resolve()) or not path.is_file():
        raise ValueError('Missing or escaping bundle file')
    return path


def inventory(directory, names=None):
    if names is not None:
        def inspect(name):
            p = payload_file(directory, name)
            info = p.stat()
            return dict(bytes=info.st_size, allocated_bytes=getattr(info, 'st_blocks', 0) * 512)
        return {name: retry_payload_io(lambda: inspect(name), directory / name) for name in sorted(names)}
    return {p.relative_to(directory).as_posix(): dict(
                bytes=p.stat().st_size,
                allocated_bytes=getattr(p.stat(), 'st_blocks', 0) * 512)
            for p in sorted(directory.rglob('*')) if p.is_file()}


def smallest_uint(values):
    maximum = int(values.max()) if values.size else 0
    for width in (1, 2, 4, 8):
        if maximum <= (1 << (8 * width)) - 1:
            return np.asarray(values, dtype=f'<u{width}')
    raise OverflowError('Unsigned value exceeds uint64')


def encode_chunk(job):
    directory, pool, begin, end, level = job
    left, right = int(pool.indptr[begin]), int(pool.indptr[end])
    arrays = [smallest_uint(pool.starts[left:right]),
              smallest_uint(pool.ends[left:right] - pool.starts[left:right]),
              smallest_uint(np.diff(pool.indptr[begin:end + 1]))]
    # Byte shuffle separates high/low bytes; zlib runs outside the Python GIL.
    shuffled = [a.view(np.uint8).reshape(-1, a.dtype.itemsize).T.tobytes() for a in arrays]
    path = directory / f'chunk-{begin:012d}.zlib'
    path.write_bytes(zlib.compress(b''.join(shuffled), level))
    return dict(file=path.name, begin=begin, end=end, blocks=right - left,
                arrays=[dict(dtype=a.dtype.str, count=len(a)) for a in arrays])


def decode_chunk(job):
    directory, chunk, *access = job
    get_file = payload_file if access and access[0] else safe_file
    descriptors = chunk['arrays']
    if len(descriptors) != 3:
        raise ValueError('Invalid chunk descriptors')
    sizes = []
    for spec in descriptors:
        if spec['dtype'] not in ('|u1', '<u2', '<u4', '<u8') or spec['count'] < 0:
            raise ValueError('Invalid unsigned chunk dtype/count')
        sizes.append(spec['count'] * np.dtype(spec['dtype']).itemsize)
    decoder = zlib.decompressobj()
    raw = decoder.decompress(get_file(directory, chunk['file']).read_bytes(), sum(sizes) + 1)
    if len(raw) != sum(sizes) or not decoder.eof or decoder.unused_data or decoder.unconsumed_tail:
        raise ValueError('Invalid compressed chunk length')
    arrays, offset = [], 0
    for spec, size in zip(descriptors, sizes):
        dtype = np.dtype(spec['dtype'])
        data = np.frombuffer(raw[offset:offset + size], dtype=np.uint8)
        arrays.append(data.reshape(dtype.itemsize, spec['count']).T.copy().reshape(-1).view(dtype))
        offset += size
    starts, lengths, counts = arrays
    starts = starts.astype(np.uint64)
    lengths = lengths.astype(np.uint64)
    if np.any(lengths > np.iinfo(np.uint64).max - starts):
        raise ValueError('Coordinate overflow')
    if (len(starts) != chunk['blocks'] or len(lengths) != len(starts)
            or len(counts) != chunk['end'] - chunk['begin']
            or np.any(counts == 0) or sum(map(int, counts)) != len(starts)):
        raise ValueError('Invalid chunk unit boundaries')
    return starts, starts + lengths, counts.astype(np.uint64)


def save_storage(pool, directory, codec, workers, chunk_units, level=6):
    directory.mkdir()
    storage = dict(codec=codec, genome_len=pool.genome_len, units=len(pool),
                   blocks=len(pool.starts), coordinate_dtype='<u4' if pool.genome_len <= 2**32 - 1 else '<u8')
    if codec == 'packed':
        jobs = ((directory, pool, begin, min(begin + chunk_units, len(pool)), level)
                for begin in range(0, len(pool), chunk_units))
        storage.update(level=level, chunk_units=chunk_units,
                       chunks=list(bounded_map(encode_chunk, jobs, workers)))
    elif codec == 'npz':
        pool.save(directory / 'pool.npz')
    elif codec == 'npy':
        def save(name):
            np.save(directory / f'{name}.npy', getattr(pool, name), allow_pickle=False)
        list(bounded_map(save, ARRAYS, workers))
    else:
        raise ValueError('Unknown storage codec')
    write_json(directory / 'storage.json', storage)
    return storage


def read_storage(directory, workers=4, *, transient_access=False):
    get_file = payload_file if transient_access else safe_file
    storage = json.loads((directory / 'storage.json').read_text())
    codec = storage['codec']
    if codec == 'npz':
        pool = BlockPool.load(get_file(directory, 'pool.npz'))
    elif codec == 'npy':
        pool = BlockPool(*(np.load(get_file(directory, f'{name}.npy'), mmap_mode='r',
                                   allow_pickle=False) for name in ARRAYS), storage['genome_len'])
    elif codec == 'packed':
        units, blocks = storage['units'], storage['blocks']
        if units < 0 or blocks < 0:
            raise ValueError('Negative geometry dimensions')
        dtype = np.dtype(storage['coordinate_dtype'])
        if dtype.str not in ('|u1', '<u2', '<u4', '<u8'):
            raise ValueError('Invalid coordinate dtype')
        starts, ends = np.empty(blocks, dtype=dtype), np.empty(blocks, dtype=dtype)
        pointers = np.empty(units + 1, dtype=np.uint64)
        pointers[0] = 0
        cursor = 0
        for chunk in storage['chunks']:
            if chunk['begin'] != cursor or not cursor < chunk['end'] <= units:
                raise ValueError('Noncontiguous chunks')
            cursor = chunk['end']
        if cursor != units:
            raise ValueError('Missing units')
        jobs = ((directory, chunk, transient_access) for chunk in storage['chunks'])
        block_cursor = 0
        for chunk, (lefts, rights, counts) in zip(storage['chunks'], bounded_map(decode_chunk, jobs, workers)):
            end = block_cursor + len(lefts)
            if end > blocks or np.any(rights > np.iinfo(dtype).max):
                raise ValueError('Chunk bounds overflow')
            starts[block_cursor:end], ends[block_cursor:end] = lefts, rights
            pointers[chunk['begin'] + 1:chunk['end'] + 1] = block_cursor + np.cumsum(counts)
            block_cursor = end
        if block_cursor != blocks:
            raise ValueError('Block count mismatch')
        pool = BlockPool(starts, ends, pointers, storage['genome_len'])
    else:
        raise ValueError('Unknown storage codec')
    if len(pool) != storage['units'] or len(pool.starts) != storage['blocks']:
        raise ValueError('Stored geometry dimensions differ')
    return pool


def validate_replay(original, loaded):
    for name in ARRAYS:
        if not np.array_equal(getattr(original, name), getattr(loaded, name)):
            raise AssertionError(f'Roundtrip differs: {name}')
    if original.genome_len != loaded.genome_len or original.total_bases != loaded.total_bases:
        raise AssertionError('Mass/layout mismatch')
    edges = window_edges(original.genome_len, min(original.genome_len, 1000))
    draws = []
    for fraction in (0., .001, .01, .1):
        target = fraction * original.total_bases
        ids = select_units(original, np.random.default_rng(401), target_bases=target)
        replay = select_units(loaded, np.random.default_rng(401), target_bases=target)
        if not np.array_equal(ids, replay):
            raise AssertionError('Selected units changed')
        values = window_sums(original, ids, edges)
        if not np.array_equal(values, window_sums(loaded, replay, edges)):
            raise AssertionError('Window sums changed')
        draws.append(dict(fraction=fraction, units=len(ids), aligned_bases=int(values.sum()),
                          ids_sha256=hashlib.sha256(ids.astype('<i8').tobytes()).hexdigest()))
    return draws


def load_reference_bundle(directory, workers=4):
    directory = Path(directory)
    if not (directory / 'status/SUCCESS').is_file() or (directory / 'status/FAILED').exists():
        raise ValueError('Reference bundle is incomplete or failed')
    if (directory / 'status/exit_code.txt').read_text().strip() != '0':
        raise ValueError('Nonzero bundle exit status')
    checksums = json.loads((directory / 'checksums.json').read_text())
    if not {'manifest.json', 'storage.json', 'report.json'} <= checksums.keys():
        raise ValueError('Incomplete checksum index')
    def verify(item):
        name, expected = item
        if sha256(safe_file(directory, name)) != expected:
            raise ValueError(f'Checksum mismatch: {name}')
    list(bounded_map(verify, sorted(checksums.items()), workers))
    manifest = json.loads((directory / 'manifest.json').read_text())
    if manifest['schema'] != SCHEMA:
        raise ValueError('Unsupported reference schema')
    layout = reference_from_layout(manifest['layout'], manifest.get('missing_header_contigs', []))
    validate_layout(layout)
    storage = json.loads((directory / 'storage.json').read_text())
    expected = storage_payload_names(storage)
    if not expected <= checksums.keys():
        raise ValueError('Unchecksummed storage payload')
    pool = read_storage(directory, workers)
    if storage.get('alignment_geometry_files'):
        from bin.alignment_geometry import Geometry
        pool.alignment_geometry = Geometry(directory/'alignment')
        replay = pool.alignment_geometry.project('read', manifest['capabilities']['coverage_alignment_policy'], manifest['policies']['mate_policy'])
        validate_replay(pool, replay)
        pool.sampling_unit = 'read'
        pool.coverage_alignment_policy = manifest['capabilities']['coverage_alignment_policy']
    if (pool.genome_len != layout.length or len(pool) != manifest['units']
            or len(pool.starts) != manifest['blocks'] or pool.total_bases != manifest['aligned_bases']):
        raise ValueError('Manifest geometry mismatch')
    return pool, layout, manifest


def build_bundle(pool, layout, outdir, sample_id, policies, sources, *, workers=4,
                 chunk_units=65536, storage_format='packed', source_check=None):
    """Core entry point; synthetic fixtures can supply a pool without pysam."""
    if workers < 1 or chunk_units < 1 or not sample_id.strip():
        raise ValueError('Require sample ID, positive workers and chunk size')
    validate_layout(layout)
    if pool.genome_len != layout.length:
        raise ValueError('Pool/reference length mismatch')
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=False)
    status = outdir / 'status'
    status.mkdir()
    (status / 'RUNNING').touch()
    started = time.perf_counter()
    try:
        options = [('packed', level) for level in (1, 6, 9)] + [('npz', 6), ('npy', 0)]
        if storage_format != 'auto':
            if storage_format not in ('packed', 'npz', 'npy'):
                raise ValueError('Unknown format')
            options = [(storage_format, 9 if storage_format == 'packed' else 6)]
        candidates = []
        # Only one sample and one candidate decompression are live at a time.
        with tempfile.TemporaryDirectory(prefix='.candidates-', dir=outdir) as temporary:
            for codec, level in options:
                directory = Path(temporary) / f'{codec}-{level}'
                tick = time.perf_counter()
                storage = save_storage(pool, directory, codec, workers, chunk_units, level)
                write_seconds = time.perf_counter() - tick
                tick = time.perf_counter()
                loaded = retry_payload_io(lambda: read_storage(directory, workers, transient_access=True), directory)
                load_seconds = time.perf_counter() - tick
                replay = validate_replay(pool, loaded)
                del loaded
                files = inventory(directory, storage_payload_names(storage) | {'storage.json'})
                candidates.append(dict(codec=codec, level=level, files=files,
                    storage_bytes=sum(f['bytes'] for f in files.values()),
                    allocated_bytes=sum(f['allocated_bytes'] for f in files.values()),
                    write_seconds=write_seconds, load_seconds=load_seconds, replay=replay,
                    directory=directory))
            best = min(candidates, key=lambda row: row['storage_bytes'])
            for name in sorted(best['files']):
                path = best['directory'] / name
                retry_payload_io(lambda: path.rename(outdir / name), path)
        if source_check is not None:
            source_check()
        geometry = getattr(pool, 'alignment_geometry', None)
        if geometry is not None:
            shutil.copytree(geometry.directory, outdir/'alignment')
            extra = sorted('alignment/'+p.name for p in (outdir/'alignment').iterdir() if p.is_file())
            storage = json.loads((outdir/'storage.json').read_text())
            storage['alignment_geometry_files'] = extra
            write_json(outdir/'storage.json', storage)
            best['files'].update(inventory(outdir, extra))
        try:
            revision = subprocess.check_output(['git', '-C', str(Path(__file__).resolve().parents[1]),
                                                 'rev-parse', 'HEAD'], text=True).strip()
        except (OSError, subprocess.CalledProcessError):
            revision = None
        manifest = dict(schema=SCHEMA, sample_id=sample_id, sources=sources, policies=policies,
            layout=layout_metadata(layout),
            reference_m5_missing=[n for n, m in zip(layout.names, layout.md5) if not m],
            missing_header_contigs=list(layout.missing_header_contigs),
            coordinates='0-based half-open; ordered concatenated contigs; CIGAR M/= /X',
            units=len(pool), blocks=len(pool.starts), aligned_bases=pool.total_bases,
            full_mean_depth=pool.total_bases / layout.length, source_revision=revision,
            implementation_sha256=sha256(Path(__file__)),
            coverage_core_sha256=sha256(Path(__file__).with_name('coverage_core.py')),
            python=sys.version, numpy=np.__version__, zlib=zlib.ZLIB_RUNTIME_VERSION,
            omitted=['SEQ', 'QUAL', 'QNAME', 'per-read MAPQ/flags'],
            retain_original_alignment=True)
        manifest['capabilities'] = (dict(geometry.metadata, coverage_alignment_policy=getattr(pool, 'coverage_alignment_policy', 'primary_only'))
            if geometry is not None else dict(geometry_schema='legacy-block-pool',
                sampling_units=[policies.get('sampling_unit', 'read')], supplementary_available=False,
                mate_relationships_available=policies.get('sampling_unit') == 'fragment', coverage_alignment_policy='primary_only'))
        if geometry is not None:
            manifest['omitted'] = ['SEQ','QUAL','QNAME']
        write_json(outdir / 'manifest.json', manifest)
        peak = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        report = dict(schema=SCHEMA, sample_id=sample_id, requested_format=storage_format,
            selected_codec=best['codec'], selected_level=best['level'], workers=workers,
            chunk_units=chunk_units, pending_chunks_limit=2 * workers,
            extraction_parallelism='serial alignment stream; threaded compression/read/hash',
            selection_objective='smallest validated data plus storage.json among measured candidates',
            candidates=[{k: v for k, v in row.items() if k != 'directory'} for row in candidates],
            numeric_pool_bytes=pool.nbytes, build_seconds=time.perf_counter() - started,
            process_peak_rss_bytes=int(peak if sys.platform == 'darwin' else peak * 1024),
            rss_scope='whole build process high-water mark, includes input and validation pools',
            replay_pass=True)
        write_json(outdir / 'report.json', report)
        files = sorted(set(best['files']) | {'manifest.json', 'report.json'})
        digests = dict(bounded_map(lambda name: (name, payload_sha256(outdir, name)), files, workers))
        write_json(outdir / 'checksums.json', digests)
        (status / 'RUNNING').unlink()
        (status / 'exit_code.txt').write_text('0\n')
        (status / 'SUCCESS').touch()
        return report
    except BaseException:
        (status / 'RUNNING').unlink(missing_ok=True)
        (status / 'FAILED').touch()
        (status / 'exit_code.txt').write_text('1\n')
        raise


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest='command', required=True)
    build = sub.add_parser('build', help='Extract and save verified sample geometry (default packed/zlib9)')
    build.add_argument('--sample-id', required=True)
    build.add_argument('--alignment', type=Path, required=True)
    reference = build.add_mutually_exclusive_group(required=True)
    reference.add_argument('--reference', type=Path)
    reference.add_argument('--reference-index', type=Path)
    build.add_argument('--outdir', type=Path, required=True)
    build.add_argument('--format', choices=['auto', 'packed', 'npz', 'npy'], default='packed',
                       help='Default packed: fixed zlib9; auto explicitly benchmarks five candidates')
    build.add_argument('--chunk-units', type=int, default=65536)
    build.add_argument('--workers', type=int, default=4)
    build.add_argument('--mate-policy', choices=['both', 'read1', 'read2'], default='both')
    build.add_argument('--sampling-unit', choices=['read', 'fragment'], default='read')
    build.add_argument('--min-mapq', type=int, default=0)
    build.add_argument('--exclude-duplicates', action='store_true')
    build.add_argument('--exclude-qcfail', action='store_true')
    inspect = sub.add_parser('inspect', help='Verify a saved reference and summarize its geometry')
    inspect.add_argument('--bundle', type=Path, required=True)
    inspect.add_argument('--workers', type=int, default=4)
    args = parser.parse_args()
    if args.workers < 1:
        parser.error('workers must be positive')
    if args.command == 'inspect':
        pool, layout, manifest = load_reference_bundle(args.bundle, args.workers)
        print(json.dumps(dict(sample_id=manifest['sample_id'], units=len(pool),
                              blocks=len(pool.starts), genome_len=layout.length,
                              aligned_bases=pool.total_bases, checksum_validation='PASS')))
        return
    if args.outdir.exists():
        raise FileExistsError(args.outdir)
    if args.chunk_units < 1 or not args.sample_id.strip():
        parser.error('Require sample ID and positive chunk size')
    inputs = [args.alignment, args.reference or args.reference_index]
    sources = list(bounded_map(lambda p: dict(path=str(p.resolve()), sha256=sha256(p),
                                             bytes=p.stat().st_size), inputs, args.workers))
    policies = {name: getattr(args, name) for name in
                ('mate_policy', 'sampling_unit', 'min_mapq', 'exclude_duplicates', 'exclude_qcfail')}
    extraction_started = time.perf_counter()
    pool, layout = extract_pool(args.alignment, args.reference,
                                reference_index=args.reference_index, **policies)
    import pysam
    sources[0].update(pysam=pysam.__version__, reference_mode='fasta' if args.reference else 'fai_alignment_fields')
    extraction_seconds = time.perf_counter() - extraction_started
    sources[0]['extraction_seconds'] = extraction_seconds
    def recheck_sources():
        for source in sources:
            if sha256(source['path']) != source['sha256']:
                raise ValueError('Source changed during build')
    report = build_bundle(pool, layout, args.outdir, args.sample_id, policies, sources,
        workers=args.workers, chunk_units=args.chunk_units, storage_format=args.format,
        source_check=recheck_sources)
    files = inventory(args.outdir)
    # External run summary avoids recursively counting a size report inside itself.
    print(json.dumps(dict(sample_id=args.sample_id, selected_codec=report['selected_codec'],
        selected_level=report['selected_level'], extraction_seconds=extraction_seconds,
        bundle_logical_bytes=sum(f['bytes'] for f in files.values()),
        bundle_allocated_bytes=sum(f['allocated_bytes'] for f in files.values()),
        files=files, units=len(pool), blocks=len(pool.starts), genome_len=layout.length,
        numeric_pool_bytes=pool.nbytes), sort_keys=True, indent=2))


if __name__ == '__main__':
    main()
