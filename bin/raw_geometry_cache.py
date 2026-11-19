"""On-demand shared, atomic cache for Raw Query decoded CRAM geometry.

Only alignment blocks and Read/Fragment relationships are persisted; no base
sequence, FASTQ/FASTA or per-draw dense coverage vectors are stored.
"""
import fcntl
import hashlib
import json
import os
from pathlib import Path
import tempfile
import zipfile

import numpy as np

from bin.coverage_core import BlockPool, Reference, merge_references

CACHE_SCHEMA = "raw-query-geometry-cache-2"
EXTRA = ("source_unit_ids", "unit_read_counts", "unit_segment_counts",
         "unit_fragment_ids", "unit_mate_sides")


def geometry_key(path, reference, policies):
    """Invalidate on changed source identity, genome sequence or read filters."""
    source = Path(path).resolve(strict=True)
    st = source.stat()
    basis = dict(schema=CACHE_SCHEMA, source=str(source), size=st.st_size,
                 mtime_ns=st.st_mtime_ns, inode=st.st_ino,
                 genome=dict(names=reference.names, lengths=reference.lengths,
                             md5=reference.md5, sort_policy=reference.sort_policy),
                 policies=policies)
    return hashlib.sha256(json.dumps(basis, sort_keys=True).encode()).hexdigest()


def _read_cache(path, reference):
    with np.load(path, allow_pickle=False) as data:
        info = json.loads(str(data["metadata"].item()))
        if info["schema"] != CACHE_SCHEMA or info["genome_len"] != reference.length:
            raise ValueError("Stale geometry cache")
        pool = BlockPool(data["starts"].copy(), data["ends"].copy(),
                         data["indptr"].copy(), reference.length)
        for field in EXTRA:
            setattr(pool, field, data[field].copy())
        if any(getattr(pool, field).ndim != 1 or
               getattr(pool, field).dtype.kind not in 'iu' or
               len(getattr(pool, field)) != len(pool) for field in EXTRA):
            raise ValueError("Incomplete geometry cache")
        if (len(np.unique(pool.source_unit_ids)) != len(pool) or
                np.any(pool.source_unit_ids < 0) or np.any(pool.unit_fragment_ids < 0) or
                np.any(pool.unit_read_counts < 1) or
                np.any(pool.unit_segment_counts < pool.unit_read_counts) or
                not np.isin(pool.unit_mate_sides, [0, 1, 2]).all()):
            raise ValueError("Invalid cached unit metadata")
        if int(pool.total_bases) != info["aligned_bases"]:
            raise ValueError("Geometry cache mass mismatch")
        pool.sampling_unit = info['sampling_unit']
        if pool.sampling_unit not in ('read', 'fragment'):
            raise ValueError('Invalid cached sampling unit')
        m5_audit = info.get("m5_audit")
        if (not isinstance(m5_audit, list) or
                any(not isinstance(r, dict) for r in m5_audit) or
                (m5_audit and [r.get("contig_id") for r in m5_audit] != list(reference.names))):
            raise ValueError("Incomplete geometry cache M5 audit")
        pool.m5_audit = m5_audit
        layout = Reference(reference.names, reference.lengths, tuple(info["md5"]),
                           tuple(info["missing"]), reference.sort_policy)
        merge_references(reference, layout)
        return pool, layout, list(info["missing"])


def _write_cache(path, pool, layout, missing):
    info = dict(schema=CACHE_SCHEMA, genome_len=pool.genome_len,
                aligned_bases=int(pool.total_bases), md5=layout.md5,
                missing=missing, sampling_unit=pool.sampling_unit,
                m5_audit=getattr(pool, "m5_audit", []))
    fd, temporary = tempfile.mkstemp(prefix=".geometry-", suffix=".npz", dir=path.parent)
    try:
        with os.fdopen(fd, "wb") as handle:
            arrays = {field: np.asarray(getattr(pool, field)) for field in EXTRA}
            np.savez(handle, metadata=np.array(json.dumps(info, sort_keys=True)),
                     starts=pool.starts, ends=pool.ends, indptr=pool.indptr, **arrays)
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, path)
    finally:
        Path(temporary).unlink(missing_ok=True)


def cached_extract(path, reference, extractor, cache_dir=None, **policies):
    """Return (pool, layout, missing, cache_hit), with per-source file locking."""
    if not cache_dir:
        pool, layout, missing = extractor(path, reference, **policies)
        return pool, layout, missing, False
    directory = Path(cache_dir).expanduser().resolve()
    directory.mkdir(parents=True, exist_ok=True)
    key = geometry_key(path, reference, policies)
    archive = directory / (key + ".npz")
    with (directory / (key + ".lock")).open("a+b") as lock:
        fcntl.flock(lock.fileno(), fcntl.LOCK_EX)
        if archive.is_file():
            try:
                pool, layout, missing = _read_cache(archive, reference)
                return pool, layout, missing, True
            except (OSError, ValueError, KeyError, EOFError, TypeError, zipfile.BadZipFile):
                archive.unlink(missing_ok=True)
        pool, layout, missing = extractor(path, reference, **policies)
        if geometry_key(path, reference, policies) != key:
            raise ValueError('Input changed during geometry extraction')
        _write_cache(archive, pool, layout, missing)
        return pool, layout, missing, False
