"""Exact aligned-block sampling and bounded-memory window coverage."""
from __future__ import annotations

from array import array
from dataclasses import dataclass, field
import gzip
import hashlib
import math
import re
import time
from pathlib import Path
from typing import Optional

import numpy as np

CONTIG_SORT_POLICY = 'natural-ascii-digit-runs-v1'

@dataclass(frozen=True)
class Reference:
    names: tuple
    lengths: tuple
    md5: tuple
    missing_header_contigs: tuple = field(default=(), compare=False)
    sort_policy: str | None = None

    @property
    def length(self):
        return sum(self.lengths)

    @property
    def offsets(self):
        out = {}
        position = 0
        for name, length in zip(self.names, self.lengths):
            out[name] = position
            position += length
        return out


def contig_sort_key(name):
    def token(index, part):
        if index % 2 == 0:
            return 0, part
        digits = part.lstrip('0')
        return 1, len(digits), digits
    tokens = tuple(token(i, part) for i, part in enumerate(re.split(r'([0-9]+)', name)))
    return tokens, name


def canonical_reference(ref):
    order = sorted(range(len(ref.names)), key=lambda i: contig_sort_key(ref.names[i]))
    missing = set(ref.missing_header_contigs)
    return Reference(tuple(ref.names[i] for i in order), tuple(ref.lengths[i] for i in order),
                     tuple(ref.md5[i] for i in order),
                     tuple(ref.names[i] for i in order if ref.names[i] in missing),
                     CONTIG_SORT_POLICY)


def layout_metadata(ref):
    data = dict(names=list(ref.names), lengths=list(ref.lengths), md5=list(ref.md5))
    if ref.sort_policy is not None:
        data.update(sort_policy=ref.sort_policy, offsets=list(ref.offsets.values()))
    return data


def reference_from_layout(data, missing=()):
    names, lengths = tuple(data['names']), tuple(data['lengths'])
    md5 = tuple(data.get('md5', [''] * len(names)))
    if (not names or any(not isinstance(n, str) or not n or any(c in n for c in '\r\n\t\0') for n in names)
            or len(names) != len(set(names)) or len(names) != len(lengths) or len(names) != len(md5)
            or any(not isinstance(d, str) for d in md5)
            or any(type(n) is not int or n <= 0 for n in lengths)):
        raise ValueError('Invalid reference layout')
    ref = Reference(names, lengths, md5, tuple(missing), data.get('sort_policy'))
    if ref.sort_policy is not None:
        if (ref.sort_policy != CONTIG_SORT_POLICY or canonical_reference(ref) != ref
                or data.get('offsets') != list(ref.offsets.values())):
            raise ValueError('Invalid natural contig order or offsets')
    elif 'offsets' in data:
        raise ValueError('Contig offsets require an explicit sort policy')
    return ref


def merge_references(left, right):
    """Validate common full geometry; retain known M5 without inventing unknowns."""
    if (left.names != right.names or left.lengths != right.lengths
            or left.sort_policy != right.sort_policy):
        raise ValueError('Reference/test layouts differ (names/order/lengths)')
    merged = []
    for name, first, second in zip(left.names, left.md5, right.md5):
        if first and second and first.lower() != second.lower():
            raise ValueError(f'Reference/test layouts differ (M5): {name}')
        merged.append((first or second).lower())
    return Reference(left.names, left.lengths, tuple(merged), sort_policy=left.sort_policy)

def read_reference(path):
    path = Path(path)
    if not path.is_file():
        raise FileNotFoundError(f"Matching reference FASTA missing: {path}")
    names, lengths, checksums = [], [], []
    name, length, digest = None, 0, None
    opener = gzip.open if path.suffix == ".gz" else open
    def finish():
        if name is not None:
            if length == 0:
                raise ValueError(f"Empty reference contig: {name}")
            names.append(name)
            lengths.append(length)
            checksums.append(digest.hexdigest())
    with opener(path, "rt") as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                finish()
                fields = line[1:].split()
                if not fields:
                    raise ValueError("Empty FASTA header")
                name, length, digest = fields[0], 0, hashlib.md5()
            else:
                if name is None:
                    raise ValueError("FASTA sequence before header")
                sequence = "".join(line.split()).upper().encode("ascii")
                length += len(sequence)
                digest.update(sequence)
    finish()
    if not names or len(names) != len(set(names)):
        raise ValueError("Reference must contain unique, nonempty contigs")
    return Reference(tuple(names), tuple(lengths), tuple(checksums))

@dataclass
class BlockPool:
    starts: np.ndarray
    ends: np.ndarray
    indptr: np.ndarray
    genome_len: int

    def __post_init__(self):
        for name in ("starts", "ends", "indptr"):
            value = np.asarray(getattr(self, name))
            if value.ndim != 1 or value.dtype.kind not in "iu":
                raise ValueError(f"{name} must be a 1-D integer array")
            setattr(self, name, value)
        if self.genome_len <= 0 or len(self.starts) != len(self.ends):
            raise ValueError("Invalid pool length")
        if not len(self.indptr) or self.indptr[0] != 0 or self.indptr[-1] != len(self.starts):
            raise ValueError("Invalid unit pointers")
        if np.any(self.indptr[1:] <= self.indptr[:-1]):
            raise ValueError("Each sampling unit must contain at least one aligned block")
        if np.any(self.starts < 0) or np.any(self.ends <= self.starts) or np.any(self.ends > self.genome_len):
            raise ValueError("Invalid aligned block coordinates")
        lengths = self.ends.astype(np.int64) - self.starts.astype(np.int64)
        self.unit_bases = (np.add.reduceat(lengths, self.indptr[:-1].astype(np.intp))
                           if len(self) else np.empty(0, dtype=np.int64))
        self.total_bases = int(self.unit_bases.sum(dtype=np.int64))
        if self.total_bases < 0:
            raise OverflowError("Aligned-base total exceeds int64")
        for value in (self.starts, self.ends, self.indptr, self.unit_bases):
            value.setflags(write=False)

    def __len__(self):
        return len(self.indptr) - 1

    @property
    def nbytes(self):
        geometry = getattr(self,'alignment_geometry',None)
        return self.starts.nbytes + self.ends.nbytes + self.indptr.nbytes + self.unit_bases.nbytes + (geometry.nbytes if geometry is not None else 0)

    @classmethod
    def from_intervals(cls, intervals, genome_len):
        values = np.asarray(intervals)
        if values.size == 0:
            values = np.empty((0, 2), dtype=np.int64)
        if values.ndim != 2 or values.shape[1] != 2 or values.dtype.kind not in "iu":
            raise ValueError("intervals must have integer shape (N,2)")
        return cls(values[:, 0].copy(), values[:, 1].copy(), np.arange(len(values) + 1, dtype=np.uint64), genome_len)

    def selected_blocks(self, ids):
        ids = np.asarray(ids)
        if ids.size == 0:
            return self.starts[:0], self.ends[:0]
        if ids.ndim != 1 or ids.dtype.kind not in "iu" or np.any(ids < 0) or np.any(ids >= len(self)):
            raise ValueError("Invalid sampling unit IDs")
        if len(np.unique(ids)) != len(ids):
            raise ValueError("Sampling unit IDs must be unique")
        ids = ids.astype(np.intp)
        counts = (self.indptr[ids + 1] - self.indptr[ids]).astype(np.intp)
        offsets = np.cumsum(counts) - counts
        indices = np.repeat(self.indptr[ids].astype(np.intp) - offsets, counts) + np.arange(int(counts.sum()))
        return self.starts[indices], self.ends[indices]

    def save(self, path):
        np.savez_compressed(path, starts=self.starts, ends=self.ends, indptr=self.indptr,
                            genome_len=np.array(self.genome_len, dtype=np.int64))

    @classmethod
    def load(cls, path):
        with np.load(path, allow_pickle=False) as archive:
            return cls(archive["starts"], archive["ends"], archive["indptr"], int(archive["genome_len"]))


def reproject_pool(pool, source, target):
    """Move each intact aligned block by name, preserving sampling-unit IDs."""
    merge_references(canonical_reference(source), target)
    if pool.genome_len != source.length or source.length > np.iinfo(np.int64).max:
        raise ValueError('Saved pool length differs from reference or exceeds int64')
    offsets = np.array(list(source.offsets.values()), dtype=np.int64)
    starts, ends = pool.starts.astype(np.int64), pool.ends.astype(np.int64)
    contigs = np.searchsorted(offsets, starts, side='right')-1
    if np.any(ends > (offsets + source.lengths)[contigs]):
        raise ValueError('Saved aligned block crosses a source contig boundary')
    target_offsets = target.offsets
    delta = (np.array([target_offsets[name] for name in source.names], dtype=np.int64)-offsets)[contigs]
    dtype = np.uint32 if pool.genome_len <= np.iinfo(np.uint32).max else np.uint64
    moved = BlockPool((starts+delta).astype(dtype), (ends+delta).astype(dtype),
                      pool.indptr, pool.genome_len)
    if moved.total_bases != pool.total_bases:
        raise AssertionError('Contig reprojection changed aligned-base mass')
    return moved

def read_reference_index(path):
    names, lengths = [], []
    with open(path) as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 2:
                raise ValueError("FAI rows require contig name and length")
            length = int(fields[1])
            if not fields[0] or length <= 0:
                raise ValueError("Invalid FAI contig")
            names.append(fields[0])
            lengths.append(length)
    if not names or len(names) != len(set(names)):
        raise ValueError("FAI must contain unique nonempty contigs")
    return Reference(tuple(names), tuple(lengths), tuple("" for _ in names))

def extract_pool(alignment, reference_path=None, *, reference_index=None, sort_contigs=False, mate_policy="both", sampling_unit="read",
                 min_mapq=0, exclude_duplicates=False, exclude_qcfail=False,
                 capture_geometry=None, coverage_alignment_policy='primary_only', geometry_namespace=''):
    """Stream primary alignments; multiple CIGAR blocks remain one sampling unit."""
    if mate_policy not in ("both", "read1", "read2") or sampling_unit not in ("read", "fragment"):
        raise ValueError("Invalid mate policy or sampling unit")
    if sampling_unit == "fragment" and mate_policy != "both":
        raise ValueError("Fragment sampling requires mate-policy both")
    if min_mapq < 0:
        raise ValueError("min_mapq must be nonnegative")
    if (reference_path is None) == (reference_index is None):
        raise ValueError("Specify matching FASTA or explicit FAI layout mode")
    ref = read_reference(reference_path) if reference_path is not None else read_reference_index(reference_index)
    if sort_contigs:
        ref = canonical_reference(ref)
    if capture_geometry is not None:
        from bin.alignment_geometry import capture
        geometry, ref = capture(alignment, ref, capture_geometry,
            namespace=geometry_namespace,
            min_mapq=min_mapq, exclude_duplicates=exclude_duplicates,
            exclude_qcfail=exclude_qcfail, reference_path=reference_path)
        return geometry.project(sampling_unit, coverage_alignment_policy, mate_policy), ref
    if coverage_alignment_policy != 'primary_only':
        raise ValueError('Supplementary coverage requires complete alignment capture')
    if not Path(alignment).is_file():
        raise FileNotFoundError(alignment)
    import pysam
    offsets, sizes = ref.offsets, dict(zip(ref.names, ref.lengths))
    digests = dict(zip(ref.names, ref.md5))
    starts, ends, pointers = array("Q"), array("Q"), array("Q", [0])
    unit_numbers, groups = array("Q"), {}
    with open(alignment, "rb") as handle:
        is_cram = handle.read(4) == b"CRAM"
    options = {"format_options": [b"required_fields=4159", b"decode_md=0"]} if is_cram else {}
    if reference_path is not None:
        options["reference_filename"] = str(reference_path)
    # Decode only the columns needed for aligned-block geometry, never SEQ/QUAL.
    with pysam.AlignmentFile(str(alignment), "rc" if is_cram else "rb", **options) as bam:
        header_refs = bam.header.to_dict().get("SQ", [])
        header_names = [row['SN'] for row in header_refs]
        if len(header_names) != len(set(header_names)):
            raise ValueError('Duplicate alignment header contig')
        if (not set(header_names) <= set(ref.names)
                or reference_path is not None and set(header_names) != set(ref.names)):
            raise ValueError("Reference contig set does not match alignment header")
        for record in bam.header.to_dict().get("SQ", []):
            name = record["SN"]
            if name not in sizes or record["LN"] != sizes[name]:
                raise ValueError(f"Reference header mismatch: {name}")
            if reference_path is not None and record.get("M5", "").lower() not in ("", digests[name]):
                raise ValueError(f"Reference sequence digest mismatch: {name}")
        if reference_index is not None:
            header_md5 = {row["SN"]: row.get("M5", "") for row in header_refs}
            ref = Reference(ref.names, ref.lengths, tuple(header_md5.get(name, '') for name in ref.names),
                            tuple(name for name in ref.names if name not in header_md5), ref.sort_policy)
        for read in bam.fetch(until_eof=True):
            if read.is_unmapped or read.is_secondary or read.is_supplementary:
                continue
            if read.mapping_quality < min_mapq or (exclude_duplicates and read.is_duplicate) or (exclude_qcfail and read.is_qcfail):
                continue
            if read.is_paired and ((mate_policy == "read1" and not read.is_read1) or (mate_policy == "read2" and not read.is_read2)):
                continue
            name = read.reference_name
            if name not in offsets:
                raise ValueError(f"Unknown aligned reference: {name}")
            blocks = read.get_blocks()
            if not blocks:
                continue
            if sampling_unit == "fragment":
                if not read.query_name:
                    raise ValueError("Fragment sampling requires query names")
                key = (read.query_name, read.get_tag("RG") if read.has_tag("RG") else "")
                if key not in groups:
                    groups[key] = len(groups)
                unit_id = groups[key]
            for start, end in blocks:
                if start < 0 or end <= start or end > sizes[name]:
                    raise ValueError("Read block out of reference bounds")
                starts.append(offsets[name] + start)
                ends.append(offsets[name] + end)
                if sampling_unit == "fragment":
                    unit_numbers.append(unit_id)
            if sampling_unit == "read":
                pointers.append(len(starts))
    dtype = np.uint32 if ref.length <= np.iinfo(np.uint32).max else np.uint64
    start_values = np.asarray(starts, dtype=dtype)
    end_values = np.asarray(ends, dtype=dtype)
    if sampling_unit == "fragment" and groups:
        units = np.asarray(unit_numbers, dtype=np.intp)
        order = np.argsort(units, kind="stable")
        start_values, end_values = start_values[order], end_values[order]
        counts = np.bincount(units, minlength=len(groups))
        pointer_values = np.concatenate(([0], np.cumsum(counts, dtype=np.uint64))).astype(np.uint64)
    else:
        pointer_values = np.asarray(pointers, dtype=np.uint64)
    return BlockPool(start_values, end_values, pointer_values, ref.length), ref

def runtime_pool(pool, saved_policies, requested_policies):
    """Select a supported unit definition; legacy relationship data is never inferred."""
    for key in ('mate_policy', 'min_mapq', 'exclude_duplicates', 'exclude_qcfail'):
        if saved_policies.get(key) != requested_policies.get(key):
            raise ValueError('Reference/query alignment filter mismatch: '+key)
    unit = requested_policies.get('sampling_unit', 'read')
    policy = requested_policies.get('coverage_alignment_policy', 'primary_only')
    geometry = getattr(pool, 'alignment_geometry', None)
    if geometry is not None:
        return geometry.project(unit, policy, requested_policies.get('mate_policy', 'both'))
    if unit != saved_policies.get('sampling_unit', 'read') or policy != 'primary_only':
        raise ValueError('Legacy geometry lacks requested mate/supplementary capability')
    pool.sampling_unit = unit
    pool.coverage_alignment_policy = policy
    return pool


def compatible_policies(saved, requested):
    keys = ('mate_policy', 'min_mapq', 'exclude_duplicates', 'exclude_qcfail')
    return all(saved.get(k) == requested.get(k) for k in keys)


def window_edges(genome_len, bin_count):
    if not isinstance(genome_len, (int, np.integer)) or not isinstance(bin_count, (int, np.integer)):
        raise ValueError("Genome length and bin count must be integers")
    if genome_len <= 0 or not 1 <= bin_count <= genome_len:
        raise ValueError("Require 1 <= bin_count <= genome length")
    if genome_len > np.iinfo(np.int64).max:
        raise ValueError("Genome length exceeds signed 64-bit coordinates")
    return np.array([(i * int(genome_len)) // int(bin_count) for i in range(bin_count + 1)], dtype=np.int64)

def lazy_depth_units(pool, rng, target_bases):
    """Uniform permutation prefix via sparse Fisher-Yates, stopping on whole units.

    Draw ranks in batches, then apply the same swap/delete operation as a full
    Fisher-Yates shuffle. Only accepted units and touched slots are retained.
    Unused ranks after stopping are discarded; this is a versioned RNG policy.
    """
    target = math.ceil(target_bases)
    remaining = len(pool)
    slots = {}
    selected = []
    mass = 0
    while remaining:
        batch = min(128, remaining)
        ranks = rng.integers(0, np.arange(remaining, remaining - batch, -1), dtype=np.int64)
        for rank in ranks:
            rank = int(rank)
            unit = slots.get(rank, rank)
            last = remaining - 1
            if rank != last:
                slots[rank] = slots.get(last, last)
            slots.pop(last, None)
            remaining -= 1
            selected.append(unit)
            mass += int(pool.unit_bases[unit])
            if mass >= target:
                return np.asarray(selected, dtype=np.intp)
    raise AssertionError("Available mass was exhausted before target")


def select_units(pool, rng, *, count=None, fraction=None, target_bases=None,
                 depth_sampler="permutation"):
    if depth_sampler not in ("permutation", "lazy"):
        raise ValueError("Unknown depth sampler")
    if sum(x is not None for x in (count, fraction, target_bases)) != 1:
        raise ValueError("Specify exactly one of count, fraction, target_bases")
    if fraction is not None:
        if not math.isfinite(fraction) or not 0 <= fraction <= 1:
            raise ValueError("fraction must be in [0,1]")
        count = math.floor(fraction * len(pool))
    if count is not None:
        if not isinstance(count, (int, np.integer)) or not 0 <= count <= len(pool):
            raise ValueError("count must be an integer in [0,N]")
        if count == 0:
            return np.empty(0, dtype=np.intp)
        if count == len(pool):
            return np.arange(len(pool), dtype=np.intp)
        return rng.choice(len(pool), size=count, replace=False, shuffle=False)
    if not math.isfinite(target_bases) or target_bases < 0 or target_bases > pool.total_bases:
        raise ValueError("Target exceeds available aligned bases or is invalid")
    if target_bases == 0:
        return np.empty(0, dtype=np.intp)
    if depth_sampler == "lazy":
        return lazy_depth_units(pool, rng, target_bases)
    # Preserve the historical RNG policy for exact seed/ID reproduction.
    order = rng.permutation(len(pool))
    mass = np.cumsum(pool.unit_bases[order], dtype=np.int64)
    stop = int(np.searchsorted(mass, math.ceil(target_bases), side="left")) + 1
    return order[:stop]

def choose_depth_sampler(pool, target_bases):
    """Calibrate once using independent probe RNGs; never consume simulation RNG."""
    if not math.isfinite(target_bases) or target_bases < 0 or target_bases > pool.total_bases:
        raise ValueError("Target exceeds available aligned bases or is invalid")
    info = {"policy": "paired-median3-margin10-v1", "probe_seeds": [91031, 91032, 91033],
            "trials_seconds": {"lazy": [], "permutation": []}}
    if target_bases == 0:
        info.update(selected="permutation", reason="zero target; calibration skipped")
        return "permutation", info
    for trial, seed in enumerate(info["probe_seeds"]):
        order = ("lazy", "permutation") if trial % 2 == 0 else ("permutation", "lazy")
        for sampler in order:
            rng = np.random.default_rng(seed)
            start = time.perf_counter()
            ids = select_units(pool, rng, target_bases=target_bases, depth_sampler=sampler)
            elapsed = time.perf_counter() - start
            info["trials_seconds"][sampler].append(elapsed)
            del ids
    medians = {key: float(np.median(values)) for key, values in info["trials_seconds"].items()}
    # Prefer historical permutation near a tie; timing is not a reproducibility promise.
    selected = "lazy" if medians["lazy"] < 0.9 * medians["permutation"] else "permutation"
    info.update(median_seconds=medians, selected=selected,
                reason="lazy requires a measured sampling-time advantage exceeding 10%")
    return selected, info


def dense_coverage(pool, ids):
    """Independent validation baseline; intentionally allocates G bases."""
    starts, ends = pool.selected_blocks(ids)
    depth = np.zeros(pool.genome_len, dtype=np.int64)
    for start, end in zip(starts, ends):
        depth[int(start):int(end)] += 1
    return depth

def window_sums(pool, ids, edges):
    """Exact O(selected blocks + bins), without a per-base coverage array."""
    edges = np.asarray(edges)
    if edges.ndim != 1 or edges.dtype.kind not in "iu" or len(edges) < 2 or edges[0] != 0 or edges[-1] != pool.genome_len or np.any(edges[1:] <= edges[:-1]):
        raise ValueError("Window edges must partition the full reference")
    starts, ends = pool.selected_blocks(ids)
    starts, ends = starts.astype(np.int64), ends.astype(np.int64)
    width = np.diff(edges).astype(np.int64)
    result = np.zeros(len(width), dtype=np.int64)
    if not len(starts):
        return result
    first = np.searchsorted(edges, starts, side="right") - 1
    last = np.searchsorted(edges, ends - 1, side="right") - 1
    same = first == last
    np.add.at(result, first[same], ends[same] - starts[same])
    different = ~same
    np.add.at(result, first[different], edges[first[different] + 1] - starts[different])
    np.add.at(result, last[different], ends[different] - edges[last[different]])
    diff = np.zeros(len(width) + 1, dtype=np.int64)
    np.add.at(diff, first[different] + 1, 1)
    np.add.at(diff, last[different], -1)
    np.cumsum(diff[:-1], out=diff[:-1])
    diff[:-1] *= width
    result += diff[:-1]
    return result


def covered_sums(pool, edges, ids=None):
    """Per-window union length of selected complete units, O(blocks log blocks).

    No per-base array and no inferred span between CIGAR blocks or mates.
    Edges partition the full virtual reference, including short final windows.
    """
    if ids is None:
        starts, ends = pool.starts, pool.ends
    else:
        starts, ends = pool.selected_blocks(ids)
    # Validate the partition even for empty selections.
    if not len(starts):
        return window_sums(pool, np.empty(0, dtype=np.int64), edges)
    order = np.argsort(starts, kind='stable')
    starts, ends = starts[order], ends[order]
    maximum = np.maximum.accumulate(ends)
    positions = np.flatnonzero(np.r_[True, starts[1:] > maximum[:-1]])
    union = BlockPool.from_intervals(
        np.column_stack((starts[positions], np.maximum.reduceat(ends, positions))), pool.genome_len)
    result = window_sums(union, np.arange(len(union)), edges)
    if np.any(result < 0) or np.any(result > np.diff(edges)) or int(result.sum()) != union.total_bases:
        raise AssertionError('Per-window union coverage is inconsistent')
    return result
