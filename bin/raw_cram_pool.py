"""Single-pass primary read/pair extraction into an ephemeral BlockPool."""
from array import array
from pathlib import Path

import numpy as np
from bin.coverage_core import BlockPool, Reference, merge_references


def extract_once(path, reference, *, sampling_unit="read", min_mapq=0,
                 coverage_alignment_policy="primary_only", strict_reference=False):
    """Keep actual mate blocks together; never use the span between mates."""
    import pysam
    if (sampling_unit not in ("read", "fragment") or min_mapq < 0
            or coverage_alignment_policy not in ("primary_only", "primary_plus_supplementary")):
        raise ValueError("Invalid sampling policy")
    with Path(path).open("rb") as handle:
        cram = handle.read(4) == b"CRAM"
    options = dict(format_options=[b"required_fields=4159", b"decode_md=0"]) if cram else {}
    groups, next_read = {}, 0
    ids, starts, ends, counts, segment_reads, fragments = (array("q") for _ in range(6))
    split_reads = {}
    read_counts, offsets = [], reference.offsets
    roles_by_read_id = {}
    sizes = dict(zip(reference.names, reference.lengths))
    with pysam.AlignmentFile(str(path), "rc" if cram else "rb", **options) as bam:
        sq = bam.header.to_dict().get("SQ", [])
        names = [r["SN"] for r in sq]
        if len(set(names)) != len(names) or not set(names) <= set(sizes):
            raise ValueError("Alignment header has duplicate/extra contigs")
        for row in sq:
            if row["LN"] != sizes[row["SN"]]:
                raise ValueError("Alignment header length differs from full FAI")
        observed_md5 = {r["SN"]: r.get("M5", "") for r in sq}
        if strict_reference:
            if set(names) != set(reference.names):
                raise ValueError("Strict FASTA reference requires complete CRAM @SQ dictionary")
            for name, expected in zip(reference.names, reference.md5):
                actual = observed_md5.get(name, "").lower()
                if not expected or actual != expected.lower():
                    raise ValueError(f"FASTA/CRAM M5 mismatch or missing: {name}")
        observed = Reference(reference.names, reference.lengths,
                             tuple(observed_md5.get(n, "") for n in reference.names),
                             tuple(n for n in reference.names if n not in names), reference.sort_policy)
        merged = merge_references(reference, observed)
        for order, read in enumerate(bam.fetch(until_eof=True)):
            role = (1 if read.is_read1 and not read.is_read2 else
                    2 if read.is_read2 and not read.is_read1 else 0 if not read.is_paired else 3)
            key = (read.query_name, read.get_tag("RG") if read.has_tag("RG") else "",
                   order if not read.query_name else None)
            group = groups.get(key)
            if group is None:
                group = [len(groups), {}, set()]
                groups[key] = group
            if role not in group[1]:
                group[1][role] = next_read
                roles_by_read_id[next_read] = role
                next_read += 1
            if read.is_secondary and read.is_supplementary:
                raise ValueError("Alignment cannot be both secondary and supplementary")
            if (read.is_unmapped or read.is_secondary
                    or (read.is_supplementary and coverage_alignment_policy == "primary_only")
                    or read.mapping_quality < min_mapq):
                continue
            blocks = read.get_blocks()
            if not blocks:
                continue
            if role == 3 or (not read.is_supplementary and role in group[2]):
                raise ValueError("Ambiguous read identity or multiple primary alignments")
            if sampling_unit == "fragment" and not read.query_name:
                raise ValueError("Pair sampling requires query names")
            if not read.is_supplementary:
                group[2].add(role)
            rid = group[1][role]
            if coverage_alignment_policy == "primary_plus_supplementary":
                from bin.alignment_geometry import query_interval
                intervals, primary = split_reads.setdefault(rid, [[], 0])
                intervals.append(query_interval(read.cigartuples, read.is_reverse))
                split_reads[rid][1] += int(not read.is_supplementary)
            segment_reads.append(rid)
            fragments.append(group[0])
            ids.append(group[1][role] if sampling_unit == "read" else group[0])
            counts.append(len(blocks))
            for left, right in blocks:
                if not 0 <= left < right <= sizes[read.reference_name]:
                    raise ValueError("Alignment block outside full FAI")
                starts.append(offsets[read.reference_name] + left)
                ends.append(offsets[read.reference_name] + right)
    del groups
    for intervals, primary in split_reads.values():
        if primary != 1 or len({row[2] for row in intervals}) != 1:
            raise ValueError("Incomplete/ambiguous split read or inconsistent query length")
        previous = -1
        for left, right, _ in sorted(intervals):
            if left < 0 or right <= left or left < previous:
                raise ValueError("Overlapping or invalid split-read query intervals")
            previous = max(previous, right)
    del split_reads
    ids = np.asarray(ids, dtype=np.int64)
    counts = np.asarray(counts, dtype=np.int64)
    order = np.argsort(ids, kind="stable")
    pointers = np.r_[0, np.cumsum(counts)]
    ordered_counts = counts[order]
    block_order = (np.repeat(pointers[order] - (np.cumsum(ordered_counts) - ordered_counts),
                             ordered_counts) + np.arange(int(ordered_counts.sum())))
    unit_ids = ids[order]
    boundaries = np.r_[0, np.flatnonzero(np.diff(unit_ids)) + 1] if len(ids) else np.empty(0, dtype=int)
    ptr = np.r_[np.r_[0, np.cumsum(ordered_counts)][boundaries], len(block_order)]
    pool = BlockPool(np.asarray(starts, dtype=np.int64)[block_order],
                     np.asarray(ends, dtype=np.int64)[block_order], ptr, reference.length)
    pool.source_unit_ids = unit_ids[boundaries]
    pool.unit_read_counts = np.diff(np.r_[boundaries, len(ids)])
    pool.unit_segment_counts = pool.unit_read_counts.copy()
    _, first_units = np.unique(ids, return_index=True)
    pool.unit_fragment_ids = np.asarray(fragments, dtype=np.int64)[first_units]
    pool.unit_mate_sides = (np.asarray([roles_by_read_id[int(uid)] for uid in pool.source_unit_ids], dtype=np.int8)
                            if sampling_unit == "read" else np.zeros(len(pool), dtype=np.int8))
    if coverage_alignment_policy == "primary_plus_supplementary":
        _, first = np.unique(np.asarray(segment_reads, dtype=np.int64), return_index=True)
        _, pool.unit_read_counts = np.unique(ids[first], return_counts=True)
    pool.sampling_unit = sampling_unit
    return pool, merged, list(observed.missing_header_contigs)

