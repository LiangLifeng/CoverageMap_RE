#!/usr/bin/env python3
"""Create paired-fragment and independent-read downsampled BAMs for CoverageMap.

The script writes two query manifests that can be passed directly to
source/main.nf: manifest_fragment.tsv and manifest_read.tsv.

Sampling units and eligible alignments follow CoverageMap's geometry semantics
at commit f7e805b. CRAM is decoded in geometry-only mode (CIGAR/coordinates,
flags, MAPQ, query name and tags); selected records are written as BAM.
For CRAM input, sequence and quality fields are intentionally omitted from the
output BAM because CoverageMap's coverage extractor only uses aligned blocks.

Requires pysam and numpy. A matching FAI is sufficient; no FASTA sequence is
needed for this geometry-only benchmark.
"""
from __future__ import annotations

import argparse
import csv
import fcntl
import hashlib
import json
import math
import os
from pathlib import Path
import re
import sqlite3
import sys
import tempfile

import numpy as np

def locate_repo_root():
    candidates = []
    env_root = os.environ.get("COVERAGEMAP_ROOT")
    if env_root:
        candidates.append(Path(env_root).expanduser().resolve())
    script = Path(__file__).resolve()
    candidates.extend(parent / "source" for parent in script.parents)
    candidates.extend(script.parents)
    cwd = Path.cwd().resolve()
    candidates.append(cwd)
    candidates.extend(cwd.parents)
    for candidate in candidates:
        if (candidate / "bin" / "reference_index.py").is_file():
            return candidate
    raise SystemExit(
        "ERROR: CoverageMap repository not found. Run from the repository root/bin, "
        "or set COVERAGEMAP_ROOT=/path/to/CoverageMap."
    )


ROOT = locate_repo_root()
sys.path.insert(0, str(ROOT))

try:
    import pysam
except ImportError as exc:  # pragma: no cover - runtime dependency message
    raise SystemExit("ERROR: pysam is required. Install it in the CoverageMap environment.") from exc

from bin.reference_index import ReferenceIndex

REQUIRED = ("Sample", "GCFID", "Depth", "CramPath")
MANIFEST_COLUMNS = (
    "Sample", "GCFID", "Depth", "CramPath", "SamplingMode", "RequestedDepth",
    "OriginalDepth", "OriginalManifestDepth", "SelectedUnits", "TotalUnits",
    "SelectedReads", "SelectedAlignedBases", "OvershootBases", "Seed",
)
AUDIT_COLUMNS = (
    "query_id", "Sample", "GCFID", "mode", "input_cram", "output_bam",
    "reference_fai", "genome_length", "requested_depth", "original_manifest_depth",
    "original_aligned_bases", "original_depth", "total_units", "total_eligible_reads",
    "selected_units", "selected_reads", "selected_aligned_bases", "actual_depth",
    "overshoot_bases", "seed", "mate_policy", "min_mapq", "exclude_duplicates",
    "exclude_qcfail", "output_bytes", "paired_manifest_included", "status", "error",
)


def sha256_file(path: Path, chunk_size: int = 8 * 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while True:
            block = handle.read(chunk_size)
            if not block:
                break
            digest.update(block)
    return digest.hexdigest()


def parse_depth(value):
    try:
        depth = float(value.strip().removesuffix("x").removesuffix("×"))
    except ValueError as exc:
        raise argparse.ArgumentTypeError("depth must be numeric, optionally followed by x or ×") from exc
    if not math.isfinite(depth) or depth <= 0:
        raise argparse.ArgumentTypeError("depth must be finite and > 0")
    return depth


def write_json_atomic(path, value):
    with tempfile.NamedTemporaryFile(mode="w", dir=path.parent, delete=False) as handle:
        temporary = Path(handle.name)
        try:
            json.dump(value, handle, indent=2, ensure_ascii=False, allow_nan=False)
            handle.write("\n")
            handle.flush()
            os.fsync(handle.fileno())
        except BaseException:
            temporary.unlink(missing_ok=True)
            raise
    os.replace(temporary, path)


def run_context(rows, manifest_path, index_path, args, policies):
    """Bind checkpoints to the executable, parameters and current input identities."""
    inputs = []
    for row in rows:
        source = row["_cram"]
        stat = source.stat()
        fai = resolve_reference_fai(row, manifest_path, args).resolve()
        row["_reference_fai"] = fai
        inputs.append(dict(query_id=row["_query_id"], cram=str(source.resolve()),
                           bytes=stat.st_size, mtime_ns=stat.st_mtime_ns,
                           fai=str(fai), fai_sha256=sha256_file(fai)))
    return dict(schema=1, script_sha256=sha256_file(Path(__file__)),
                manifest=str(manifest_path), manifest_sha256=sha256_file(manifest_path),
                index=str(index_path), index_sha256=sha256_file(index_path),
                target_depth=args.target_depth, seed=args.seed, unit_tag=args.unit_tag,
                policies=policies, inputs=inputs)


def checkpoint_result(path, context, policies, threads):
    checkpoint = json.loads(path.read_text())
    if checkpoint["context"] != context:
        raise ValueError(f"Checkpoint inputs or parameters changed: {path}")
    result = checkpoint["result"]
    bam = Path(result["output_bam"])
    if sha256_file(bam) != checkpoint["bam_sha256"]:
        raise ValueError(f"Checkpoint BAM checksum mismatch: {bam}")
    reads, bases = count_output_bases(bam, policies, threads)
    if (reads, bases) != (result["selected_reads"], result["selected_aligned_bases"]):
        raise ValueError(f"Checkpoint BAM accounting mismatch: {bam}")
    return result


def safe_name(value: str) -> str:
    cleaned = re.sub(r"[^A-Za-z0-9_.-]+", "_", value).strip("._")
    return cleaned[:100] or "sample"


def read_manifest(path: Path):
    rows = []
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames or len(set(reader.fieldnames)) != len(reader.fieldnames):
            raise ValueError("Input manifest must have unique TSV headers")
        missing = set(REQUIRED) - set(reader.fieldnames)
        if missing:
            raise ValueError("Input manifest missing required columns: " + ", ".join(sorted(missing)))
        for number, row in enumerate(reader, 1):
            for key in REQUIRED:
                if row.get(key) is None or not row[key].strip():
                    raise ValueError(f"Input manifest row {number}: empty required field {key}")
            try:
                manifest_depth = float(row["Depth"])
            except ValueError as exc:
                raise ValueError(f"Input manifest row {number}: Depth is not numeric") from exc
            if not math.isfinite(manifest_depth) or manifest_depth < 0:
                raise ValueError(f"Input manifest row {number}: Depth must be finite and >= 0")
            cram = Path(row["CramPath"])
            if not cram.is_absolute():
                cram = (path.parent / cram).resolve()
            if not cram.is_file() or cram.stat().st_size == 0:
                raise ValueError(f"Input manifest row {number}: CRAM/BAM missing or empty: {cram}")
            row["_query_id"] = f"q{number:08d}"
            row["_manifest_depth"] = manifest_depth
            row["_cram"] = cram
            rows.append(row)
    if not rows:
        raise ValueError("Input manifest contains no query rows")
    return rows, reader.fieldnames


def load_index_context(index_path: Path):
    index = ReferenceIndex(index_path)
    try:
        if index.metadata.get("storage_mode") != "native":
            raise ValueError("The supplied SQLite index must be a native CoverageMap index")
        policies = dict(index.metadata["config"]["policies"])
        if policies.get("mate_policy") != "both":
            raise ValueError(
                "This paired-vs-independent benchmark requires the index mate_policy=both; "
                f"observed {policies.get('mate_policy')!r}"
            )
        for key, default in (("min_mapq", 0), ("exclude_duplicates", False), ("exclude_qcfail", False)):
            policies.setdefault(key, default)
        return policies
    finally:
        index.close()


def get_layout(index_path: Path, gcf: str):
    index = ReferenceIndex(index_path)
    try:
        layout = index.group(gcf)
        names = tuple(layout["names"])
        lengths = tuple(int(x) for x in layout["lengths"])
        if not names or len(names) != len(lengths):
            raise ValueError(f"Invalid indexed layout for {gcf}")
        return dict(zip(names, lengths)), sum(lengths)
    finally:
        index.close()


def resolve_reference_fai(row, manifest_path: Path, args) -> Path:
    value = args.reference_fai
    if not value:
        value = row.get(args.reference_fai_column, "")
    if not value:
        raise ValueError(
            f"{row['_query_id']}: provide --reference-fai or a "
            f"{args.reference_fai_column!r} column in the input manifest"
        )
    path = Path(value)
    if not path.is_absolute():
        path = (manifest_path.parent / path).resolve()
    if not path.is_file():
        raise FileNotFoundError(f"Reference FAI not found: {path}")
    return path


def read_fai(path: Path) -> dict[str, int]:
    result = {}
    with path.open(encoding="utf-8") as handle:
        for number, line in enumerate(handle, 1):
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 2 or not fields[0]:
                raise ValueError(f"Invalid FAI row {number}: {path}")
            try:
                length = int(fields[1])
            except ValueError as exc:
                raise ValueError(f"Invalid FAI length at row {number}: {path}") from exc
            if length <= 0 or fields[0] in result:
                raise ValueError(f"Invalid or duplicate FAI contig at row {number}: {path}")
            result[fields[0]] = length
    if not result:
        raise ValueError(f"Empty FAI: {path}")
    return result


def open_alignment(path: Path, threads: int):
    with path.open("rb") as handle:
        is_cram = handle.read(4) == b"CRAM"
    mode = "rc" if is_cram else "rb"
    options = {}
    if is_cram:
        # Same geometry-only decode policy used by CoverageMap.extract_pool.
        options["format_options"] = [b"required_fields=4159", b"decode_md=0"]
    return pysam.AlignmentFile(str(path), mode, threads=threads, **options)


def validate_header(alignment, layout_lengths, reference_fai: Path):
    fai_lengths = read_fai(reference_fai)
    if fai_lengths != layout_lengths:
        missing = sorted(set(layout_lengths) - set(fai_lengths))[:5]
        extra = sorted(set(fai_lengths) - set(layout_lengths))[:5]
        mismatched = [name for name in set(fai_lengths) & set(layout_lengths)
                      if fai_lengths[name] != layout_lengths[name]][:5]
        raise ValueError(
            "FAI does not match the indexed GCF layout "
            f"(missing={missing}, extra={extra}, length_mismatch={mismatched})"
        )
    sq = alignment.header.to_dict().get("SQ", [])
    names = [item.get("SN") for item in sq]
    if len(names) != len(set(names)):
        raise ValueError("Duplicate reference contig names in alignment header")
    if not set(names) <= set(layout_lengths):
        unknown = sorted(set(names) - set(layout_lengths))
        raise ValueError("Alignment header contains contigs absent from indexed GCF: " + ", ".join(unknown[:10]))
    for item in sq:
        name = item["SN"]
        length = int(item["LN"])
        if length != layout_lengths[name] or length != fai_lengths[name]:
            raise ValueError(
                f"Alignment/FAI/index contig length mismatch for {name}: "
                f"header={length}, FAI={fai_lengths[name]}, index={layout_lengths[name]}"
            )


def eligible_read(read, policies):
    if read.is_unmapped or read.is_secondary or read.is_supplementary:
        return False
    if read.mapping_quality < int(policies["min_mapq"]):
        return False
    if policies["exclude_duplicates"] and read.is_duplicate:
        return False
    if policies["exclude_qcfail"] and read.is_qcfail:
        return False
    # The paired-vs-independent experiment deliberately retains both mates.
    return True


def aligned_block_bases(read) -> int:
    blocks = read.get_blocks()
    return sum(int(end) - int(start) for start, end in blocks)


def build_tagged_pool(
    input_cram: Path,
    reference_fai: Path,
    layout_lengths: dict[str, int],
    policies: dict,
    mode: str,
    temp_bam: Path,
    unit_db_path: Path,
    unit_tag: str,
    threads: int,
):
    """Write eligible alignments to a temporary BAM, tagging each sampling unit."""
    if mode not in ("read", "fragment"):
        raise ValueError(f"Unknown sampling mode: {mode}")
    conn = sqlite3.connect(unit_db_path)
    conn.execute("PRAGMA journal_mode=OFF")
    conn.execute("PRAGMA synchronous=OFF")
    conn.execute("PRAGMA temp_store=MEMORY")
    conn.execute("PRAGMA cache_size=-100000")
    conn.execute("CREATE TABLE units(unit_id INTEGER PRIMARY KEY, aligned_bases INTEGER NOT NULL)")
    if mode == "fragment":
        conn.execute(
            "CREATE TABLE unit_keys(qname TEXT NOT NULL, rg TEXT NOT NULL, unit_id INTEGER NOT NULL, "
            "PRIMARY KEY(qname, rg)) WITHOUT ROWID"
        )

    next_unit = 0
    eligible_reads = 0
    total_bases = 0
    tag_collisions = 0
    try:
        with open_alignment(input_cram, threads) as source:
            validate_header(source, layout_lengths, reference_fai)
            header = source.header
            with pysam.AlignmentFile(str(temp_bam), "wb", header=header, threads=threads) as tagged:
                for read in source.fetch(until_eof=True):
                    if not eligible_read(read, policies):
                        continue
                    bases = aligned_block_bases(read)
                    if bases <= 0:
                        continue
                    if read.has_tag(unit_tag):
                        tag_collisions += 1
                        raise ValueError(
                            f"Input alignment already contains temporary tag {unit_tag}; "
                            "rerun with a different --unit-tag"
                        )

                    if mode == "read":
                        unit_id = next_unit
                        next_unit += 1
                        conn.execute("INSERT INTO units(unit_id, aligned_bases) VALUES (?, ?)", (unit_id, bases))
                    else:
                        qname = read.query_name
                        if not qname:
                            raise ValueError("Fragment sampling requires nonempty query names")
                        rg = read.get_tag("RG") if read.has_tag("RG") else ""
                        found = conn.execute(
                            "SELECT unit_id FROM unit_keys WHERE qname=? AND rg=?", (qname, rg)
                        ).fetchone()
                        if found is None:
                            unit_id = next_unit
                            next_unit += 1
                            conn.execute(
                                "INSERT INTO unit_keys(qname, rg, unit_id) VALUES (?, ?, ?)",
                                (qname, rg, unit_id),
                            )
                            conn.execute(
                                "INSERT INTO units(unit_id, aligned_bases) VALUES (?, ?)", (unit_id, bases)
                            )
                        else:
                            unit_id = int(found[0])
                            conn.execute(
                                "UPDATE units SET aligned_bases=aligned_bases+? WHERE unit_id=?",
                                (bases, unit_id),
                            )

                    if unit_id > 2_147_483_647:
                        raise OverflowError("Sampling-unit count exceeds BAM integer-tag capacity")
                    read.set_tag(unit_tag, int(unit_id), value_type="i")
                    tagged.write(read)
                    eligible_reads += 1
                    total_bases += bases
                    if eligible_reads % 100_000 == 0:
                        conn.commit()
        conn.commit()
        unit_count = int(conn.execute("SELECT COUNT(*) FROM units").fetchone()[0])
        if unit_count == 0 or total_bases == 0:
            raise ValueError("No eligible aligned reads found after applying index filters")
        unit_bases = np.fromiter(
            (int(row[0]) for row in conn.execute("SELECT aligned_bases FROM units ORDER BY unit_id")),
            dtype=np.int64,
            count=unit_count,
        )
        if len(unit_bases) != unit_count or int(unit_bases.sum(dtype=np.int64)) != total_bases:
            raise AssertionError("Sampling-unit aligned-base accounting mismatch")
        return conn, unit_bases, eligible_reads, total_bases, unit_count
    except BaseException:
        conn.close()
        raise


def select_units(unit_bases: np.ndarray, target_bases: float, seed: int):
    """Match CoverageMap's permutation-prefix whole-unit selection policy."""
    if not math.isfinite(target_bases) or target_bases <= 0:
        raise ValueError("Target aligned bases must be finite and positive")
    total_bases = int(unit_bases.sum(dtype=np.int64))
    required = math.ceil(target_bases)
    if required > total_bases:
        raise ValueError(
            f"Requested target requires {required:,} aligned bases, but only "
            f"{total_bases:,} are available"
        )
    rng = np.random.default_rng(seed)
    order = rng.permutation(len(unit_bases))
    cumulative = np.cumsum(unit_bases[order], dtype=np.int64)
    stop = int(np.searchsorted(cumulative, required, side="left")) + 1
    selected_mask = np.zeros(len(unit_bases), dtype=np.bool_)
    selected_mask[order[:stop]] = True
    selected_bases = int(cumulative[stop - 1])
    return selected_mask, stop, selected_bases


def count_output_bases(path: Path, policies: dict, threads: int):
    reads = 0
    bases = 0
    with open_alignment(path, threads) as alignment:
        for read in alignment.fetch(until_eof=True):
            if eligible_read(read, policies):
                value = aligned_block_bases(read)
                if value > 0:
                    reads += 1
                    bases += value
    return reads, bases


def write_sampled_bam(
    temp_bam: Path,
    output_bam: Path,
    selected_mask: np.ndarray,
    unit_tag: str,
    threads: int,
):
    selected_reads = 0
    with pysam.AlignmentFile(str(temp_bam), "rb", threads=threads) as source:
        with pysam.AlignmentFile(str(output_bam), "wb", header=source.header, threads=threads) as destination:
            for read in source.fetch(until_eof=True):
                unit_id = int(read.get_tag(unit_tag))
                if selected_mask[unit_id]:
                    read.set_tag(unit_tag, None)
                    destination.write(read)
                    selected_reads += 1
    if not output_bam.is_file() or output_bam.stat().st_size == 0:
        raise IOError(f"Failed to write nonempty output BAM: {output_bam}")
    return selected_reads


def sample_one_mode(
    row, mode, output_dir, layout_lengths, genome_length, policies, target_depth,
    seed, unit_tag, threads,
):
    query_id = row["_query_id"]
    sample = row["Sample"]
    gcf = row["GCFID"]
    input_cram = row["_cram"]
    reference_fai = row["_reference_fai"]
    mode_dir = output_dir / mode
    mode_dir.mkdir(parents=True, exist_ok=True)
    output_bam = mode_dir / f"{query_id}_{safe_name(sample)}_{safe_name(gcf)}_{mode}.bam"
    if output_bam.exists():
        raise FileExistsError(f"Refusing to overwrite existing output: {output_bam}")

    with tempfile.TemporaryDirectory(prefix=f"{query_id}_{mode}_", dir=output_dir / ".tmp") as tmp:
        tmp = Path(tmp)
        temp_bam = tmp / "eligible_tagged.bam"
        unit_db = tmp / "units.sqlite"
        conn, unit_bases, eligible_reads, original_bases, total_units = build_tagged_pool(
            input_cram, reference_fai, layout_lengths, policies, mode,
            temp_bam, unit_db, unit_tag, threads,
        )
        try:
            target_bases = target_depth * genome_length
            selected_mask, selected_units, selected_bases = select_units(unit_bases, target_bases, seed)
        finally:
            conn.close()

        pending_bam = tmp / "sampled.bam"
        selected_reads = write_sampled_bam(
            temp_bam, pending_bam, selected_mask, unit_tag, threads,
        )
        # Publish only fully closed, independently verified BAMs.
        verified_reads, verified_bases = count_output_bases(pending_bam, policies, threads)
        if verified_bases != selected_bases or verified_reads != selected_reads:
            raise AssertionError(
                f"Output verification mismatch: selected={selected_reads} reads/{selected_bases} bases, "
                f"observed={verified_reads} reads/{verified_bases} bases"
            )
        os.replace(pending_bam, output_bam)

    actual_depth = verified_bases / genome_length
    overshoot = verified_bases - math.ceil(target_depth * genome_length)
    return dict(
        query_id=query_id, Sample=sample, GCFID=gcf, mode=mode,
        input_cram=str(input_cram), output_bam=str(output_bam.resolve()),
        reference_fai=str(reference_fai), genome_length=genome_length,
        requested_depth=target_depth, original_manifest_depth=row["_manifest_depth"],
        original_aligned_bases=original_bases, original_depth=original_bases / genome_length,
        total_units=total_units, total_eligible_reads=eligible_reads,
        selected_units=selected_units, selected_reads=verified_reads,
        selected_aligned_bases=verified_bases, actual_depth=actual_depth,
        overshoot_bases=overshoot, seed=seed, mate_policy=policies["mate_policy"],
        min_mapq=policies["min_mapq"], exclude_duplicates=policies["exclude_duplicates"],
        exclude_qcfail=policies["exclude_qcfail"], output_bytes=output_bam.stat().st_size,
        status="COMPLETE", error="",
    )


def write_tsv(path: Path, columns, rows):
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, required=True,
                        help="Original TSV with Sample, GCFID, Depth, CramPath")
    parser.add_argument("--index", type=Path, required=True,
                        help="The same native SQLite reference index used by source/main.nf")
    parser.add_argument("--outdir", type=Path, required=True,
                        help="New output directory; it must be empty or not yet exist")
    parser.add_argument("--target-depth", type=parse_depth, required=True,
                        help="Target mean aligned depth; both modes are sampled to this target")
    parser.add_argument("--reference-fai", type=Path,
                        help="One matching FAI for all rows; otherwise use a per-row ReferenceFai column")
    parser.add_argument("--reference-fai-column", default="ReferenceFai",
                        help="Input manifest column containing the matching FAI path (default: ReferenceFai)")
    parser.add_argument("--seed", type=int, default=401)
    parser.add_argument("--threads", type=int, default=2)
    parser.add_argument("--resume", action="store_true",
                        help="Reuse verified mode checkpoints with identical inputs and parameters")
    parser.add_argument("--unit-tag", default="ZQ",
                        help="Temporary BAM tag used internally; choose an unused two-character SAM tag")
    args = parser.parse_args()

    if not math.isfinite(args.target_depth) or args.target_depth <= 0:
        parser.error("--target-depth must be finite and > 0")
    if args.seed < 0 or args.threads < 1:
        parser.error("--seed must be >= 0 and --threads must be >= 1")
    if not re.fullmatch(r"[A-Za-z][A-Za-z0-9]", args.unit_tag):
        parser.error("--unit-tag must be a valid two-character SAM tag, e.g. ZQ")
    manifest_path = args.manifest.resolve()
    index_path = args.index.resolve()
    if not manifest_path.is_file() or not index_path.is_file():
        parser.error("--manifest and --index must exist")
    policies = load_index_context(index_path)
    rows, input_columns = read_manifest(manifest_path)

    # Check global FAI accessibility before creating output or processing any query.
    context = run_context(rows, manifest_path, index_path, args, policies)
    context_id = hashlib.sha256(json.dumps(context, sort_keys=True).encode()).hexdigest()
    outdir = args.outdir.resolve()
    outdir.mkdir(parents=True, exist_ok=True)
    lock = (outdir / ".run.lock").open("a")
    try:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
    except BlockingIOError:
        parser.error(f"Another sampling process is using {outdir}")
    existing = [path for path in outdir.iterdir() if path.name != ".run.lock"]
    config_path = outdir / "run_context.json"
    if existing:
        if not args.resume:
            parser.error(f"--outdir must be empty; use --resume for verified checkpoints: {outdir}")
        if not config_path.is_file() or json.loads(config_path.read_text()) != context:
            parser.error("Resume requires the same recorded inputs, code and parameters")
    else:
        write_json_atomic(config_path, context)
    (outdir / ".tmp").mkdir(exist_ok=True)
    checkpoints = outdir / ".checkpoints"
    checkpoints.mkdir(exist_ok=True)

    manifest_rows = {"read": [], "fragment": []}
    audit_rows = []
    failures = []
    layout_cache = {}

    for row in rows:
        query_id = row["_query_id"]
        row_results = {}
        try:
            gcf = row["GCFID"]
            if gcf not in layout_cache:
                layout_cache[gcf] = get_layout(index_path, gcf)
            layout_lengths, genome_length = layout_cache[gcf]
            reference_fai = resolve_reference_fai(row, manifest_path, args)
            row["_reference_fai"] = reference_fai
            for mode in ("read", "fragment"):
                try:
                    checkpoint = checkpoints / f"{query_id}_{mode}.json"
                    if args.resume and checkpoint.is_file():
                        result = checkpoint_result(checkpoint, context_id, policies, args.threads)
                        print(f"[{query_id}] {mode}: resumed verified BAM", flush=True)
                    else:
                        result = sample_one_mode(
                            row, mode, outdir, layout_lengths, genome_length, policies,
                            args.target_depth, args.seed, args.unit_tag, args.threads,
                        )
                        write_json_atomic(checkpoint, dict(
                            context=context_id, result=result,
                            bam_sha256=sha256_file(Path(result["output_bam"])),
                        ))
                    result["paired_manifest_included"] = False
                    row_results[mode] = result
                    audit_rows.append(result)
                    print(
                        f"[{query_id}] {mode}: depth={result['actual_depth']:.6g}x "
                        f"units={result['selected_units']}/{result['total_units']} "
                        f"reads={result['selected_reads']} output={result['output_bam']}",
                        flush=True,
                    )
                except Exception as error:
                    failures.append(dict(query_id=query_id, Sample=row["Sample"], GCFID=gcf,
                                         mode=mode, error=f"{type(error).__name__}: {error}"))
                    audit_rows.append(dict(query_id=query_id, Sample=row["Sample"], GCFID=gcf,
                                           mode=mode, input_cram=str(row["_cram"]),
                                           output_bam="", reference_fai=str(reference_fai),
                                           genome_length=genome_length, requested_depth=args.target_depth,
                                           original_manifest_depth=row["_manifest_depth"],
                                           paired_manifest_included=False, status="FAILED",
                                           error=f"{type(error).__name__}: {error}"))
                    print(f"[{query_id}] {mode} FAILED: {type(error).__name__}: {error}", file=sys.stderr, flush=True)
            # Keep exactly the same query cohort in both manifests. If either mode
            # fails for a source sample, omit that sample from both benchmark runs.
            if set(row_results) == {"read", "fragment"}:
                for mode in ("read", "fragment"):
                    result = row_results[mode]
                    result["paired_manifest_included"] = True
                    manifest_rows[mode].append(dict(
                        Sample=row["Sample"], GCFID=gcf,
                        # Depth is the measured output depth, because CoverageMap uses
                        # manifest Depth to route the low/high branch.
                        Depth=f"{result['actual_depth']:.12g}",
                        CramPath=result["output_bam"], SamplingMode=mode,
                        RequestedDepth=f"{args.target_depth:.12g}",
                        OriginalDepth=f"{result['original_depth']:.12g}",
                        OriginalManifestDepth=row["_manifest_depth"],
                        SelectedUnits=result["selected_units"], TotalUnits=result["total_units"],
                        SelectedReads=result["selected_reads"],
                        SelectedAlignedBases=result["selected_aligned_bases"],
                        OvershootBases=result["overshoot_bases"], Seed=args.seed,
                    ))
            else:
                print(f"[{query_id}] excluded from BOTH manifests because one sampling mode failed", file=sys.stderr, flush=True)
        except Exception as error:
            failures.append(dict(query_id=query_id, Sample=row.get("Sample", ""),
                                 GCFID=row.get("GCFID", ""), mode="both",
                                 error=f"{type(error).__name__}: {error}"))
            print(f"[{query_id}] FAILED: {type(error).__name__}: {error}", file=sys.stderr, flush=True)

    for mode in ("read", "fragment"):
        write_tsv(outdir / f"manifest_{mode}.tsv", MANIFEST_COLUMNS, manifest_rows[mode])
    write_tsv(outdir / "sampling_audit.tsv", AUDIT_COLUMNS, audit_rows)
    write_tsv(outdir / "failures.tsv", ("query_id", "Sample", "GCFID", "mode", "error"), failures)
    summary = dict(
        input_manifest=str(manifest_path), reference_index=str(index_path),
        target_depth=args.target_depth, seed=args.seed, policies=policies,
        input_queries=len(rows), read_manifest_rows=len(manifest_rows["read"]),
        fragment_manifest_rows=len(manifest_rows["fragment"]), failures=len(failures),
        manifests={mode: str((outdir / f"manifest_{mode}.tsv").resolve()) for mode in ("read", "fragment")},
        audit=str((outdir / "sampling_audit.tsv").resolve()),
    )
    write_json_atomic(outdir / "summary.json", summary)
    print(json.dumps(summary, indent=2, ensure_ascii=False), flush=True)
    return 1 if failures else 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (ValueError, OSError, KeyError, sqlite3.Error) as error:
        print(f"ERROR: {type(error).__name__}: {error}", file=sys.stderr)
        raise SystemExit(2)
