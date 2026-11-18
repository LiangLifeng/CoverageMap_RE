#!/usr/bin/env python3
"""Simulate lower-depth window coverage from complete aligned sampling units."""
from __future__ import annotations
import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import sys
import time

import numpy as np
from coverage_core import extract_pool, select_units, choose_depth_sampler, window_edges, window_sums

def sha256(path):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()

def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--alignment", type=Path, required=True)
    reference = parser.add_mutually_exclusive_group(required=True)
    reference.add_argument("--reference", type=Path, help="Matching FASTA sequence")
    reference.add_argument("--reference-index", type=Path, help="Explicit FAI layout + alignment-only CRAM decoding (no SEQ/QUAL)")
    parser.add_argument("--outdir", type=Path, required=True)
    target = parser.add_mutually_exclusive_group(required=True)
    target.add_argument("--reads", type=int, help="Number of sampling units retained; fragment mode counts fragments")
    target.add_argument("--fraction", type=float, help="Fixed floor(fraction*N) units, not Bernoulli thinning")
    target.add_argument("--depth", type=float, help="Target mean aligned-base depth; last unit retained whole")
    parser.add_argument("--bins", type=int, default=1000)
    parser.add_argument("--replicates", type=int, default=100)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--mate-policy", choices=["both","read1","read2"], default="both")
    parser.add_argument("--sampling-unit", choices=["read","fragment"], default="read")
    parser.add_argument("--min-mapq", type=int, default=0)
    parser.add_argument("--exclude-duplicates", action="store_true")
    parser.add_argument("--exclude-qcfail", action="store_true")
    parser.add_argument("--save-pool", action="store_true")
    parser.add_argument("--depth-sampler", choices=["permutation", "lazy", "auto"], default="permutation",
                        help="Depth RNG policy; auto calibrates once, lazy generates only a prefix")
    parser.add_argument("--output-format", choices=["tsv", "npy"], default="tsv",
                        help="npy stores int64 masses and shared edges; mean depth = mass/window width")
    return parser.parse_args()

class WindowOutput:
    """Stream either legacy text rows or disk-backed fixed-shape int64 arrays."""
    def __init__(self, outdir, edges, replicates, output_format):
        self.outdir = outdir
        self.edges = edges
        self.widths = np.diff(edges)
        self.format = output_format
        self.handle = None
        self.matrix = None
        if output_format == "tsv":
            self.handle = (outdir / "simulations.tsv").open("w", newline="")
            self.writer = csv.writer(self.handle, delimiter="\t")
            self.writer.writerow(["replicate", "window", "start", "end", "aligned_bases", "mean_depth"])
        else:
            np.save(outdir / "edges.npy", np.asarray(edges, dtype="<i8"), allow_pickle=False)
            self.matrix = np.lib.format.open_memmap(outdir / "simulations.npy", mode="w+",
                dtype="<i8", shape=(replicates, len(edges) - 1))

    def rows(self, label, values):
        return ([label, i, int(self.edges[i]), int(self.edges[i+1]), int(mass),
                 format(float(mass) / self.widths[i], ".12g")] for i, mass in enumerate(values))

    def full(self, values):
        if self.format == "tsv":
            with (self.outdir / "full.tsv").open("w", newline="") as handle:
                writer = csv.writer(handle, delimiter="\t")
                writer.writerow(["replicate", "window", "start", "end", "aligned_bases", "mean_depth"])
                writer.writerows(self.rows("full", values))
        else:
            np.save(self.outdir / "full.npy", np.asarray(values, dtype="<i8"), allow_pickle=False)

    def replicate(self, index, values):
        if self.format == "tsv":
            self.writer.writerows(self.rows(index, values))
        else:
            self.matrix[index] = values

    def close(self):
        if self.handle is not None:
            self.handle.close()
        if self.matrix is not None:
            self.matrix.flush()
            self.matrix._mmap.close()
            self.matrix = None


def main():
    args = parse_args()
    if args.outdir.exists():
        raise FileExistsError(f"Output exists: {args.outdir}")
    if args.replicates < 1 or args.seed < 0:
        raise ValueError("replicates must be positive and seed nonnegative")
    started = time.perf_counter()
    stages = {key: 0.0 for key in ("extract_pool", "select_units", "sampler_calibration", "full_coverage",
              "simulation_coverage", "hash_inputs", "write_results", "save_pool")}
    t = time.perf_counter()
    pool, ref = extract_pool(args.alignment, args.reference, reference_index=args.reference_index,
        mate_policy=args.mate_policy, sampling_unit=args.sampling_unit, min_mapq=args.min_mapq,
        exclude_duplicates=args.exclude_duplicates, exclude_qcfail=args.exclude_qcfail)
    stages["extract_pool"] += time.perf_counter() - t
    edges = window_edges(ref.length, args.bins)
    target = {"count": args.reads} if args.reads is not None else (
        {"fraction": args.fraction} if args.fraction is not None else {"target_bases": args.depth * ref.length})
    sampler = args.depth_sampler
    calibration = None
    if sampler == "auto":
        t = time.perf_counter()
        if "target_bases" in target:
            sampler, calibration = choose_depth_sampler(pool, target["target_bases"])
        else:
            sampler = "permutation"
            calibration = {"selected": sampler, "reason": "count/fraction mode; depth policy unused"}
        stages["sampler_calibration"] = time.perf_counter() - t
    rng = np.random.default_rng(args.seed)
    t = time.perf_counter()
    first_ids = select_units(pool, rng, **target, depth_sampler=sampler)
    stages["select_units"] += time.perf_counter() - t
    t = time.perf_counter()
    identity = {"alignment_sha256": sha256(args.alignment),
        "reference_input_sha256": sha256(args.reference or args.reference_index),
        "reference_mode": "fasta" if args.reference is not None else "fai_alignment_fields",
        "reference_header_md5": list(ref.md5), "reference_names": list(ref.names),
        "missing_header_contigs": list(ref.missing_header_contigs),
        "reference_lengths": list(ref.lengths)}
    stages["hash_inputs"] += time.perf_counter() - t
    args.outdir.mkdir(parents=True)
    t = time.perf_counter()
    full = window_sums(pool, np.arange(len(pool)), edges)
    stages["full_coverage"] += time.perf_counter() - t
    t = time.perf_counter()
    output = WindowOutput(args.outdir, edges, args.replicates, args.output_format)
    summaries = []
    try:
        output.full(full)
        stages["write_results"] += time.perf_counter() - t
        for replicate in range(args.replicates):
            if replicate == 0:
                ids = first_ids
            else:
                t = time.perf_counter()
                ids = select_units(pool, rng, **target, depth_sampler=sampler)
                stages["select_units"] += time.perf_counter() - t
            t = time.perf_counter()
            sums = window_sums(pool, ids, edges)
            stages["simulation_coverage"] += time.perf_counter() - t
            mass = int(pool.unit_bases[ids].sum(dtype=np.int64))
            if int(sums.sum(dtype=np.int64)) != mass:
                raise AssertionError("Aligned-base mass was not conserved")
            summaries.append({"replicate": replicate, "units": len(ids), "aligned_bases": mass,
                "mean_depth": mass/ref.length,
                "selected_ids_sha256": hashlib.sha256(np.asarray(ids, dtype="<i8").tobytes()).hexdigest()})
            t = time.perf_counter()
            output.replicate(replicate, sums)
            stages["write_results"] += time.perf_counter() - t
    finally:
        t = time.perf_counter()
        output.close()
        stages["write_results"] += time.perf_counter() - t
    if args.save_pool:
        t = time.perf_counter()
        pool.save(args.outdir / "pool.npz")
        stages["save_pool"] += time.perf_counter() - t
    elapsed = time.perf_counter() - started
    stages["other_inside_main"] = elapsed - sum(stages.values())
    provenance = {"schema": "coverage-map-2", "inputs": identity, "sampling_unit": args.sampling_unit,
        "mate_policy": args.mate_policy,
        "coverage": "M/= /X aligned read bases; pair overlaps counted for each read",
        "target": target, "seed": args.seed, "depth_sampler": sampler, "depth_sampler_requested": args.depth_sampler,
        "sampler_calibration": calibration,
        "depth_sampler_version": "numpy-permutation-v1" if sampler == "permutation" else "lazy-fisher-yates-batch128-v1",
        "output_format": args.output_format, "units_total": len(pool), "blocks_total": len(pool.starts),
        "aligned_bases_total": pool.total_bases, "mean_depth_full": pool.total_bases/ref.length,
        "bins": args.bins, "pool_numeric_bytes": pool.nbytes, "replicates": summaries,
        "simulation_matrix_policy": "disk-backed int64 NPY" if args.output_format == "npy" else "streamed TSV",
        "stages_seconds": stages, "elapsed_seconds": elapsed}
    # Provenance serialization happens after the reported stage clock; process wall includes it.
    (args.outdir / "provenance.json").write_text(json.dumps(provenance, indent=2, allow_nan=False) + "\n")
    print(json.dumps({key: provenance[key] for key in
          ("units_total", "blocks_total", "mean_depth_full", "pool_numeric_bytes", "elapsed_seconds")}), file=sys.stderr)

if __name__ == "__main__":
    try:
        main()
    except (ValueError, FileNotFoundError, FileExistsError) as error:
        print(f"ERROR: {error}",file=sys.stderr)
        sys.exit(2)
