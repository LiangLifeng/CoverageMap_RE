#!/usr/bin/env python3
"""Plan, score and collect one-query raw-CRAM tasks from a metadata SQLite."""
import argparse
import base64
from collections import defaultdict
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import sqlite3
import shutil
import sys

for name in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[name] = "1"
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import numpy as np
from bin.coverage_core import canonical_reference, read_reference_index, read_reference, window_sums, covered_sums, select_units
from bin.raw_query_sampling import query_arms, compatible_reference_ids, select_nearest_units
from bin.raw_geometry_cache import cached_extract
from bin.raw_reference_replacement import resolve_reference_choices, InsufficientReferenceCapacity
from bin.gcf_windows import GCF_CONCAT, edges_for, coordinate_rows, COORDINATE_COLUMNS
from bin.query_manifest_stats import similarity, aggregate, reference_row_means
from bin.raw_cram_pool import extract_once
from bin.jaccard_breadth import thresholds
from bin.query_decision_methods import parse_methods, validate_options, COLUMNS as METHOD_COLUMNS
from bin.query_manifest import DISTRIBUTION_COLUMNS
from bin.raw_query_outputs import legacy_outputs, unscored_outputs, read_rows

VERSION = "1.3.0"
SUMMARY = ("query_id", "Sample", "GCFID", "status", "error", "arm", "manifest_depth",
           "actual_query_depth", "available_count", "eligible_count", "requested_draws",
           "valid_draws", "unique_references", "duplicate_draws", "decoded_reference_files",
           "alpha", "jaccard", "jaccard_p", "jaccard_status", "jaccard_decision",
           "pearson", "pearson_p", "pearson_status", "pearson_decision", "detection_status",
           "sampling_tolerance_bp", "tolerated_undershoot_draws", "tolerated_overshoot_draws",
           "out_of_tolerance_draws", "geometry_cache_hits", "query_read1_bases", "query_read2_bases",
           "insufficient_reference_policy", "replaced_draws", "replacement_attempts")
SUMMARY = tuple(dict.fromkeys((*SUMMARY, *DISTRIBUTION_COLUMNS, "query_exit_code", "result_source")))
DRAWS = ("draw_id", "reference_id", "reference_sample", "reference_path", "status", "error",
         "rng_sha256", "selected_ids_sha256", "selected_units", "selected_reads",
         "actual_bases", "target_bases", "actual_depth", "jaccard", "pearson",
         "intersection", "union", "jaccard_status", "pearson_status",
         "mate_side", "sampling_status", "sampling_delta_bases", "available_bases",
         "tolerance_bp", "tolerance_accepted", "geometry_cache_hit",
         "original_reference_id", "replacement_attempts")


class UnscoredQuery(ValueError):
    def __init__(self, status, reason):
        super().__init__(reason)
        self.status = status


def rng(seed, *parts):
    token = hashlib.sha256(json.dumps([seed, *parts], ensure_ascii=False, separators=(",", ":")).encode()).digest()
    return np.random.default_rng(np.frombuffer(token, dtype="<u4")), token.hex()


def dump(path, value):
    Path(path).write_text(json.dumps(value, ensure_ascii=False, indent=2, allow_nan=False) + "\n")


def table(path, columns, rows):
    with Path(path).open("w", newline="") as handle:
        writer = csv.DictWriter(handle, columns, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def stamp(path):
    path = Path(path)
    if not path.is_absolute():
        raise ValueError(f"Require absolute input path: {path}")
    stat = path.stat()
    if not path.is_file() or stat.st_size == 0:
        raise ValueError(f"Require nonempty input file: {path}")
    return dict(size=stat.st_size, mtime_ns=stat.st_mtime_ns)


def check_stamps(spec):
    for path, expected in spec["stamps"].items():
        if stamp(path) != expected:
            raise ValueError(f"Input changed after planning: {path}")


def fai_mapping(path):
    if path is None:
        return {}
    with Path(path).open(encoding="utf-8-sig", newline="") as handle:
        first = handle.readline()
        handle.seek(0)
        rows = csv.DictReader(handle, delimiter="\t" if "\t" in first else ",")
        headers = rows.fieldnames or []
        fields = ("GCFID", "FaiPath") if {"GCFID", "FaiPath"} <= set(headers) else ("GCF_ID", "GCF_path")
        if not set(fields) <= set(headers):
            raise ValueError("FAI map requires GCFID/FaiPath or GCF_ID/GCF_path")
        result = {}
        for row in rows:
            gcf, fai = row[fields[0]], row[fields[1]]
            if not gcf or not fai or gcf in result:
                raise ValueError("Empty or duplicate GCF in FAI map")
            result[gcf] = fai
        return result


def candidates(connection, query, low):
    # Source CSV depth is a planning value; actual insufficient mass is never upsampled.
    rows = connection.execute(
        'SELECT id,Sample,GCFID,CramFilePath,Depth FROM reference_records '
        'WHERE GCFID=? ORDER BY Sample,id', (query["GCFID"],)).fetchall()
    selected, seen = [], set()
    available = 0
    for row in rows:
        row = dict(row)
        depth = float(row["Depth"])
        if not math.isfinite(depth) or depth < 0:
            raise ValueError("Invalid reference CSV Depth")
        path = str(Path(row["CramFilePath"]).resolve())
        if depth <= 0 or row["Sample"] == query["Sample"] or path == str(Path(query["CramPath"]).resolve()):
            continue
        # A CRAM is a source, not a repeated CSV row. Conflicting ownership is rejected.
        if path in seen:
            continue
        seen.add(path)
        available += 1
        if low and depth < float(query["Depth"]):
            continue
        row["CramFilePath"] = path
        selected.append(row)
    return available, selected


def plan(args):
    parse_methods(getattr(args, "jaccard_methods", "normal"))
    validate_options(args.alpha, getattr(args, "bootstrap_replicates", 1000),
                     getattr(args, "bootstrap_confidence", .95))
    thresholds(args.jaccard_coverage_threshold, getattr(args, "jaccard_coverage_thresholds", None))
    mapping = fai_mapping(args.fai_map)
    config = dict(reference_draws=args.reference_draws, alpha=args.alpha, seed=args.seed,
                  sampling_tolerance_bp=getattr(args, "sampling_tolerance_bp", 150),
                  insufficient_reference_policy=getattr(args, "insufficient_reference_policy", "unscored"),
                  geometry_cache_dir=str(Path(args.geometry_cache_dir).expanduser().resolve()) if getattr(args, "geometry_cache_dir", None) else None,
                  depth_threshold=args.depth_threshold, bin_count=args.bin_count,
                  window_size=args.window_size, sampling_unit=args.sampling_unit,
                  min_mapq=args.min_mapq, jaccard_coverage_threshold=args.jaccard_coverage_threshold,
                  threads=args.threads, min_samples=args.min_samples,
                  jaccard_methods=getattr(args, "jaccard_methods", "normal"),
                  bootstrap_replicates=getattr(args, "bootstrap_replicates", 1000),
                  bootstrap_confidence=getattr(args, "bootstrap_confidence", .95),
                  jaccard_coverage_thresholds=getattr(args, "jaccard_coverage_thresholds", None),
                  coverage_unit=getattr(args, "coverage_unit", None),
                  coverage_alignment_policy=getattr(args, "coverage_alignment_policy", "primary_only"))
    with Path(args.manifest).open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if (not reader.fieldnames or len(set(reader.fieldnames)) != len(reader.fieldnames)
                or not {"Sample", "GCFID", "Depth", "CramPath"} <= set(reader.fieldnames)):
            raise ValueError("Manifest requires unique Sample,GCFID,Depth,CramPath columns")
        queries = list(reader)
    if not queries:
        raise ValueError("Manifest is empty")
    output = Path(args.output)
    output.mkdir()
    records, roster = [], []
    with sqlite3.connect(Path(args.index).resolve().as_uri() + "?mode=ro", uri=True) as connection:
        connection.row_factory = sqlite3.Row
        schema = connection.execute("SELECT value FROM metadata WHERE key='schema'").fetchone()
        if schema is None or json.loads(schema[0]) != "coverage-raw-cram-index-1":
            raise ValueError("Require coverage-raw-cram-index-1 SQLite")
        has_catalog = connection.execute(
            "SELECT 1 FROM sqlite_master WHERE type='table' AND name='reference_genomes'").fetchone() is not None
        for number, query in enumerate(queries, 1):
            qid = f"q{number:08d}"
            spec = dict(schema="raw-cram-query-1", query_id=qid, query=query, config=config,
                        references=[], choices=[], stamps={}, reference_errors={}, error="")
            dependencies = []
            try:
                if None in query or any(v is None for v in query.values()):
                    raise ValueError("Malformed manifest record")
                if any(not query[k].strip() for k in ("Sample", "GCFID", "Depth", "CramPath")):
                    raise ValueError("Empty required query field")
                depth = float(query["Depth"])
                if not math.isfinite(depth) or depth < 0:
                    raise ValueError("Query Depth must be finite and nonnegative")
                if query.get("Coverage"):
                    unit = config["coverage_unit"]
                    if unit not in ("fraction", "percent"):
                        raise ValueError("Nonempty Coverage requires --coverage-unit fraction|percent")
                    value = float(query["Coverage"])
                    if not math.isfinite(value) or not 0 <= value <= (1 if unit == "fraction" else 100):
                        raise ValueError("Query Coverage outside the declared unit range")
                low = depth <= args.depth_threshold
                spec["arm"] = "low" if low else "high"
                optional = [query[key] for key in ("CraiPath", "BedGraphPath") if query.get(key)]
                for path in dict.fromkeys([query["CramPath"], *optional]):
                    spec["stamps"][path] = stamp(path)
                    dependencies.append(path)
                available, refs = candidates(connection, query, low)
                spec.update(available_count=available, eligible_count=len(refs))
                if len({r["Sample"] for r in refs}) <= args.min_samples:
                    status = ("INSUFFICIENT_REFERENCE_DEPTH" if low and available > args.min_samples
                              else "INSUFFICIENT_REFERENCE_SAMPLES")
                    raise UnscoredQuery(status, "Requires eligible distinct reference samples > min-samples")
                fna = None
                if has_catalog:
                    genome = connection.execute(
                        "SELECT FnaFilePath,FaiFilePath FROM reference_genomes WHERE GCFID=?",
                        (query["GCFID"],)).fetchone()
                    if genome is None or not genome["FnaFilePath"] or not genome["FaiFilePath"]:
                        raise ValueError("GCF lacks a complete SQL FNA/FAI catalog entry")
                    fna, fai = genome["FnaFilePath"], genome["FaiFilePath"]
                else:
                    fai = (mapping.get(query["GCFID"]) if args.fai_map else
                           str(Path(args.fai_dir) / (query["GCFID"] + ".fai")) if args.fai_dir else None)
                    if not fai:
                        raise ValueError("Legacy CRAM-only SQL requires --fai-map or --fai-dir")
                spec["fna"] = str(Path(fna).absolute()) if fna else None
                spec["fai"] = str(Path(fai).absolute())
                generator, token = rng(args.seed, query["GCFID"], query["Sample"], qid, "reference-samples")
                choices = generator.integers(0, len(refs), size=args.reference_draws)
                spec["references"] = (refs if low and config["insufficient_reference_policy"] == "replace"
                                      else [refs[i] for i in sorted(set(int(i) for i in choices))])
                spec["choices"] = [refs[int(i)]["id"] for i in choices]
                spec["source_rng_sha256"] = token
                for genome_file in (spec["fna"], spec["fai"]):
                    if genome_file:
                        spec["stamps"][genome_file] = stamp(genome_file)
                        dependencies.append(genome_file)
                for ref in spec["references"]:
                    path = ref["CramFilePath"]
                    try:
                        spec["stamps"][path] = stamp(path)
                        dependencies.append(path)
                    except (OSError, ValueError) as error:
                        spec["reference_errors"][str(ref["id"])] = str(error)
            except UnscoredQuery as error:
                spec.update(unscored_status=error.status, unscored_reason=str(error))
            except (OSError, ValueError, TypeError) as error:
                spec["error"] = str(error)
            if any(any(c in p for c in "|\n\r\t") for p in dependencies):
                spec["error"] = "Input paths cannot contain pipe, newline or tab"
                dependencies = []
            payload = base64.b64encode(json.dumps(spec, sort_keys=True, ensure_ascii=False).encode()).decode()
            records.append(dict(query_id=qid, payload=payload, inputs="|".join(dict.fromkeys(dependencies))))
            roster.append(dict(query_id=qid, Sample=query.get("Sample", ""), GCFID=query.get("GCFID", ""),
                               manifest_depth=query.get("Depth", "")))
        table(output / "tasks.tsv", ("query_id", "payload", "inputs"), records)
        dump(output / "roster.json", roster)
    print(f"Planned {len(records)} query processes", flush=True)


def calculate(spec, output, *, extractor=extract_once):
    query, config, qid = spec["query"], spec["config"], spec["query_id"]
    count = config["reference_draws"]
    summary = dict(query_id=qid, Sample=query.get("Sample", ""), GCFID=query.get("GCFID", ""),
                   status="ERROR", error=spec["error"], arm=spec.get("arm", ""),
                   manifest_depth=query.get("Depth", ""), requested_draws=count, valid_draws=0,
                   manifest_coverage=query.get("Coverage", ""),
                   decoded_reference_files=0, alpha=config["alpha"], detection_status="NOT_EVALUATED",
                   available_count=spec.get("available_count", 0), eligible_count=spec.get("eligible_count", 0))
    if spec["error"]:
        return summary
    check_stamps(spec)
    if spec.get("unscored_status"):
        summary.update(status=spec["unscored_status"], error=spec["unscored_reason"])
        return summary
    strict = bool(spec.get("fna"))
    layout = canonical_reference(read_reference(spec["fna"]) if strict else read_reference_index(spec["fai"]))
    if strict:
        indexed = canonical_reference(read_reference_index(spec["fai"]))
        if (layout.names, layout.lengths) != (indexed.names, indexed.lengths):
            raise ValueError("SQL Genome FASTA/FAI contig dictionary mismatch")
    policies = dict(sampling_unit=config["sampling_unit"], min_mapq=config["min_mapq"],
                    coverage_alignment_policy=config.get("coverage_alignment_policy", "primary_only"),
                    strict_reference=strict)
    cache_dir = config.get("geometry_cache_dir") if strict else None
    query_pool, layout, missing, query_hit = cached_extract(
        query["CramPath"], layout, extractor, cache_dir=cache_dir, **policies)
    summary["geometry_cache_hits"] = int(query_hit)
    edges = edges_for(layout, GCF_CONCAT, count=config["bin_count"], size=config["window_size"])
    lengths = np.diff(edges)
    target = int(query_pool.total_bases)
    vector = window_sums(query_pool, np.arange(len(query_pool)), edges)
    breadth = covered_sums(query_pool, edges)
    summary.update(actual_query_depth=target / layout.length, error="",
                   actual_query_coverage=int(breadth.sum()) / layout.length,
                   gcf_reference_length=layout.length, windows=len(lengths),
                   requested_bin_count=config["bin_count"], nominal_window_size=config["window_size"],
                   window_policy=GCF_CONCAT, jaccard_coverage_threshold=config["jaccard_coverage_threshold"])
    np.save(output / "query_windows.npy", vector, allow_pickle=False)
    np.save(output / "query_covered_bases.npy", breadth, allow_pickle=False)
    # Mate-stratified sampling is activated by the new SQL FNA/FAI catalog.
    # Old CRAM-only catalogs retain their 202e289 sampling/statistics contract.
    sampling_arms = (query_arms(query_pool, config["sampling_unit"], count)
                     if strict and spec["arm"] == "low" else [])
    side_vectors = {a.name: window_sums(query_pool, a.query_ids, edges) for a in sampling_arms}
    side_breadths = {a.name: covered_sums(query_pool, edges, a.query_ids) for a in sampling_arms}
    if len(sampling_arms) == 2:
        summary.update(query_read1_bases=sampling_arms[0].target_bases,
                       query_read2_bases=sampling_arms[1].target_bases)
        np.save(output / "query_mate_covered_bases.npy",
                np.stack([side_breadths[a.name] for a in sampling_arms]), allow_pickle=False)
    np.save(output / "edges.npy", edges, allow_pickle=False)
    np.save(output / "window_lengths.npy", lengths, allow_pickle=False)
    table(output / "windows.tsv", COORDINATE_COLUMNS,
          coordinate_rows(layout, edges, config["window_size"]))
    dump(output / "layout.json", dict(names=layout.names, lengths=layout.lengths, md5=layout.md5,
                                      missing_query_header_contigs=missing))
    del query_pool
    if target == 0:
        summary.update(status="ZERO_QUERY_DEPTH", error="No eligible aligned query bases")
        return summary
    matrix = np.zeros((count, len(lengths)), dtype=np.int64)
    covered = np.zeros_like(matrix)
    references = {ref["id"]: ref for ref in spec["references"]}
    draw_arms = ([arm for arm in sampling_arms for _ in range(arm.slots)]
                 if sampling_arms else [None] * count)
    if len(draw_arms) != count:
        raise AssertionError("Mate arm draw partition differs from reference_draws")
    choices = list(spec["choices"])
    replacements = []
    preflight = {}
    source_policy = config.get('insufficient_reference_policy', 'unscored')
    summary.update(insufficient_reference_policy=source_policy, replaced_draws=0, replacement_attempts=0)
    if source_policy == 'replace' and spec['arm'] == 'low':
        if not strict or not cache_dir:
            raise ValueError('Replacement requires SQL FASTA catalog and shared geometry cache')
        def capacity(ref, arm):
            refid = ref['id']
            if refid not in preflight:
                if str(refid) in spec['reference_errors']:
                    raise ValueError(spec['reference_errors'][str(refid)])
                if os.path.samefile(query['CramPath'], ref['CramFilePath']):
                    raise ValueError('Reference aliases query source')
                candidate, _, _, hit = cached_extract(
                    ref['CramFilePath'], layout, extractor, cache_dir=cache_dir, **policies)
                summary['decoded_reference_files'] += int(not hit)
                summary['geometry_cache_hits'] += int(hit)
                info = dict(reference_id=refid, cache_hit=hit, capacities={})
                for candidate_arm in sampling_arms:
                    eligible = compatible_reference_ids(candidate, candidate_arm.mate_side)
                    info['capacities'][candidate_arm.name] = (
                        int(candidate.unit_bases[eligible].sum()), len(eligible))
                preflight[refid] = info
                del candidate
            return preflight[refid]['capacities'][arm.name]
        try:
            choices, replacements = resolve_reference_choices(spec, sampling_arms, capacity, rng)
        except InsufficientReferenceCapacity as error:
            replacements = error.events
            dump(output/'replacement_audit.json', replacements)
            table(output/'replacement_audit.tsv', tuple(replacements[0]) if replacements else ('draw_id','reason'), replacements)
            dump(output/'replacement_capacity_audit.json', list(preflight.values()))
            summary.update(status='INSUFFICIENT_REFERENCE_CAPACITY', error=str(error))
            return summary
        summary.update(replaced_draws=sum(a != b for a,b in zip(spec['choices'], choices)),
                       replacement_attempts=len(replacements))
    dump(output/'replacement_audit.json', replacements)
    table(output/'replacement_audit.tsv', tuple(replacements[0]) if replacements else ('draw_id','reason'), replacements)
    dump(output/'replacement_capacity_audit.json', list(preflight.values()))
    replacement_counts = defaultdict(int)
    for event in replacements:
        replacement_counts[event['draw_id']] += 1
    groups, draws = defaultdict(list), []
    for number, refid in enumerate(choices):
        ref = references[refid]
        _, token = rng(config["seed"], qid, query["Sample"], query["GCFID"], "read-units", number)
        draws.append(dict(draw_id=number, reference_id=refid, reference_sample=ref["Sample"],
                          query_id=qid, Sample=query["Sample"], GCFID=query["GCFID"],
                          sampling_unit=config["sampling_unit"],
                          coverage_alignment_policy=config.get("coverage_alignment_policy", "primary_only"),
                          reference_path=ref["CramFilePath"], rng_sha256=token,
                          original_reference_id=spec["choices"][number],
                          replacement_attempts=replacement_counts[number],
                          mate_side=draw_arms[number].mate_side if draw_arms[number] else "both",
                          status="ERROR", error="",
                          target_bases=draw_arms[number].target_bases if draw_arms[number] else target,
                          tolerance_bp=config.get("sampling_tolerance_bp", 150)))
        groups[refid].append(number)
    decode_records = []
    for refid, slots in groups.items():
        ref, pool = references[refid], None
        audit = dict(reference_id=refid, path=ref["CramFilePath"], draw_slots=slots,
                     decode_attempts=0, full_coverage_calculations=0, query_similarity_calculations=0)
        try:
            if str(refid) in spec["reference_errors"]:
                raise ValueError(spec["reference_errors"][str(refid)])
            if os.path.samefile(query["CramPath"], ref["CramFilePath"]):
                raise ValueError("Reference aliases query source")
            pool, layout, omitted, cache_hit = cached_extract(
                ref["CramFilePath"], layout, extractor, cache_dir=cache_dir, **policies)
            audit["decode_attempts"] = int(not cache_hit)
            summary["decoded_reference_files"] += int(not cache_hit)
            summary["geometry_cache_hits"] += int(cache_hit)
            audit.update(units=len(pool), actual_bases=int(pool.total_bases),
                         geometry_cache_hit=cache_hit, missing_header_contigs=omitted)
            full = None
            if spec["arm"] == "high":
                full_ids = np.arange(len(pool))
                full = (window_sums(pool, full_ids, edges), covered_sums(pool, edges))
                audit["full_coverage_calculations"] = 1
                full_score = similarity(vector, full[0], covered_x=breadth, covered_y=full[1],
                                        lengths=lengths, threshold=config["jaccard_coverage_threshold"])
                audit["query_similarity_calculations"] = 1
            for number in slots:
                draw = draws[number]
                draw['geometry_cache_hit'] = cache_hit
                if full is None:
                    generator, _ = rng(config["seed"], qid, query["Sample"], query["GCFID"], "read-units", number)
                    arm = draw_arms[number]
                    if arm is not None:
                        eligible = compatible_reference_ids(pool, arm.mate_side)
                        selected = select_nearest_units(
                            pool, eligible, generator, arm.target_bases,
                            config.get("sampling_tolerance_bp", 150))
                        ids = selected.indices
                        draw.update(sampling_status=selected.sampling_status,
                                    sampling_delta_bases=selected.delta_bases,
                                    available_bases=selected.available_bases,
                                    tolerance_accepted=selected.within_tolerance)
                        if not len(ids):
                            draw.update(status="ERROR", error="No eligible reference units for requested mate-side")
                            continue
                        comparison_vector = side_vectors[arm.name]
                        comparison_breadth = side_breadths[arm.name]
                        draw_status = "READY" if selected.within_tolerance else "OUT_OF_TOLERANCE"
                    else:
                        ids = select_units(pool, generator, target_bases=target, depth_sampler="permutation")
                        comparison_vector, comparison_breadth = vector, breadth
                        draw_status = "READY"
                        draw.update(sampling_status="LEGACY", available_bases=int(pool.total_bases),
                                    tolerance_accepted=True)
                    data, cov = window_sums(pool, ids, edges), covered_sums(pool, edges, ids)
                    scores = similarity(comparison_vector, data, covered_x=comparison_breadth,
                                        covered_y=cov, lengths=lengths,
                                        threshold=config["jaccard_coverage_threshold"])
                    audit["query_similarity_calculations"] += 1
                    identity = pool.source_unit_ids[ids]
                    draw["selected_ids_sha256"] = hashlib.sha256(identity.astype("<i8").tobytes()).hexdigest()
                else:
                    ids = full_ids
                    data, cov = full
                    scores = full_score
                    draw_status = "READY"
                    draw.update(sampling_status="FULL_REFERENCE", available_bases=int(pool.total_bases),
                                tolerance_accepted=True)
                matrix[number], covered[number] = data, cov
                draw.update(scores, status=draw_status, selected_units=len(ids),
                            selected_reads=int(pool.unit_read_counts[ids].sum()),
                            selected_fragments=len(np.unique(pool.unit_fragment_ids[ids])),
                            selected_segments=int(pool.unit_segment_counts[ids].sum()),
                            selected_blocks=int((pool.indptr[ids+1] - pool.indptr[ids]).sum()),
                            sampling_delta_bases=int(data.sum())-int(draw["target_bases"]),
                            overshoot_bases=max(0, int(data.sum())-int(draw["target_bases"])) if full is None else None,
                            target_depth=int(draw["target_bases"])/layout.length if full is None else None,
                            actual_bases=int(data.sum()), actual_depth=int(data.sum()) / layout.length)
            audit["status"] = ("READY" if all(draws[i]["status"] == "READY" for i in slots) else "PARTIAL")
        except (OSError, ValueError, OverflowError) as error:
            audit.update(status="ERROR", error=str(error))
            for number in slots:
                draws[number].update(status="ERROR", error=str(error))
        finally:
            pool = None
        decode_records.append(audit)
    check_stamps(spec)
    summary.update(valid_draws=sum(d["status"] == "READY" for d in draws),
                   tolerated_undershoot_draws=sum(d["status"] == "READY" and
                                                    d.get("sampling_status") == "UNDERSHOOT" for d in draws),
                   tolerated_overshoot_draws=sum(d["status"] == "READY" and
                                                   d.get("sampling_status") == "OVERSHOOT" for d in draws),
                   out_of_tolerance_draws=sum(d["status"] == "OUT_OF_TOLERANCE" for d in draws),
                   sampling_tolerance_bp=config.get("sampling_tolerance_bp", 150),
                   unique_references=len(groups), duplicate_draws=count-len(groups))
    complete = summary["valid_draws"] == count
    nulls = reference_row_means(matrix, config["threads"], covered=covered, lengths=lengths,
                               threshold=config["jaccard_coverage_threshold"]) if complete else {}
    # Preserve the v2 distribution mean, population SD, lower normal tail and AND decisions.
    stats = aggregate(draws, nulls, "low", count, stable_tail=True)
    for metric, result in stats.items():
        for source, suffix in (("value", ""), ("mean", "_mean"), ("max", "_max"), ("p", "_p"),
                               ("status", "_status"), ("null_count", "_null_count"),
                               ("null_mean", "_null_mean"), ("null_sd", "_null_sd")):
            summary[metric + suffix] = result[source]
        summary[metric + "_decision"] = ("NOT_EVALUATED" if result["p"] is None else
                                         "PASS" if result["p"] >= config["alpha"] else "FAIL")
    decisions = [summary[k + "_decision"] for k in ("jaccard", "pearson")]
    if "NOT_EVALUATED" not in decisions:
        summary["detection_status"] = "DETECTED" if decisions == ["PASS", "PASS"] else "NOT_DETECTED"
    hard_error = any(d["status"] == "ERROR" for d in draws)
    summary.update(status="COMPLETE" if complete else "ERROR" if hard_error else "COMPLETE_WITH_GAPS",
                   error="" if complete else "Reference draws not all within tolerance or failed; inspect draw audit")
    np.save(output / "reference_windows.npy", matrix, allow_pickle=False)
    np.save(output / "reference_covered_bases.npy", covered, allow_pickle=False)
    table(output / "draws.tsv", DRAWS, draws)
    table(output / "null.tsv", ("metric", "left_draw", "value"),
          (dict(metric=k, left_draw=i, value=v) for k, values in nulls.items() for i, v in enumerate(values)))
    dump(output / "statistics.json", stats)
    dump(output / "decode_audit.json", decode_records)
    query_side_covered = (np.stack([side_breadths[a.name] for a in sampling_arms])
                          if len(sampling_arms) == 2 else breadth)
    legacy_outputs(output, summary, draws, query, vector, matrix, query_side_covered, covered,
                   edges, config, nulls)
    if summary["status"] == "COMPLETE" and summary["detection_status"] == "NOT_EVALUATED":
        summary["status"] = "COMPLETE_WITH_GAPS"
    return summary


def score(args):
    output = Path(args.output)
    output.mkdir()
    spec = json.loads(Path(args.spec).read_text())
    if args.stub:
        summary = dict(query_id=spec["query_id"], Sample=spec["query"].get("Sample", ""),
                       GCFID=spec["query"].get("GCFID", ""), status="STUB",
                       detection_status="NOT_EVALUATED")
    else:
        try:
            summary = calculate(spec, output)
        except (OSError, ValueError, OverflowError) as error:
            summary = dict(query_id=spec["query_id"], status="ERROR", error=str(error))
    if not (output / "decision_methods.tsv").exists():
        unscored_outputs(output, summary, spec["config"])
    dump(output / "summary.json", summary)
    table(output / "summary.tsv", SUMMARY, [summary])
    if summary["status"] == "ERROR":
        audit = getattr(args, 'failure_audit_dir', None)
        if audit is not None:
            audit = Path(audit)
            audit.mkdir(parents=True, exist_ok=True)
            for name in ('summary.json', 'summary.tsv', 'draws.tsv', 'decode_audit.json', 'statistics.json', 'replacement_audit.json', 'replacement_audit.tsv', 'replacement_capacity_audit.json'):
                source = output / name
                if source.is_file():
                    shutil.copyfile(source, audit / name)
        print(summary.get("error", "Query failed"), file=sys.stderr)
        return 1
    return 0


def collect(args):
    methods = parse_methods(getattr(args, "jaccard_methods", "normal"))
    roster = json.loads(Path(args.roster).read_text())
    rows = []
    for query in roster:
        summary = Path(args.results_root) / query["query_id"] / "summary.json"
        rows.append(json.loads(summary.read_text()) if summary.is_file() else
                    dict(query, status="ERROR", error="See Nextflow task log/work directory",
                         detection_status="NOT_EVALUATED"))
    output = Path(args.output)
    output.mkdir()
    for row in rows:
        failed = row["status"] in ("TASK_FAILED", "ERROR")
        row.update(query_exit_code=1 if failed else 0, result_source="NO_RECEIPT" if failed else "SCORE")
    table(output / "summary.tsv", SUMMARY, rows)
    dump(output / "summary.json", rows)
    dump(output / "counts.json", dict(queries=len(rows), missing=sum(r["result_source"] == "NO_RECEIPT" for r in rows)))
    method_rows, sensitivity_rows = [], []
    for row in rows:
        directory = Path(args.results_root) / row["query_id"]
        saved = directory / "decision_methods.tsv"
        method_rows.extend(read_rows(saved) if saved.is_file() else [
            dict(query_id=row["query_id"], Sample=row.get("Sample", ""), GCFID=row.get("GCFID", ""),
                 status=row["status"], method=method, method_status="missing_query_method_result",
                 jaccard_decision="NOT_EVALUATED", pearson_decision="NOT_EVALUATED",
                 detection_status="NOT_EVALUATED") for method in methods])
        saved = directory / "sensitivity.tsv"
        sensitivity_rows.extend(read_rows(saved) if saved.is_file() else [
            dict(row, jaccard_coverage_threshold=tau) for tau in thresholds(
                getattr(args, "jaccard_coverage_threshold", 0),
                getattr(args, "jaccard_coverage_thresholds", None))])
    table(output / "decision_methods.tsv", METHOD_COLUMNS, method_rows)
    table(output / "sensitivity.tsv", SUMMARY, sensitivity_rows)
    errors = [r for r in rows if r["status"] in ("TASK_FAILED", "ERROR")]
    for name, subset in (
            ("errors", errors),
            ("negative", [r for r in rows if r not in errors and r["detection_status"] == "NOT_DETECTED"]),
            ("unscored", [r for r in rows if r not in errors and r["detection_status"] == "NOT_EVALUATED"])):
        table(output / (name + ".tsv"), SUMMARY, subset)
    counts = {state: sum(r["detection_status"] == state and r not in errors for r in rows)
              for state in ("DETECTED", "NOT_DETECTED", "NOT_EVALUATED")}
    dump(output / "completion.json", dict(
        queries=len(rows), receipts=len(rows)-len(errors), missing_receipts=len(errors),
        completion_status="PARTIAL_FAILURE" if errors else "COMPLETE",
        qc_status="EXECUTION_ERRORS" if errors else "HAS_UNRESOLVED" if counts["NOT_EVALUATED"] else "COMPLETE",
        ERROR=len(errors), **counts))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--version", action="version", version=VERSION)
    sub = parser.add_subparsers(dest="command", required=True)
    p = sub.add_parser("plan")
    for name in ("manifest", "index", "output"):
        p.add_argument("--" + name, type=Path, required=True)
    group = p.add_mutually_exclusive_group()
    group.add_argument("--fai-map", type=Path)
    group.add_argument("--fai-dir", type=Path)
    p.add_argument("--reference-draws", type=int, default=100)
    p.add_argument("--sampling-tolerance-bp", type=int, default=150)
    p.add_argument("--geometry-cache-dir", type=Path)
    p.add_argument("--insufficient-reference-policy", choices=("unscored", "replace"), default="unscored")
    p.add_argument("--alpha", type=float, required=True)
    p.add_argument("--seed", type=int, default=401)
    p.add_argument("--depth-threshold", type=float, default=5)
    group = p.add_mutually_exclusive_group(required=True)
    group.add_argument("--bin-count", type=int)
    group.add_argument("--window-size", type=int)
    p.add_argument("--sampling-unit", choices=("read", "fragment"), default="read")
    p.add_argument("--min-mapq", type=int, default=0)
    p.add_argument("--min-samples", type=int, default=10)
    p.add_argument("--threads", type=int, default=1)
    p.add_argument("--jaccard-coverage-threshold", type=float, default=0)
    p.add_argument("--coverage-unit", choices=("fraction", "percent"))
    p.add_argument("--coverage-alignment-policy", choices=("primary_only", "primary_plus_supplementary"),
                   default="primary_only")
    p.add_argument("--bootstrap-replicates", type=int, default=1000)
    p.add_argument("--bootstrap-confidence", type=float, default=.95)
    s = sub.add_parser("score")
    s.add_argument("--spec", type=Path, required=True)
    s.add_argument("--output", type=Path, required=True)
    s.add_argument("--stub", action="store_true")
    s.add_argument("--failure-audit-dir", type=Path)
    c = sub.add_parser("collect")
    c.add_argument("--roster", type=Path, required=True)
    c.add_argument("--results-root", type=Path, required=True)
    c.add_argument("--output", type=Path, required=True)
    for stage in (p, c):
        stage.add_argument("--jaccard-methods", default="normal")
        stage.add_argument("--jaccard-coverage-thresholds")
    c.add_argument("--jaccard-coverage-threshold", type=float, default=0)
    args = parser.parse_args()
    try:
        if args.command == "plan":
            if (args.reference_draws < 2 or args.sampling_tolerance_bp < 0 or args.min_samples < 0 or args.threads < 1 or args.min_mapq < 0
                    or not math.isfinite(args.alpha) or not 0 < args.alpha < 1
                    or not math.isfinite(args.depth_threshold) or args.depth_threshold < 0
                    or not math.isfinite(args.jaccard_coverage_threshold)
                    or not 0 <= args.jaccard_coverage_threshold <= 1
                    or (args.bin_count is not None and args.bin_count < 1)
                    or (args.window_size is not None and args.window_size < 1)):
                raise ValueError("Invalid draw, depth, window, threshold or CPU parameter")
            plan(args)
        elif args.command == "score":
            return score(args)
        else:
            collect(args)
    except (OSError, ValueError, sqlite3.Error) as error:
        parser.exit(1, f"ERROR: {error}\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())

