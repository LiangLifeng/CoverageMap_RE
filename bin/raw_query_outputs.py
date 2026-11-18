"""Legacy query feature adapters over the same in-memory raw-CRAM draws."""
import csv
from pathlib import Path

import numpy as np
from bin.jaccard_sensitivity import write_sensitivity
from bin.query_decision_methods import COLUMNS as METHOD_COLUMNS, evaluate
from bin.query_manifest import DRAW_COLUMNS, NULL_COLUMNS, DISTRIBUTION_COLUMNS, table


RAW_DRAW_COLUMNS = (*DRAW_COLUMNS, "mate_side", "sampling_status", "sampling_delta_bases",
                    "available_bases", "tolerance_bp", "tolerance_accepted", "geometry_cache_hit")


def legacy_outputs(output, summary, draws, query, test, matrix, breadth, covered,
                   edges, config, primary_nulls):
    """Threshold grids/methods consume frozen draws; never re-decode or re-sample."""
    canonical_draws = [{**{key: row.get(key) for key in DRAW_COLUMNS},
                        **{key: row.get(key) for key in RAW_DRAW_COLUMNS[len(DRAW_COLUMNS):]}}
                       for row in draws]
    canonical_summary = {key: summary.get(key) for key in DISTRIBUTION_COLUMNS}
    numeric_query = dict(query, query_id=summary["query_id"])
    extra, *_ = write_sensitivity(output, canonical_summary, canonical_draws,
                                  test, matrix, breadth, covered, edges, config,
                                  numeric_query, primary_nulls, draw_columns=RAW_DRAW_COLUMNS)
    summary.update(extra)
    null_rows = [dict(metric=key, left_draw=i, right_draw="", value=value,
                      left_reference_id=draws[i]["reference_id"],
                      left_reference_sample=draws[i]["reference_sample"])
                 for key, values in primary_nulls.items() for i, value in enumerate(values)]
    table(output / "draws.tsv.gz", RAW_DRAW_COLUMNS, canonical_draws)
    table(output / "null.tsv.gz", NULL_COLUMNS, null_rows)
    methods = evaluate(summary, canonical_draws, null_rows,
                       methods=config.get("jaccard_methods", "normal"),
                       alpha=config["alpha"],
                       bootstrap_replicates=config.get("bootstrap_replicates", 1000),
                       bootstrap_confidence=config.get("bootstrap_confidence", .95),
                       seed=config["seed"])
    table(output / "decision_methods.tsv", METHOD_COLUMNS, methods)


def unscored_outputs(output, summary, config):
    """Reference support/zero depth are successful NOT_EVALUATED results."""
    from bin.jaccard_breadth import thresholds
    methods = evaluate(summary, [], [], methods=config.get("jaccard_methods", "normal"),
                       alpha=config["alpha"],
                       bootstrap_replicates=config.get("bootstrap_replicates", 1000),
                       bootstrap_confidence=config.get("bootstrap_confidence", .95), seed=config["seed"])
    table(output / "decision_methods.tsv", METHOD_COLUMNS, methods)
    rows = [dict({key: summary.get(key) for key in DISTRIBUTION_COLUMNS},
                 jaccard_coverage_threshold=tau) for tau in thresholds(
                     config["jaccard_coverage_threshold"], config.get("jaccard_coverage_thresholds"))]
    table(output / "sensitivity.tsv", DISTRIBUTION_COLUMNS, rows)


def read_rows(path):
    with Path(path).open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))

