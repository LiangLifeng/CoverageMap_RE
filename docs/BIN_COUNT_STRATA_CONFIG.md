# Bin Count stratum configuration and runtime diagnostics

The `plan` and `run` commands accept `--strata-config PATH`. A config
defines an ordered, contiguous partition and an explicit list of positive
integer nominal bin widths for each layer.

```json
{
  "strata": [
    {
      "min_genome_bp": 1000,
      "max_genome_bp_exclusive": 10000,
      "nominal_bin_widths_bp": [110, 125, 150, 175]
    },
    {
      "min_genome_bp": 10000,
      "max_genome_bp_exclusive": 100000,
      "nominal_bin_widths_bp": [110, 125, 150, 175, 200]
    },
    {
      "min_genome_bp": 100000,
      "max_genome_bp_exclusive": null,
      "nominal_bin_widths_bp": [110, 125, 150, 175, 200, 300]
    }
  ]
}
```

Lower bounds are inclusive. The final upper bound must be null. A config
replaces the legacy stratum and width arguments for that command. Mixing
the two interfaces is an error. Legacy invocations and existing plans and
matrices remain readable. The resolved config is recorded in the plan.
Bin counts still use floor(shortest genome length / nominal width). Widths
that yield the same count are recorded as aliases. The row order, bin edge
policy, normalization and scoring are unchanged.

For cutoff scans, `--ordination none` runs scoring without building a full
sample pairwise distance matrix. `--jobs N` defaults to one and is accepted
by `build`, `evaluate`, and `run`.

- Build processes up to N sample–GCF records concurrently. Each worker
  verifies and reads one geometry at a time and owns one row in every
  candidate matrix. The number of pending rows is bounded, and workers
  finish before the mappings close. Choose N to fit memory and shared
  storage throughput. The disk limit still covers matrix payload bytes.
- Evaluation computes GCF centers and within-group distances in parallel,
  with deterministic ordering and aggregation. Its working-memory guard
  reserves additional scratch space for concurrently scored groups.
- For ordination with more than 64 rows, N greater than one enables tiled
  pairwise distances for all four metrics. Tiles use the serial argument
  direction, reflect the upper triangle and preserve a zero diagonal.
  Pearson centers profiles once before applying cosine distance, matching
  SciPy's correlation kernel.

Build caches the integer boundary vectors once per GCF and candidate.
Stderr reports per-stratum row-loop wall time, summed worker time for
geometry and binning, and explicit flush time. Worker times can overlap.
Evaluation separates pairwise compute and write time by candidate and
distance. Result tables remain in deterministic order.
