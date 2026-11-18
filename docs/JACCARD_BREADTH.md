# Jaccard coverage breadth

Issue #8 adds breadth scoring to both manifest query CLIs and their Nextflow
entry point. Each runtime selected complete Read/Fragment contributes its exact
alignment blocks. The per-window union C counts each reference position once,
including overlapping mates, reads and CIGAR blocks. No fragment span is invented.
C is bounded by the actual window length L, including short tails and windows
crossing virtual contig boundaries. All window counts sum to the full union.

`tau=0` means C>0. Positive tau means C/L>=tau, with no rounded percentages.
Both query/reference comparisons and every off-diagonal reference comparison
use the same state function and tau. Default tau remains zero. Pearson uses
unchanged aligned sums. Sampling choices, unit stopping, seeds, routing, alpha
and the AND decision are unchanged. The original Tools.py and historical
concordance replay interfaces remain frozen historical implementations; their
old outputs cannot be interpreted as positive-threshold results.

## Use

```bash
query_manifest_v2.py --manifest /absolute/queries.tsv --index /absolute/index.sqlite \
  --outdir /absolute/new-output --bin-count 60 --coverage-unit percent --alpha .05 \
  --jaccard-coverage-threshold 0 --jaccard-coverage-thresholds default
```

The grid is 0,.01,.05,.10,.15,.20,.30. A comma-separated custom grid is also
accepted, with finite fractions in [0,1]. The primary scalar threshold is always
included first; duplicates are removed without sorting. Omission of the grid
scores only the primary threshold. Nextflow uses `--jaccard_coverage_threshold`
and `--jaccard_coverage_thresholds`. Thresholds affect scoring tasks; the
independent reference builder is unchanged.

Every query is decoded and every draw sampled once. All thresholds reuse those
same matrices and draw IDs. Default single-threshold files retain one primary
summary row. `sensitivity.tsv` contains every query/threshold, including skipped
and failed queries. Nextflow receipts/collector preserve these rows; missing
receipts remain errors for all requested thresholds.

## Outputs and invalid comparisons

- `query_windows.npy`, `reference_windows.npy`: original aligned sums for Pearson.
- `query_covered_bases.npy`, `reference_covered_bases.npy`, `window_lengths.npy`:
  exact int64 breadth inputs, stored once and shared across thresholds.
- `breadth.tsv.gz`: streaming per-bin covered bases, length and breadth, with
  query/draw ID. `windows.tsv` and its SHA identify coordinates.
- `jaccard_thresholds/tNN/`: threshold summary, query draw I/U and scores,
  reference row means in `null.tsv.gz`, and `zero_unions.tsv.gz` with both slot
  identities and reference IDs for every undefined pair. The primary threshold
  also has the conventional top-level draws/null files.
- Summary zero-union counts and fractions use ready query pairs and enumerated
  reference pairs as denominators. Missing pair populations yield blank ratios.
- Input hashes, source identities, runtime policies, selection tokens and per-draw
  selected-ID/RNG hashes remain in run/audit/draw outputs and recursive checksums.

Union=0 produces blank Jaccard, `zero_union`, and invalidates the corresponding
reference row mean. Incomplete query draws or null rows are non-estimable;
no rows are dropped to create a smaller valid denominator. V2 yields
NOT_EVALUATED when either metric is non-estimable. Constant/null-SD handling and
normal-tail calculation otherwise remain unchanged. Only zero-union handling
intentionally differs at tau0; all finite tau0 scores must match historical output.

Use `jaccard_sensitivity_report.py --tp /absolute/tp/sensitivity.tsv
--fp /absolute/fp/sensitivity.tsv --outdir /absolute/new-report` for a table,
SVG and HTML comparison. It checks exact cohort identity and Pearson invariance
across thresholds, retains all rows in denominators, and does not select a best
threshold. These are sample–GCF query counts, not unique patient counts.

## Resources and acceptance

For M selected blocks, union computation is O(M log M) time/O(M) scratch, no
per-genome array. Full high references are cached per unique source/window/policy;
low draws union only selected blocks. Batching and threads retain existing CPU
budgets. Numeric window budget now reserves five aligned-matrix equivalents
for aligned/covered vectors, states and statistics scratch. Geometry/cache RSS
is separate and observed on the target cluster. Seven thresholds store matrices
once and score 7*N*(N-1)/2 pairs; zero-pair audits are streamed.

Current validation design: user-confirmed TP344 nonreference FeverPatients and
FP334 InHouseCtrlSample, same viral G1000_S10 index, B60/min-samples9,
alpha.05/depth5/seed401,100 draws. B60/N100 requires 48,000 bytes per matrix;
678 queries have 23,492,700 reference pair scores across seven thresholds.
Real I/O/RSS/time and full performance acceptance require cluster measurements;
synthetic tests do not establish cohort throughput or clinical calibration.
