# Query-depth coverage-profile comparison

This is the default native SQLite query policy as of 2026-10-23. Existing
alignment/bundle/legacy-depth modes are unchanged. `--comparison-policy nearest`
explicitly replays the previous native nearest-reference workflow and its depth
arguments. The reference library and SQLite are never written by queries.

## Query metric and windows

For each Sample/CramFilePath/GCFID row in `--query-manifest`, extract the query
once using the index's policies and complete FAI geometry. `--depth-metric mean`
(currently the only supported metric) is eligible primary M/=/X aligned-base mass
divided by complete FAI length, including uncovered/missing-header contigs.
Do not average unequal-length window means to obtain whole-genome mean depth.
A later median metric will require an explicit definition and index support.

On `gcf-concat-v1` indexes, `--window-size 1000` uses the whole GCF and `--bin-count B` creates exactly B balanced windows. References and queries are rebinned from verified geometry to the same edges. Statistics use aligned-base mass divided by actual window width. Legacy indexes retain fixed per-contig windows. See the [window contract](GCF_CONCAT_WINDOWS.md).

## High-depth route

Query original depth `>= --low-depth-threshold` (default5X) uses ALL selected,
positive-depth same-GCF references in the inclusive interval query mean depth
`+/- --depth-tolerance` (default1X). Index filtering, `--reference-samples` ID
filter occur before comparison. Native depth-band queries include same-sample ID
and identical-source-SHA references by default (user correction 2026-10-24).
`--exclude-self` explicitly excludes either match; `--include-self` is optional.
Nearest/legacy modes retain their prior exclusion default for replay compatibility.
No reference-limit truncation or query downsampling within the depth band.

If the band is empty, compare the nearest available positive-depth reference(s)
within the same GCF, keeping all exact ties in absolute depth difference. Explicit
index selection, ID filters and --exclude-self still apply; no excluded row is
restored. Output selection_mode=nearest_depth_fallback and depth_band_matched=false,
with band_candidate_count=0. Such outside-band comparisons retain descriptive
similarities and null p/acceptance (outside_requested_depth_band).

Only a GCF absent from the index is SKIPPED/NO_REFERENCE_GCF. A GCF present but
removed by index filtering is processed as COMPLETE_WITH_GAPS/NO_SELECTED_REFERENCE_GCF,
without restoring filtered rows. Zero/unsupported depth or an empty explicitly
filtered source set is also COMPLETE_WITH_GAPS, not a skipped query. Corrupt input
or reference files remain failures. No usable statistics are fabricated.

Follow `concordance_test_high_cov` in original Tools.py (import56c9df6): for each
metric, take the MAXIMUM query-to-reference similarity and apply the lower-tail
normal z comparison. Original Tools.py accepts externally prepared reference
similarity pools; this implementation explicitly constructs its pool from ALL
unordered, distinct reference-to-reference pairs within the selected band. It
uses population mean/std, no diagonal/self pairs. Moments are streamed with O(R)
scratch rather than retaining an R-by-R matrix. Full reference vectors occupy
O(R*W) numeric memory; cache budgets do not cap this working set.

Jaccard is positive-window intersection/union; Pearson compares numerical depth
profiles. There is no similarity-based reference selection. A single reference
can provide descriptive similarity but no reference distribution; constant
profiles and zero-variance distributions produce null p/acceptance and explicit
non_estimable status. These preserve accepted corrected semantics, not the
original epsilon/zero-variance corner cases. No statistical calibration is claimed.

## Low-depth route

Query original mean `<5X` uses multiple requested simulation target layers:

- `--simulation-depths 0.1,0.5,1` specifies X, not replicate indices. Default is
  the query's original mean, preserving the full query profile in that layer.
- `--replicates 100` generates N independent complete-unit draws PER layer PER
  selected reference sample. At least two repeats are required.
- `--simulation-references all`, `20` or `20%` selects source samples within the
  query's GCF. Counts are capped at the available count. Percentages use
  ceil(eligible count * percentage/100). Draw without replacement with fixed
  seed and GCF-keyed RNG, independent of row order, worker scheduling and layers.
  The same selected sources are reused across layers. Default is all.
- Source selection includes same-sample/source records and has NO depth restriction by default. Optional
  `--simulation-source-depth-range 5:100` restricts ORIGINAL source means to an
  inclusive range; `--simulation-match-depth-band` additionally applies the
  query's +/- depth-tolerance band. These filters precede the count/percent draw.
- `--reference-samples a,b` remains an optional explicit ID filter, independent
  of the count/percentage control.

References with insufficient mass for a layer receive UNSUPPORTED_REFERENCE_DEPTH;
layers above query mass receive UNSUPPORTED_QUERY_DEPTH. They are never upsampled,
silently clipped to a different target or replaced by another source. Use layers
at or below query depth. Zero query mass is ZERO_QUERY_DEPTH, not similarity evidence.

For supported layers below the original query depth, make one independently
seeded query draw at that depth and reuse it across the selected references.
Each reference is simulated separately; no biological sample pooling. Follow
Tools.py's LOW-coverage structure: mean query-to-simulation similarity versus
the off-diagonal simulation-to-simulation row-mean distribution and lower-tail
normal z calculation. At the exact full reference mass all whole-unit draws
are identical; retain N audit rows but mark degenerate statistics non_estimable.
Complete-unit draws may overshoot the target by a unit; record actual mass/depth.

## Outputs and integrity

`summary.tsv` records reference count/candidate count, depth metric, window size,
query/reference/target/actual depths, route, similarities, p/acceptance/status.
Each query has selection.json with candidates count, sampled IDs, original source
depths, root-resolved manifests/SHA, source-depth bounds and self-consistency scope;
query_mean_depth.npy, edges.npy, provenance/timing and terminal status.
High route writes high-depth/concordance.json with selection_mode/depth_band_matched, selected reference IDs and
pair-distribution definition/moments. Low route writes layer-*/test.npy/test_ids.npy
and per-reference concordance.json with all N replicate ID hashes and actual depths.
Simulations remain in bounded retained-cache storage; temporary N-by-W matrices
and statistical working sets are additional RAM. No persistent reference writes.

Index identity, query bytes and loaded reference file stamps/status are checked
before success. Payload hashes/registered metadata and layout/M5/policies are
validated on load. Native simulation cache keys include the actual window edges,
so changing profile windows cannot reuse an incompatible matrix.
Unknown M5 remains unverified, not evidence of sequence identity.

## Example

```bash
python -B /absolute/source/bin/coverage_concordance.py \
  --reference-sqlite /absolute/selected.sqlite \
  --query-manifest "$PWD/query_input.csv" \
  --window-size 1000 --depth-tolerance 1 \
  --simulation-depths 0.1,0.5,1 --replicates 100 \
  --simulation-references 20% --seed 401 \
  --query-jobs 4 --query-scheduling gcf \
  --outdir "$PWD/new_concordance_result"
```

The layer flags affect low-depth queries only. Replace 20% with all or a fixed
count. Use an absolute source/index path and a new output directory; source roots
registered in the SQLite must remain accessible. Indexed mode uses CPU/pysam,
`--workers 1`; query concurrency uses `--query-jobs`.

Acceptance: independently verify original Tools.py maximum/mean formulas on
identical nondegenerate vectors, inclusive +/-1X boundaries, exactly5X route,
source exclusion/filters/count/percentage/layers, unsupported mass, zero/missing/
degenerate data, window rebins, multi-root selection and serial/process equality.
Real-input cohort throughput and inferential calibration are outside this change.
