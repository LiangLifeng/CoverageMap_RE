# Manifest query v2

Reference pool size N is configurable: Python `--reference-draws N`; Nextflow
`--reference_draws N`. Default **100**; integer **N >= 2**. Both arms use N.
This controls reference draws, separately from `--bootstrap_replicates`.
N enters the run fingerprint and Nextflow scoring command; changing it requires
new scoring and invalidates incompatible resume results. Default and explicit
100 preserve the previous sampling/statistical formulas.

Use the independent `runs/native-resume-reference-draws-20261113/code/main.nf`
deployment for configurable N. The older `runs/native-resume-20261113/code`
does not support this Nextflow parameter. Merely adding an unknown parameter
to an older Nextflow entry does not enable the feature. Startup logs now print
the source directory, N and the separate bootstrap count. The task wrapper
rejects a successful backend result whose `requested_draws` differs from N or
whose COMPLETE `valid_draws` differs from N, writes an ERROR receipt and exits
nonzero so the mismatch cannot be cached as a successful query.

`bin/query_manifest_v2.py` compares each query with its own pool of N same-GCF reference coverage maps. Use the existing native SQLite index (including verified legacy per-contig geometry, rebinned read-only to GCF-concat windows) and TSV columns `Sample`, `GCFID`, `Depth`, and `CramPath`. Optional nonempty `Coverage` requires an explicit `--coverage-unit`. Choose `--window-size` or `--bin-count` using the [existing manifest window convention](QUERY_MANIFEST_V1.md#original-jaccardpearson-mode).

## Reference pool

- Remove candidates with the same Sample or the same input-file SHA-256 as the query.
- When the manifest depth is at or below `--depth-threshold` (default 5×), draw N sources with replacement from references with enough aligned bases to match the query. Independently downsample each draw to the query's measured aligned-base mass, retaining the final complete sampling unit.
- Above the depth threshold, draw N eligible sources with replacement and use their full coverage maps. A reference can occupy more than one draw slot in either branch. The eligible source count must be strictly greater than `--min-samples` (default 10).

## Similarity and detection

For both Jaccard similarity (`1 − Jaccard distance` on covered windows) and signed Pearson R on window coverage, the query statistic is the mean of its N pairwise similarities. Each reference draw contributes a mean against the other N−1 draw slots. Omit its own slot and retain comparisons with other slots drawn from the same source.

For each metric, the lower-tail normal approximation evaluates the query mean against the mean and population standard deviation of the N reference means. Pass `--alpha` explicitly. Both p-values must be at least alpha to produce `DETECTED`. If both p-values are available and either is below alpha, the result is `NOT_DETECTED`. An undefined metric, zero reference variance, zero query depth, insufficient references, or an incomplete pool yields `NOT_EVALUATED`. This exploratory test uses a shared, dependent reference distribution. The p-values and calls require independent calibration for interpretation as detection performance.

`summary.tsv` contains each metric's query mean, p-value and decision as well as the combined call. Each complete draw pool has `draws.tsv.gz` for the query-to-reference comparisons and `null.tsv.gz` for the N reference means per metric, with draw and reference IDs. The original window arrays, checksums, run identity and resume behavior remain available through the shared runner.

Use the Python environment that runs the existing manifest CLI. Set `MANIFEST`, `INDEX`, `OUTDIR`, `BINS`, and `ALPHA` for the intended evaluation.

```bash
python -B /absolute/path/to/bin/query_manifest_v2.py \
  --manifest "$MANIFEST" --index "$INDEX" --outdir "$OUTDIR" \
  --bin-count "$BINS" --alpha "$ALPHA"
```

## CPU parallelism

`--threads` controls reference loading, coverage binning, simulations and pair scoring within each query. It defaults to two and requires at least two. Each distinct reference loads once per query in batches bounded by the thread count. Repeated low-depth draws simulate independently, and repeated high-depth draws reuse the loaded coverage map. `--jobs` defaults to one and runs different GCFs in separate processes in batch mode. Allocate `jobs × threads` CPU slots. NumPy BLAS threads are fixed to one.

The reference calculation prepares packed coverage masks and Pearson vectors once, scores the 4,950 unordered pairs, and adds scores in a fixed order. The output arrays and decisions use the same draw order and arithmetic at every thread count.

### One query per task

Pass `--query-id q00000003 --jobs 1 --threads 4` to run the third data row of the original manifest with four CPU slots. IDs count data rows after the header and retain leading zeros. The original manifest row number is used for the random seed, so a task and the corresponding row of a batch run select the same references and sampling units. Give each task its own output directory. A Nextflow process can use the original manifest and one query ID per input tuple, and pass its allocated CPUs as `--threads`.

## Nextflow batch workflow

Use Nextflow 26.04 or later with the Conda executable available on `PATH`. The workflow enables Conda globally and activates the installed `/home/lianglifeng/01.WorkDir/env/conda-envs/coverage-map` environment for every process. Use `--coverage_conda_env /absolute/path/to/existing/environment` on systems with a different prefix. The environment and input files must be readable from the executor. The manifest's alignment paths are absolute. Keep the index and its registered reference payloads immutable during the workflow.

```bash
NXF_SYNTAX_PARSER=v2 nextflow run /absolute/path/to/main.nf -profile local \
  --manifest "$MANIFEST" --index "$INDEX" --outdir "$OUTDIR" \
  --bin_count "$BINS" --alpha "$ALPHA" --query_cpus 4
```

The `slurm` profile submits each task to `cpu,tai`, allowing Slurm to select one of those partitions. Set any site-specific account requirements in the launch config. `--query_memory`, `--query_time`, `--cache_mib` and `--work_mib` control query resources. Use `--window_size` instead of `--bin_count` for fixed windows and supply `--coverage_unit` when needed.

Use `-profile slurm,nhb` on NHB. Its Slurm partitions advertise 1 MB of memory per node, so the `nhb` profile omits scheduler memory requests. The program still applies `--cache_mib` and `--work_mib` to its cache and numeric work.

The planner creates one stable query ID per original data row, including rows with missing alignment files. Each score task receives its own ID and uses `task.cpus` threads. The collector matches compact task receipts to the original roster and writes `collected/summary.tsv` in input order.

| Output under `collected/` | Rows |
| --- | --- |
| `summary.tsv` | Every original query |
| `negative.tsv` | `NOT_DETECTED`, with both metrics evaluated |
| `errors.tsv` | Query errors and tasks with no receipt after retry |
| `unscored.tsv` | `NOT_EVALUATED` without an execution error |
| `completion.json` | Row counts by outcome and number of missing receipts |

Task logs and per-query scoring files stay in the Nextflow work directory. Use
the trace in `pipeline_info/` to locate the task by query ID. The score wrapper
preserves the original nonzero exit (signals become 128 + signal), writes its
diagnostic receipt and returns nonzero. A failed task's receipt is not emitted
as a successful Nextflow output. Exit 75, 137 or 143 receives at most one retry;
deterministic input/schema errors are not retried automatically. Independent
queries finish and the collector represents absent receipts as `NO_RECEIPT`,
`ERROR` and `NOT_EVALUATED` for every requested method. An all-failed cohort
still publishes its error ledger. `completion.json` explicitly records
`completion_status=PARTIAL_FAILURE` and `qc_status=EXECUTION_ERRORS` when needed.
`workflow.failOnIgnore=true` makes the workflow return nonzero **after** partial
collection, rather than canceling healthy queries. Optional diagnostic-rendering
failures also return nonzero, with the original code in the task trace and
`diagnostics/plot_exit_code.txt` in the failed task's work directory.

After fixing an external execution fault, use the same launch/work directories
and `-resume`: successful score tasks can be reused while failed tasks rerun.
Receipts are sorted by query ID before the collector's process boundary. The
planner deliberately runs each time to recheck the index and alignment-file
stat identities; its newly created roster may therefore rerun collection even
when score tasks are cached. This is separate from receipt-arrival ordering.
Changing a declared input or command parameter invalidates its dependent tasks.

Upgrade limitation: historical wrapper versions could cache `status=ERROR` as
a successful task. The changed wrapper/command is a one-time cache boundary;
do not delete or rewrite old work/cache to conceal it. Before a production
upgrade, inspect generated commands, traces and receipts to confirm failed
queries were actually rerun. Query IDs and RNG still depend on original row
numbers. Reordering a manifest changes that replay identity; this release does
not promise reorder-invariant draws. The complete manifest and staged Python
code collection still enter score-task hashes, so append/metadata/unrelated-code
changes can invalidate more tasks than scientifically necessary. Refining those
boundaries requires a separate interface change preserving the original row IDs.

## Coverage breadth threshold

See [Jaccard breadth](JACCARD_BREADTH.md). Default tau0 remains C>0; positive
thresholds and sensitivity mode require exact selected-unit covered-bases vectors.


## Query vNext sampling and exploratory Jaccard decisions (2026-11-12)

New Geometry vNext reference bundles keep Read, Fragment/Mate, Supplementary,
alignment segments and blocks independent of runtime sampling. Use
`--sampling-unit read` (independent complete Read units) or
`--sampling-unit fragment` (entire observed Fragment; R1/R2 together).
Older coordinate-only bundles must not claim dynamic Fragment capability.
The query CRAM uses the same sampling policy for its coverage, but it is **not
downsampled**. Low-depth Reference draws are sampled without replacement
**within each draw** by complete Read/Fragment to the measured query aligned-base
mass; the final unit is retained, even with depth overshoot. The selection of
N Reference sources is still **with replacement**, with independent per-draw
unit RNG. High-depth Reference draws remain **full-reference, unsampled**.
The manifest Depth threshold, Pearson, geometry reference, coverage alignment
policy and existing official AND decision are unchanged.

For a 1-bp window, use `--jaccard-coverage-threshold 0`: each window is a
single covered/not-covered base. Positive thresholds up to 1 generate exactly
the same binary states, so a breadth threshold scan adds no information.

Additional *exploratory* Jaccard decisions can be requested, e.g.:

```bash
python -B bin/query_manifest_v2.py --manifest "$MANIFEST" --index "$INDEX" \
  --outdir "$OUTDIR" --window-size 1 --alpha .05 \
  --sampling-unit fragment --jaccard-coverage-threshold 0 \
  --jaccard-methods normal,empirical,minimum,bootstrap \
  --bootstrap-replicates 1000 --bootstrap-confidence .95
```

The default `--jaccard-methods normal` preserves the official
`summary.tsv` P-value and the original Pearson/AND decision. A new
`decision_methods.tsv` in the run root and in each query directory records:

- `normal`: existing lower-tail normal P-value; SD=0 remains nonestimable.
- `empirical`: compare the **query mean Jaccard against all N draws**
  to the noninterpolated `ceil(alpha * N)`-th ascending **reference row mean**.
- `minimum`: compare against the minimum of the N reference row means.
  This is NOT the minimum of the 4,950 pairwise reference comparisons.
- `bootstrap`: resample distinct Reference **source groups**, retaining the
  saved draw slots assigned to each source, to measure variability in the
  empirical cutoff. Report point cutoff, confidence limits and bootstrap pass
  fraction. PASS requires query >= upper bootstrap cutoff; FAIL requires query
  < lower cutoff; otherwise UNCERTAIN/NOT_EVALUATED. This is a *conditional
  stability diagnostic*, not a calibrated independent-sample prediction limit.

All four methods re-use the identical N Reference draws and Pearson P-value.
Incomplete, zero-union or missing reference row means are never dropped to make
the other methods appear estimable. A constant null can be evaluated using the
empirical and minimum *ordering* even when the normal P-value is undefined.
The N row means are dependent and may repeat Reference source identities;
do not infer nominal frequentist error rates from these bounds without
independent TP/FP validation. Compare TPR and FPR jointly.

To re-score **existing completed Query outputs**, without CRAM decoding or
new random draws, use:

```bash
python -B -m bin.query_decision_methods \
  --run-dir "/absolute/previous-query-run" \
  --out "/absolute/new/method_comparison.tsv" \
  --methods normal,empirical,minimum,bootstrap --alpha .05
```

A completed run retains its original results and its original checksums.
Offline output must use a new file path, not overwrite the run's frozen tables.
The Nextflow interface supports `--jaccard_methods`,
`--bootstrap_replicates`, and `--bootstrap_confidence`. Its collector
writes `collected/decision_methods.tsv`; missing task receipts are represented
as NOT_EVALUATED, not silently excluded.

The entry workflow validates method names/duplicates, bootstrap option bounds
and explicit `--diagnostics true|false` before planning tasks. With `--window_size
1`, use scalar `--jaccard_coverage_threshold 0` without a redundant breadth grid.
If G is the number of bases/windows, an N-draw int64 matrix uses `8*N*G` bytes.
Both aligned and covered-base matrices are retained, and statistics, query
vectors, geometry cache, staging and export require additional memory/storage.
The work guard's estimate is not a measured process-memory cap. Validate the
largest GCF, full query entrypoint and concurrent I/O before choosing queue size.

Native wrapper regression (synthetic scorer, no BAM decoding or calibration):

```bash
COVERAGEMAP_TEST_NEXTFLOW=/absolute/path/to/nextflow \
  python -B -m pytest -q tests/test_query_flow_nextflow.py
```

This checks actual `.command.sh`, `.exitcode` and traces for deterministic errors,
signals, retry exhaustion, missing output, all-failed/partial collection,
diagnostic failures and a three-run failure/recovery/cache sequence. The normal
Python suite skips these tests unless the executable is explicitly supplied.


After obtaining separate TP and FP decision tables from the same frozen
experimental design, calculate full-denominator sensitivity and false-positive
rates (including NOT_EVALUATED in denominators):

```bash
python -B -m bin.query_method_report \
  --tp /absolute/tp/decision_methods.tsv \
  --fp /absolute/fp/decision_methods.tsv \
  --out /absolute/new/tp_fp_methods.tsv
```

The report rejects different method rosters within either cohort and unequal
method sets across TP/FP. It shows DETECTED, NOT_DETECTED, NOT_EVALUATED,
TPR/FPR over all queries and evaluated-only rates separately, and changes
relative to Normal. It does not claim a 99% target has been achieved unless
the independent cohort data show it, nor does it treat samples as patients.
