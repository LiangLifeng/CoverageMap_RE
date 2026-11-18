# CoverageMap_RE Raw Query vNext 1.2.1

Development baseline: `db11a9f21700178b8571566970480a1f8d7dc3e4`.
Feature branch: `codex/raw-query-vnext-integration-20261118`.
Only CoverageMap_RE receives new development; CoverageMap remains historical.
This is a local engineering handoff. Native Nextflow/NHB/data acceptance is pending.

## Integrated behavior

The FNA/FAI SQL catalog selects strict FASTA-backed validation. FASTA and FAI
contig names/lengths must agree and alignment @SQ must match FASTA SN/LN/M5.
The older CRAM-only SQL plus explicit FAI retains its historical sampler; it does
not silently enable the new algorithm. Use the catalog SQL for vNext.

Low-depth read mode partitions a fully mate-labelled query into Read1/Read2 arms,
with floor(N/2) and remainder slots. Each arm uses its own actual aligned bases;
query depth is not blindly halved. Single/mixed-mate queries remain unstratified.
Fragment mode preserves observed whole units and mate blocks, without an invented
span or missing mate. High-depth mode reuses full reference coverage.

Each low-depth slot selects the closest nonempty permutation prefix, without
replacement or unit truncation. Under/overshoots within the default 150 bp
absolute tolerance are retained and audited. Out-of-tolerance slots remain in the
QC arrays but cannot yield an accepted statistical distribution. No redraw,
upsampling, missing-slot removal or smaller denominator is introduced.

All N slots retain the existing pooled off-diagonal reference row-mean null,
query mean and Normal/Pearson AND. Mate arms do not introduce a new stratified
null model. Optional normal/empirical/minimum/bootstrap and threshold grids reuse
the same frozen arrays. Statistical calibration has not been established.

The output adapter now preserves mate/tolerance fields in both main and threshold
draw tables. Diagnostics show full query coverage plus slot-weighted query mate
coverage, and can render nearly constant score distributions without a histogram
error. Undefined scores/variance retain NOT_EVALUATED.

Shared cache schema2 stores only geometry/relationships and sampling-unit metadata.
Keys bind source path/stat identity, reference layout/M5 and extraction policies.
Per-key file locks and atomic replace protect concurrent users. Invalid or corrupt
entries are rebuilt; a source change during extraction prevents publication.
Cold/warm results must have identical scientific arrays, statistics and decisions.
Sources must remain immutable; this is stat-based identity, not full CRAM hashing.
The cache is optional/disposable and separate from Nextflow's query resume cache.

## Local verification

Python3.12.14, NumPy2.3.5, pysam0.24.1, pytest9.2.0.
Matplotlib3.11.2/SciPy1.18.1 used for plotting/oracles.
Focused suite: 94 passed; one optional native Nextflow test skipped locally.
Includes BAM/CRAM, Read/Fragment, low/high, odd N, 1bp dense selected-ID coverage,
undershoot/tolerance rejection, strict M5 mismatch, concurrent/corrupt/changed-source
cache, cold/warm equivalence, methods/grid, and four-method rendered Diagnostics.
Compileall, diff whitespace and launcher bash syntax pass. A synthetic minimum
PNG was visually inspected.

A broad historical suite initially returned 479 passed, 14 failed, 24 skipped,
35 subtests passed. Twelve failures were missing scikit-learn; after installing
the declared local dependency, all 31 Bin Count tests and 8 subtests passed.
Two unrelated historical tests also fail on pristine db11a9f: the legacy native
Slurm-shell mock and the REFERENCE_CHUNK cache=false assertion. They are not fixed
or hidden here. No complete-repository green claim is made.

## Prepared NHB defaults and real command

`raw_query/nextflow.config` defaults already point to the prepared environment:

- coverage_conda_env: `/home/lianglifeng/01.WorkDir/env/conda-envs/coverage-map`
- diagnostics_bin_dir: `/home/lianglifeng/01.WorkDir/env/conda-envs/coverage-map/bin`

The supplied launcher uses the user's SQL and 344-row TP manifest, verified
read-only for schema/catalog presence and manifest columns/count. Individual
CRAM/FASTA/FAI availability, sequence consistency and node mounts are not yet
verified. No clinical records or environment contents were returned.

Source directory:
`/home/lianglifeng/01.WorkDir/01.Project/03.CoverageMap/runs/coveragemap-re-vnext-20261118/code`

From the user's prepared NHB login environment:

```bash
bash /home/lianglifeng/01.WorkDir/01.Project/03.CoverageMap/runs/coveragemap-re-vnext-20261118/code/examples/nhb/run_query_tp.sh -bg
```

Recovery uses the same launch/work/output paths:

```bash
bash /home/lianglifeng/01.WorkDir/01.Project/03.CoverageMap/runs/coveragemap-re-vnext-20261118/code/examples/nhb/run_query_tp.sh -resume -bg
```

The expanded command is in `examples/nhb/run_query_tp.sh`. It retains the prior
launcher settings: Read/N100/B60/alpha.05/min_samples9/depth_threshold5/seed401,
4 query CPUs, minimum method, primary_only, percent Coverage, Diagnostics on.
`minimum` does not invoke bootstrap; summary.tsv still reports the Normal result,
and minimum decisions live in collected/decision_methods.tsv and Diagnostics.
The shell accepts trailing native Nextflow options/parameter overrides.

NHB node scratch/TMPDIR is explicitly /disk in examples/nhb/raw_query.config.
Slurm uses cpu,tai, queue size100 and no explicit memory/time requests.
Shared geometry persists under the independent run root, outside scratch/work.
No old source, reference library, SQL, user job or output directory is replaced.

## User-run engineering acceptance before real data

Run the focused Python suite in the prepared NHB environment. Native DAG test
uses only synthetic inputs and checks methods/Diagnostics, query cache resume and
Fragment re-execution:

```bash
cd /home/lianglifeng/01.WorkDir/01.Project/03.CoverageMap/runs/coveragemap-re-vnext-20261118/code
COVERAGEMAP_TEST_NEXTFLOW=/home/lianglifeng/01.WorkDir/env/app/nextflow-official-26.04.6/bin/nextflow /home/lianglifeng/01.WorkDir/env/conda-envs/coverage-map/bin/python -m pytest -q tests/test_raw_query_nextflow.py --basetemp=/home/lianglifeng/01.WorkDir/01.Project/03.CoverageMap/runs/coveragemap-re-vnext-20261118/synthetic_acceptance
```

Use that basetemp only for this disposable synthetic test; pytest clears its own
basetemp on rerun. Keep real data/results/work directories separate.

Check real outputs: collected/completion.json, errors.tsv, unscored.tsv,
decision_methods.tsv and diagnostics/plot_index.tsv. Per-query draws include
mate_side, sampling_status, sampling_delta_bases, tolerance_accepted and
geometry_cache_hit. Source errors are nonzero task failures; out-of-tolerance
or undefined statistics remain explicit unscored results. Keep launch directory
and work for -resume. The agent has not run formal NHB queries or Slurm jobs.
