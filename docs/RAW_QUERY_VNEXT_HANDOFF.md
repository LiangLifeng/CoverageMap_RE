> Current v1.3.2: optional `--m5_policy geometry_only` accepts header M5
> mismatches only when contig SN/LN and the authoritative FASTA layout agree;
> source-level M5 audit is written for every decoded Query and Reference.
>
> Current v1.3.1: Reference source Depth QC is enforced at query planning with
> Nextflow --min_reference_depth (default 5x). It is separate from --depth_threshold.
>
> Current v1.3.0: rebuild with --mate-audit to exclude unidentifiable reference mates.
> Nextflow --insufficient_reference_policy unscored|replace selects capacity behavior; this run uses replace.
> Details below supersede the historical no-redraw source policy for replace mode only.

# CoverageMap_RE Raw Query vNext 1.2.1

Development baseline: `db11a9f21700178b8571566970480a1f8d7dc3e4`.
Feature branch: `codex/raw-query-vnext-integration-20261118`.
Only CoverageMap_RE receives new development; CoverageMap remains historical.
This is a local engineering handoff. Native Nextflow/NHB/data acceptance is pending.

## Integrated behavior

The FNA/FAI SQL catalog selects strict FASTA-backed validation. FASTA and FAI
contig names/lengths must agree. Alignment @SQ may be a subset of FASTA; present
contigs must match FASTA SN/LN/M5, omitted contigs contribute zero coverage
with the full FASTA/FAI denominator and are recorded in the missing-contig audit.
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


## Failure-audit retention repair — 2026-11-19

RAW_QUERY_SCORE passes --failure-audit-dir pointing to the persistent output failure_audits/query_id
directory configured from launchDir and outdir. On ERROR, score copies only summary.json/tsv, draws.tsv,
decode_audit.json and statistics.json before returning nonzero. Dense arrays and
raw alignments are not copied. Successful output/scoring/cache semantics remain
unchanged. This prevents scratch cleanup from removing the reason for failed
reference slots. Local focused97tests PASS; NHB native failure-path validation
pending. This does not resolve or alter the pending mate-label sampling decision.


## v1.3.0 labelled catalog and source replacement

The complete347993-row primary-mate inventory preserves rowID,Sample,GCF,path,
observed labels, scan completeness, observed primary count and errors. Presence
of both labels can stop scanning early; negative claims require a complete scan.
Counts from that inventory are observed lower bounds, not full distributions.
Original CSV/FNA/FAI and CRAMs remain unchanged.

`build_cram_index.py --mate-audit mate_audit.tsv` (Nextflow --mate_audit)
requires complete original-row coverage and exact Sample/GCF/path matches.
References with no identifiable mate labels or observed invalid mate flags are
excluded from reference_records; identifiable one-sided sources remain labelled,
with unsupported-arm capacity0 handled by source replacement. Original input
rowIDs remain stable despite exclusions. reference_mate_audit retains provenance
and included flags separately. A read error is not proof of missing mate labels
and prevents final filtering until resolved. Standalone excluded_references.tsv
retains all original CSV fields plus FNA/FAI paths and full audit evidence.

`--insufficient_reference_policy unscored` retains the existing source draws and
reports out-of-tolerance results. `replace` makes all SQL-eligible sources tracked
inputs, validates actual eligible bases separately for each Query arm, and replaces
only sources that cannot meet target within the unchanged tolerance. Complete
units are never truncated or upsampled. Each replacement uses independently keyed
source RNGs; read-unit draw RNGs and the1000-slot denominator remain unchanged.
Sources invalid for one arm remain usable for the other. Actual distinct-sample
support must still exceed min_samples. Decode/identity/M5 errors are not silently
converted into low capacity. Exhausted support is explicit NOT_EVALUATED.

replacement_audit.json/tsv retain each rejection, replacement source ID, reason,
target, actual available mass, tolerance and RNG token. Draw tables retain the
original reference ID and replacement attempt count. Geometry capacity checks
reuse the policy-keyed shared cache, so some archives are loaded for preflight
and again for scoring while raw decoding is cached. No speed claim is made.
All scientific statistics and optional methods/grid consume only the final frozen
slot arrays. No threshold, tolerance, clinical label or manifest row is altered.

## v1.3.1 Reference minimum-depth quality gate (2026-11-19)

- New Nextflow flag: `--min_reference_depth 5` (Python planner:
  `--min-reference-depth 5`). This is a **minimum reference-source
  manifest/CSV Depth**, inclusive, for **both** low and high Query branches.
  It is not a restriction on query depth or a re-estimate from CRAM alignment.
  `--depth_threshold 5` continues to choose the Query sampling branch only.
- Candidate rules: same GCF, exclude same sample/path and collapse duplicate
  source paths, positive source Depth, `ref.Depth >= min_reference_depth`
  for **all** queries; additionally `ref.Depth >= query.Depth` for a low
  query (`query.Depth <= depth_threshold`). High queries still compare full
  source coverage; no full-depth matching or upsampling is added.
- Exclusions happen **before** seeded source selection. The replacement pool
  is built from that same filtered set, so `replace` cannot reintroduce a
  below-floor source. `--min_reference_depth 0` restores the old selection
  floor (zero-depth records remain excluded, as previously).
- Source CSV Depth is **metadata only**. Actual eligible aligned bases and
  Read1/Read2 capacities continue to be checked under the existing
  150-bp tolerance and replacement policy. A declared >=5x source can still
  have insufficient actual usable capacity.
- Per-query `reference_depth_audit.json` records available distinct sources
  before depth filtering, distinct sample support, the number excluded below
  the reference floor, the number excluded by the low-Query depth check, and
  the final eligible pool. The same exclusion counts and configured floor
  appear in summary.json/summary.tsv, including unscored queries. Fatal SCORE
  failure audits retain the depth-audit JSON.
- The raw CRAM SQLite index, mate-label filter, random seed, per-slot
  sampling, null model, AND decision and scientific threshold methods are not
  modified. A change in the eligible source roster may change RNG-selected
  identities and resulting scores: this is an intentional source-QC
  policy change, not a numerical parity promise.
- For the NHB 344-query cohort, compare baseline vs filtered
  `draws.tsv` joined to SQLite `reference_records` by reference_id,
  count sub-5x selected slots separately by low/high arms and check the
  changed eligible-support/unscored rates before approval.

## v1.3.2 Configurable CRAM/BAM header M5 policy (2026-11-19)

Use Nextflow `--m5_policy strict|geometry_only`, passed to Python as
`raw_cram_query.py plan --m5-policy`. The default is **strict**.

- `strict`: each observed contig's alignment-header `@SQ M5` must match
  the digest computed from the exact SQL-cataloged FNA sequence. A missing
  M5 is rejected; missing entire contigs remain allowed as before.
- `geometry_only`: keep the authoritative FNA/FAI contig ID, length,
  ordering, offsets and computed M5 for all coverage geometry. Allow a
  differing or missing alignment-header M5 **without replacing the FNA
  checksum or modifying the alignment**, and record the discrepancy.
- A SQL FNA/FAI catalog is mandatory for `geometry_only`. Legacy CRAM-only
  indexes without a full FASTA are rejected for this mode.
- SN duplicates/unknown contigs and LN disagreements are rejected under
  **both** policies. A contig absent from alignment @SQ contributes zero
  coverage under the full FNA denominator. `skip` is intentionally unsupported.
- If the underlying CRAM requires a source reference for decoding and
  HTSlib cannot decode it, `geometry_only` still fails; no reference
  reconstruction, sequence bypass, or M5 rewriting is attempted.

Each scored Query writes `m5_audit.json`. It includes Query and each
decoded Reference source identifier, path, policy, cache hit, per-source
mismatch/missing-M5 counts, and per-contig `contig_id`, `length`,
`fna_m5`, `cram_m5`, `status` (`MATCH`, `MISMATCH`,
`MISSING_M5`, `MISSING_CONTIG`). The source status is `PASS`,
`MISMATCH_ALLOWED`, `NOT_VALIDATED` (legacy FAI-only), or `ERROR`.
Strict failures include per-contig M5 evidence in the persisted failure
audit. The cache serializes the audit and includes `m5_policy` in its
key, so policy changes cannot silently reuse a differently validated
geometry cache. The numeric scoring contract and reference replacement
rules are not changed.

**Interpretation limitation:** SN/LN identity does not prove that two
different FASTA sequences are coordinate-equivalent. Use
`geometry_only` only when provenance establishes that the CRAM
coordinates are comparable; follow up with the upstream generation team
for any unexpected digest disagreement. Do not interpret a permitted M5
mismatch as evidence that the underlying sequences are identical.

**Codex NHB validation:** test Query and Reference mismatches separately,
strict error audit and failure persistence, missing M5, missing/unknown
contigs, SN/LN mismatch rejection, `replace` path, Read/Fragment,
cold/warm geometry cache and policy-specific cache keys; then test
reproducibility and output comparisons in native Nextflow with both
`--m5_policy` modes.
