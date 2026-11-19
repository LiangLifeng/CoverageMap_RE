> CoverageMap_RE vNext 1.3.2: see [integration status and NHB test handoff](docs/RAW_QUERY_VNEXT_HANDOFF.md). Formal NHB acceptance remains pending.

# CoverageMap

## Automatic whole-GCF reference batches

Supply the full CSV to `build_reference.nf` as before. Preprocessing checks the
global roster and file metadata, preserves complete GCF groups, and balances
batches using total source bytes. Full source SHA256 and building run inside
each independently scheduled batch. `--chunk_size` is a target average record
count, so a large GCF can exceed it. `--plan_cpus` defaults to 64;
`--cpus_per_chunk` bounds batch workers and `-qs` bounds queued tasks.

Start this version in a fresh outdir. Inspect `plan/batches.tsv` in the PLAN
work directory, then retain `outdir/batches/*/database`: the one final
`outdir/reference.sqlite` points to these complete batch databases, without a
second global payload copy. `outdir/database` holds aggregation/input receipts,
not a standalone native database. Full checksum/mutation checks and verified
record checkpoint reuse remain required. See [the contract](docs/REFERENCE_GCF_BATCHES.md).

## Same-sample query correction (2026-10-24)

Native depth-band queries include same-sample IDs and identical CRAM content by default; `--exclude-self` is opt-in. An empty +/-1X band falls back to same-GCF nearest available depth (all exact ties), labelled outside-band with descriptive scores and null p/acceptance. A one-reference pool still emits similarity; only an absent indexed GCF is skipped. Existing filtered/unsupported cases are reported explicitly as COMPLETE_WITH_GAPS and corrupt inputs still fail. Full contract: [DEPTH_BAND_CONCORDANCE.md](docs/DEPTH_BAND_CONCORDANCE.md). This supersedes the earlier native self-exclusion/no-fallback defaults; legacy/nearest replay defaults are unchanged.

## Query-depth band and multi-layer simulations (2026-10-23)

The default native SQLite query now follows [DEPTH_BAND_CONCORDANCE.md](docs/DEPTH_BAND_CONCORDANCE.md): full-FAI mean depth, versioned whole-GCF windows on new indexes, inclusive +/-1X high-depth reference pool and Tools high-coverage maximum; queries below5X use multi-depth/N-repeat simulations with all/count/percent source selection and optional source-depth restrictions. This supersedes nearest-per-reference defaults for native queries; explicit `--comparison-policy nearest` retains the previous workflow. Read the contract for input/output, reference pair-distribution definition, non-estimable cases, integrity and synthetic acceptance.

**Current reference construction stores original coverage with 1 kb whole-GCF windows by default.** See [native reference build/query and files](docs/NATIVE_REFERENCE_V2.md) and editable [local](examples/local.sh) / [Slurm](examples/slurm.sh) wrappers. Builder simulation parameters and compatibility aliases have been removed; query sampling parameters remain available.
Exact complete-read or complete-fragment downsampling of real aligned coverage, with direct window aggregation.
The original Tools.py is preserved byte-identically in the initial main commit; current public names have documented compatibility changes.

## Inputs and interpretation
Use a matching FASTA (--reference) or explicit FAI layout mode (--reference-index).
For CRAM only required alignment fields are decoded; sequences/qualities and MD/NM reconstruction are omitted. FAI mode checks names and lengths against CRAM SQ and records M5; it cannot support DNA reconstruction.
Reads are sampled without replacement. Fragment mode groups query name plus read group and counts both aligned ends; overlap counts twice as read-base coverage. read1-only remains an explicit mate policy, not overlap deduplication.
Only CIGAR M/= /X blocks contribute. Original coverage preferences are retained.

## CLI
See bin/coverage_map.py --help. Required arguments: alignment, one reference mode, new outdir, one of reads/fraction/depth.
Optional: bins (1000), replicates (100), seed (0), mate-policy both/read1/read2, sampling-unit read/fragment, MAPQ/duplicate/QC-fail filters.
Output full.tsv and simulations.tsv retain exact coordinates, aligned-base sums and mean depth; provenance.json records targets, selected-unit fingerprints, source hashes and achieved depth. --save-pool writes a compact lossless compressed numeric pool.
--reads counts sampling units; --fraction uses floor(f*N), not Bernoulli thinning. --depth uses a random-order complete-unit stopping rule and reports overshoot.

## Staged optimization
1. Retain mate and CIGAR semantics, complete sampling units, independent seeded repeats.
2. Store block coordinates and unit pointers, without Python tuple pools or sequences.
3. Directly accumulate boundary contributions and full-window difference counts. Per-repeat coverage computation is O(selected blocks + windows); no per-base or replicates-by-genome matrix is retained.
The sparse selection ID check includes uniqueness validation. Depth-target selection explicitly requires a full random ordering of the pool; fixed low-count sampling avoids this path.
Use bin/validate_coverage.py on NH for identical-ID comparison with an independent dense baseline and a separate synthetic 5Mb scale case.

## Tests and limits
Run python -m pytest tests. Real synthetic CRAM integration needs pysam; local core tests do not.
Concordance functions retain the prior normal-approximation model with corrected degeneracy handling. Their p-values are not calibrated biological acceptance evidence.
See docs/IO_CONTRACT.md for deliberate API changes and docs/ACCEPTANCE.md for validation criteria.


## Reusable per-sample reference bundles

`bin/reference_store.py build` extracts one alignment without merging samples,
compares lossless disk formats, saves the smallest validated candidate and
records source identity, filtering, reference layout and storage measurements.
Compression/read/hash work uses bounded threads; extraction streams one sample.

```bash
python -B bin/reference_store.py build --sample-id SAMPLE_ID \
  --alignment /absolute/sample.cram --reference-index /absolute/reference.fna.fai \
  --workers 4 --format packed --outdir /absolute/new-reference-bundle
python -B bin/reference_store.py inspect --bundle /absolute/new-reference-bundle --workers 4
```

Downstream Python code can use `from bin.reference_store import load_reference_bundle`
and `pool, layout, manifest = load_reference_bundle(bundle_path, workers=4)`;
`select_units` and `window_sums` then operate on the returned pool unchanged.
The concordance CLI accepts the bundle directly. Keep original CRAMs for any
future policy changes. See [storage contract](docs/REFERENCE_STORAGE.md).

```bash
python -B bin/coverage_concordance.py \
  --reference-bundle /absolute/reference-bundle/storage.json \
  --test-alignment /absolute/new-sample.cram \
  --depths 0.01,1,10,100 --replicates 100 --bins 1000 \
  --workers 4 --executor thread --depth-sampler permutation \
  --outdir /absolute/new-results
```

Saved policies/layout are inherited and verified; do not copy storage.json away
from its bundle files. An additional FASTA/FAI is optional in this geometry-only
mode. See [concordance contract](docs/CONCORDANCE.md).


## GCF original reference database

`bin/ref_depth_database.py` builds original reference geometry and fixed-length coverage windows only. See [build contract](docs/REF_DEPTH_DATABASE.md), [native storage/query details](docs/NATIVE_REFERENCE_V2.md) and editable [local](examples/local.sh)/[Slurm](examples/slurm.sh) scripts. No preset depth, simulation-count, seed or fixed-bin parameters are accepted by the builder.

```bash
python -B /absolute/source/bin/ref_depth_database.py run --mode local \
  --records /absolute/records.csv --references /absolute/fai_path.csv \
  --outdir /absolute/new-database --window-size 1000 --cpus 4 --jobs 2
python -B /absolute/source/bin/reference_index.py build \
  --database /absolute/new-database --out /absolute/new-database.sqlite --verify-payloads
```

For Slurm use `run --mode slurm --jobs N`, default partition cpu, memory omitted unless explicitly set. `--dry-run` writes plans/scripts only, without worker/scheduler calls. Normal run prepares, executes/submits, monitors and collects. Script wrappers also build SQLite after successful collection; both default DRY_RUN=1. Recovery commands: execute(local), submit then status --watch --scheduler --refill(slurm), collect. One concurrency option: --jobs (default1). New source/config requires a new plan directory; never overwrite existing results.

Each GCF/sample attempt stores geometry/, edges.npy, native_coverage.npz, result.json, timing.json, checksums.json and status/. Original windows preserve contig/local boundaries, short-window lengths, aligned-base sums, mean depths and union breadth. Root catalogs are GCF_INDEX.tsv, SAMPLE_INDEX.tsv, NATIVE_INDEX.tsv, RESOURCE_INDEX.tsv and SUMMARY.json. Downstream native queries read the separately built SQLite. No depth matrices or sampled-ID archives are built. Full FAI including zero-covered missing-header contigs remains the coordinate/depth denominator; optional CRAI support and all required file checks remain.

Native query extraction occurs once; indexed nearest original-depth selection and matched sampling run at query time. Reference cache is bounded per worker, same-GCF queries reuse it. See native guide for explicit target depths, nearest-reference limit, descriptive-only full-depth scores and non-estimable statistics. Existing normal approximation is not calibrated; no speedup or large-cohort throughput claimed.

The SQLite section below documents existing depth-library query support; these query arguments do not re-enable historical builder parameters.

## SQLite reference index and concurrent query batches

Keep the current depth/geometry storage unchanged. Build a separate metadata index
from a completed database (`--verify-payloads` adds full product checksum checks;
default checks metadata and records pinned payload hashes):

```bash
python -B bin/reference_index.py build --database /absolute/depth-database --out /absolute/references.sqlite --verify-payloads
python -B bin/reference_index.py inspect --index /absolute/references.sqlite
```

Then compare one query against all compatible indexed references:

```bash
python -B bin/coverage_concordance.py --reference-sqlite /absolute/references.sqlite \
  --test-alignment /absolute/query.cram --query-sample SAMPLE --gcf GCF_000913595.1 \
  --depths 0.1,1 --seed 123456 --outdir /absolute/new-query-output
```

For multiple queries, CSV columns are `Sample,CramFilePath,GCFID`; relative CRAM
paths resolve against the CSV parent. Each sample–GCF must be unique:

```bash
python -B bin/coverage_concordance.py --reference-sqlite /absolute/references.sqlite \
  --query-manifest /absolute/queries.csv --depths 0.1,1 --seed 123456 \
  --query-jobs 4 --reference-cache-mb 256 --reference-batch-size 32 \
  --outdir /absolute/new-batch-output
```

SQLite mode extracts each query CRAM once, reuses each depth/sampler test draw,
and compares cached matrices without loading reference geometry or raw reference
CRAMs. Full reference layout/filter/known-M5 compatibility and exact-depth/matrix
integrity checks remain. Default excludes same-sample and identical-source-hash
references; `--include-self` explicitly enables self-consistency. Optional
`--reference-samples ID1,ID2` limits references. Stored bins/replicates/policies are
authoritative. Unsupported query/reference depth and absent exact cached depths
are explicit summary rows, with no interpolation; zero usable comparisons fail.

`--query-jobs` controls CPU query processes, one native numerical thread each.
Keep `--workers 1` in this mode. Cache is per process, not shared decoded RAM;
its budget accounts arrays/JSON, not full RSS. Metadata is fetched in bounded batches;
matrices are loaded one at a time and retained by a bounded LRU when they fit.
Versioned reference files must remain immutable; cache misses verify checksums,
hits check file stats/status, parent verifies index hash before/after. Use
`--database-root /relocated/database` if payload root moved; index stores relative
payload locations. Do not edit a published SQLite file in place. New mode is CPU-only;
legacy alignment/bundle/depth-manifest modes and GPU options remain available.

Outputs: root `summary.tsv`, `provenance.json`, `queries.json`, `status/`;
per-query `summary.tsv`, `provenance.json`, `progress.json`, shared `edges.npy`,
`draws/*/test.npy` and `test_ids.npy`, and per-reference/depth `concordance.json`.
A failed query does not prevent remaining queries; final status is nonzero for
errors/no usable reference. Per-query extraction count, cache hits/misses, wall/CPU
and worker cumulative peak RSS are recorded. Statistical model/calibration unchanged.
See [SQLITE_BATCH_CONCORDANCE.md](docs/SQLITE_BATCH_CONCORDANCE.md) for the contract.

### GCF affinity and resource stages

SQLite queries default to `--query-scheduling gcf`: all queries for the same GCF share one worker and its bounded reference cache. `--query-jobs` limits concurrently active GCF groups. A single GCF runs without spawning a child process. Use `--query-scheduling query` to retain independent query scheduling when same-GCF latency is more important than cache reuse. Neither policy guarantees faster execution; compare on the intended workload.

Read output `stage_timing.tsv` or per-query `provenance.json` → `stage_timing` for exclusive wall/user/system measurements, plus `worker_pid` and cache hits. New reference builds also write checksummed per-sample `attempt-*/timing.json`; existing libraries need no rebuilding for query timing or GCF scheduling. Stages preserve payload formats, seeds, IDs and statistics. CLI startup is outside query timers; CPU from concurrent worker threads can overlap callback/write measurement windows.

Queries match the reference database by GCF, not query sample name. An absent GCF is reported as `SKIPPED / NO_REFERENCE_GCF` without CRAM extraction; self/filter-only reference sets as `SKIPPED / NO_ELIGIBLE_REFERENCE`; no usable depths as `SKIPPED / NO_COMPARABLE_DEPTH`. Per-depth summary rows keep reasons and blank similarities. Expected skips do not abort other queries or make processing exit nonzero: root state is `COMPLETE_WITH_GAPS`, with explicit ready/skipped query counts (including all-skipped runs). A skipped comparison is not a similarity pass. Corrupt files and calculation errors still fail with nonzero exit.

### CRAM-only input and required payload integrity

New reference builds accept CRAM without CRAI. An adjacent `<full-cram-filename>.crai`, if present, is validated as readable/nonempty, pinned by stat/SHA and archived; an absent CRAI is recorded as null and is neither generated nor required. Optional recorded files still must remain unchanged while a plan executes. Existing frozen libraries/indexes remain readable.

All other required input and payload integrity checks remain: readable/nonempty inputs, mandatory geometry/chunks/arrays, coordinates, READY matrices/reference metadata, sampled-ID audits, new-build timing metadata, checksum registration and full SHA validation where enabled. Declared CRAI archives remain checked for libraries that contain them. Empty status marker files are intentional. Use a new output plan for new builds after code changes.

### Slurm memory defaults and dry run

Omit `--memory` (or set `--memory default`) to use scheduler-default memory; child scripts and submit commands omit all memory flags. Explicit `--memory 8G` remains supported. `--partition cpu` applies to child jobs; `--cpus 10` controls their CPU request, independently of an outer controller allocation. For omitted memory, child submission clears inherited `SBATCH_MEM_PER_NODE`, `SBATCH_MEM_PER_CPU`, `SBATCH_MEM_PER_GPU`.

Add `--dry-run` to `run --mode slurm` to validate inputs and generate plan/scripts plus `dry_run.json`, without running workers or any scheduler command. Existing plans support `execute --mode slurm --dry-run` and `submit --dry-run`. This is an offline preview, not scheduler feasibility validation. A preview directory can later use `submit --database DIR` followed by `status --database DIR --watch --scheduler`; a fresh `run` requires a new output directory. Existing immutable plans retain their recorded memory; generate a new plan to change resource settings.

Slurm child jobs now default to `-p cpu` (`--partition` overrides). `run --mode slurm --jobs 20` caps this database's tracked nonterminal jobs and each round's new submission attempts, automatically refilling slots until all GCFs finish. This is not an account-wide cap or a memory reservation. `dry_run.json` includes the limit and first-round bound; `submission_limits.json` records submitted/deferred counts. Failures do not retry automatically; uncertain submissions still require reconciliation. For a prepared/dry-run plan, use `submit --database DIR --jobs 20`, then `status --database DIR --watch --scheduler --refill --jobs 20` to continue explicitly; ordinary status/watch without `--refill` never submits. No explicit memory request remains the default.

Concurrency now uses one public option, `--jobs N`, consistently across run/execute/submit/status. In local mode it limits simultaneous GCF processes; in Slurm it limits tracked queued/running GCF jobs and each round of submission attempts. Default1 in either mode. The historical concurrency alias is no longer accepted. For Slurm previews use `run --mode slurm --jobs 5 --dry-run`; for explicit refill use `status --watch --scheduler --refill --jobs 5`.

## Build progress

Batch builds now emit flushed timestamped progress and append `<database>/build.log`; per-GCF sample logs remain in `worker.log`. Observe an existing batch read-only with `status --database /absolute/database --watch --interval 10 --no-collect`. See [progress/logging semantics](docs/BUILD_PROGRESS.md), including active-batch source protection and sample–GCF counter definitions.

Production geometry defaults to fixed packed/zlib9 after the full IH format survey; see [format policy and evidence](docs/GEOMETRY_FORMAT.md). This evaluates one validated candidate per sample. Explicit auto remains a diagnostic benchmark only.

## Native sample scheduling

Native `run --mode local` (including an outer Slurm allocation) now uses persistent sample workers across GCFs, bounded by `--jobs`; `--cpus` remains the per-worker geometry-thread budget. `--collect-jobs 4` independently bounds final integrity-validation threads and in-flight records. Sample logs: `<database>/<GCF>/samples/<sample-folder>/worker.log`; aggregate build.log/progress remains. Built-in Slurm keeps GCF jobs. See [execution/integrity contract](docs/SAMPLE_EXECUTION.md). New code requires new plans; existing libraries are not rebuilt or rewritten.


Multi-database SQLite and configurable sample-depth/GCF-support filtering: see
[the schema3 guide](docs/MULTI_DATABASE_INDEX.md) and examples/build_reference_sqlite.sh.


To index valid samples from a partially failed library explicitly, add
`reference_index.py build --skip-invalid-samples`. Omitted samples are recorded
in SQLite and `report` exports skipped_samples.tsv; original libraries stay unchanged.

Use `--skipped-samples-tsv /absolute/new/skipped_samples.tsv` with build/filter
to keep invalid-sample details only outside the final SQLite. Later filters inherit
the external audit; report verifies and exports it.


## Read-only Sample–GCF export (2026-10-24)

`bin/reference_export.py` supports exact repeatable Sample/GCF IDs, inclusive depth bounds and existing selected-only flags. See [query/export contract](docs/REFERENCE_EXPORT.md). Standard-library-only schema-3 metadata queries, deterministic TSV, read transaction, no alignment/payload access, and exclusive atomic output. Acceptance covers combined filters, zero matches, provenance, SQL parameters, invalid bounds, unknown schemas, source SHA preservation and publication failures.

Builder checksum reliability and bounded transient I/O recovery: [docs/PAYLOAD_RETRY.md](docs/PAYLOAD_RETRY.md).

Multi-sample TSV queries with replacement: [query_manifest.py contract](docs/QUERY_MANIFEST_V1.md). This separate CLI implements the 2026-10-25 user flowcharts. High p uses the user-defined distinct-Sample draw-pair distribution by default, preserving repeated draw weights and retaining at most1000 pairs. Degenerate reference distributions remain explicitly non-estimable.

Nextflow v2.2 and the native builder use natural contig order and virtual GCF coordinates. `--bin_count B` in Nextflow produces exactly B balanced windows. `--window_size S` produces fixed-bp windows with one possible short GCF tail. Use a rebuilt or reprojected native index and re-score older results. See the [window and migration contract](docs/GCF_CONCAT_WINDOWS.md).

Manifest query memory budgets: `query_manifest.py --cache-mib 0 --work-mib 0` disables the cache retention and numeric matrix budgets independently. Defaults are 131072 MiB (128 GiB) each; positive values divide across jobs, negatives fail. These are application budgets, not OS/RSS limits. See `docs/QUERY_MANIFEST_V1.md`.
# Query input / reference file-access diagnosis

Run inside the original query execution environment:

```bash
python -B bin/audit_query_files.py --query-run /absolute/query/output --outdir /absolute/new/audit/output
```
The standalone standard-library tool traverses every query input and every indexed
native reference dependency, including unselected records, and records filesystem
errors separately from historical draw failures. See [scope, outputs and exit
codes](docs/QUERY_FILE_AUDIT.md). This is a diagnostic audit; it does not rerun queries.

## Stratified Bin Count evaluation

See [CLI and scientific contract](docs/BIN_COUNT_EVALUATION.md) for frozen preflight + SQLite inputs, geometry-only reconstruction and independent distance/ordination choices.

Pearson is available as `--distance correlation` (distance = 1 - r).
Use `evaluate --matrices /absolute/existing/matrices --outdir /absolute/new/output`
to add this distance without rebuilding coverage. See the [Pearson contract](docs/BIN_COUNT_PEARSON.md)
for undefined-profile diagnostics and saved-matrix validation.
