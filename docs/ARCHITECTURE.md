# Architecture

## Jaccard coverage breadth

[Jaccard breadth](JACCARD_BREADTH.md) defines shared query/reference union counts,
threshold0 regression, seven-threshold shared-draw scoring, zero-union audit,
output contracts, fixed Pearson/sampling semantics and paired TP/FP acceptance.

## Complete alignment geometry and reference-build workflow

[Geometry vNext](GEOMETRY_VNEXT.md) defines the Fragment/Read/Segment/Block schema,
runtime sampling capabilities, legacy limits, coverage policy, independent
`build_reference.nf` DAG, bounded record workers and exact receipt/SQLite acceptance.
Complete capture supersedes primary-only information extraction for new builds;
primary-only remains the default coverage calculation. NHB requests CPUs and
concurrency without explicit Slurm memory or time limits.

## BAM sampling benchmark — 2026-11-11

The existing content-aware sequential alignment extractor accepts BAM in CramPath. Sampling adds atomic validated BAM publication, per-mode checkpoints and an output lock. Nextflow stages the separate saved-score plotter after collection; rendering does not change query scoring. See [the benchmark contract](BAM_SAMPLING_BENCHMARK.md).

## GCF virtual windows

The [GCF window contract](GCF_CONCAT_WINDOWS.md) supersedes the historical per-contig geometry below for new native libraries and Nextflow v2.2. It specifies the version, layout, payload, migration and acceptance checks.


## Multi-CPU reference index build — 2026-11-03

`reference_index.py build --jobs N` adds bounded spawn-process read/validation with ordered single-writer SQLite, default1. Existing strict/skip/audit/config/source-change and atomic read-only publication guarantees remain. Serial/parallel native, legacy, multi-root, corruption and CLI acceptance plus queue bounds are required. See [index contract](MULTI_DATABASE_INDEX.md#multi-cpu-build).

## Same-sample query correction (2026-10-24)

Native depth-band queries include same-sample IDs and identical CRAM content by default; `--exclude-self` is opt-in. An empty +/-1X band falls back to same-GCF nearest available depth (all exact ties), labelled outside-band with descriptive scores and null p/acceptance. A one-reference pool still emits similarity; only an absent indexed GCF is skipped. Existing filtered/unsupported cases are reported explicitly as COMPLETE_WITH_GAPS and corrupt inputs still fail. Full contract: [DEPTH_BAND_CONCORDANCE.md](DEPTH_BAND_CONCORDANCE.md). This supersedes the earlier native self-exclusion/no-fallback defaults; legacy/nearest replay defaults are unchanged.

## Query-depth band and multi-layer simulations (2026-10-23)

The default native SQLite query now follows [DEPTH_BAND_CONCORDANCE.md](DEPTH_BAND_CONCORDANCE.md): full-FAI mean depth, 1000bp per-contig windows, inclusive +/-1X high-depth reference pool and Tools high-coverage maximum; queries below5X use multi-depth/N-repeat simulations with all/count/percent source selection and optional source-depth restrictions. This supersedes nearest-per-reference defaults for native queries; explicit `--comparison-policy nearest` retains the previous workflow. Read the contract for input/output, reference pair-distribution definition, non-estimable cases, integrity and synthetic acceptance.
bin/coverage_core.py: reference layout, validated ragged block pool, CRAM stream extraction, selection, independent dense baseline and vectorized direct window aggregation.
bin/coverage_map.py: argparse CLI and deterministic outputs.
Tools.py: original API compatibility; preserve simple ndarray interval use and add full-block pool support. Declare behavior changes and deprecations explicitly.
tests/: independent geometry and sampling tests, synthetic indexed/CRAM integration when pysam is available.
Reference concatenation follows supplied FASTA order. Contig boundaries and widths are retained. Whole-genome windows may span contig boundaries without inventing bridging read coverage.
Use compressed numeric block arrays plus per-unit pointers; no sequence copies and no n×G retained simulation arrays.

Optional concordance acceleration: bin/gpu_core.py holds read-block geometry on one CUDA device, accumulates batches of CPU-selected IDs via int64 boundary/difference kernels, and returns compact window counts. The thread executor shares one serialized GPU stream; process GPU coverage is rejected. Statistics independently selects blocked FP64 GPU computation. bin/pool_bridge.py explicitly exchanges numeric geometry from a compatible pysam Python when the GPU interpreter uses a different ABI. Defaults and simulation RNG remain on CPU.


## Separate sample reference storage

The authorized storage builder, formats, thread limits, downstream reader, source
identity and lossless acceptance contract are specified in [REFERENCE_STORAGE.md](REFERENCE_STORAGE.md).


## Concordance saved-reference input

`coverage_concordance.py --reference-bundle` consumes the verified per-sample
bundle directory or its storage.json alongside a new test alignment.
`bin/concordance_inputs.py` handles immutable input identity, policy inheritance,
layout derivation and M5 compatibility without changing sampling/statistics.
Backward compatibility, negative checks and CRAM/bundle equivalence acceptance
are specified in [CONCORDANCE.md](CONCORDANCE.md).


## GCF reference depth database v1

The user-approved directory hierarchy, explicit-depth CLI, local/Slurm modes,
source identity, failure/availability semantics and acceptance requirements are
fixed in [REF_DEPTH_DATABASE.md](REF_DEPTH_DATABASE.md). Real build deferred.


## SQLite reference index / batched query extension

User-approved interfaces, immutable index, one-extraction query concurrency, unchanged storage/kernel semantics and acceptance are specified in [SQLITE_BATCH_CONCORDANCE.md](SQLITE_BATCH_CONCORDANCE.md).


## Native reference v2 supersession (2026-10-19)

Default reference builder now preserves originals only with 1000-bp per-contig windows. No mandatory depth grid/replicate persistence; historical builder parameters and aliases are removed. See [v2 contract](NATIVE_REFERENCE_V2.md) for indexed nearest actual-depth selection, on-demand matched sampling, non-estimable statistics, file integrity and acceptance.

## Build progress

ref_depth_database.ProgressLog reads atomic per-GCF counters at bounded heartbeat intervals, emits flushed stderr lines and appends build.log for execution. Workers emit sample lifecycle events and refresh progress after each outcome. Status --no-collect is a read-only observer; no payload verification in the display polling path. See BUILD_PROGRESS.md.

## Production geometry policy (2026-10-21)

[GEOMETRY_FORMAT.md](GEOMETRY_FORMAT.md) fixes production defaults to packed/zlib9 and one verified candidate; explicit auto remains diagnostic. All existing formats/levels remain readable and existing plans/payloads are preserved. Default/explicit packed round-trip, selected level, candidate count and independent dense/same-ID invariants must pass before deployment.

## Native local pool and verification

bin/native_execution.py owns bounded spawned persistent sample workers and controller-only GCF progress. ref_depth_worker.execute_sample shares an exclusive per-sample attempt lock across worker modes. native_reference verifies bounded samples concurrently and publishes deterministic catalogs from one controller, with pre/post per-GCF FAI identity and common-layout/M5 validation. See SAMPLE_EXECUTION.md.


## Multi-database metadata index and depth selection (2026-10-21)

[MULTI_DATABASE_INDEX.md](MULTI_DATABASE_INDEX.md) defines schema3 multi-root payload resolution,
audited overlapping-copy collapse/conflict rejection, per-sample full-genome depth,
post-depth distinct-sample support, immutable metadata-only reselection and reports.
Original library payloads and old schema1/2 readers remain preserved.


## Explicit partial SQLite indexing (2026-10-22)

`reference_index.py build --skip-invalid-samples` may omit invalid source samples
with transactional rollback and a persistent error audit, while preserving valid
sample validation, strict defaults and fatal global conflicts/source changes.
Partial index provenance/counts and report semantics are specified in
[MULTI_DATABASE_INDEX.md](MULTI_DATABASE_INDEX.md#explicit-omission-of-invalid-samples).


## External-only skipped-sample audit (2026-10-22)

`build` and `filter --skipped-samples-tsv NEW_PATH` retain detailed omission
records in an exclusive external TSV, purge their table/deleted bytes from the new
SQLite, preserve count/provenance and keep later filtering external. Report validates
and copies that TSV. Defaults and valid-reference validation remain unchanged.
See [MULTI_DATABASE_INDEX.md](MULTI_DATABASE_INDEX.md#keep-omission-details-only-in-an-external-tsv).


## Read-only Sample–GCF export (2026-10-24)

`bin/reference_export.py` supports exact repeatable Sample/GCF IDs, inclusive depth bounds and existing selected-only flags. See [query/export contract](REFERENCE_EXPORT.md). Standard-library-only schema-3 metadata queries, deterministic TSV, read transaction, no alignment/payload access, and exclusive atomic output. Acceptance covers combined filters, zero matches, provenance, SQL parameters, invalid bounds, unknown schemas, source SHA preservation and publication failures.

## Required payload I/O recovery (2026-10-27)

Builders enumerate mandatory payloads from storage declarations and apply bounded transient filesystem retries. Strict integrity failures and old frozen plans remain protected. See [PAYLOAD_RETRY.md](PAYLOAD_RETRY.md) for the exact errno, delay, logging and acceptance contract.

For current native local builds, final verified FAILED sample–GCF pairs get at most two additional rounds with unchanged resource bounds, preserved healthy outputs, independent attempts and sample_retries.json audit. Persistent failures remain INCOMPLETE; interruptions are not retried. This does not automatically resubmit raw Slurm or legacy jobs.

## Query manifest v1 — 2026-10-25

New user-authorized opt-in TSV query workflow: see [QUERY_MANIFEST_V1.md](QUERY_MANIFEST_V1.md) for the controlling ARCHITECTURE additions, replacement SAMPLE draws, strict thresholds, per-contig Tools-style windows, explicit p-value limits, output/identity contracts and acceptance. Existing builder and old CLI behavior are preserved.


### High-depth null completed — user definition, 2026-10-25

The user selected different-Sample pairs among this query's drawn references, with repeated draw multiplicities and at most1000 unordered pairs. This is now the default high-null policy; exact pair identities/I/U and population/sample counts are retained, and normal-CDF p is computed where estimable. See QUERY_MANIFEST_V1.md for the controlling definition. Tests must cover replacement weights, the1000-pair cap, degenerate null status, independent p reconstruction and NH serial/process equality while retaining all prior results.


### GCF-global bin-count correction — 2026-10-25

The user's bin count is a divisor of total GCF reference length: shared W=ceil(sum(all contig lengths)/B), then apply W within each contig and retain actual short contigs/tails. It is NOT a target count independently applied per contig. --bin-count is the accurate CLI name; --windows-per-contig is a deprecated warning-emitting alias with corrected semantics. Fixed-bp windows remain unchanged. See QUERY_MANIFEST_V1.md and global-width/dense/CLI tests. Previous N-mode results are preserved but superseded; no database format change/rebuild.

## Manifest query unlimited memory budgets — 2026-10-29

User-authorized `query_manifest.py`: `--cache-mib 0` and `--work-mib 0` independently mean unlimited; defaults are131072 MiB (128 GiB) each, positive budgets divide across jobs, negative values fail. Derived unlimited bytes are JSON null; NativeCache uses None so other internal zero-cache callers stay bounded. Stamp verification, numerical results and resume identity checks remain required. See [contract](QUERY_MANIFEST_V1.md). Validate cache retention/invalidation, work guard boundary, CLI negative/default cases and serial/spawn numerical equality for bounded, mixed and unlimited budgets. Deploy after local/NH synthetic regression; do not rerun TP queries.

## Automatic whole-GCF batches — 2026-11-11

The user provides the full CSV; preprocessing groups intact GCFs into deterministic
byte-balanced batches. Metadata checks stay global, full source SHA256 moves into
each independently published batch, and the final registry merges every complete
batch without another full payload copy. Existing scientific semantics, strict
validation, and production jobs remain protected. Scope, CLI, output-layout change
and acceptance: [contract](REFERENCE_GCF_BATCHES.md).
## 2026-11-16 Raw CRAM SQL index scope

The independent build_cram_index.nf entry point imports BuildInput CSV metadata
into a GCF-indexed SQLite database without reading alignments or saving geometry.
See RAW_CRAM_SQL_INDEX.md for the schema, validation, output and resume contract.
This stage does not implement the future raw-CRAM similarity query backend.
Acceptance requires CSV roundtrip/index lookup and error protection tests plus
an actual synthetic Nextflow run and resume; formal scientific runs remain user-owned.

## 2026-11-16 Raw CRAM query supersession

The user authorizes a new raw_query/main.nf: one query per process, raw SQLite GCF
lookup, replacement source draws, single decode per unique selected reference,
independent low-depth read/pair simulations, high-depth vector/similarity reuse,
original v2 Normal Jaccard/Pearson AND decisions and native query-level resume.
See RAW_CRAM_QUERY.md for the controlling minimal contract and acceptance.
This supersedes the earlier proposed per-reference/GCF process and pre-scan designs.

## 2026-11-16 Raw query legacy feature parity

User requires all query features from the pinned17f8fd8 main.nf, waiving only
cache_mib/work_mib. RAW_QUERY_PARITY.md supersedes the earlier limited raw query
interface: original methods/bootstrap/sensitivity/coverage policies/diagnostics
and collected ledgers are restored atop SQL, one-query processes and decode reuse.
Acceptance includes original-function method/threshold/plot and split-read oracle
checks, option forwarding, real synthetic Nextflow diagnostics and native resume.



## CoverageMap_RE vNext integration — 2026-11-18

Continue migrated db11a9f on a local feature branch. Preserve the existing pooled
N-slot row-mean null, Normal/Pearson AND and optional methods. Connect mate-aware
whole-unit tolerance sampling and shared geometry through threshold outputs and
Diagnostics. Cache is disposable, policy/source/reference keyed, and must not
change frozen scientific arrays. Validate locally with synthetic BAM/CRAM, dense
selected-unit oracles, all methods/grid, cold/warm/corrupt cache, low/high and
read/fragment. Formal NHB Nextflow/resume/data tests remain user-owned.


## Reference mate audit and configurable replacement — 2026-11-19

User supplies BuildInput_all_virus_347993.csv plus c50_20261114 fna.csv/fai_path.csv
for a new independent SQL catalog and read-label inventory. Keep every CSV row;
never repair flags or silently filter source data. The mate-label CLI scans
mapped primary reads at MAPQ0; it may stop after proving both labels present.
Missing/one-sided/invalid/empty/error states require a complete scan and remain
explicit; observed counts are not full read-count estimates. Full347993-file
audit runs only in bounded NHB compute jobs; no full CRAM hashes or recopies.
User supersedes earlier no-redraw choice: insufficient-reference behavior becomes
external unscored|replace; this rerun selects replace, with per-slot replacement
audit, deterministic RNG, complete-unit sampling, unchanged150bp tolerance and
1000-slot denominator. Unlabelled-reference compatibility remains a separate
user choice. Deliver prior49-query/152capacity-slot detailed tables as well.

User follow-up locks final exclusion of unidentifiable mate references, with
standalone provenance. Optional --mate-audit catalog filtering preserves original
row IDs and rejects incomplete/error inventories. Capacity replacement is the
only new sampling-policy parameter; no fallback guessing of unlabelled mates.
