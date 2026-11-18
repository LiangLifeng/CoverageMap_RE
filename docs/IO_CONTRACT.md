# I/O contract

## Capture scratch and reference-free field decoding

[Reference capture I/O](REFERENCE_CAPTURE_IO.md) defines the geometry-only CRAM
field mask, execution-only scratch environment settings, reserve checks and
temporary cleanup guarantees. Geometry, checkpoint and final-output contracts
remain unchanged.

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

CramPath contains absolute BAM or CRAM paths, including geometry-only BAM with no sequence/quality or BAI. The sampler accepts numeric depths with optional x/×. Resume binds code, index, manifest, FAI, source file stamps and sampling parameters; existing unknown outputs are rejected. See [the benchmark contract](BAM_SAMPLING_BENCHMARK.md).

## GCF virtual windows

The [GCF window contract](GCF_CONCAT_WINDOWS.md) supersedes the historical per-contig geometry below for new native libraries and Nextflow v2.2. It specifies the version, layout, payload, migration and acceptance checks.


## Multi-CPU reference index build — 2026-11-03

`reference_index.py build --jobs N` adds bounded spawn-process read/validation with ordered single-writer SQLite, default1. Existing strict/skip/audit/config/source-change and atomic read-only publication guarantees remain. Serial/parallel native, legacy, multi-root, corruption and CLI acceptance plus queue bounds are required. See [index contract](MULTI_DATABASE_INDEX.md#multi-cpu-build).

## Same-sample query correction (2026-10-24)

Native depth-band queries include same-sample IDs and identical CRAM content by default; `--exclude-self` is opt-in. An empty +/-1X band falls back to same-GCF nearest available depth (all exact ties), labelled outside-band with descriptive scores and null p/acceptance. A one-reference pool still emits similarity; only an absent indexed GCF is skipped. Existing filtered/unsupported cases are reported explicitly as COMPLETE_WITH_GAPS and corrupt inputs still fail. Full contract: [DEPTH_BAND_CONCORDANCE.md](DEPTH_BAND_CONCORDANCE.md). This supersedes the earlier native self-exclusion/no-fallback defaults; legacy/nearest replay defaults are unchanged.

## Query-depth band and multi-layer simulations (2026-10-23)

The default native SQLite query now follows [DEPTH_BAND_CONCORDANCE.md](DEPTH_BAND_CONCORDANCE.md): full-FAI mean depth, 1000bp per-contig windows, inclusive +/-1X high-depth reference pool and Tools high-coverage maximum; queries below5X use multi-depth/N-repeat simulations with all/count/percent source selection and optional source-depth restrictions. This supersedes nearest-per-reference defaults for native queries; explicit `--comparison-policy nearest` retains the previous workflow. Read the contract for input/output, reference pair-distribution definition, non-estimable cases, integrity and synthetic acceptance.
Coordinates: 0-based, half-open. Coverage includes CIGAR M,=,X, excludes D,N,I,S,H,P; gaps remain gaps.
Read sampling: each eligible primary mapped alignment is one unit; paired mates remain eligible under default both policy. read1/read2-only are optional.
Fragment sampling: group eligible primary alignments by query_name plus read group, select group together; coverage sums both aligned ends (overlap counts twice as read-base coverage). It is not insert-size span coverage or overlap deduplication.
Optional read filters: min MAPQ (default 0), exclude duplicate/QC-fail (defaults false); default unmapped/secondary/supplementary excluded.
Fixed count: uniform without replacement over units. Fraction: floor(f*N) fixed count, clearly distinct from Bernoulli thinning. Seed defines reproducible independent draws.
Target mean depth: choose complete units in uniform random order until target aligned bases reached; last unit is retained whole, overshoot reported, demand above parent total fails. Variable-length stopping is not fixed count sampling; it is explicitly an alternate stopping model.
Do not truncate units or recycle the pool. Empty pool/zero count produces zero descriptive coverage; constant correlation remains non-estimable.
Window edges: integer partition of [0,G] into exactly B nonempty windows, require 1<=B<=G. Retain actual widths. sum is aligned-base mass, mean=sum/width.
Unmatched reference names/lengths, malformed positions/pointers, invalid count/depth/fraction fail early.
CRAM may be streamed until EOF without CRAI. A missing matching external FASTA is an explicit input gap; embedded reference behavior must be checked on NH.
Explicit FAI mode: use --reference-index with CRAM required_fields=4159 (QNAME/FLAG/RNAME/POS/MAPQ/CIGAR/RGAUX) and decode_md=0. No SEQ/QUAL decoding or reconstruction is requested. Validate the exact contig set and lengths against the CRAM header and retain its M5 identities. FAI is layout metadata, not a substitute DNA sequence. This mode has been demonstrated on the supplied NH CRAM; full sequence analysis remains unsupported without matching FASTA.
Compatibility: extract_read_intervals/get_cached_reads_pool now return BlockPool rather than an (N,2) span array, because multi-block units cannot be represented faithfully by one interval. SubsamplingState still accepts legacy integer (N,2) arrays. Positive target_bases keeps the last unit whole (may overshoot) and requests above the original mass fail. Old unverified .npy caches are not reused. Constant Pearson is NaN; concordance with non-estimable evidence raises ValueError. The normal-approximation concordance model remains uncalibrated.
Cache identity includes source content hash, FASTA identity, layout/order and filtering/grouping version. Legacy unproven caches must not silently be reused.
No claiming concordance p-values calibrated; statistical model change outside scope.

## Explicit performance options

Defaults remain `--depth-sampler permutation --output-format tsv`, preserving historical selected IDs for the same seed and workload. `--depth-sampler lazy` uses sparse Fisher-Yates permutation prefixes with batches of at most 128 uniform ranks. It keeps whole units, uniform order without replacement and the same target-base stopping rule. Unused ranks in the final batch are discarded. It has a distinct seed/ID policy (`lazy-fisher-yates-batch128-v1`), recorded in provenance; do not require its IDs to equal the historical permutation policy. Count/fraction modes retain their existing selection regardless of this depth-only flag.

`--output-format npy` writes little-endian int64 `edges.npy` (B+1), `full.npy` (B), and disk-backed `simulations.npy` (R,B). Counts are aligned-base masses; mean depth is counts / diff(edges), computed when read. Row index is replicate index and column index is window index. No N×G coverage matrix is retained. The R×B output is a memory-mapped disk file; its mapped pages may contribute to RSS. Binary arrays and TSV with the same sampler/seed have exactly equal counts, coordinates and selected-ID hashes. Pool files and provenance are retained as before.

Schema `coverage-map-2` records sampler identity, output format and stage clocks for extraction, selection, full/simulation coverage, input hashing, output writes/flush, optional pool compression and remaining main-function work. Imports and final provenance serialization are included in external process wall but excluded from these main-stage totals. Failed partial output is not an accepted run; output-directory reuse remains prohibited.

Auto sampling: explicit `--depth-sampler auto` calibrates depth stopping once on the loaded pool, using three paired trials per sampler with independent fixed probe seeds and alternating order. It selects lazy only when its median sampling time is more than 10% below permutation; near ties choose permutation. Probe time is recorded separately and included in process wall. Probes do not advance the simulation RNG, and all replicates use the resolved policy. Zero targets skip calibration; count/fraction modes keep their unchanged selection and do not calibrate. Defaults stay permutation. Record requested/resolved policy, probe timings/seeds and version in provenance. Auto decisions can vary with host load; exact seed replay uses the recorded explicit policy. Validate one-time resolution, RNG isolation, and exact auto/explicit resolved IDs and counts on local fixtures and NH input. Auto minimizes measured sampler time, not a guarantee of globally optimal workflow time.


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


## User-approved FAI subset headers — 2026-10-18

FAI mode accepts alignment headers that are strict subsets of the full FAI.
Full FAI order, lengths, offsets and genome length remain authoritative; absent
contigs contribute zero coverage and remain in all window/depth denominators.
Record missing_header_contigs and their count/bases, separate from missing M5.
Reject extra header names, duplicate header names and shared-name length conflicts.
Present known M5 conflicts remain errors. Across GCF samples and reference/test
inputs, permit unknown M5 and merge known identities for subsequent consistency
checks; never invent identities for absent contigs. Retain incomplete-identity
provenance. FASTA mode keeps exact contig-set validation. Persist absence metadata
through reference bundles and extraction bridges. No changes to read filters,
complete-unit sampling, RNG or statistical kernels. Validate synthetic subset CRAM,
full-coordinate offsets, zero windows, full-genome depth, multiple samples missing
different contigs and downstream cached/reference input compatibility; real rebuild
is performed by the user after sync. Source/config fingerprint changes require a
new database output plan; preserve previous results.


## SQLite reference index / batched query extension

User-approved interfaces, immutable index, one-extraction query concurrency, unchanged storage/kernel semantics and acceptance are specified in [SQLITE_BATCH_CONCORDANCE.md](SQLITE_BATCH_CONCORDANCE.md).


## Native reference v2 supersession (2026-10-19)

Default reference builder now preserves originals only with 1000-bp per-contig windows. No mandatory depth grid/replicate persistence; historical builder parameters and aliases are removed. See [v2 contract](NATIVE_REFERENCE_V2.md) for indexed nearest actual-depth selection, on-demand matched sampling, non-estimable statistics, file integrity and acceptance.

## Progress output

Additional timestamped stderr/build.log lines and worker.log lifecycle records; existing JSON/status schemas retained with an additive status.progress aggregate. Processed sample–GCF percentage includes successes and known failures, with success/failure/reuse separate. No read percentage/ETA. --no-collect disables terminal watch catalog writes and rejects --refill. Existing frozen plans are not rewritten. See BUILD_PROGRESS.md.

## Production geometry policy (2026-10-21)

[GEOMETRY_FORMAT.md](GEOMETRY_FORMAT.md) fixes production defaults to packed/zlib9 and one verified candidate; explicit auto remains diagnostic. All existing formats/levels remain readable and existing plans/payloads are preserved. Default/explicit packed round-trip, selected level, candidate count and independent dense/same-ID invariants must pass before deployment.

## Sample concurrency

Native local --jobs bounds sample tasks, built-in Slurm --jobs still bounds GCF jobs. --collect-jobs defaults4 on prepare/run/execute/collect/status; stored scientific payload formats are unchanged. Native summaries add collection_jobs/collection_elapsed_seconds; local execution adds scheduling/task/worker metadata and local_sample_exit_codes.json. Sample progress/logs live outside attempt payloads. Source fingerprints include native_execution.py; old frozen plans are not rewritten. See SAMPLE_EXECUTION.md.


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

New user-authorized opt-in TSV query workflow: see [QUERY_MANIFEST_V1.md](QUERY_MANIFEST_V1.md) for the controlling IO_CONTRACT additions, replacement SAMPLE draws, strict thresholds, per-contig Tools-style windows, explicit p-value limits, output/identity contracts and acceptance. Existing builder and old CLI behavior are preserved.


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


## Reference draw delivery correction (2026-11-14)

Continue the accepted native-resume branch. Keep reference sampling and statistical formulas unchanged. The workflow must log its source directory and requested draw count; score must reject missing/mismatched requested_draws and inconsistent COMPLETE valid_draws, write an ERROR receipt and return nonzero. Preserve valid COMPLETE_WITH_GAPS handling. Validate 1000 requested versus stale100, correct1000, missing count and Nextflow forwarding with synthetic inputs. Deploy only to the independent reference-draws checkout, preserving the old in-use checkout and real results. Provide exact corrected launch commands; no real cohort run or GitHub publication. Allowed changes: main.nf, bin/coverage_query_flow.py, focused tests and relevant query documentation.
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
