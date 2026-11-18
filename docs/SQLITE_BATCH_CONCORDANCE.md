# SQLite reference index and batched concordance — 2026-10-18

## Query-depth band and multi-layer simulations (2026-10-23)

The default native SQLite query now follows [DEPTH_BAND_CONCORDANCE.md](DEPTH_BAND_CONCORDANCE.md): full-FAI mean depth, 1000bp per-contig windows, inclusive +/-1X high-depth reference pool and Tools high-coverage maximum; queries below5X use multi-depth/N-repeat simulations with all/count/percent source selection and optional source-depth restrictions. This supersedes nearest-per-reference defaults for native queries; explicit `--comparison-policy nearest` retains the previous workflow. Read the contract for input/output, reference pair-distribution definition, non-estimable cases, integrity and synthetic acceptance.

User authorizes building a SQLite metadata index from the completed depth database,
batch reference comparison with one extraction per query sample–GCF CRAM, concurrent
query samples, current-library real tests on NH, then delivery on the same feature
branch. Existing geometry/depth files and scientific kernels remain unchanged.
No reference database rebuild, package installation, Slurm or GitHub/main action.

CLI: reference_index.py build --database DIR --out FILE [--verify-payloads]; inspect
--index FILE. New coverage_concordance.py --reference-sqlite FILE with either
--test-alignment CRAM --query-sample ID --gcf GCF or --query-manifest CSV
(Sample,CramFilePath,GCFID), explicit --depth/--depths/--depth-grid, new --outdir,
--query-jobs (default1), --reference-cache-mb (default256 per process),
--reference-batch-size (default32), optional --reference-samples comma list,
--include-self. Existing source/input modes remain compatible. New mode is CPU-only,
uses one native thread per process, --workers=1; query-jobs controls processes.

Index: finalized DB only, normalized reference samples and exact float.hex depth
keys, composite/PK lookups, relative paths rooted at recorded database location
(optional query --database-root relocation). Pin sample manifest/checksum registry
and selected matrix/edges/reference metadata hashes; include both supported and
unsupported depth states. Default index construction verifies metadata and payload
existence/stat, --verify-payloads verifies complete successful sample payloads.
Atomic no-overwrite publication, deterministic records, SQLite integrity check;
query read-only connections. No raw reference CRAMs/FAI or geometry decoded by new
cached-matrix path. Index binds to immutable database plan/version; raw reference
sources retained as provenance hashes, not query dependencies.

Single extraction per query CRAM into common full-FAI coordinates; validate per-GCF
geometry and inherited policies. Known M5 conflicts still reject; missing M5 recorded.
Each depth/resolved sampler test draw reuses task_rng(seed,depth,role1,rep0), exactly
matching old cached-reference CLI. Default exclude same sample or same source hash,
explicit include-self labels self consistency. Unsupported/missing ref/query depth
and no-reference states reported, not interpolated; expected no-match cases skip as defined below.
Per-reference errors remain explicit and fail run while other queries continue.
No new clinical decisions/calibration. Same evaluate_windows kernel/statistics.

Bounded process concurrency, at most2*jobs pending tasks, read-only index per worker,
per-process bounded LRU of verified matrices+metadata with stat-based invalidation;
checksum verified on load, files expected immutable throughout run. Cache accounts
array+JSON bytes (not total RSS), one current uncached matrix may exceed budget;
no global reference list/matrix materialization. Query geometry/draws local to a
query. Matrix cache can serve later queries in same process; no claim of cross-
process decoded-array sharing. Results stream to deterministic per-query and global
TSVs; per-comparison JSON, shared per-query test arrays, status/provenance/resources.
Index hash/stat verified before/after parent run. Keep old files unmodified.

Files: bin/reference_index.py, bin/concordance_batch.py, CLI routing changes,
tests/test_reference_index.py, usage/spec docs. Synthetic acceptance includes SQLite
join/negative cases/query plan, read-only/index collision, metadata/payload corruption,
relocation/no raw-reference dependencies, exact old/new statistics/IDs, extraction
call counts, bounded cache hits/eviction, serial/concurrent equality, failed-query
continuation/self-exclusion/depth gaps. Real NH acceptance: existing s24-B031/GCF
000913595.1 against5other references at0.1/1x, plus second-query serial/process
comparison; baseline old cached CLI statistics equal, record wall/CPU/RSS and
metadata/extraction/cache counts without unmeasured speed claims. Local synthetic
then committed candidate to NH for exact-SHA full tests and real-library acceptance.

## GCF affinity and exclusive stage timing (2026-10-18)

Default --query-scheduling gcf runs all queries of a GCF sequentially in one task/worker, with at most query-jobs GCFs active. A single GCF uses the calling process without spawn; --query-scheduling query retains prior independent query scheduling for comparisons. GCF tasks are ordered by descending query count, ties by GCF; summary order remains deterministic. Cache stays bounded per process, no shared decoded arrays. Cache eviction can still prevent reuse when a GCF working set exceeds the budget.

Root stage_timing.tsv summarizes all query stages. Query provenance adds worker_pid, stage_timing (exclusive wall/user/system/calls; nesting charged only to innermost stage), and scheduling metadata. Timers distinguish CRAM SHA, extraction, SQL lookup, reference metadata/checksum validation, payload SHA, NPZ/NPY load/decompression, reference compatibility, test sampling/counting, statistics, and output writes. Remaining setup/loop time is other. Loading/decompression includes filesystem read wait; timers do not prove disk versus codec causality. Output is not compressed in the query path. Snapshot excludes its own final report/status publication and CLI startup/imports; parent timing scope remains explicit. Failed query reports also retain timing.

New reference builds add checksummed timing.json with the same timer schema, separating CRAM hashes/extraction, geometry format selection/write/readback, sampling/counting, ID audit writes, matrix compression/write/readback, metadata and final payload checksums. No reference payload encoding or RNG changes. Geometry phase remains aggregate unless a future codec profile is implemented. Existing libraries/indexes remain readable; user reruns real workloads. Acceptance: unchanged IDs/statistics and failures; same GCF worker/cache reuse; multi-GCF bounded scheduling; finite nonnegative timing and nested attribution.

Stage CPU is process-wide over each wall interval; concurrent simulation threads can execute while the main thread writes results. Exclusive wall partitioning avoids double counting but cannot attribute concurrent CPU instructions to a particular function. Use function/thread profiling for that question.

## Queries without matching references (2026-10-18)

A new query sample need not occur in the reference database. Match by GCF. Missing GCF is an expected SKIPPED outcome (reason NO_REFERENCE_GCF) with no CRAM extraction. Known GCF with all references excluded by self/source/sample filters is SKIPPED/NO_ELIGIBLE_REFERENCE; with eligible references but no supported exact depth is SKIPPED/NO_COMPARABLE_DEPTH. Emit explicit summary rows per requested depth when there is no reference row, with empty similarities and a diagnostic reason. Existing per-reference depth-gap rows remain.

Expected skips use status/SKIPPED and exit_code0, not SUCCESS; processing errors/corrupt references remain FAILED/exit1. Any skips/depth gaps make root COMPLETE_WITH_GAPS; even all-skipped returns0 for successful processing, with ready_query_count0 and explicit skipped_query_count/outcomes. Never imply a skipped query passed similarity. Root summary/provenance must retain all queries and reasons. Existing index/payload formats and RNG/statistics stay unchanged. Tests cover mixed present/absent GCFs, all absent, only-self, filter removal, unsupported depths, serial and process scheduling plus true failure continuation. No real rerun required before delivery.
