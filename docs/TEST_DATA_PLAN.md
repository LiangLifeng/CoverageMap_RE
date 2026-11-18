# Test data
Mac: tiny synthetic arrays and synthetic alignment objects only. Use independent per-position summation rather than reproducing direct-window algorithm.
NH: current 03.CoverageMap inputs, inspected before any processing. Run from a new external run directory and preserve inputs.
For performance comparisons, run direct and dense algorithms on identical pools, unit IDs and bin edges. Real current genome results cannot prove large-bacterial-genome scale performance; use declared synthetic stress data for that separate case.


## GCF reference depth database v1

The user-approved directory hierarchy, explicit-depth CLI, local/Slurm modes,
source identity, failure/availability semantics and acceptance requirements are
fixed in [REF_DEPTH_DATABASE.md](REF_DEPTH_DATABASE.md). Real build deferred.


## SQLite reference index / batched query extension

User-approved interfaces, immutable index, one-extraction query concurrency, unchanged storage/kernel semantics and acceptance are specified in [SQLITE_BATCH_CONCORDANCE.md](SQLITE_BATCH_CONCORDANCE.md).


## Current builder CLI cleanup (2026-10-19)

ref_depth_database.py accepts original-only construction with --window-size (1000bp) and one --jobs concurrency budget. Historical depth-grid/simulation/fixed-bin options and hidden concurrency alias are removed, not deprecated aliases. New native configs/results omit unused simulation parameters. [Current builder contract](REF_DEPTH_DATABASE.md) supersedes earlier depth-builder commands; coverage_concordance query sampling controls remain separate.

## Automatic whole-GCF batches — 2026-11-11

The user provides the full CSV; preprocessing groups intact GCFs into deterministic
byte-balanced batches. Metadata checks stay global, full source SHA256 moves into
each independently published batch, and the final registry merges every complete
batch without another full payload copy. Existing scientific semantics, strict
validation, and production jobs remain protected. Scope, CLI, output-layout change
and acceptance: [contract](REFERENCE_GCF_BATCHES.md).
