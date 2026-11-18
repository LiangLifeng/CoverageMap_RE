# MODULE_PLAN

Orchestrator ref_depth_database.py; catalog/manifest/scheduler module ref_depth_catalog.py; per-sample worker/storage reader ref_depth_worker.py. Reuse existing extraction, complete-unit RNG and bounded simulate_grid. See REF_DEPTH_DATABASE.md.


## Current builder CLI cleanup (2026-10-19)

ref_depth_database.py accepts original-only construction with --window-size (1000bp) and one --jobs concurrency budget. Historical depth-grid/simulation/fixed-bin options and hidden concurrency alias are removed, not deprecated aliases. New native configs/results omit unused simulation parameters. [Current builder contract](REF_DEPTH_DATABASE.md) supersedes earlier depth-builder commands; coverage_concordance query sampling controls remain separate.

## Automatic whole-GCF batches — 2026-11-11

The user provides the full CSV; preprocessing groups intact GCFs into deterministic
byte-balanced batches. Metadata checks stay global, full source SHA256 moves into
each independently published batch, and the final registry merges every complete
batch without another full payload copy. Existing scientific semantics, strict
validation, and production jobs remain protected. Scope, CLI, output-layout change
and acceptance: [contract](REFERENCE_GCF_BATCHES.md).
