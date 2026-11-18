# Automatic whole-GCF reference batches

## Approved scope (2026-11-11)

The user continues supplying the complete records CSV and references CSV to
`build_reference.nf`. Starting from b00fe99, preprocessing must keep every GCF
together, balance batches using total source-file bytes, build isolated batch
databases, and combine all successful batches into one metadata index. Do not
stop, update, or restart an existing production workflow. Validate locally and
in a new NHB checkout, then deliver a new code directory; no GitHub publication.

The problem is a global full-CRAM SHA256 pass blocking every build task. This
change does not change read/fragment geometry, coverage windows, thresholds,
clinical definitions, source data, or source/index integrity requirements. No
production speed ratio or completion-time estimate follows from toy tests.

## Input and preprocessing

- Existing required CSV columns, duplicate Sample/GCF and reference-ID checks,
  assembly IDs, reference mappings, FAI layout checks, optional CRAI checks,
  source readability and source size/mtime snapshots remain strict.
- PLAN reads file metadata to obtain byte sizes, but does not hash full CRAM
  content in a production run. The source identity becomes complete inside its
  batch, before checkpoint lookup; source mutation against the PLAN stat fails.
- For N records and G nonempty GCF groups, create
  `min(G, ceil(N / chunk_size))` batches. `chunk_size` is a target average record
  count, not a hard cap; an individual GCF is never split. Assign descending
  GCF byte totals to the lightest batch, breaking ties by record count and
  stable batch/GCF IDs. Resource worker counts do not affect membership.
- Export a deterministic `batches.tsv` with ID, record/GCF counts and bytes.
  Source size is an estimate of I/O workload, not decoded memory or runtime.
- Explicit smoke scope remains visibly partial and retains its existing
  full-roster metadata checks and selected-source hashing contract.

## Execution and publication

`PLAN -> independent GCF batches -> completeness merge -> one SQLite index`.

Each REFERENCE_CHUNK runs bounded record workers, completes full source SHA256
and existing pre/post-build mutation checks, and publishes its own complete
database under `outdir/batches/<batch-id>/database`. Checkpoints and work
directories remain isolated per batch. Content-addressed record checkpoints
use one-record process-pool messages, avoiding repeated serialization of an
entire large GCF batch for each record. Existing record identities are retained.
The checkpoints
remain reusable only after current source-content validation. Nextflow chunk
wrapper caching is disabled because the lightweight PLAN no longer carries
full source hashes; `-resume` still reuses verified record checkpoints.

Every invoked reference script and its Python helpers must be staged as task
inputs and placed first on PATH. This fixes inherited PATH resolving a script
from a different checkout and makes helper changes visible to Nextflow.

The global collector rejects missing, duplicate, foreign, failed, changed or
incomplete batch/record receipts. It verifies exact global record and GCF
membership, batch database identities, complete summaries and byte/count
accounting. It publishes `outdir/database/BATCHES.json`, the global immutable
plan, input snapshots and summary, without copying all batch payloads again.
This directory is an aggregation manifest, not a standalone native database.

The final existing registry builder fully validates every batch payload and
creates `outdir/reference.sqlite`, with global depth/sample selection only after
all batches have been included. Its source roots point to the batch databases;
those directories must remain available. No per-batch support filtering or
SQLite concatenation is allowed. Final SUCCESS follows validated indexing.
Existing eager-plan single-database receipts remain readable by the Python
collector/index path; new batches require a fresh outdir. No existing completed
database is silently replaced.

## CLI and outputs

The existing Nextflow command remains valid. `--plan_cpus` defaults to 64;
`--cpus_per_chunk` bounds each batch's workers and `-qs` bounds queued tasks.
`ref_build_pipeline.py chunk` accepts optional `--database` for its independent
batch publication; absence retains legacy record-only behavior for frozen plans.
Batch-mode chunks require this argument. New outputs are the batch databases,
`batches.tsv`, aggregation receipts, one global index and `reference_build.json`.

## Correctness and acceptance

Before business implementation, the allowed changes are the builder/batch
helpers, reference workflow/modules/config, focused tests and these contracts.
Required acceptance:

1. Exact disjoint/exhaustive Sample/GCF and complete-GCF partition; byte totals
   and deterministic membership across manifest order and CPU counts.
2. No full source hashing in production PLAN; full source SHA in each batch;
   missing/changed sources and optional-index corruption still fail strictly.
3. Concurrent worker bounds, checkpoint replay, same-size/mtime content changes,
   interrupted/failed batches and no incomplete final publication.
4. Rejection of missing/duplicate/foreign/tampered receipts, keys and batch roots;
   immutable publication/replay and no overwrite of conflicting outputs.
5. Whole-GCF batch index versus eager build: identical record/GCF membership,
   source hashes, geometry, aligned/covered arrays and depths, with global support
   semantics intact. Filesystem paths/provenance may differ.
6. Existing reference/index/geometry regression tests; portable script, version,
   PATH and staged-helper checks; strict Nextflow parser/stub DAG and a real
   synthetic end-to-end run in a new NHB run directory. No real cohort build.

Command records and exact local/transport/HPC commit identities belong in the
source-external project control directory and STATUS. Current production jobs
and older code directories remain untouched.
