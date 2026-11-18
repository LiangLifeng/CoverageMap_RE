# Original reference database builder — current contract (2026-10-19)

`ref_depth_database.py` constructs one original reference per GCF/sample. The prior multi-depth builder CLI is removed, without hidden aliases. [Native storage and query contract](NATIVE_REFERENCE_V2.md) defines numerical layout and output fields.

## Accepted CLI

- prepare: required --records, --references, --outdir; --window-size (1000 bp default), --cpus (4), optional --memory (omitted default), --partition (cpu), --python, --geometry-format (auto/packed/npy/npz); explicit CSV-column mappings and read/fragment/filter policies.
- run: prepare arguments plus required --mode local|slurm, --jobs (1), --interval (10 seconds), --dry-run.
- execute: --database, --mode, --jobs, --dry-run.
- submit: --database, --jobs, --retry, --dry-run.
- status: --database, --jobs, --scheduler, --watch, --refill, --interval. Refill requires watch+scheduler; otherwise monitor never submits.
- collect: --database. Internal worker: --database, --gcf.

Removed parameters are rejected as unknown: --depths, --depth-grid, --replicates, --bins, --seed, --batch-size, --sampler, --matrix-format, --insufficient-depth, --max-inflight. They are not persisted as unused fields in new native configs/results. The corresponding query sampling controls in coverage_concordance.py are separate and remain valid.

## Safety and execution

CSV Sample/CramFilePath/GCFID plus GCF_ID/GCF_path; nonempty inputs, duplicate identities/ambiguous references/invalid paths fail before output. Immutable plan pins code, input manifests and source stats. Full FAI defines reference coordinates including zero-covered omitted header contigs; known length/M5 conflicts fail. No source modification. Optional adjacent CRAI validated/archived when supplied; absent CRAI remains supported by sequential extraction and retains HTSlib reminder.

prepare/dry-run only write plan, scripts and pending catalogs. run prepares -> bounded local workers or Slurm submission/refill -> monitor -> verified final collection. One GCF worker processes samples sequentially, threaded geometry compression/readback/hash uses cpus; jobs bounds GCF concurrency (Slurm queued/running plus per-round attempts). Local memory is advisory; Slurm uses scheduler memory defaults unless explicitly set, clearing inherited SBATCH_MEM flags when omitted. Default child partition cpu, explicit -p cpu on submission. No time/account/node directives.

Per-GCF exclusive locks, numbered sample attempts, verified successful reuse, corruption rebuild into new attempt, ambiguous submission records, no automatic failed/unknown retry. A lost scheduler entry never means success. Three successful missing-job polls terminate with NEEDS_ATTENTION; raw scheduler query errors terminate. No automatic lock stealing. Changed code/config requires a new output plan. Do not resume old depth-grid construction with this CLI.

## Products and integrity

Store original geometry, contig-fixed-bp coverage NPZ and edges, full source/policy/layout/native-depth metadata, timing/checksum/status. No reference replicate matrices or sampled-ID archives. No new simulation seed/grid fields. All required geometry payloads, coverage, edges, result, timing and any declared optional index must be readable/nonempty and registered; verify digests before reuse/collect/index when full verification enabled. Original CRAM/FAI remain needed for collection/recovery checks, query can use relocated saved payloads without raw reference input.

Collect GCF_INDEX, SAMPLE_INDEX, NATIVE_INDEX, RESOURCE_INDEX and SUMMARY. Failure/pending samples stay visible. Peak RSS is cumulative per process, not additive across samples. Compression auto compares measured lossless candidates, no global-optimality claim. SQLite is published separately by reference_index.py build; provided local/slurm wrappers invoke it after complete build. Wrapper DRY_RUN=1 never indexes or starts jobs.

## Acceptance

Reject every removed option before output/worker/scheduler effects; inspect help/config to ensure none survive. Verify native local/generated-Slurm-script numerical identity, synthetic CRAM both with/without CRAI, full FAI subsets, required file and source integrity, deterministic windows/union breadth against independent dense baseline, empty/short references, concurrency/dry-run/refill/failure/retry invariants. Full NH regression at identical local/bare/checkout SHA. Older saved depth libraries remain readable by downstream readers; their direct-library test fixtures do not expose a builder compatibility CLI. No real construction or Slurm submission authorized during this cleanup.
