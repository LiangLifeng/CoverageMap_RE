# CLI

Geometry capture uses node-local temporary scratch by default. Optional
`COVERAGEMAP_SPOOL_DIR` selects an existing absolute directory and
`COVERAGEMAP_SPOOL_MIN_FREE_MIB` sets its nonnegative free-space reserve
(default1024). These are worker environment settings, not scientific CLI
parameters. See [capture I/O](REFERENCE_CAPTURE_IO.md) for field decoding,
scratch selection, errors, cleanup and capacity limits.

python bin/coverage_map.py --alignment PATH --reference PATH --outdir NEW_DIRECTORY [--reads N | --fraction F | --depth D] --bins B --replicates R --seed S [--mate-policy both|read1|read2] [--sampling-unit read|fragment] [--min-mapq N] [--exclude-duplicates] [--exclude-qcfail] [--save-pool].
All scripts support --help. Output collisions fail. Emit full.tsv, simulations.tsv and provenance.json; pool.npz is optional. Human targets are explicit; no default numerical sampling target.
--reference-index PATH is mutually exclusive with --reference and selects alignment-field-only decoding plus FAI layout validation. This does not decode or reconstruct read/reference DNA.


## GCF reference depth database v1

The user-approved directory hierarchy, explicit-depth CLI, local/Slurm modes,
source identity, failure/availability semantics and acceptance requirements are
fixed in [REF_DEPTH_DATABASE.md](REF_DEPTH_DATABASE.md). Real build deferred.


Depth database entry point: `python -B bin/ref_depth_database.py`.

- `prepare --records CSV --references CSV --outdir NEW_DIRECTORY (--depths LIST | --depth-grid START:STOP:STEP) [--cpus 4 --memory 4G --replicates 100 --bins 1000 --seed 401]`: validate/count all selected records and generate only a plan and scripts.
- `execute --database DIR --mode local --jobs 1`: bounded GCF processes, each using configured CPU threads; retries preserve prior attempts and reuse verified successes.
- `execute --database DIR --mode slurm`: regenerate scripts only, never submit.
- `submit --database DIR [--retry]`: explicit sbatch call for an authorized operator. Unknown submissions are never duplicated automatically.
- `status --database DIR [--watch --interval 10 --scheduler]`: filesystem status, optional squeue annotation; watch collects after all workers have terminal records.
- `collect --database DIR`: verify source stats/reference hashes and product checksums, then write complete and usable catalogs; incomplete returns nonzero.
- `worker --database DIR --gcf ASSEMBLY`: internal shared worker used by both modes.

No default depths. CPU/memory are per GCF job; local total CPU budget is jobs*cpus, memory is advisory locally. Default compressed depth matrices remain lossless int64. Missing target depth is explicit UNSUPPORTED_DEPTH (or --insufficient-depth error).

coverage_concordance.py also accepts --reference-depth-manifest PATH (mutually exclusive with --alignment/--reference-bundle). Select sample_manifest from READY_DEPTH_INDEX.tsv. Requested depths must be available exactly; cached bins/replicates and policies are authoritative. The test remains an independent draw; --seed controls it, while cached_reference_seed records the saved reference RNG.


## One-command construction (2026-10-18)

`ref_depth_database.py run --mode local|slurm` accepts the same required manifests,
new outdir and explicit depth options as prepare, plus jobs/interval. It performs
prepare -> local execution or Slurm submit -> wait -> verified final collection.
Normal success requires no intervening command or confirmation. Existing individual
subcommands remain available for inspection/recovery. pipeline.json records phases,
terminal state and exit code. Errors remain nonzero, incomplete catalogs remain
visible, and jobs absent from squeue for three consecutive successful polls without
terminal worker files cause NEEDS_ATTENTION rather than an endless wait or false
success. sbatch ambiguity is never duplicated. A controller interruption preserves
submitted jobs; component commands allow recovery. No automatic retry of data or
resource failures. Tests use synthetic data and mocked schedulers only; no real run.


## SQLite reference index / batched query extension

User-approved interfaces, immutable index, one-extraction query concurrency, unchanged storage/kernel semantics and acceptance are specified in [SQLITE_BATCH_CONCORDANCE.md](SQLITE_BATCH_CONCORDANCE.md).


## Current builder CLI cleanup (2026-10-19)

ref_depth_database.py accepts original-only construction with --window-size (1000bp) and one --jobs concurrency budget. Historical depth-grid/simulation/fixed-bin options and hidden concurrency alias are removed, not deprecated aliases. New native configs/results omit unused simulation parameters. [Current builder contract](REF_DEPTH_DATABASE.md) supersedes earlier depth-builder commands; coverage_concordance query sampling controls remain separate.

## Native local sample pool (2026-10-21)

Native local --jobs now bounds sample workers, with persistent spawn processes; built-in Slurm continues bounded GCF submission. --collect-jobs INT defaults4 for native final verification, independent of compute jobs/cpus. Collection emits separate verified/total/failed progress. See SAMPLE_EXECUTION.md for locks, logs, deterministic publication and new-plan requirements.

## Automatic whole-GCF batches — 2026-11-11

The user provides the full CSV; preprocessing groups intact GCFs into deterministic
byte-balanced batches. Metadata checks stay global, full source SHA256 moves into
each independently published batch, and the final registry merges every complete
batch without another full payload copy. Existing scientific semantics, strict
validation, and production jobs remain protected. Scope, CLI, output-layout change
and acceptance: [contract](REFERENCE_GCF_BATCHES.md).
