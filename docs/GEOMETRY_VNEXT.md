# Geometry vNext (Issue #7)

The confirmed base is 4b75a61. Geometry capture is independent of runtime sampling.
Primary-only coverage, statistics, thresholds and complete-unit stopping retain
their existing definitions. Read mode selects a read with its eligible segments;
fragment mode selects all observed eligible mates. Fragment coverage is aligned
read-base mass, including mate overlap, never insert-span coverage.

## Resource envelope and execution

The user estimates 100,000 records. Read-only inventory of the supplied vOTU CSV
on 2026-11-11 observed 137,595 records, 11,196 samples, 7,534 GCFs, up to 1,368
records per GCF, and no duplicate Sample/GCF pair. References has 11,185 rows.
An evenly spaced 64-record stat sample spans 2,242–43,717,684 CRAM bytes; this is
not the full maximum and does not bound decoded alignment count. Full-file and
worst-case read-count performance remain unverified.

The later supplied `BuildInput_all_486938.csv` was independently counted at
486,938 rows (96,981,989 bytes, same five input columns). At target size100 this
requests up to 4,870 whole-GCF batches, capped by the number of GCFs. Its per-file extrema and full-build resources have not been
measured; the earlier vOTU sample is not a bound on this larger input.

The DAG is metadata validation/plan -> whole-GCF byte-balanced batches -> bounded
record source hashing/build -> independent complete batch databases ->
completeness collector -> one multi-root SQLite index. See
[automatic GCF batches](REFERENCE_GCF_BATCHES.md). Batch workers do not write
global catalogs. The Slurm profile uses native `executor.queueSize=100`
per workflow, overridable with `-qs`. Per-process chunk maxForks is removed;
at four CPUs/chunk the maximum requested chunk budget is 400 CPUs per workflow.
The user manually budgets concurrent workflows and other jobs. Each
record uses one decoding/geometry worker, with BLAS/OMP/compression nesting capped.
The collector has its own CPU budget. No production-scale completion-time or
storage-expansion estimate is justified before measurements.

Capture uses an indexed on-disk spool in worker-local scratch, then compact memory-mapped numeric tables
and a UTF-8 string dictionary. It avoids a cohort-wide or per-alignment Python
object graph. Numeric projection is O(segments + blocks), grouping O(S log S),
with O(S+B) numeric working memory; disk spool/cache remains bounded. Temporary
storage includes spool, final tables, projected legacy-compatible block payload,
validation reads and eventual collector copies. Shared filesystem I/O, checksum
passes and small-file counts must be measured across the full entrypoint.
The worker default no longer places the disposable spool in the shared output
directory; [capture I/O](REFERENCE_CAPTURE_IO.md) defines explicit scratch and
reserve settings. TLEN is not consumed and is excluded from the CRAM field mask
to avoid unnecessary external-reference loading. Final geometry is unchanged.
Unknown maximum decoded file size is retained explicitly; CPU concurrency must
be chosen from measured per-record RSS/I/O before full-table production.

NHB uses no explicit Slurm memory or time request and no Nextflow memory/time
limit, as instructed by the user. CPU and concurrent chunk bounds still apply.
Generic profiles may request memory; estimates are separate from requests.

## Planning resources and progress

`--plan_cpus` (default 64) controls metadata-validation I/O threads and the CPU
request. The planner checks global keys/mappings, each required FAI once, source
stats/readability and optional CRAI. Full CRAM SHA256 moves into each batch's
bounded record workers. At most twice the metadata worker count of futures is
pending. Canonical GCF/Sample ordering makes partitioning deterministic.

The planner reuses file metadata checks, omits the unused four-byte magic read
for files already selected by the `.cram` suffix, and streams the final plan JSON
once. Full source/index hashes and mutation checks remain required in the batch.
`chunk_size` sets a target average record count; descending source-byte GCF groups
are assigned to the lightest batch. It is not a hard cap, and a large individual
GCF is never split. `batches.tsv` reports actual counts and bytes. Explicit
smoke builds still validate metadata for the entire manifest and hash only the
selected canonical subset. Progress reports reference count, completed records,
worker count and elapsed time. `planning.json` records zero full source bytes
hashed for a production PLAN, total source bytes and the per-batch validation
stage; worker count does not enter scientific chunk identity.

The immutable plan still retains O(records + references) metadata. Each digest
worker uses a bounded read buffer; shared filesystem throughput and metadata
latency may limit gains. The default allows 64 concurrent file checks; override
`--plan_cpus` to tune file access for the target storage. An explicit command-line
value such as `--plan_cpus 2` overrides the default, so omit it to use 64.
Neither a full-table speedup nor a completion-time estimate is established by
synthetic tests.

Unchanged inputs produce byte-identical batch rosters across metadata worker
counts. Batch-wrapper caching is disabled: every replay validates source content
before reusing record checkpoints. Start this output-layout upgrade in a fresh
outdir; subsequent retries may use its original launch/outdir/work with `-resume`
after stopping its old controller. Keep `outdir/batches/*/database`, which the
final SQLite references without another complete payload copy. Existing
single-database outputs are not silently replaced.

## Storage and identities

Schema `coverage-geometry-2` stores Fragment, Read, Segment and Block tables in
numeric `.npy` files, plus a UTF-8 dictionary. All coordinates are 0-based,
half-open in the explicitly recorded contig order. Fragment/read integer IDs
are deterministic first-observation indices; fragment identity additionally
stores a full SHA256 of namespace, QNAME and RG with collision detection during
capture. Missing QNAME uses a per-segment orphan identity and is explicitly
incomplete. Missing RG is an explicit audit flag. Multiple primary records for
one identity remain ambiguous; no mate is invented.

Coverage policy defaults to `primary_only`. Supplementary records, CIGAR, query
interval, flags, MAPQ, mate declaration and observed mates remain saved even
when not counted. Secondary metadata remains distinguishable; secondary coverage
is disabled. Supplementary coverage uses a strict non-overlapping original-read
query-interval policy: overlap, missing primary, or inconsistent query lengths
fails explicitly, preserving the captured evidence for review. Reverse-strand
query intervals are transformed into original-read orientation, including hard
clipping in the full query coordinate system.

Capabilities are explicit. Legacy bundles retain their original read/fragment
capability only, with supplementary_available=false; no mate relation is guessed.
New geometry carries both sampling capabilities; runtime configuration records
sampling unit and coverage policy rather than treating unit as a build constraint.

## Acceptance boundary

Required tests include primary-only aligned/covered/window mass against the old
extractor, read vs fragment whole-unit selection, split/cross-contig relations,
unmapped/incomplete/ambiguous mates, CIGAR query coordinates, checksums/read-back,
serial/parallel equality, receipt conservation, injected failure recovery,
native Nextflow resume and changed-input invalidation. Exact committed code must
pass a small real CRAM NHB Slurm end-to-end run before runtime acceptance. These
checks do not establish full 137,595-record production performance.

## Table fields

Table row numbers are zero-based IDs. `-1` denotes unavailable integer metadata;
AS/NM additionally have explicit presence bits so a real negative score remains
distinct from absence. QNAME is kept only in the temporary collision-check spool.

| File | Fields |
| --- | --- |
| fragments.npy | identity (SHA256 hex), rg dictionary ID, nonexclusive audit flags, observed read1/read2 row IDs (or -1) |
| reads.npy | fragment ID, role (0 unpaired, 1 R1, 2 R2, 3 unknown), primary flags when observed, rg ID, declared mate contig/start/unmapped/reverse, audit flags |
| segments.npy | read ID, type (0 primary, 1 supplementary, 2 secondary), contig ID, reference start/end, original-read query start/end/length, reverse, MAPQ, raw flags, CIGAR dictionary ID, declared mate fields, SA dictionary ID, AS/NM and presence flags, filter eligibility, block begin/end pointers |
| blocks.npy | segment ID, explicit contig ID, contig-local start/end; row number is block ID |
| strings.npy / strings.utf8 | dictionary byte offset/length and UTF-8 data |
| geometry.json | schema, contig order/layout, namespaces, coordinate/overlap rules, capabilities and audit flag legend |

The bundle checksum ledger covers all tables and dictionaries. The default
primary block projection remains readable by existing bundle consumers; rich
capabilities are attached by the current reader and never inferred from pointers.

## Execution and recovery

```bash
nextflow run /absolute/source/build_reference.nf -profile slurm,nhb -qs 100 \
  --records /home/share/lianglifeng/CoverageMap/Data/InputTable_share/BuildInput_vOTU.csv \
  --references /home/share/lianglifeng/CoverageMap/Data/InputTable_share/fai_path.csv \
  --outdir /absolute/run/geometry-vnext --chunk_size 100 \
  --cpus_per_chunk 4 --collect_cpus 2
```

Ordinary `nextflow run` is supported without a custom sbatch wrapper, shared
account lock, environment marker or custom submission retry. The native queue
size counts this workflow's in-flight tasks, including submitted pending tasks;
it does not query or reserve other workflows' account slots. The user observes
the cluster's total100 convention and can choose, for example, `-qs 50` for each
of two simultaneous workflows. `--max_chunks` is retired/ignored with a warning,
including old commands that specified two. CPU count/chunk size remain configurable.
See [Nextflow executor configuration](https://docs.seqera.io/nextflow/reference/config/executor/#executorqueuesize).

NHB omits memory/time directives. Generic
profiles support `--memory_per_chunk` and `--reference_time`. `--partition` sets
the Slurm queue. New native CLI builds use `--alignment-capture complete` by
default; explicit `legacy` retains the old immutable-unit representation.

The final native database is `outdir/database`; its single-writer metadata index
is `outdir/reference.sqlite`. Checkpoints stay under `outdir/checkpoints` and
Nextflow reports/traces under `outdir/pipeline_info`. No chunk writes global
catalogs or the global SQLite index. A failed chunk exits nonzero after saving
record receipts. Successful record checkpoints are content-checked before reuse.
The collector rejects missing, duplicated, changed, failed or skipped results.

Use the same launch/work directory with `-resume`. Planning deliberately runs
again to detect CRAM/CRAI/FAI changes hidden behind CSV paths; unchanged chunk
rosters remain deep-cache inputs. Only the affected roster changes when one
source changes. Changing chunk size/CPU allocation may rerun task wrappers,
while valid per-record scientific checkpoints remain reusable. Control stages
revalidate publication; the global SQLite is built once and verified on replay.
Never remove a stale lock until its owner/job has been checked.

`--smoke_records N` explicitly validates the full input roster but builds only
the first N records in canonical GCF/Sample order. Full and selected counts are
recorded separately, and the summary labels this as a smoke subset. This is not
a full-table build. Production defaults to all records.

Query CLI uses `--sampling-unit read|fragment` and optional
`--coverage-alignment-policy primary_only|primary_plus_supplementary`; Nextflow
Query uses the corresponding underscore parameter names. The same runtime
choice controls query units, reference draws and reference-reference nulls.
Draw audit records target/actual mass/depth, overshoot, selected units, observed
reads/fragments, eligible coverage segments/blocks and existing RNG provenance.
Missing legacy counts remain unavailable, never guessed.

At this cohort size, per-record payload/status files and checkpoint plus final
publication can create millions of filesystem entries and multiple payload
copies. Global checksum validation and collector copies remain material I/O
costs. Retention and full-scale throughput require measurement before production;
small smoke runs cannot establish an expansion factor or completion estimate.
Queue size is execution policy; it does not change payloads, sampling,
scientific thresholds or native task inputs.
Bundled entrypoint scripts use native Nextflow tracking. Python helper-only
changes require an explicit affected-task rerun/new run directory; recursive
helper tracking is not claimed and no source SHA is passed as a workflow input.
