# Computed CRAM index (2026-11-21)

User inputs: BuildInput_all_virus_347993.csv, fna.csv, fai_path.csv. New
`build_cram_index.nf` defaults to computed indexing, with maximum chunk size100, 4096 fixed identity buckets and four
worker processes per chunk. Formal real-data runs remain user-owned.

Preserve `coverage-raw-cram-index-1`, the exact input column names/order/TEXT
schema in `reference_records`, original ordinal `id`, `reference_gcf_sample`,
`reference_genomes`, `gcf_catalog_audit`, `reference_cram_catalog` and
`reference_mate_audit`. Existing queries continue reading `Depth` and `Coverage`.
Their values are now calculated; source numeric strings are never parsed or
used to choose records. `source_records` retains every original cell and ID,
including excluded unlabelled references. Additional input columns stay intact.

Use the unchanged raw query `extract_once` primary-only, MAPQ0 default,
read-unit geometry with full FASTA/FAI layout. Depth is sum of M/=/X block lengths
/ full reference length; overlapping mates both contribute to depth. Coverage
is union covered bp / full length *100 (percent). D/N gaps do not contribute.
Full FNA versus FAI layouts are checked. Present-header M5 audit is preserved;
default geometry_only matches the user's latest runs, and strict is explicit.
Missing header contigs contribute zero within the full denominator.

`reference_metrics`, keyed by the unchanged record id, stores genome_length,
aligned_bases, covered_bases, read_count, fragment_count, R1/R2/unlabelled counts
and aligned bases, per-side depth and covered bp/coverage percent, layout and M5
evidence. Side depth uses the full genome denominator, not half a genome. These
are full scans, not the former early-exit mate presence lower bounds. Scalar
capacity queries need no CRAM decoding; actual read sampling still requires
geometry or CRAM. No persisted per-read geometry is added by this index.

Retain the existing labelled-reference rule: at least one R1/R2 label. Excluded
unlabelled/empty inputs remain in source_records, reference_metrics and the full
mate audit; no flag repair or inferred mate identity. Decode, ambiguous identity,
reference-layout or strict M5 errors fail that chunk and prevent final publication;
there is no silently partial final database. Successfully computed chunks remain
available for Nextflow -resume. Source-file metadata dependencies and tracked
calculation helpers invalidate affected tasks; no per-file SHA or CRAM copies.

DAG: CSV rows enter the Nextflow channel, combine with exact FNA/FAI mappings,
then group by a fixed deterministic hash bucket of (Sample,CramFilePath,GCFID).
Each bucket is sorted, duplicate identities are collapsed, and subdivided into
chunks of at most100. Per-chunk ProcessPool workers decode each unique identity
once. Serial SQLite merge reads the current original CSV, restores its ordinal
IDs/duplicate rows and every original cell, and copies scalar metrics by identity.
There is no whole-CSV PLAN task. No shared concurrent SQLite writer. No new
Slurm memory/time requests; executor queue size bounds cluster concurrency.
Cache invalidation uses ordinary Nextflow file fingerprints; it does not detect
malicious changes that preserve all file metadata. Use a fresh output directory;
never replace the current in-use database.

Parameters: --index_chunk_size100 --index_chunk_buckets4096 --index_chunk_cpus4 --index_min_mapq0
--index_m5_policy geometry_only. Old --mate_audit is rejected by the computed
entry because flags are fully scanned here. The old metadata-only Python API is
retained for existing callers and explicitly does not recalculate values.

Acceptance: known CIGAR/flags/overlap/zero/subset-header BAM and CRAM metrics
versus independent dense counts; correct R1/R2 identities and sums; invalid
source Depth/Coverage ignored; unchanged SQL field types/order/views/index and
old query candidate lookup; serial/parallel equality; disjoint/exhaustive chunks,
missing/duplicate/policy-mismatched shards rejected; failure/no overwrite;
strict native Nextflow initial/resume/changed-source and stub acceptance.

## SQL examples (no CRAM access)

```sql
SELECT r.id,r.Sample,r.GCFID,r.Depth,r.Coverage,
       m.read1_count,m.read2_count,m.read1_bases,m.read2_bases,
       m.read1_depth,m.read2_depth,m.read1_coverage,m.read2_coverage
FROM reference_records r JOIN reference_metrics m ON r.id=m.id
WHERE r.GCFID=? ORDER BY r.Sample,r.id;
```

The unchanged `reference_records` fields remain TEXT, including numeric strings
for recalculated Depth/Coverage; metrics in the auxiliary table have numeric
SQLite types. Coverage values are percent; counts and capacities are integers.
`source_records` is provenance and must not be used as the candidate pool.

On NHB add `-c <cpu_only.config>,<code>/examples/nhb/computed_cram_index.config`
for explicit /disk scratch. Use `-profile slurm,nhb -qs 100` and omit memory/time
requests. Formal runs are user-submitted; never reuse/overwrite the labelled
production SQL. After a failed run use the same source, inputs, work directory,
outdir and `-resume`; valid chunks are cached. A deliberate completed-dataset
revision requires a new output directory because published indexes are immutable.

## Resume after source changes or input deletions

Keep the bucket count and chunk size unchanged. A CSV deletion/reorder does not
change unrelated buckets; only affected bucket parts are recomputed. Deleting a
duplicate identity may require no metric recomputation. Each CRAM and associated
FNA/FAI is a tracked file dependency, so ordinary content/mtime changes invalidate
its chunk without invalidating unrelated chunks. The final SQL always reflects
the current input membership and ordinal IDs; deleted records cannot persist.
Use the same work/cache and -resume, but a new --outdir after revising a completed
dataset. This keeps the prior complete SQL intact. This design does not reuse
old whole-CSV PLAN caches. Chunk source identities exclude upstream numeric
Depth/Coverage, so correcting only those provenance values reruns merge rather
than CRAM calculations. Fixed buckets intentionally trade small scheduler overhead
for predictable cache reuse; no measured cluster speed ratio is claimed.

## Query scalar-capacity reuse

The existing candidate API still reads recalculated reference_records.Depth.
For replace preflight, the planner attaches only compact reference_metrics scalar
capacities when MAPQ, primary coverage and M5 policies match. The scorer checks
CRAM/FNA/FAI size and mtime before using capacities. Stale or absent metrics and
legacy databases fall back to actual extraction. Audit rows report capacity_source
sqlite or geometry. No SQL-derived geometry or random IDs are invented: retained
references still undergo the original sampling/extraction and M5 checks. Selected
IDs, coverage arrays, scores and replacement RNG remain unchanged. The synthetic
rejection test proves an insufficient reference is not decoded solely to check
capacity, while downstream scientific arrays match the original path exactly.
