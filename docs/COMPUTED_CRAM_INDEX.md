# Computed CRAM index (2026-11-21)

User inputs: BuildInput_all_virus_347993.csv, fna.csv, fai_path.csv. New
`build_cram_index.nf` defaults to computed indexing, with chunk size100 and four
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

DAG: plan validates CSV associations and emits fixed-row chunks; per-chunk
ProcessPool workers decode each row once; serial SQLite merge copies only small
scalar tables and adds catalogs. No shared concurrent SQLite writer. No new
Slurm memory/time requests; executor queue size bounds cluster concurrency.
Cache invalidation uses ordinary Nextflow file fingerprints; it does not detect
malicious changes that preserve all file metadata. Use a fresh output directory;
never replace the current in-use database.

Parameters: --index_chunk_size100 --index_chunk_cpus4 --index_min_mapq0
--index_m5_policy geometry_only. Old --mate_audit is rejected by the computed
entry because flags are fully scanned here. The old metadata-only Python API is
retained for existing callers and explicitly does not recalculate values.

Acceptance: known CIGAR/flags/overlap/zero/subset-header BAM and CRAM metrics
versus independent dense counts; correct R1/R2 identities and sums; invalid
source Depth/Coverage ignored; unchanged SQL field types/order/views/index and
old query candidate lookup; serial/parallel equality; disjoint/exhaustive chunks,
missing/duplicate/policy-mismatched shards rejected; failure/no overwrite;
strict native Nextflow initial/resume/changed-source and stub acceptance.
