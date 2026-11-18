# Explicit taxonomy partition

`bin/split_reference_taxonomy.py --index /absolute/all.sqlite --taxonomy /absolute/RefFaiManifest.csv --outdir /absolute/new_directory`

Requires a complete native schema-3 index with zero skipped samples. Exact assembly_accession matches with domain bacteria or viral enter separate physical SQLite indexes. Unmatched, archaea, fungi and other domains remain only in the unchanged parent and appear in excluded_gcfs.tsv. No domain inference from ID prefixes or source libraries. Conflicting classifications and unknown SQL objects fail; identical duplicate classifications are accepted.

The new directory contains two `<parent_stem>_{bacteria,viral}.sqlite`, full gcf_classification.tsv, excluded_gcfs.tsv, split_audit.json and a terminal SUCCESS marker. Existing output directories are refused. All files are staged and validated before exclusive publication; an interrupted publication lacks SUCCESS and must be investigated rather than overwritten.

IDs, numerical metadata, JSON, paths, all source roots, depth-selection flags/thresholds and query indexes are retained. Each retained table is compared cell-for-cell against ordered parent rows. Scoped aggregate counts are recomputed; parent aggregates and input hashes remain in taxonomy_partition metadata. Sources and their full-build counts/fingerprints retain original provenance, explicitly labelled as such. The original empty external skipped-audit dependency remains and is hash-checked. No CRAM, geometry, coverage or source pool is copied or modified. Payload validation is inherited, not repeated. No acceleration claim is made.

Acceptance: synthetic partition/other-domain exclusion, duplicate reference provenance, unselected rows, preserved IDs/JSON/indexes/thresholds, source byte immutability, existing output refusal, changed-input detection, conflicts/partial/unknown-schema rejection; exact committed cluster fixture test before production. Every output receives integrity_check, foreign_key_check, row comparison and SHA-256; both input hashes are rechecked before publication.

## Scoped integrity and progress repair (2026-11-03)

Parent integrity/foreign keys are checked once. Each output checks only `main.integrity_check` and `main.foreign_key_check`; attaching the parent must not trigger redundant full-parent scans. Full retained-row comparisons and all hash checks remain. Stderr now reports timestamped START/DONE/FAIL, elapsed seconds, per-table copies/readback (every50,000rows), indexes, integrity, hashes, publication and completion. CLI flags, filtering and output schema stay unchanged. Regression attaches a corrupt unrelated database: unscoped integrity fails, while the repaired splitter correctly validates only its output. No measured production speed ratio or ETA is implied. Editing the file does not accelerate a Python process that has already loaded it; do not stop/restart user jobs automatically.
