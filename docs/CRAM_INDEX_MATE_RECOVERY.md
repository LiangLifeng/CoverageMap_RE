# Strict R1/R2 index recovery

Use the separate `build_cram_index_mates.nf` entry with the original launch
directory, work directory, session UUID, bucket parameters and original project
checkout path. The original bundled executable remains byte-identical. The canonical CRAM_METRICS_CHUNK
process and its inputs retain the old task hashes. Failed tasks route to a
record-aware worker, which continues after a recognized mate-identity problem.
Already known bad chunks can go directly to that worker through the comma
separated `--mate_recovery_buckets` argument. Use the NHB recovery configuration
so a legacy chunk failure enters this path.

A file qualifies when the eligible mapped primary reads include at least one
R1 and one R2, every eligible read has an identifiable mate role, and the
extractor accepts its read identities. Ambiguous identities, missing labels,
and no eligible primary records are explicit exclusions. Read failures and
reference or file integrity errors still require attention. The final merge
requires every input identity to have measured metrics or an explicit exclusion.
Read counts from successful legacy shards enforce the same gate without opening
the CRAM during merge. Same-side rules apply at read level and do not infer pairs.

The final SQL retains the established schema and original input row IDs for
retained files. Excluded files appear only in the separate
`cram_index/excluded_crams.tsv`, with source row number, identity and reason.
