# Diagnostic coverage breadth (v1.2.0)

With `--diagnostics true`, each query PNG additionally reports Query coverage
breadth and Ref draw mean, minimum–maximum and READY draw n/N. Breadth is
unique covered bases / full GCF reference length, displayed as a percentage;
this is independent of the Jaccard occupied-bin threshold. Reference draws are
weighted equally, including repeated source selections. Low uses saved
matched-depth draws; high uses saved original-reference draws. The existing
route label identifies which applies. Breadth annotations do not change scoring or decisions. N is the requested reference draw count (`--reference_draws`, default 100); all plot denominators and slot checks follow N.

Query uses `summary.tsv:actual_query_coverage` (fraction 0–1). Reference uses
`reference_covered_bases.npy`, `window_lengths.npy`, `draws.tsv.gz` adjacent to
`null.tsv.gz`, and `summary.tsv:gcf_reference_length`. Validate identity, unique
contiguous draw slots, integer covered counts within actual window lengths,
full-length agreement and READY status. Sum covered bases before dividing by
full length, including uncovered positions; never average per-window fractions
or substitute `reference_covered_bins / windows`. Non-READY draws are excluded
and incomplete n is explicit. Corrupt or inconsistent inputs are errors. Missing
legacy coverage inputs display NA and PLOTTED_WITH_GAPS; never substitute zero.

The workflow stages existing per-query result directories by symlink so the
plotter sees the saved scores and coverage inputs together. No large arrays are
copied or rebuilt by the renderer. Arrays are memory-mapped, scanned one draw
at a time, and mappings explicitly closed. Standalone `--null-root` and
`--work-dir --trace` retain their interfaces; old flat `--staged-nulls` works
with NA if its coverage siblings are absent. PNG Description metadata includes
the numeric coverage values, n and definition for audit. Similarity panels and
normal p-values remain unchanged.
