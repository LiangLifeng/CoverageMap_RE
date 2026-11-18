# Read and fragment BAM benchmark

`sample_cram_modes.py` is imported unchanged from the user-provided NHB script
(SHA256 `66c13c6ffc7968c80601045aa5fce305e6077f3f706254da8be81776c0f50e04`)
before the scoped startup and checkpoint changes. Sampling and eligibility rules,
permutation RNG, whole-unit overshoot and paired-cohort inclusion are preserved.

The manifest column `CramPath` accepts BAM or CRAM. Detection uses file content,
and query reads stream sequentially without BAI/CRAI. Geometry-only BAMs can omit
SEQ/QUAL. The full indexed layout supplies the depth denominator. Nextflow's
`main.nf` now stages the saved-score plotting utility and enables per-query PNG
diagnostics: red query mean and blue reference row-mean distribution for each
metric. Numerical scoring remains at the f7e805bf baseline.
Coverage breadth annotations and their saved-input requirements are specified
in [QUERY_DIAGNOSTIC_COVERAGE.md](QUERY_DIAGNOSTIC_COVERAGE.md).

The sampler accepts `0.5`, `0.5x` or `0.5×`. A sibling `source/` checkout is found
automatically, or use `COVERAGEMAP_ROOT`. All FAIs are checked for accessibility
before output creation. On NHB, copy an authoritative FAI to the shared project
if its original group mount is absent on compute nodes; verify identical SHA256.

`--resume` requires the same script, manifest and index hashes, input alignment
paths/sizes/nanosecond mtimes, FAI paths/hashes, seed, depth, tag and policies.
Completed per-mode checkpoints require matching BAM SHA256 and independently
recounted aligned bases and reads. An advisory output-directory lock prevents
concurrent writers. Publish BAMs only after validation, then save atomic mode
checkpoints. Uncheckpointed BAMs fail explicitly instead of being overwritten.
Only queries with both modes completed enter both final manifests. An older
failed directory without checkpoints must be retained separately before starting
the new run; it is never silently imported.

This batch has user-confirmed TP labels only. FP counts/rates are not applicable.
Compare detected TP / all paired input queries, both metric NOT_EVALUATED counts,
zero null-SD counts and the distribution of nonmissing null SD, and measured
mean depth. Preserve missing/error/unevaluated rows and report their denominators.
Use the locked B=60, min_samples=9, alpha=.05, depth threshold=5 and seed=401.

## Legacy native index compatibility

V2 queries retain the f7e805bf GCF-concat window definition when reading an
older per-contig native index. Validate the original registered metadata,
checksums and geometry under its original index policy. Rebase aligned block
coordinates into natural contig order, preserving sampling-unit order,
pointers and per-unit aligned-base masses, then bin those blocks into the
query's exact GCF-concat edges. Never combine, interpolate or relabel old
window totals. Reuse native totals only when contig coordinates and all edges
match exactly. Original indexes and payloads remain read-only. Each query audit
records both window policies and whether geometry was rebinned.

Acceptance compares independent legacy and GCF-concat native libraries built
from the same multi-contig BAMs (reversed header order, unequal contig lengths,
query windows crossing boundaries). Low/high queries must have byte-identical
scores, query/reference matrices, null distributions and window tables;
reference draws must retain identical unit IDs and seeds. Also reject source
blocks crossing a contig boundary and verify all original inputs remain hashed
unchanged after execution.
