# Pearson correlation distance — 2026-11-04

User authorizes a fourth selectable distance on the accepted Bin Count feature
line: `--distance correlation`, defined as `1 - Pearson r`. All existing
distances/defaults, frozen records, strata, candidate Bin Counts, normalization,
mean-profile centers, equal-GCF within weighting and unordered-center-pair
between weighting are unchanged. Scores use the original normalized vectors.
Pearson r is similarity; the final ratio R remains mean between distance divided
by mean within distance. No r-squared, absolute-r conversion or significance
test is introduced.

`evaluate --matrices EXISTING/matrices --distance correlation` reuses sealed
matrices and verifies their hashes, without reading SQLite, geometry, FASTA or
CRAM. It writes only to a new exclusive output directory. The previous three
distance results remain intact. Use `--distance euclidean cosine jensenshannon
correlation` only if intentionally recomputing all four distances. Both t-SNE
and PCoA use the selected dissimilarity as precomputed input; scores do not
depend on the embedding. Correlation dissimilarity is not guaranteed Euclidean;
the existing PCoA negative-spectrum reporting remains required.

## Undefined and near-constant profiles

Fewer than two bins or an exactly constant row/center makes Pearson undefined.
Fail the evaluation with a nonzero exit and no SUCCESS, retaining a
`pearson_diagnostics.tsv` and `FAILURE.json` in the candidate's correlation
directory. Report sample/GCF row identifiers and center GCF identifiers; never
drop records, impute zero correlation, add pseudocounts, change binning or fall
back to another metric. This numerical guard is not a biological filter.

Near-constant profiles with nonzero range are retained and explicitly warned
when `norm(x-mean(x)) < 1e-13 * abs(mean(x))`, following the SciPy pearsonr
documentation's numerical warning criterion. They appear in the same diagnostic
table. All distances must be finite; correlation roundoff is clipped to [0,2].
The diagnostic table always exists for a correlation candidate, even if empty.
The ordinary single-GCF/zero-within null-ratio status remains unchanged.

Correlation within-row and between-center tables additionally export
`pearson_r = 1 - distance`. Score JSON also records macro-averaged within r and
pair-averaged between r, separately from ratio R. Arithmetic mean-profile
centers are retained, not redefined as optimizing correlation barycenters.

## Acceptance and deployment

Verify against independent scipy.stats.pearsonr calculations, including positive,
negative and zero correlations, scale/offset invariance, bounded symmetric
pairwise distances, unequal GCF sizes, constant rows, constant centers,
near-constant warnings and one-bin rejection. Exercise a saved-matrix CLI with
all four distances and both embeddings; prove evaluation works after synthetic
SQLite/geometry are removed, matrix hashes/row order stay identical, repeated
scores agree and existing outputs cannot be overwritten. Preserve the three
legacy formulas and regression tests. Local tests precede a clean exact-SHA NHB
synthetic acceptance. No real evaluation or scheduler submission by the agent.

References:
- https://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.distance.correlation.html
- https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.pearsonr.html
