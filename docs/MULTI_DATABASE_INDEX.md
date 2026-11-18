# Multiple reference databases and auditable depth selection

## Same-sample query correction (2026-10-24)

Native depth-band queries include same-sample IDs and identical CRAM content by default; `--exclude-self` is opt-in. An empty +/-1X band falls back to same-GCF nearest available depth (all exact ties), labelled outside-band with descriptive scores and null p/acceptance. A one-reference pool still emits similarity; only an absent indexed GCF is skipped. Existing filtered/unsupported cases are reported explicitly as COMPLETE_WITH_GAPS and corrupt inputs still fail. Full contract: [DEPTH_BAND_CONCORDANCE.md](DEPTH_BAND_CONCORDANCE.md). This supersedes the earlier native self-exclusion/no-fallback defaults; legacy/nearest replay defaults are unchanged.

`reference_index.py build` now accepts repeated `--database` arguments. Schema3
registers libraries without moving/copying their payloads; completed libraries are
required by default. Explicit sample omission is described below.
Native and legacy modes may each be indexed, but cannot be mixed. Scientific
configuration, same-GCF names/order/lengths and known M5 must agree, including
references later excluded. Operational CPU, Python paths and geometry codec may differ.

Each reference resolves its own registered absolute root in serial/process query
workers. Old schema1/2 indexes remain readable. Single-root `--database-root`
relocation remains supported; multi-root relocation via that ambiguous flag fails
explicitly (rebuild metadata using the new roots). Sources must stay accessible.

## Multi-CPU build

Add `build --jobs N` to use N spawned worker processes for per-reference file
reads, checksum/decoder validation, JSON parsing and SQL-row preparation. Default
`--jobs 1` keeps serial operation; zero, negative, boolean and noninteger API values
are rejected. This is CPU concurrency, not the number of libraries or SQLite
writers. Request corresponding CPUs in the external scheduler allocation.

Workers never open or share the output SQLite connection. One controller writes
in original plan order, preserving row IDs, source precedence, counts, selections,
omission details and schema. Up to `2*N` batches of16 references are submitted at
once; results are consumed in submission order, with pending work cancelled and
running workers joined on failure/interruption. Sources are indexed sequentially,
each using the same worker budget, so multiple roots do not multiply concurrency.
Processes receive root/config only once; no complete cohort plan is replicated.
SQLite aggregation, integrity checks and exclusive publication remain serial.

All verification defaults and `--verify-payloads`, `--skip-invalid-samples`, external
TSV audits, corruption handling, source-change checks and exclusive outputs retain
their existing semantics. Serial and parallel builds on identical sources produce
identical logical tables and metadata; synthetic acceptance also checks complete
SQLite file SHA equality. Existing indexes and reference payloads are untouched.
`filter`, `report` and `inspect` do not gain a jobs option.

Example: append `--jobs 32` to an existing `reference_index.py build` command.
More CPUs do not imply a measured speedup: small libraries can be slower due to
process startup, and shared storage/SQLite can limit scaling. No production
throughput or optimal worker count is claimed without a controlled measurement.

## Counts and overlapping libraries

Sample identity is `(GCF, sample ID)`. Repeated/aliased input roots fail. Copies
across distinct roots are collapsed only when source content SHA, scientific sample
metadata, full layout/M5, registered numerical payload SHA maps and geometry
decoding descriptor SHA match. Optional archived CRAI presence does not change
sample identity; its provenance remains in the original manifest. Different
encodings conservatively conflict; equal depths alone never establish identity.
Conflicts fail without publishing. Both storage locations and manifest/checksum
registry hashes remain in `reference_sources`; lexically first resolved root is the
canonical payload location, with no silent fallback. `duplicate_copies` and
`input_references` distinguish copies from unique references. Same source with
*different* sample IDs is not silently identified as one patient; self-exclusion
by source SHA at query time remains in effect. Samples are not independent patients.

Full verification (`--verify-payloads`) verifies files once per source reference,
including copies, before publication. Default metadata verification does not claim
full payload validation. Temporary SQLite files, per-row streaming and indexed
aggregation avoid loading geometry into memory. INDEX progress prints each input
root and every100 samples. No real database rebuild, alignment extraction or job
submission is part of indexing.

## Explicit omission of invalid samples

`build --skip-invalid-samples` allows an incomplete source SUMMARY and omits only
individual samples that lack a successful validated attempt or fail sample-level
metadata/file/checksum/decoder validation. The default remains strict. This does
not repair a failed sample or change SUCCESS/checksum files. `--verify-payloads`
keeps full hashing for samples actually indexed; required-file/checksum registry
checks also remain active when full hashing is not requested.

Each sample uses an SQLite savepoint: reference rows, depth rows and counters are
all rolled back on omission, even if the error occurs late during row insertion.
SQLite/systemic output failures, interruptions, root plan/summary parse errors,
source changes, scientific configuration/layout/M5 or cross-root duplicate conflicts
remain fatal. Per-sample read/validation errors are printed as INDEX SKIP with
GCF/sample/type/reason. Progress counts processed/indexed/skipped separately.

The SQLite `skipped_samples` table and report `skipped_samples.tsv` retain source
root, GCF, sample, sample directory, exception type and reason. `skipped_samples`
counts SOURCE occurrences, not unique biological samples. A healthy copy from
another source can still be indexed normally. Metadata records the explicit flag,
planned_input_references, successful input_references, skipped_samples,
partial_index and each source's original summary completeness/matching state.
A successfully published partial index exits0 but clearly prints INDEX PARTIAL;
that is not a claim that every source sample completed successfully. If every
sample is invalid, the index remains inspectable with zero usable references.

Depth distributions/counts use only successfully validated UNIQUE references.
`gcf_depths.tsv` additionally reports skipped_copies per GCF (source occurrences).
Refiltering preserves omission audit/provenance and cannot restore invalid samples;
a new build is required after their source has been repaired. Existing schema3
indexes without an omission table remain readable/reportable with zero audit rows.
The wrapper supports SKIP_INVALID_SAMPLES=1 (default0).

## Keep omission details only in an external TSV

Add `--skipped-samples-tsv /absolute/new/skipped_samples.tsv` to `build` alongside
`--skip-invalid-samples`. The final SQLite has no `skipped_samples` detail table;
only total/source/GCF counts and external audit path/SHA/row count remain. Staging
indexes may contain errors temporarily; final VACUUM removes deleted detail bytes
before publication. This affects invalid-sample audits only, not valid reference
metadata or ordinary depth-filter selection flags/reasons.

The TSV is generated during build, not delayed until report. Both destinations
must be new and outside source databases. Its filename is user-selected, and the
TSV is published before SQLite; a failed SQLite publication removes only the TSV
created by that attempt. An interrupted process may leave an orphan TSV; no index
is published before its audit. Pre-existing user files are never replaced.

`filter` without the new option inherits the external TSV reference and keeps
errors outside SQLite. `filter --skipped-samples-tsv NEW_PATH` can also export/remove
errors from an existing in-SQLite audit, or copy an existing external audit to a
new path. It creates a new SQLite, leaving the original index unchanged. `report`
verifies the external audit SHA/count/header and copies its full rows into the
report's skipped_samples.tsv; a missing/changed TSV is an explicit report error.
Keep the external TSV for later audit/report; valid reference queries do not load
it. The wrapper accepts `SKIPPED_SAMPLES_TSV=/absolute/new/skipped_samples.tsv`.

## Filtering order and depth definition

1. Keep each sample if `full_mean_depth >= --min-depth` (inclusive).
2. Keep a GCF if its distinct eligible sample count is `>= --min-samples`.
3. Query only retained samples in retained GCFs. Self/source exclusions remain
   query-specific and can reduce usable references further.

`full_mean_depth = aligned_bases / full FAI genome length`, following stored
primary/read policy, with missing-header contigs contributing zero. This is NOT
median covered-region depth or depth at each genome position. Zero depth is kept
in metadata/selection at min-depth0, but nearest queries require positive depth.
Defaults0/1 introduce no empirical quality cutoff; users specify their thresholds.

Every original reference remains indexed with `eligible_depth`, `selected` and
`exclusion_reason`: BELOW_MIN_DEPTH or INSUFFICIENT_SAMPLES_AFTER_DEPTH_FILTER.
`counts` always describes unique indexed products; `selected_counts` describes the
current selection. Zero selected GCFs is a valid inspectable selection, not a usable
reference guarantee. Query reports missing/excluded GCF before alignment extraction.

## Reports and changing thresholds

`report --index INDEX --outdir NEW_DIRECTORY` exports deterministic ordered TSVs:

- `gcf_depths.tsv`: total/eligible/selected counts, selected state/reason, and all
  versus depth-eligible min/q25/median/q75/max/mean. Quantiles interpolate linearly
  at `(n-1)*q`; empty distributions are blank, not zero.
- `sample_depths.tsv`: full mean depth, selection flags/reason, canonical root,
  source SHA, manifest and number of copies.
- `reference_sources.tsv`: all canonical/replica locations and their identities.
- `selection.json`: input roots, plan/summary SHA, scientific config, thresholds,
  verification provenance and counts.

`filter --index ALL.sqlite --out NEW.sqlite --min-depth VALUE --min-samples COUNT`
reuses audited indexed metadata, with no geometry reads or payload hash pass. It
re-evaluates ALL retained metadata rows, so thresholds may be tightened or relaxed.
The source SQLite hash and inherited verification provenance are recorded. It does
not newly verify live payloads; query loaders still verify the selected payloads.
Older schema1/2 indexes require one rebuild to obtain these metadata tables.
All output SQLite targets are exclusive/atomic and read-only after publication;
existing SQLite/library files remain unchanged. Report directory must be new.

## Example (thresholds supplied by user)

```bash
python -B /home/lianglifeng/01.WorkDir/01.Project/03.CoverageMap/source/bin/reference_index.py build \
  --database /home/share/fandcj/DataShare/CoverageFeatureMap/RefPool/RefPoolCheck \
  --database /home/share/fandcj/DataShare/CoverageFeatureMap/RefPool/RefPoolIH \
  --out /home/share/fandcj/DataShare/CoverageFeatureMap/RefPool/RefPool_all.sqlite \
  --verify-payloads

python -B /home/lianglifeng/01.WorkDir/01.Project/03.CoverageMap/source/bin/reference_index.py filter \
  --index /home/share/fandcj/DataShare/CoverageFeatureMap/RefPool/RefPool_all.sqlite \
  --out /home/share/fandcj/DataShare/CoverageFeatureMap/RefPool/RefPool_selected.sqlite \
  --min-depth "$MIN_DEPTH" --min-samples "$MIN_SAMPLES"

python -B /home/lianglifeng/01.WorkDir/01.Project/03.CoverageMap/source/bin/reference_index.py report \
  --index /home/share/fandcj/DataShare/CoverageFeatureMap/RefPool/RefPool_selected.sqlite \
  --outdir /home/share/fandcj/DataShare/CoverageFeatureMap/RefPool/RefPool_selected_report
```

The wrapper example supports original DATABASE/SQLITE_INDEX/VERIFY_PAYLOADS plus
DATABASES_FILE (one absolute database path per line; exclusive with DATABASE),
MIN_DEPTH and MIN_SAMPLES; outputs the report beside the SQLite.

Acceptance: multi-root native/legacy serial/process queries, actual per-root
payload resolution, source/self exclusions, exact boundary depth then sample count,
empty/relaxed selections, known-M5/config/duplicate conflicts, audited copy collapse,
old-schema reading, deterministic reports, original hashes unchanged and exclusive
publication. Synthetic tests only; production library indexing is user-run.
