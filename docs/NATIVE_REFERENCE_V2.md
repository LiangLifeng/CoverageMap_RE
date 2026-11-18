# Original reference database v2

## Same-sample query correction (2026-10-24)

Native depth-band queries include same-sample IDs and identical CRAM content by default; `--exclude-self` is opt-in. An empty +/-1X band falls back to same-GCF nearest available depth (all exact ties), labelled outside-band with descriptive scores and null p/acceptance. A one-reference pool still emits similarity; only an absent indexed GCF is skipped. Existing filtered/unsupported cases are reported explicitly as COMPLETE_WITH_GAPS and corrupt inputs still fail. Full contract: [DEPTH_BAND_CONCORDANCE.md](DEPTH_BAND_CONCORDANCE.md). This supersedes the earlier native self-exclusion/no-fallback defaults; legacy/nearest replay defaults are unchanged.

## Query-depth band and multi-layer simulations (2026-10-23)

The default native SQLite query now follows [DEPTH_BAND_CONCORDANCE.md](DEPTH_BAND_CONCORDANCE.md): full-FAI mean depth, versioned native windows, inclusive +/-1X high-depth reference pool and Tools high-coverage maximum; queries below5X use multi-depth/N-repeat simulations with all/count/percent source selection and optional source-depth restrictions. This supersedes nearest-per-reference defaults for native queries; explicit `--comparison-policy nearest` retains the previous workflow. Read the contract for input/output, reference pair-distribution definition, non-estimable cases, integrity and synthetic acceptance.

Default builder has no depth arguments and stores originals only. Builder parameters for depth grids, simulation counts, seeds, batch sampling, sampler selection, simulation matrix format, insufficient depth and the old concurrency alias have been removed entirely; unknown arguments fail before side effects. Query-time sampling parameters remain valid in coverage_concordance.py. Use a new output directory for v2. A legacy library is not converted by rebuilding SQLite.

## Build and index

To migrate a per-contig database, use a fresh output and re-extract or add `--reuse-geometry-from /absolute/old-database` to the build command. Build a fresh SQLite for the new output.

Edit examples/local.sh or examples/slurm.sh: CPUS, MEMORY (empty means scheduler default), DATABASE, WINDOW_SIZE=1000, JOBS, paths. DRY_RUN=1 generates offline plans/scripts only. DRY_RUN=0 builds and collects all GCFs, then creates DATABASE.sqlite with full checksum verification. Slurm defaults -p cpu and bounds queued/running GCF jobs by JOBS. No extra memory request by default.

```bash
python -B /absolute/source/bin/ref_depth_database.py run --mode local \
  --records /absolute/records.csv --references /absolute/fai_path.csv \
  --outdir /absolute/new-native-db --window-size 1000 --cpus 4 --jobs 4
python -B /absolute/source/bin/reference_index.py build \
  --database /absolute/new-native-db --out /absolute/new-native-db.sqlite --verify-payloads
```

## Stored files

Per GCF/sample/attempt: geometry/ (lossless original aligned blocks and sampling-unit pointers), edges.npy, native_coverage.npz, result.json, timing.json, checksums.json, status/. Optional source.cram.crai only when supplied. Original CRAM retained by its owner; absent CRAI supported with existing HTSlib reminder. No coverage/depth directories, replicate matrices or sampled-ID archives in native build. Geometry validation may perform small replay checks, not build a reference simulation distribution.

New native databases use [GCF virtual windows](GCF_CONCAT_WINDOWS.md). The default is fixed 1000 bp over the complete GCF; `--bin-count` is an alternative. The NPZ carries global bin coordinates, a multi-contig segment mapping and coverage measurements. Original mean depth uses eligible CIGAR M/=/X mass over the full FAI length. Fragment policy counts aligned mate blocks, and the geometry preserves sampling-unit boundaries.

Result records sample/GCF, policies, layout/M5/missing header contigs, source hashes/stats, original aligned bases/mean depth/unit/block counts, breadth, window count/size, unweighted window-depth quantiles, CPU/wall/RSS. Root GCF_INDEX.tsv, SAMPLE_INDEX.tsv, NATIVE_INDEX.tsv, RESOURCE_INDEX.tsv, SUMMARY.json are build catalogs. SQLite is the downstream lookup file, not a coverage payload; preserve it with the database directory. --database-root supports moving payloads. Numeric index (gcf,full_mean_depth,sample) enables two range seeks for nearest candidates. Entire sample metadata remains in record_json. Zero-mass originals recorded, excluded from nearest usable references.

## Query

```bash
python -B /absolute/source/bin/coverage_concordance.py \
  --reference-sqlite /absolute/new-native-db.sqlite --query-manifest /absolute/queries.csv \
  --outdir /absolute/new-query-output --query-jobs 2 --reference-cache-mb 256 \
  --reference-limit 1 --replicates 100 --seed 123456
```

Without depth arguments, target is query original mean depth. --depths "0.5,5" can request explicit targets; targets above query mass are recorded unsupported. Select same-GCF nearest original mean depths by absolute difference, deterministic sample tie break. --reference-limit defaults 1, positive values allow more; filters and default same-source/sample exclusion occur BEFORE nearest selection. No similarity-based candidate selection or invented depth tolerance. Include-self is explicit self-consistency only.

For each selected reference, common target=min(requested query depth,reference original mean). Downsample query to it. If reference original depth is higher, generate N independent whole-unit reference draws on demand, then compare at matching target; record actual overshoot and selected-ID hashes. Simulations are only in bounded per-worker RAM, not persisted into the reference library. Cache covers original pools and simulated matrices, reused within same-GCF worker; no persistent cross-run simulation cache yet. Original query extraction occurs once; --query-jobs/--query-scheduling control query concurrency, independent of builder JOBS. Native sampler explicit permutation(default) or lazy; auto rejected for reproducibility.

If reference is full at common depth, compare its original vector directly, mark descriptive_only (no manufactured full-depth replicate distribution/p-value). Dynamic comparisons use window mean depths for Pearson to account for terminal short windows; Jaccard uses positive windows. Existing normal approximation remains uncalibrated. Zero-variance reference distribution or constant/single-window Pearson is non_estimable with null p/acceptance; direct descriptive similarity can remain available. No valid statistical significance claim from a one-window genome.

summary.tsv/concordance.json record original reference/query means, requested/common/actual achieved depths, absolute difference/ratio, selected reference, comparison mode, scores/statistical status, audit hashes, self/source scope. stage_timing.tsv and provenance.json contain exclusive stages, per-worker cache stats, wall/user/system CPU and cumulative process peak RSS. Cache budget accounts retained numeric payload/metadata, not total RSS/temporary working sets. No runtime/resource speedup claimed without a matched benchmark.
