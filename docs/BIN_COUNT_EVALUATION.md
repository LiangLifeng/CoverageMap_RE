# Stratified Bin Count evaluation CLI

## Scope and authority

User-authorized local implementation and synthetic validation only. Continue existing CoverageMap feature branch; no deployment, remote scientific run, Slurm or publication. Input frozen preflight directory and original native SQLite; SQLite resolves source roots and layouts, geometry supplies aligned intervals. No FASTA/CRAM, source mutation, sample inference/subsampling or clinical outcomes. Keep all accepted sample–GCF pairs (not global sample filtering). Default genome strata [1000,10000), [10000,100000), [100000,infinity) bp. Defaults bin sizes [1,50,100] bp; per stratum B=floor(actual shortest genome/size), deduplicate B without losing size aliases. Range alternative uses inclusive upper limit only when naturally reached by start+k*step, never forces endpoints. All GCFs in the stratum have exactly B windows: w=L//B, edges [0,w,...,(B-1)w,L], remainder entirely in last bin. Contigs concatenated in manifest order without gaps or invented aligned blocks. Reject B<1 or B>shortest length.

## Interfaces and modules

`bin/bin_count_evaluate.py` argparse entry: `plan`, `build`, `evaluate`, `run` (all in sequence).

- plan: --preflight-dir --index --outdir; optional --genome-size-breaks-bp 1000 10000 100000; --bin-sizes-bp 1 50 100 OR --bin-size-start-bp/--bin-size-stop-bp/--bin-size-step-bp; --source-root-map JSON (explicit source-id to absolute root, only for known relocated geometry). Hash frozen preflight payloads and original index; validate exact record/index/layout identities, unique Sample+GCF key, counts/thresholds/membership and source pointers. Store resolved rows, plan, stratum/candidate summary and projected disk/working memory. Planning requires no geometry decompression. Empty strata explicitly skipped. No automatic sample/GCF cap.
- build: --plan PLAN_DIR --outdir; --max-matrix-gib default64 bounds total raw+normalized on-disk numeric payload before work. Read verified geometry once per record, reuse pool for candidate bins. Geometry inner checksums/status validated by existing load_reference_bundle; result/checksum SHA match SQLite. Validate layout names/order/length/M5 and source record identity/mass/depth. Save raw mean depths and row-normalized sums=1e6 as float64 disk-backed NPY, ordered row keys and per-GCF edges. Chunk/row processing; no per-base production allocation; npy matrices source-external.
- evaluate: --matrices BUILD_DIR --outdir; --distance euclidean [cosine jensenshannon correlation] (default euclidean); --ordination tsne [pcoa none] (default tsne); --seeds default20261103; --perplexity default30; --iterations default1000; --learning-rate defaultauto; --jobs default1; --max-working-gib default8. Explicit limits reject work, never silently subsample or replace methods. No automatic perplexity adjustment; n<=perplexity reports ordination unavailable but retains distance scores. PCA-like warm reduction is NOT silently inserted. t-SNE uses chosen metric through precomputed distances, random init and recorded seed/parameters. PCoA uses the same distance matrix, records full spectrum and negative-eigenvalue magnitude; 2D plots are approximations. PCoA seed-independent, computed once. none scores only.
- run accepts combined options and writes plan/, matrices/, evaluation/ into a new exclusive directory.

Core in `bin/bin_count_eval_io.py` (frozen input/planning/build/persistence), `bin/bin_count_eval_math.py` (windows, metrics, centers/scores and embeddings), `bin/bin_count_eval_render.py` (scatter/score curves). All reused code stays unchanged.

## Scoring definitions (accepted implementation proposal)

Fourth distance `correlation` = 1 - Pearson r is defined in the
[Pearson extension contract](BIN_COUNT_PEARSON.md), including constant-profile
errors and near-constant diagnostics. Mean-profile centers and aggregation remain
the same for all four distances. Existing matrices can be evaluated directly.

Evaluate in original normalized coverage space, NEVER embedding coordinates. Each GCF center is arithmetic mean of its normalized rows; for cosine/JS this is explicitly the mean-profile center, not a claimed minimizing metric barycenter. Euclidean/cosine/JS used consistently for sample-to-center and center-to-center. JS converts profiles to probability sums1 and uses square root JS divergence (natural log); cosine clips floating errors and rejects zero norms. Full nonzero finite/nonnegative rows required; zero-mass row is an error rather than imputing composition. Per-GCF within mean then equal-weight mean across GCFs; between is mean over unordered distinct center pairs. R=between/within. Single group/zero within return null with explicit reason; never infinity or fake best. Winners per stratum and distance only, retaining ties. No comparison of scores across metrics as a universal ranking. GCF labels never used to fit ordination. Same-sample rows across GCFs are not declared independent patients.

## Outputs, safety and validation

New output dirs only, SUCCESS after complete publication; checksum ledger for consumed inputs, source code and all produced numerical artifacts. No overwrite/resume. Plan persists identities/row key JSON tuple (Sample,GCF), source roots/layouts/candidates and source SHA. Build persists input SHA, raw+normalized matrices, bin edges, row order and geometry provenance. Evaluation verifies saved matrix hashes, scores, per-row and per-GCF within distances, centers, all center-pair distances, embeddings, PDF+PNG figures, versions, parameters, ordination status and warnings. No 317-color mandatory legend: larger groups use a color key table and plot legend omitted above20 GCFs. Plots never imply all GCFs should separate.

Tests: hand-computed remainder windows and independent dense interval coverage incl cross-contig windows; arbitrary multi-contig geometry packed+npz/npy reuse; mass/mean/normalization; identities/hash corruption/missing geometry/source-map/duplicate rejection; exact shortest stratification boundaries/empty and B aliases; Euclidean/cosine/JS against hand/SciPy references; unequal GCF sample sizes prove macro weighting; single group/zero within; PCoA Euclidean reconstruction and negative spectrum; t-SNE determinism and invalid perplexity; end-to-end CLI plan/build/evaluate/run with frozen synthetic preflight and multi-root SQLite, re-evaluation reproducible scores, no overwrite and memory guard; PDF rendering reviewed. Scope local synthetic, no performance/speedup or real scientific-result claims.

## Required input closure, verified read-only

User-selected preflight: `/home/lianglifeng/share/lianglifeng/CoverageMap/Data/step.3.query/01.vrius_bin_count/viral_bin_count_preflight_V20261103_G1000_S10`.
User-selected index: `/home/share/lianglifeng/CoverageMap/Data/step.2.ref_index/taxonomy_split_V20261103/RefPool_IH_MAG_vOTU_V20261103_c50_all_viral.sqlite`.
SQLite sources resolve geometry beneath these three NHB roots (each root and representative geometry metadata found through NH's g8220100 shared mount; not all payloads revalidated yet):

- `/home/share/lianglifeng/CoverageMap/Data/step.1.ref_depth_database/RefPool/RefPoolIH_V20261026`
- `/home/share/lianglifeng/CoverageMap/Data/step.1.ref_depth_database/RefPool/RefPoolMAG_V20261026`
- `/home/share/lianglifeng/CoverageMap/Data/step.1.ref_depth_database/RefPool/RefPoolvOTU_V20261026`

No new biological input is required. The filtered directory alone is not a self-contained geometry archive. Run the CLI on the host with these stored paths, or explicitly supply --source-root-map and the actual index location; never silently substitute inputs. No remote readiness claim for new SciPy/scikit-learn/Matplotlib dependencies; installation/testing in this task is local only.

## Example staged commands

Use full checked-out CoverageMap source (entry point imports sibling modules and the existing geometry reader), not an isolated copy of only the CLI file.

```bash
python /absolute/source/bin/bin_count_evaluate.py plan \
  --preflight-dir /absolute/preflight \
  --index /absolute/original_viral.sqlite \
  --genome-size-breaks-bp 1000 10000 100000 \
  --bin-sizes-bp 1 50 100 \
  --outdir /absolute/new_plan
```

```bash
python /absolute/source/bin/bin_count_evaluate.py build \
  --plan /absolute/new_plan --max-matrix-gib 64 --outdir /absolute/new_matrices
```

```bash
python /absolute/source/bin/bin_count_evaluate.py evaluate \
  --matrices /absolute/new_matrices \
  --distance euclidean cosine jensenshannon \
  --ordination tsne pcoa --seeds 20261103 \
  --perplexity 30 --iterations 1000 --jobs 1 \
  --max-working-gib 8 --outdir /absolute/new_evaluation
```

`run` combines the same options, placing plan/matrices/evaluation under one new --outdir. `evaluate --ordination none` calculates scores and score-curve figures without an embedding. A second evaluate command against the same matrices reuses them and does not read the original geometry. Centroids remain arithmetic mean profiles for all four distances. The score ratio is descriptive and not a calibrated statistical test. For each stratum/metric, candidate_maxima.json retains all tied maximum-ratio candidate counts; ordination changes never alter this choice.

Pairwise matrices and PCoA eigendecomposition scale quadratically/cubically; max-working-gib is a conservative admission estimate, not a hard process RSS cap. max-matrix-gib limits projected raw+normalized disk payload, not total output storage or geometry RSS. First deployment still requires exact-source HPC acceptance and runtime checks; do not launch production from this local-only result.

## 2026-11-04 Shared-filesystem mapping lifecycle

NHB synthetic assertions passed but TemporaryDirectory cleanup failed errno39 under matrices/stratum_02/B31; directory later empty, consistent with delayed NFS mapped-file cleanup. Explicitly close all raw/normalized NPY mappings on normal and exceptional exits, not merely flush/drop Python references. Apply locally via context manager/ExitStack; test holding mapping references across normal/error context exit. Numerical formulas/inputs/CLI unchanged. Rerun local and exact-commit NHB9+regression tests on shared TMPDIR before publishing runnable command.
