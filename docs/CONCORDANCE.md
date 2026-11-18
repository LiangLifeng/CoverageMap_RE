# Matching-depth concordance CLI

## Same-sample query correction (2026-10-24)

Native depth-band queries include same-sample IDs and identical CRAM content by default; `--exclude-self` is opt-in. An empty +/-1X band falls back to same-GCF nearest available depth (all exact ties), labelled outside-band with descriptive scores and null p/acceptance. A one-reference pool still emits similarity; only an absent indexed GCF is skipped. Existing filtered/unsupported cases are reported explicitly as COMPLETE_WITH_GAPS and corrupt inputs still fail. Full contract: [DEPTH_BAND_CONCORDANCE.md](DEPTH_BAND_CONCORDANCE.md). This supersedes the earlier native self-exclusion/no-fallback defaults; legacy/nearest replay defaults are unchanged.

## Query-depth band and multi-layer simulations (2026-10-23)

The default native SQLite query now follows [DEPTH_BAND_CONCORDANCE.md](DEPTH_BAND_CONCORDANCE.md): full-FAI mean depth, 1000bp per-contig windows, inclusive +/-1X high-depth reference pool and Tools high-coverage maximum; queries below5X use multi-depth/N-repeat simulations with all/count/percent source selection and optional source-depth restrictions. This supersedes nearest-per-reference defaults for native queries; explicit `--comparison-policy nearest` retains the previous workflow. Read the contract for input/output, reference pair-distribution definition, non-estimable cases, integrity and synthetic acceptance.

## Saved reference input

`--reference-bundle /absolute/bundle/storage.json` (or the bundle directory)
replaces `--alignment`. The index must remain beside its manifest, checksums,
payload chunks and successful status; an isolated copied index is insufficient.
`--test-alignment` remains the new sample CRAM/BAM. The original reference CRAM
is not opened in bundle mode. Loading/checksums use `--workers` bounded threads.

The saved five extraction policies are inherited. Explicit flags must match the
saved policies; conflicting requests fail instead of silently refiltering a
geometry-only pool. Without a bundle the prior defaults remain unchanged.
The saved ordered contig names/lengths can provide the test extraction layout
when neither FASTA nor FAI is given: a private temporary two-column layout table
is passed to the accepted field-only extractor and removed afterwards. This is
coordinate metadata, not a fabricated FASTA index or DNA reference. Exact layout
and CRAM-header M5 agreement is still required; missing M5 proves no sequence
identity. Explicit FASTA/FAI input remains supported.

Seeds, matching-depth N reference + 1 independent test draws, CPU/thread/process
executors, optional GPU calculation and statistics are unchanged. Input hashes
include all reference bundle payload/metadata and its checksum index, plus the
new alignment and any explicit reference. Verify these again before SUCCESS.
Provenance records bundle sample/schema/manifest identity, inherited policies,
layout derivation, load time and content-based self-consistency scope. Different
input sources are not automatically independent biological samples.

Acceptance: compare bundle and original CRAM CLI outputs for fixed samplers,
including matrices, every selected-ID hash, test IDs and both statistics across
serial/thread/process CPU executors. Reject incomplete/corrupt bundles, policy
and layout/M5 conflicts, unsupported depths and colliding outputs. Run the full
NH CPU suite and real saved-bundle concordance; do not claim new GPU hardware
or scientific p-value calibration from CPU acceptance.

```bash
python -B bin/coverage_concordance.py \
  --reference-bundle /absolute/reference-bundle/storage.json \
  --test-alignment /absolute/new-sample.cram \
  --depths 0.01,1,10,100 --replicates 100 --bins 1000 \
  --depth-sampler permutation --executor thread --workers 4 \
  --outdir /absolute/new-result-directory
```

## Original alignment input

`python bin/coverage_concordance.py --alignment reference.cram --test-alignment test.cram --reference-index reference.fna.fai --depths 0.01,1,10,100 --replicates 100 --bins 1000 --depth-sampler permutation --workers 4 --executor process --outdir NEW_DIRECTORY`

Both pools use the same reference layout and read/fragment/filter policy. Each depth independently draws N reference replicates and one test draw, using complete eligible units without replacement. Depth must be supported by both pools. The test draw's seed role differs from every reference replicate. Using the same input twice is explicitly marked `same-source-self-consistency`, reuses immutable pool storage, and does not constitute independent-sample evidence.

Depths accept a single `--depth`, comma-separated `--depths`, or decimal `--depth-grid start:stop:step`. A grid includes its stop even if the step does not land on it; `0.01:100:0.1` has 1001 points. At least two reference replicates are required. Default N is 100. Each depth runs all N replicates, not one replicate per depth.

CPU execution uses a single bounded executor across depth/replicate batches. Options are serial, thread, and process; process workers mmap immutable numeric geometry rather than receiving pool copies per task. Worker unit-base summaries remain private. At most twice the worker count of batches are pending. Native numerical thread counts are fixed to one per worker. Thread speedup depends on the fraction of Python versus NumPy work; no speedup is guaranteed. Process startup and shared-pool files are charged to workflow time/storage. One worker always executes serially.

The new frontend seed policy `depth-value-role-replicate-v1` hashes the numeric depth and uses independent role/replicate seed streams. Worker count, batch size and depth ordering do not alter fixed-policy draws. This policy deliberately differs from the older CLI's sequential RNG; the older CLI remains unchanged. Auto probes each reference depth in the parent, separately from formal RNG streams, and records time/choice. Replay uses the recorded explicit sampler; auto timing decisions can vary across runs.

The statistics retain the accepted Tools.py low-coverage structure: mean test-versus-reference similarity, reference off-diagonal row-mean similarity distribution, and lower-tail normal z calculation. Jaccard means presence above zero; Pearson uses centered counts. Computation proceeds in row blocks, retaining N-by-B count/working arrays plus block-by-N similarity arrays rather than full N-by-N similarity matrices. Constant coverage makes Pearson non-estimable with null fields; zero coverage is rejected. No normal-approximation calibration is claimed. This is not the original high-coverage maximum-similarity method.

Original Tools.py had epsilon-regularized Pearson, subtraction of a presumed unit diagonal, and a different zero-variance z rule. Accepted corrected Tools.py uses explicit non-estimability, exact off-diagonal exclusion and lower-tail zero-variance behavior. The new CLI retains those accepted semantics. Original-relative native workflows also differ in reads/CIGAR/windows/clipping/cursor policy; do not interpret differing p values or timings as same-output comparisons.

Outputs include per-depth simulations.npy (N-by-B int64), independent test.npy/test_ids.npy, sampled-ID audit archives, concordance.json (both tests, p values, acceptance, sampler, actual masses and timing), ordered summary.json, input/version/policy provenance and terminal status. Output reuse is refused. Failures retain partial artifacts with FAILED; only SUCCESS confirms complete execution. Process shared-pool files are retained and count toward output storage.

Acceptance requires statistical agreement with accepted Tools.py on nondegenerate fixed windows, exact serial/thread/process draws across grid/batch order changes, and independent dense replay. NH must additionally validate exact SHA/input hashes and matching original comparisons, report five fresh-process timing repetitions with medians/ranges, CPU/RSS/storage, and both gains and slowdowns before making efficiency claims. GPU backends are optional; hardware/library availability, exact-count agreement and synchronized transfer-inclusive performance must pass NH acceptance before a speed claim.

## Optional GPU backends

`--coverage-backend gpu --executor thread --workers 4 --batch-size 16` keeps CPU complete-unit sampling and immutable seed streams, sends bounded batches of selected IDs to a single GPU, and accumulates exact int64 window counts. Geometry and block-window indices remain resident; boundary atomic additions plus difference-prefix sums handle arbitrarily long blocks without per-base matrices. A lock serializes CUDA work shared by CPU sampling threads. Multi-process GPU coverage is rejected. `--statistics-backend gpu` independently selects blocked FP64 GEMM for Jaccard intersections and centered Pearson; original z calculations remain on CPU. Defaults stay CPU; unavailable GPU raises an error without silent fallback.

Use the user-specified NH python-hpc-latest CuPy environment. GPU counts and selected IDs must match CPU exactly; nondegenerate FP64 similarities and p-values use explicit acceptance tolerances, and non-estimable states must agree. Record synchronized transfer-inclusive time, first initialization/compile, device memory, CPU time/RSS and storage, and run five fresh-process repetitions against original Tools.py and optimized CPU. Small workloads may regress. GPU sampling is not introduced, and normal-test calibration remains out of scope.

When GPU Python has no compatible pysam wheel, pass `--extraction-python /home/lianglifeng/01.WorkDir/env/conda-envs/coverage-map/bin/python`. This explicitly delegates extraction to the same committed source using that interpreter, exchanges numeric arrays via unique temporary files next to the output, records interpreter/library identities and temporary file sizes, then removes only its own temporary directory. Hashing, subprocess startup, export/import and allocation all count toward workflow time; GPU timing is not measured on a pre-extracted pool when reporting full CLI time. CPU controls can use the same bridge for attribution.

The internal GPU path consumes IDs directly from validated complete-unit selection and avoids sorting them again to prove uniqueness. Public arbitrary-ID calls still enforce uniqueness; both paths retain shape/range checks, and CLI batches verify exact aligned mass. This avoids a redundant CPU O(K log K) check before GPU transfers without changing sample order or results.
