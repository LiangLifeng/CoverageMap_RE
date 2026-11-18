# GCF virtual window contract

New native reference builds and Nextflow v2.2 use `gcf-concat-v1`. The same GCF layout and query edges define every coverage vector in a comparison, including simulated references. The Jaccard and Pearson statistics and their detection thresholds are unchanged.

## Contig coordinates

`natural-ascii-digit-runs-v1` compares each name as alternating case-sensitive text and ASCII digit runs. Text runs use Unicode code-point order and digit runs compare as nonnegative integers. If token sequences compare equal, the original name breaks ties by code-point order. For example, `contig_1`, `contig_2`, `contig_10` appear in that order. Leading zeroes remain part of the original name.

The full FAI defines the name set and lengths. FASTA and alignment header entries are matched by name, including lengths and known sequence digests where available. Missing alignment header contigs in FAI mode occupy their full length and contribute zero coverage. Each contig has virtual offset equal to the sum of preceding lengths. Layout metadata records names, lengths, offsets, sort policy and known M5 values. Coordinates are 0-based and half-open. Alignments are mapped within their original contigs, and FASTA and CRAM remain separate and unchanged.

## Window boundaries and payloads

- `--window-size S` or `--bin-size S` gives boundaries `0, S, 2S, …, L`. The final window may be shorter. Nextflow uses `--window_size S`.
- `--bin-count B` gives exactly B nonempty windows with edges `floor(i*L/B)` for `i=0..B`. Nextflow uses `--bin_count B`. B must be between 1 and L inclusive. L is the sum of contig lengths.
- For lengths 1300, 800, 700, the fixed 1000 bp edges are `[0, 1000, 2000, 2800]` and the three balanced-bin edges are `[0, 933, 1866, 2800]`.
- `edges.npy` contains the global edges. The native NPZ stores `global_start`, `global_end`, `length`, and `segment_ptr` (one pointer per window edge). Slice the `segment_contig_index`, `segment_local_start`, and `segment_local_end` arrays using adjacent pointers. A segment's contig name comes from the layout.
- Query `windows.tsv` has one row per contig segment with a shared zero-based `bin_index` and global bounds for segments of the same window. Query `audit.json` has the layout and window policy. Native vectors store aligned-base sums, covered-base interval unions, mean depth and covered fraction. Both fractions use the true window length.

The native index binds each sample to its declared policy, geometry, parameters and payload checksums. A query at a different window size or count reconstructs reference vectors from verified geometry on the query's exact edges. An exact-edge match can use the validated native vector directly. Repeated simulations use those same edges. A window crossing a contig boundary aggregates separate aligned blocks from the contigs.

## Migration and legacy results

Build to a fresh directory and create a fresh SQLite. With a verified saved geometry bundle, use `--reuse-geometry-from /absolute/old-database` on `ref_depth_database.py prepare` or `run`. This verifies the old result and source identity, checks layout/read policies, converts each old global block to its contig-local interval and applies the new offset by name. Sampling-unit order and aligned-base mass are preserved. Missing, inconsistent or unreadable reusable geometry fails explicitly. Run without the reuse option to extract from the source alignments.

Keep the old database, SQLite and result artifacts for historical replay. Re-indexing an old database does not convert its payloads. Recompute every coverage vector, similarity and reference distribution used by the new workflow. Nextflow and `query_manifest_v2.py` reject old indexes, which remain readable by the earlier version and the legacy paths.

- The native depth-band and nearest CLI can read the original per-contig policy. For new indexes, depth-band uses the chosen GCF geometry and nearest uses the stored geometry.
- The standalone legacy manifest grid-distance rules and old Bin Count matrix plans retain their original geometries. The grid-distance thresholds require a legacy index. New Bin Count plans from a `gcf-concat-v1` index use balanced bins, and shard merges validate the window policy.
- `Tools.py`, standalone FAI-order alignment/bundle replay and retained legacy depth matrices keep their historical coordinate contracts. A standalone comparison using a naturally sorted geometry bundle reads the query into that same layout.

## Acceptance

`python -m pytest tests` includes direct checks of ordering, offsets, size/count boundaries, windows spanning three contigs, overlapping and empty coverage, invalid inputs, and incompatible policy/layout/edge failures. The synthetic native smoke builds from FAI and BAM in different contig orders, migrates an old database without modifying it, indexes the new databases, and scores both low- and high-depth query arms.
