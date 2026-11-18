# Raw query legacy feature parity — v1.1.0

Authoritative feature baseline: native-resume-reference-draws-20261113/code/main.nf,
commit17f8fd89431024ae01c825eaecda0c61a65f3ebb. The raw SQL/CRAM execution substrate
remains one-query-one-process with one decode per selected reference.

## Supported query options

The old top-level query params (except the two explicitly waived budget controls)
are declared and wired: manifest,index,outdir,alpha,window_size,bin_count,
reference_draws,seed,depth_threshold,min_samples,sampling_unit,
coverage_alignment_policy,coverage_unit,jaccard_coverage_threshold,
jaccard_coverage_thresholds,jaccard_methods,bootstrap_replicates,
bootstrap_confidence,query_cpus,diagnostics,diagnostics_cpus,
diagnostics_bin_dir,query_memory,query_time,coverage_conda_env.
cache_mib/work_mib are accepted compatibility placeholders and intentionally ignored.
Profiles local/slurm/nhb exist; Slurm/NHB omit explicit memory/time per project policy.
Local query_memory/query_time remain configurable.

- jaccard_methods: distinct normal,empirical,minimum,bootstrap, including comma lists.
  Reuses the original query_decision_methods.evaluate without formula changes.
- bootstrap_replicates:10..10000; bootstrap_confidence strictly0..1, original
  source-cluster bootstrap with repeated source multiplicities and deterministic seed.
  Bootstrap runs when bootstrap is selected; minimum does not use bootstrap.
- The original Normal+Pearson output remains summary.tsv. Requested method-specific
  Jaccard+Pearson calls are in decision_methods.tsv, exactly as in the old entry.
- Scalar/grid Jaccard thresholds use the original write_sensitivity and the same
  frozen draws/covered-bases arrays. No re-decoding or extra sampling for thresholds.
  default grid:0,.01,.05,.10,.15,.20,.30; primary scalar remains first.
- coverage_unit fraction|percent validates nonempty manifest Coverage;
  manifest_coverage retains the original declared value, as in the legacy output;
  actual coverage still comes from CRAM unique covered bases/full FAI length.
  Optional CRAI/BedGraph remain validated/tracked inputs when supplied.
- coverage_alignment_policy primary_only|primary_plus_supplementary. Single-pass
  in-memory extraction groups all eligible split segments into the same read/pair,
  requiring exactly one eligible primary, equal original query length, non-overlap
  in original-read query coordinates. Reuses original query_interval validation.
  Secondary records remain excluded. No persistent alignment geometry is built.
- diagnostics true|false: QUERY_DIAGNOSTICS uses frozen nulls and summaries, with decision labels
  and thresholds from decision_methods.tsv. A single method retains the original
  PNG names; multiple methods write separate method directories indexed by plot_index.tsv.
  Supports diagnostics_cpus and optional diagnostics_bin_dir, produces per-query
  PNGs, plot_index.tsv, exit code/stderr files. A plot failure makes workflow fail
  while preserving scored-query outputs/caches.

## Compatible outputs

Per query: summary.tsv/json, draws.tsv.gz, null.tsv.gz, query/reference windows,
query/reference covered-bases arrays, window_lengths.npy, original breadth.tsv.gz,
sensitivity.tsv and jaccard_thresholds/t*/zero_unions/draws/null/summary,
decision_methods.tsv plus raw decode audit. Plain draws/null TSV remain for the
previous raw-query interface.

Collected: summary.tsv/json, sensitivity.tsv, decision_methods.tsv,
negative.tsv,errors.tsv,unscored.tsv,completion.json (plus raw counts.json).
Execution failures use the legacy ERROR/NO_RECEIPT status in summary and errors.tsv
and completion ERROR counts; source errors are not silently rescored/replaced.
Healthy query caches survive another query failure and resume.
Reference-support/depth insufficiency and ZERO_QUERY_DEPTH remain cacheable
NOT_EVALUATED results in unscored.tsv, not execution errors. Invalid inputs,
reference decode failures or changed files remain nonzero ERROR tasks.

The raw SQL candidate/depth selection contract in RAW_CRAM_QUERY.md remains
explicit; feature parity does not claim old geometry candidate-set or bitwise
whole-cohort equivalence. Formal clinical runs remain user-owned.

