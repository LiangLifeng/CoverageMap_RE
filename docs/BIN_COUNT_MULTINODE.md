# Multi-node Bin Count runs

The existing run command uses one node by default. Under a Slurm allocation,
`run --nodes N --jobs T` uses one worker process on each of N distinct nodes
and up to T threads per worker. Request the same number of nodes and tasks,
one task per node, and at least T CPUs per task. All ranks need the same
source, Python environment, input data and shared output path. The launcher
checks the allocation before creating output and uses srun's bad-exit
handling. It does not submit an allocation.

## Output and scheduling

- One frozen plan feeds every build rank. Rank r handles rows r, r+N, … in
  each stratum, retaining per-record geometry reuse across its candidates.
  Ranks write exclusive shards. The coordinator verifies their seals,
  parent-plan SHA, rank identities, candidate metadata and exact row keys,
  then restores the frozen row order in standard matrices. Allow space for
  both the temporary matrix shards and merged matrices. Their combined
  numeric payload is checked against the matrix limit.
- Evaluation assigns each (stratum, B, distance) once using a stable
  longest-estimated-work-first heuristic. Workers validate their selected
  matrix files. The coordinator verifies rank results, configuration and
  embedding completeness, then emits the ordinary score tables and plots
  in canonical order. Each task stays on one node, including its pairwise
  calculation and any embedding. Scores-only mode skips full pairwise
  calculation and ordination as in the single-node path.
- The run only seals after complete validation. It retains dispatch, host
  and timing receipts under `_multinode/`. Successfully merged temporary
  shard directories are removed. Failure leaves diagnostics and an
  unsealed output. Retry uses a new output directory.

Use `--nodes 2 --jobs 64` within a two-node, 64-CPU-per-task allocation.
Standalone `build` and `evaluate` also accept `--nodes N --jobs T`.
For independent distance jobs, run `plan` and `build` once. Then invoke
one `evaluate` per distance with the same sealed `--matrices` and a
distinct output directory. With Slurm, use an afterok dependency on the
build job. The staged worker receipts live alongside each output in an
exclusive `.multinode` directory. `evaluate --min-bin-count 2` explicitly
excludes the one-bin endpoint when computing Pearson on a shared matrix
that includes it. The selection is recorded in the evaluation manifest,
and the original undefined-profile handling remains in effect for included
candidates.

Keep the original matrix, scoring and ordination options when comparing
one versus multiple nodes. Report job queue wait separately from in-job
wall time. A synthetic acceptance timing is specific to the fixture and
storage, and does not predict a full cohort's runtime.
