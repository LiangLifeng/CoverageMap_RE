# Build progress and logs

`ref_depth_database.py run` and `execute --mode local` print timestamped progress
to stderr and append to `<database>/build.log`. No extra dependency is required.
Output is flushed and line oriented, including when redirected or run by Slurm.
Local execution reports immediately, every 10 seconds, and on exit. Routine GCF
start/successful-exit lines are omitted; failed exits retain IDs, exit codes and
worker log paths. Slurm `run` reports
at each `--interval` poll (default 10 seconds). Phase messages identify execution,
collection/integrity verification, completion, failure or interruption.

The bar measures processed sample–GCF records (`successful + failed`) / planned
records. Success, failure and reused-success counts are separate. A 100% bar
does not mean all records succeeded or final integrity collection passed.
GCF counters separately report terminal/total, successful, failed/interrupted,
running and pending. Reused samples are included in successes, never added twice.
Current sample IDs (up to five) and elapsed time since this invocation are shown.
No read-level percentage or ETA is inferred. A long single CRAM can leave the
sample percentage unchanged; heartbeat lines continue. Counters come from atomic
worker progress files, without rehashing payloads on every refresh.

`<database>/<GCF>/worker.log` contains flushed sample start/success/failure,
exception messages, final counts and elapsed time in both local and Slurm modes.
Per-sample progress JSON updates after each sample, including failures.
Existing status JSON fields remain and `status` adds aggregate `progress`.

## Observe an existing batch without changing it

```bash
python -B /absolute/source/bin/ref_depth_database.py status \
  --database /absolute/existing-database --watch --interval 10 --no-collect
```

This reads old plans regardless of their source fingerprint, prints progress to
stderr and existing JSON snapshots to stdout. It does not write `build.log`,
collect catalogs, spawn workers, or query/submit scheduler jobs. For an external
Slurm cluster, add `--scheduler` only when using its own login environment.
`--no-collect` rejects `--refill`. Without `--no-collect`, historical terminal
watch collection behavior remains. Stop the viewer with Ctrl-C; the build runs
independently.

Already running processes keep their loaded code. Do not replace source files
under an active batch: later workers check the frozen code fingerprint. Use a
separate updated source copy for read-only status, and keep the old source and
plan for the active batch. New builds with this version require a new plan/output
directory; never patch a frozen plan to force compatibility. Build data formats,
sampling/coverage calculations, integrity checks and scheduler policy are unchanged.
