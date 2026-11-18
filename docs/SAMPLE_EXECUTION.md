# Native sample scheduling and bounded collection

The user-approved 2026-10-21 scope optimizes the observed deployment: an outer
Slurm allocation running `run --mode local`. Native local execution schedules
individual sample–GCF records across long-lived spawned worker processes.
`--jobs` bounds total active sample tasks, including tasks sharing one GCF.
Largest source files are queued first (ties by GCF/sample); this is a scheduling
heuristic, not a change to sample selection or an expected-speedup guarantee.
Workers load the immutable plan once, reuse their interpreter, and keep native
math-library threads at one; geometry compression still uses the explicit cpus
budget. A sample-level exclusive lock protects attempt allocation/reuse.
The controller owns GCF progress files; workers own sample progress/logs.

Storage remains GCF/sample/attempt with fixed packed9 and unchanged numeric
payloads. Existing successful samples are verified before reuse; corrupt/failed
attempts are preserved and rebuilt in a new attempt. Completion must reject
incompatible same-GCF known M5 identities, changed sources, corrupt/missing
required payloads and failed/interrupted tasks. Failure cannot become success
merely because all tasks have returned.

`--collect-jobs` defaults to 4 and bounds final native validation threads and
in-flight sample metadata. All integrity checks remain; a per-GCF FAI identity
is checked once before validation and rechecked after the group. Output sizes
come from one directory walk per validated sample. Only the controller writes
catalogs, in plan order, regardless of completion order. Collection reports its
own processed/total counters and does not present the build's 100% as validation
completion. Serial collection (`--collect-jobs 1`) must produce identical tables.

Built-in `--mode slurm` retains one GCF job per scheduler submission and its
existing bounded jobs/refill policy; it benefits from bounded final collection.
No scheduler submissions are authorized for acceptance. Legacy depth builds
retain GCF-local execution. New code requires new native plans; do not rewrite
old plans, source copies used by active workers, existing libraries, or IH data.

Acceptance: one heavy GCF demonstrably uses multiple persistent sample workers
while respecting jobs; independent dense windows and matched sample payloads
equal serial GCF builds; deterministic serial/parallel catalogs; sample locks,
retry/reuse/corruption/source changes/cross-sample M5 conflicts remain enforced;
failure and interruption stay non-successful. Test synthetic local execution,
mock scheduler regression and full NH suite on the exact committed SHA.
No performance claim relative to Tools.py or cohort throughput without matched
resource measurements. The completed IH library is not rebuilt for this change.
