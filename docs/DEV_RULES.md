# DEV_RULES

Use existing feature branch; local synthetic fixtures only, then exact-SHA NH synthetic tests. Do not run demo CRAMs or submit Slurm in this task. No source edits on NH. Existing global hpc-rule/loop-engineering files were not found in prior searched Codex roots; explicit project rules apply.

## Automatic whole-GCF batches — 2026-11-11

The user provides the full CSV; preprocessing groups intact GCFs into deterministic
byte-balanced batches. Metadata checks stay global, full source SHA256 moves into
each independently published batch, and the final registry merges every complete
batch without another full payload copy. Existing scientific semantics, strict
validation, and production jobs remain protected. Scope, CLI, output-layout change
and acceptance: [contract](REFERENCE_GCF_BATCHES.md).
