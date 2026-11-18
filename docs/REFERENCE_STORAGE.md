# Per-sample reference storage v1

Authorized on 2026-10-17: continue `codex/coverage-downsampling` from
`5b4daae7a2b93441803af9b8ec162767b22edb43`. Keep samples separate. Build and
measure one reusable geometry bundle before deciding on a 20,000-sample rollout.

`bin/reference_store.py build` takes explicit sample ID, CRAM/BAM, matching
FASTA or FAI, extraction policies, output directory, worker count and chunk size.
It reuses the accepted extractor without altering read order, complete-unit
boundaries or filtering. Extraction is a serial stream; compression, input hashes
and reading independent chunks use bounded threads. One sample pool is in memory;
this is not an out-of-core CRAM extractor or a cohort launcher.

The bundle stores starts, aligned lengths, unit block counts and ordered reference
names/lengths/M5. Those arrays reconstruct the exact BlockPool. It records sample
identity, source hashes, policy, schema, source revision, versions and counts.
No SEQ/QUAL, query names, per-read flags or MAPQ are retained; changing filters
requires the original alignment. Missing M5 is explicit; FAI does not prove DNA
sequence identity. Original CRAM files must be retained.

Formats compared by `--format auto`: original compressed NPZ; raw NPY for mmap;
and per-unit-chunk byte-shuffled minimal-width unsigned arrays compressed by zlib.
For compact chunks, widths are selected without truncation for starts, lengths
and counts separately. Unit and block order are never sorted. Compression levels
1, 6 and 9 are measured. Each candidate is read back, arrays checked and the
same seeded complete-unit selections/window sums replayed. Candidates are
temporary; only the smallest validated candidate (data + storage-index bytes)
is retained. Auto optimizes measured bytes among these candidates, not global
optimality or query speed. Load/write seconds and logical/allocated bytes are
reported separately. Peak RSS is the complete build process high-water mark,
not exclusive compressed-buffer memory or total future process-tree memory.

Publication: unique destination, RUNNING then SUCCESS/exit 0 only after checks,
source hashes before/after, payload checksum manifest and final report. On error,
retain FAILED/exit 1 and partial artifacts for diagnosis; readers reject them.
Workers and pending chunk jobs are bounded (at most twice worker count), no
cohort-size number of open files/futures. Output collisions fail. No deletion or
overwrite of input/previous outputs. Temporary candidate directories are owned
by this invocation and removed after selection. Capacity measurements include
metadata and checksums; no saved simulation grid is needed for resampling.

`bin.reference_store.load_reference_bundle(path, workers=4)` returns
`(BlockPool, Reference, manifest)` after checksum/schema/geometry checks.
`inspect` validates and prints metadata without reconstructing SEQ/QUAL.
Raw-NPY mode maps geometry read-only, but BlockPool still computes unit masses.
The concordance CLI accepts this directory or its storage.json through
`--reference-bundle`; see [CONCORDANCE.md](CONCORDANCE.md) for policy/layout
checks and examples. Reading uses NumPy `allow_pickle=False` and safe relative paths.

Tests must include gapped/multi-block/fragment-like units, empty pools, widths
above uint32, compression/NPZ/NPY equivalence, thread-count independence, corrupted
payloads/incomplete state and source/policy preservation. NH acceptance uses
the same committed SHA and real test input, checks saved/reloaded IDs and window
counts plus independently computed dense counts. The small reference is not a
whole-sample capacity estimate; 20,000-times arithmetic is a labelled scenario.

## Production geometry policy (2026-10-21)

[GEOMETRY_FORMAT.md](GEOMETRY_FORMAT.md) fixes production defaults to packed/zlib9 and one verified candidate; explicit auto remains diagnostic. All existing formats/levels remain readable and existing plans/payloads are preserved. Default/explicit packed round-trip, selected level, candidate count and independent dense/same-ID invariants must pass before deployment.

## Required payload I/O recovery (2026-10-27)

Builders enumerate mandatory payloads from storage declarations and apply bounded transient filesystem retries. Strict integrity failures and old frozen plans remain protected. See [PAYLOAD_RETRY.md](PAYLOAD_RETRY.md) for the exact errno, delay, logging and acceptance contract.
