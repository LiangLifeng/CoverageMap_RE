# Reference capture: fields and temporary I/O

Reference capture reads QNAME, flags, coordinates, CIGAR, mate coordinates, RG
and stored auxiliary tags. It does not consume TLEN, SEQ or QUAL. CRAM decoding
therefore requests `required_fields=6399,decode_md=0`. In HTSlib 1.24 requesting
unused TLEN can load an external reference for multi-reference slices. Removing
that request does not remove mate relationships or change aligned-base coverage.
It does not promise every possible CRAM encoding is reference-independent,
nor enable sequence reconstruction or generated MD/NM without a reference.

Only the disposable SQLite spool moves to temporary scratch. Final geometry
arrays, dictionaries, checksums, receipts and recoverable checkpoints stay in the
configured output directory. Fragment/Read/Segment/Block schema and IDs are
unchanged. Successful record checkpoint reuse still requires original full
source-content validation; failed attempts remain failed and can be retried.

Scratch directory selection, evaluated inside each worker:

1. `COVERAGEMAP_SPOOL_DIR`, when nonempty.
2. `SLURM_TMPDIR`, when nonempty.
3. Python's system temporary directory (`TMPDIR`/`TMP`/`TEMP`, normally `/tmp`).

The selected directory must exist, be absolute and permit creating a private
temporary subdirectory. Invalid settings fail; there is no silent fallback to
the output/shared filesystem. Verify the selected filesystem is node-local
before a large run; locality cannot be inferred merely from a directory name.

`COVERAGEMAP_SPOOL_MIN_FREE_MIB` is a nonnegative integer, default 1024. Free
space is checked before capture and every 10,000 decoded records. Insufficient
space raises ENOSPC with the scratch path and configuration guidance. This is
a reserve guard, not a reservation or quota check; concurrent workers and other
users can still exhaust space between checks. Size local scratch and concurrent
workers from representative spool measurements, including the largest records.

The temporary SQLite keeps its bounded 16MiB cache and disk-backed working
tables. Journaling and synchronous durability are disabled only for this private,
rebuildable spool. Persistent databases, final integrity checks, complete-state
publication and error propagation are unchanged. Python exceptions clean the
private temporary directory; a hard kill or node failure may leave disposable
scratch for the cluster's job cleanup policy. Do not delete another live job's
scratch or alter original inputs to recover.

`CAPTURE START` logs the actual scratch path. `CAPTURE SPOOL_DONE` reports decoded
records, final SQLite file bytes and spool-stage elapsed seconds. This is not
the whole build runtime or a measurement of total temporary/persistent storage.
Environment settings are execution configuration, not scientific plan identity
inputs. Existing Nextflow source tracking may rerun a wrapper when code changes;
record checkpoint validation remains independent.

Relation audits use one segment stream ordered by fragment, read role and query
interval. The existing unique(fragment,role) constraint bounds each fragment
to at most four read states; segment lists are not retained. Read/fragment flag
updates are flushed in batches of approximately10,000 reads. The persistent
tables and audit-bit definitions are unchanged; this replaces repeated per-read
SQLite lookups while preserving supplementary/unmapped/ambiguous-mate evidence.

First-observation fragment/read IDs use private-spool counters instead of a
SQL MAX lookup for every new identity. The string dictionary retains a bounded
LRU of at most4096 values and1MiB charged Python storage; large tags bypass it.
Eviction falls back to the indexed dictionary and never changes assigned IDs.

## Native Nextflow resume

Reference planning, chunk execution, collection and indexing use Nextflow's
default task cache. No reference-build selector forces `cache = false`. With
`-resume`, Nextflow decides whether a successful task can be reused; cache hits
do not invoke Python record validation. Tasks that Nextflow executes retain the
existing record checkpoints and payload checks for partial-failure recovery.
The caller manages the external input lifecycle under this cache policy.

Keep the original launch directory, Nextflow cache, work directory, output
directory and batching parameters to retain the available recovery state.
Changing source/configuration can invalidate task hashes; enabling caching
does not retroactively turn a running or previously uncached task into a hit.
