# Query file-access audit

`bin/audit_query_files.py` is a standalone Python standard-library CLI. It does
not import query code or change sampling, source libraries, indexes or queries.
Run it with the same account, node, container, mounts and Python as the failing
query. A later audit cannot reconstruct an earlier process's filesystem state.

Inputs: repeatable `--query-run` (reads frozen `run.json`), `--manifest` (query
TSV, independent of filename extension), `--index` (SQLite); required new
absolute `--outdir`; `--jobs` bounded threads (default 8). At least one input is
required. `--rechecks` (default 0; enable for repeat observations), `--recheck-delay` (default .2 seconds) retain
initial failures and subsequent observations. `--hash-payloads` optionally reads
every payload fully; default probes stat plus a one-byte read, and hashes small
metadata against registered identities. No path substitution is performed.

Scope: every manifest row's CramPath and explicit CraiPath/BedGraphPath; every
SQLite reference, including unselected rows (selection annotated); indexed source
roots and database.json/SUMMARY.json; each reference's result/checksum metadata,
success/exit markers, mandatory native/geometry declarations, storage chunks and
all checksum-ledger members. Original reference CRAM paths recorded in result
sources are also inspected as provenance, labelled query_required=0: their absence
does not establish failure of a query using stored geometry. Undeclared/unrelated
files elsewhere in the library are outside this audit's scope.

Historical draw logs are streamed for all queries listed in summary.tsv; every
`Missing file: ` row is preserved and freshly probed, with query/draw/reference
identity. Missing or malformed metadata/logs are explicit incomplete-enumeration
issues; never silently report a complete audit of unknown descendants.

Outputs: `files.tsv.gz` (one row per dependency/observation, ordered by input,
reference ID, path, recheck); `issues.tsv.gz` (non-OK rows); `historical_missing.tsv.gz`
(historical failure versus current access); `environment.json` (UTC, UID/GID/groups,
host, Python, cwd, mount namespace and mountinfo; no full environment dump);
`summary.json` with counts, traversal completeness and interpretation limits.
Rows contain errno, first inaccessible component, symlinks along failed paths,
owner/mode, size, actual is_file behaviour, initial/recheck index and UTC time.
Empty SUCCESS files are valid. Metadata parse failures, unsafe stored paths,
SHA mismatch and inaccessible payloads remain distinct. Output directory is
exclusive; partial runs have no successful completion receipt.

Exit 0: complete enumeration and no required-file issues; 2: complete enumeration
with required-file issues; 1: incomplete enumeration or fatal error. Historical
failures alone do not cause exit 2 when currently accessible. Recovered initial
failures remain issues. Provenance-only failures remain visible but do not cause
exit 2. No automatic repair, permission change, retry of query calculations or
Slurm submission.

Validation: synthetic single/federated indexes, unselected references, missing
payload and manifest, malformed ledger and input, permission errno, broken and
escaping symlinks, empty markers, transient recovery, frozen-run/historical
mapping, input-byte preservation, and serial/thread output equivalence (excluding
observation timestamps). Local/NH tests must use the identical committed source.
