# Raw CRAM + reference catalog SQL index — v1.1.0

Run `build_cram_index.nf` with `--records BuildInput.csv --fna fna.csv --fai fai_path.csv --outdir NEW_OUTPUT`.
This separate Nextflow entry point imports CSV metadata only; it never opens CRAM,
CRAI or FAI files and writes no geometry. Use a new output directory.
The generated file is `NEW_OUTPUT/cram_index/reference_cram.sqlite`.
Use the same command and work directory with `-resume` to reuse a successful task.

Required columns: Sample,CramFilePath,GCFID,Coverage,Depth.
All columns and cell strings are retained. IDs are 1-based input record ordinals.
Duplicate GCF and Sample–GCF rows are retained; the table does not invent a unique
biological key. Empty input, malformed rows, duplicate/case-colliding column names,
the reserved id column, padded/empty identifiers, relative CRAM paths and invalid
numeric Depth/Coverage fail. Existing CLI output files are never overwritten.
A failed CLI build removes only the new file it created.

Schema: `coverage-raw-cram-index-1` in the metadata table (JSON values).
Table: reference_records(id INTEGER PRIMARY KEY, original CSV columns TEXT NOT NULL).
Index: reference_gcf_sample(GCFID, Sample, id).
Example parameterized lookup:

```sql
SELECT * FROM reference_records WHERE GCFID=? ORDER BY Sample,id;
```

Depth/Coverage remain source metadata, not exact values recalculated under a
sampling policy. Paths are not checked individually during indexing. These checks,
authoritative FAI/layout binding, self-source exclusions and actual reference
decoding belong to the future raw-CRAM query loader. This database deliberately has
a distinct schema from the existing geometry reference index. It cannot yet be
passed to the existing query main.nf as --index; similarity is not changed or
validated by this indexing task. The query manifest remains unchanged.

The process has one CPU, no explicit memory/time requests, a tracked script input,
and a structural stub (its empty SQLite placeholder is not a usable database).
Nextflow supplies task caching and standard pipeline_info outputs.


## Reference catalog extension (2026-11-18)

Nextflow now requires both reference CSVs as tracked file inputs, so editing
any of the three input files invalidates the build cache. The Python API/CLI
retains CRAM-only mode for legacy callers; supplying only one catalog fails.
No CRAM rows are filtered. All reference IDs from either catalog are stored,
including IDs without CRAMs. Match full GCFID/GCF_ID strings including accession
prefix and version exactly; never substitute a different version or GCA/GCF.

FNA columns must be `GCF_ID,GCF_path`; FAI columns `GCF_ID,FAI_path`.
Paths must be absolute and the filename must be ID.fna, ID.fai or ID.fna.fai
as appropriate. Identical duplicate ID/path mappings are collapsed, with source
row counts retained in metadata; conflicting paths fail. Malformed/empty catalogs
fail. Any CRAM ID missing either FNA or FAI fails the build, reporting counts and
up to 20 missing IDs. No incomplete new database is retained on failure.

Additive tables/views retain the existing `coverage-raw-cram-index-1` schema
and original `reference_records` table/index for existing raw query readers:

- `reference_genomes(GCFID PRIMARY KEY, FnaFilePath, FaiFilePath)` stores the
  union of reference IDs. NULL means that catalog lacks an unused reference ID.
- `reference_cram_catalog` joins every CRAM row with its FNA and FAI paths.
- `gcf_catalog_audit` lists every reference ID, paths, `cram_records`, and
  `MATCHED`/`MISSING_FNA`/`MISSING_FAI`. Missing both cannot occur in the union.
- JSON metadata includes `catalog_schema=coverage-reference-catalog-1`, source
  paths/counts, `missing_cram_records=0`, `all_crams_catalog_matched=true`, and
  `reference_files_checked=false`.

```sql
SELECT Sample,GCFID,CramFilePath,FnaFilePath,FaiFilePath
FROM reference_cram_catalog WHERE GCFID=? ORDER BY Sample,id;

SELECT * FROM gcf_catalog_audit ORDER BY GCFID;
```

This proves CSV association and filename identity only. No referenced file is
opened/statted: availability, CRAM header SN/LN/M5 and actual FNA/FAI sequence
consistency remain unverified. Existing query code can read the CRAM table;
automatically consuming the new reference paths in queries is a separate change.
Use a fresh outdir/workdir for this upgrade; preserve existing SQL and work caches.
