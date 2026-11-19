#!/usr/bin/env python3
"""Index a BuildInput CSV without opening its alignment files."""
import argparse
import csv
import json
import math
from pathlib import Path
import sqlite3

__version__ = "1.2.0"
REQUIRED = {"Sample", "CramFilePath", "GCFID", "Coverage", "Depth"}


def quoted(name):
    return '"' + name.replace('"', '""') + '"'


def load_catalog(path, path_column):
    """Validate an explicit ID/path mapping without touching referenced files."""
    mapping = {}
    rows = 0
    with Path(path).open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames != ["GCF_ID", path_column]:
            raise ValueError(f"{path}: expected columns GCF_ID,{path_column}")
        for rows, row in enumerate(reader, 1):
            if None in row or any(v is None for v in row.values()):
                raise ValueError(f"{path} record {rows}: wrong field count")
            gcf, filename = row["GCF_ID"], row[path_column]
            if any(not v.strip() or v != v.strip() for v in (gcf, filename)):
                raise ValueError(f"{path} record {rows}: empty or padded ID/path")
            suffixes = (".fna",) if path_column == "GCF_path" else (".fai", ".fna.fai")
            if not Path(filename).is_absolute() or Path(filename).name not in [gcf + x for x in suffixes]:
                raise ValueError(f"{path} record {rows}: absolute path basename must match full ID {gcf}")
            if gcf in mapping and mapping[gcf] != filename:
                raise ValueError(f"{path} record {rows}: conflicting paths for {gcf}")
            mapping[gcf] = filename
    if not rows:
        raise ValueError(f"{path}: no data records")
    return mapping, rows


def add_catalog(connection, fna, fai):
    fna_map, fna_rows = load_catalog(fna, "GCF_path")
    fai_map, fai_rows = load_catalog(fai, "FAI_path")
    connection.execute("""CREATE TABLE reference_genomes (
        GCFID TEXT PRIMARY KEY, FnaFilePath TEXT, FaiFilePath TEXT)""")
    connection.executemany("INSERT INTO reference_genomes VALUES (?,?,?)",
                          ((g, fna_map.get(g), fai_map.get(g)) for g in sorted(fna_map.keys() | fai_map.keys())))
    connection.execute("""CREATE VIEW gcf_catalog_audit AS
        SELECT g.*, COALESCE(r.cram_records,0) AS cram_records,
               CASE WHEN g.FnaFilePath IS NULL THEN 'MISSING_FNA'
                    WHEN g.FaiFilePath IS NULL THEN 'MISSING_FAI' ELSE 'MATCHED' END AS status
        FROM reference_genomes g LEFT JOIN
          (SELECT GCFID, COUNT(*) AS cram_records FROM reference_records GROUP BY GCFID) r
        ON g.GCFID=r.GCFID""")
    missing = connection.execute("""SELECT r.GCFID,COUNT(*) FROM reference_records r
        LEFT JOIN reference_genomes g ON r.GCFID=g.GCFID
        WHERE g.FnaFilePath IS NULL OR g.FaiFilePath IS NULL
        GROUP BY r.GCFID ORDER BY r.GCFID""").fetchall()
    if missing:
        raise ValueError(f"Missing FNA/FAI mapping for {len(missing)} GCFs / "
                         f"{sum(n for _, n in missing)} CRAM records: {missing[:20]}")
    connection.execute("""CREATE VIEW reference_cram_catalog AS
        SELECT r.*,g.FnaFilePath,g.FaiFilePath FROM reference_records r
        JOIN reference_genomes g ON r.GCFID=g.GCFID""")
    return dict(catalog_schema="coverage-reference-catalog-1",
                fna_csv=str(Path(fna).resolve()), fai_csv=str(Path(fai).resolve()),
                fna_rows=fna_rows, fai_rows=fai_rows,
                catalog_gcfs=len(fna_map.keys() | fai_map.keys()),
                missing_cram_records=0, all_crams_catalog_matched=True,
                reference_files_checked=False,
                matching_basis="exact full GCF_ID/GCFID including version; path basename checked")


def import_mate_audit(connection, path):
    """Store provenance separately; unidentifiable sources are not references."""
    connection.execute('''CREATE TABLE reference_mate_audit (
        id INTEGER PRIMARY KEY, Sample TEXT NOT NULL, GCFID TEXT NOT NULL,
        CramFilePath TEXT NOT NULL, status TEXT NOT NULL, included INTEGER NOT NULL,
        read1_seen INTEGER NOT NULL, read2_seen INTEGER NOT NULL,
        unlabelled_seen INTEGER NOT NULL, invalid_flags_seen INTEGER NOT NULL,
        scan_complete INTEGER NOT NULL)''')
    count = 0
    boolean = {'True': 1, 'False': 0}
    with Path(path).open(newline='') as handle:
        reader = csv.DictReader(handle, delimiter='\t')
        required = {'id','Sample','GCFID','CramFilePath','status','read1_seen','read2_seen',
                    'unlabelled_seen','invalid_flags_seen','scan_complete','error'}
        if not required <= set(reader.fieldnames or []):
            raise ValueError('Incomplete mate-audit columns')
        batch = []
        for row in reader:
            if row['status'] == 'ERROR' or row['error']:
                raise ValueError('Mate audit has input errors; resolve them before filtering')
            try:
                flags = [boolean[row[k]] for k in ('read1_seen','read2_seen','unlabelled_seen',
                                                   'invalid_flags_seen','scan_complete')]
            except KeyError as error:
                raise ValueError('Invalid boolean in mate audit') from error
            status = row['status']
            if status not in ('BOTH_LABELS_PRESENT','ONLY_READ1','ONLY_READ2',
                              'NO_MATE_LABELS','INVALID_MATE_FLAGS','NO_ELIGIBLE_PRIMARY'):
                raise ValueError('Unknown mate-audit status')
            included = int(bool(flags[0] or flags[1]) and not flags[3])
            if not included and not flags[4] and not flags[3]:
                raise ValueError('Excluded source requires a complete mate scan')
            count += 1
            batch.append((int(row['id']),row['Sample'],row['GCFID'],row['CramFilePath'],
                          status,included,*flags))
            if len(batch) == 10000:
                connection.executemany('INSERT INTO reference_mate_audit VALUES (?,?,?,?,?,?,?,?,?,?,?)',batch)
                batch.clear()
        if batch:
            connection.executemany('INSERT INTO reference_mate_audit VALUES (?,?,?,?,?,?,?,?,?,?,?)',batch)
    bounds = connection.execute('SELECT MIN(id),MAX(id) FROM reference_mate_audit').fetchone()
    if not count or bounds != (1,count):
        raise ValueError('Mate audit must cover all original contiguous input row IDs')
    return count


def build_index(records, output, fna=None, fai=None, mate_audit=None):
    if (fna is None) != (fai is None):
        raise ValueError("Supply both --fna and --fai, or neither for legacy CRAM-only mode")
    records, output = Path(records), Path(output)
    with records.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        columns = reader.fieldnames
        if (not columns or len(columns) != len(set(columns))
                or any(not c.strip() for c in columns)
                or len({c.casefold() for c in columns}) != len(columns)
                or "id" in {c.casefold() for c in columns}
                or not REQUIRED <= set(columns)):
            raise ValueError("CSV requires unique columns, no id column, and " + ", ".join(sorted(REQUIRED)))
        # Exclusive creation protects existing files, including an input/output collision.
        with output.open("xb"):
            pass
        connection = None
        try:
            connection = sqlite3.connect(output)
            connection.execute("BEGIN")
            definition = ", ".join(quoted(c) + " TEXT NOT NULL" for c in columns)
            connection.execute("CREATE TABLE reference_records (id INTEGER PRIMARY KEY, " + definition + ")")
            connection.execute("CREATE TABLE metadata (key TEXT PRIMARY KEY, value TEXT NOT NULL)")
            audit_count = import_mate_audit(connection, mate_audit) if mate_audit else None
            insert = ("INSERT INTO reference_records VALUES (" +
                      ",".join("?" for _ in range(len(columns) + 1)) + ")")
            batch, count, indexed = [], 0, 0
            for count, row in enumerate(reader, 1):
                if None in row or any(value is None for value in row.values()):
                    raise ValueError(f"CSV record {count}: wrong field count")
                for name in ("Sample", "GCFID", "CramFilePath"):
                    if not row[name].strip() or row[name] != row[name].strip():
                        raise ValueError(f"CSV record {count}: empty or padded {name}")
                if not Path(row["CramFilePath"]).is_absolute():
                    raise ValueError(f"CSV record {count}: CramFilePath must be absolute")
                for name in ("Depth", "Coverage"):
                    value = float(row[name])
                    if not math.isfinite(value) or value < 0:
                        raise ValueError(f"CSV record {count}: {name} must be finite and nonnegative")
                if mate_audit:
                    observed = connection.execute(
                        'SELECT Sample,GCFID,CramFilePath,included FROM reference_mate_audit WHERE id=?',
                        (count,)).fetchone()
                    if observed is None or observed[:3] != (row['Sample'],row['GCFID'],row['CramFilePath']):
                        raise ValueError(f'Mate audit identity differs from original input row {count}')
                    if not observed[3]:
                        continue
                indexed += 1
                batch.append((count, *(row[c] for c in columns)))
                if len(batch) == 10000:
                    connection.executemany(insert, batch)
                    batch.clear()
                    print(f"Indexed {count} records", flush=True)
            if not count:
                raise ValueError("CSV has no data records")
            if mate_audit and audit_count != count:
                raise ValueError('Mate audit row count differs from original CSV')
            if not indexed:
                raise ValueError('No identifiable mate-labelled references remain')
            if batch:
                connection.executemany(insert, batch)
            connection.execute("CREATE INDEX reference_gcf_sample ON reference_records(GCFID, Sample, id)")
            gcfs = connection.execute("SELECT COUNT(DISTINCT GCFID) FROM reference_records").fetchone()[0]
            metadata = dict(schema="coverage-raw-cram-index-1", version=__version__,
                            source_csv=str(records.resolve()), columns=columns,
                            records=indexed, gcfs=gcfs,
                            row_order="1-based original CSV data-record order" if mate_audit else "1-based CSV data-record order",
                            coverage_depth_basis="source CSV, not recalculated",
                            alignment_files_checked=False)
            if mate_audit:
                metadata.update(source_records=count, excluded_mate_records=count-indexed,
                                mate_audit_file=str(Path(mate_audit).resolve()),
                                mate_filter='At least one identifiable mate label; no observed invalid flags',
                                mate_presence_audited=True)
            if fna is not None:
                metadata.update(add_catalog(connection, fna, fai))
            connection.executemany("INSERT INTO metadata VALUES (?,?)",
                                   ((key, json.dumps(value, ensure_ascii=False)) for key, value in metadata.items()))
            connection.commit()
            if fna is not None:
                print(f"Catalog: {metadata['catalog_gcfs']} IDs; all {indexed} indexed CRAM records matched FNA + FAI", flush=True)
            print(f"Complete: {indexed} indexed / {count} input records, {gcfs} GCFs -> {output}", flush=True)
            return metadata
        except BaseException:
            if connection is not None:
                connection.close()
                connection = None
            output.unlink()
            raise
        finally:
            if connection is not None:
                connection.close()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--records", type=Path, required=True, help="BuildInput CSV")
    parser.add_argument("--output", type=Path, required=True, help="New SQLite file; parent must exist")
    parser.add_argument("--fna", type=Path, help="FNA mapping CSV: GCF_ID,GCF_path")
    parser.add_argument("--fai", type=Path, help="FAI mapping CSV: GCF_ID,FAI_path")
    parser.add_argument("--mate-audit", type=Path, help="Complete mate presence TSV; omit unidentifiable sources from reference_records")
    parser.add_argument("--version", action="version", version=__version__)
    args = parser.parse_args()
    try:
        build_index(args.records, args.output, args.fna, args.fai, args.mate_audit)
    except (OSError, ValueError, sqlite3.Error) as error:
        parser.exit(1, f"ERROR: {error}\n")


if __name__ == "__main__":
    main()

