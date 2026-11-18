"""External skipped-sample audit, with no per-sample error details in final SQLite."""
import csv
import os
from pathlib import Path
import tempfile
import json

from bin.ref_depth_catalog import digest

COLUMNS = ['database_root','gcf','sample','sample_path','error_type','error']


def has_table(connection, name):
    return connection.execute("SELECT 1 FROM sqlite_master WHERE type='table' AND name=?", (name,)).fetchone() is not None


def audit_metadata(connection):
    row = connection.execute("SELECT value FROM metadata WHERE key='skipped_samples_audit'").fetchone()
    return json.loads(row[0]) if row else None


def validate_external(connection):
    audit = audit_metadata(connection)
    if audit and audit['storage'] == 'external-tsv':
        path = Path(audit['path'])
        if digest(path) != audit['sha256']:
            raise ValueError('External skipped-sample audit checksum mismatch')
    return audit


def iter_skipped(connection):
    if has_table(connection, 'skipped_samples'):
        yield from connection.execute('SELECT * FROM skipped_samples ORDER BY database_root,gcf,sample')
        return
    audit = validate_external(connection)
    if not audit:
        return  # Older schema3 indexes without any omission audit.
    count = 0
    with Path(audit['path']).open(newline='') as handle:
        reader = csv.reader(handle, delimiter='\t')
        if next(reader, None) != COLUMNS:
            raise ValueError('Invalid external skipped-sample audit header')
        for row in reader:
            if len(row) != len(COLUMNS):
                raise ValueError('Invalid external skipped-sample audit row')
            count += 1
            yield row
    if count != audit['rows']:
        raise ValueError('External skipped-sample audit row count mismatch')
    validate_external(connection)


def prepare_external(connection, destination):
    """Stage TSV and remove detail table; caller publishes TSV before SQLite.

    The caller must VACUUM before publishing SQLite to remove freed row contents.
    No original source/index is edited.
    """
    destination = Path(destination)
    descriptor, name = tempfile.mkstemp(prefix='.skipped-samples-', dir=destination.parent)
    temporary = Path(name)
    try:
        count = 0
        with os.fdopen(descriptor, 'w', newline='') as handle:
            writer = csv.writer(handle, delimiter='\t'); writer.writerow(COLUMNS)
            for row in iter_skipped(connection):
                writer.writerow(row); count += 1
        expected = connection.execute("SELECT value FROM metadata WHERE key='skipped_samples'").fetchone()
        if expected and count != json.loads(expected[0]):
            raise ValueError('Skipped-sample audit count differs from index metadata')
        if has_table(connection, 'skipped_samples'):
            connection.execute('CREATE TABLE IF NOT EXISTS skipped_sample_counts(gcf TEXT PRIMARY KEY,copies INTEGER NOT NULL)')
            connection.execute('DELETE FROM skipped_sample_counts')
            connection.execute('INSERT INTO skipped_sample_counts SELECT gcf,COUNT(*) FROM skipped_samples GROUP BY gcf')
            connection.execute('DROP TABLE skipped_samples')
        audit = dict(storage='external-tsv', path=str(destination), sha256=digest(temporary), rows=count)
        connection.execute('INSERT OR REPLACE INTO metadata VALUES (?,?)', ('skipped_samples_audit',json.dumps(audit,sort_keys=True)))
        return temporary
    except BaseException:
        temporary.unlink(missing_ok=True)
        raise


def preflight_destination(output, skipped_tsv, roots):
    if skipped_tsv is None:
        return None
    target = Path(skipped_tsv)
    if target.is_symlink() or target.exists():
        raise FileExistsError(target)
    target = target.resolve()
    if target == output or any(target.is_relative_to(Path(root)) for root in roots):
        raise ValueError('Skipped-sample TSV must be separate from SQLite and outside source databases')
    target.parent.mkdir(parents=True, exist_ok=True)
    return target
