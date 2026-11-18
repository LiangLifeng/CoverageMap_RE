import csv
import json
import sqlite3

import pytest

from bin.build_cram_index import build_index


def write_input(path, rows):
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["Sample", "CramFilePath", "GCFID", "Coverage", "Depth", "Note"])
        writer.writerows(rows)


def test_roundtrip_and_index(tmp_path):
    source, output = tmp_path / "input.csv", tmp_path / "index.sqlite"
    rows = [
        ["s2", "/missing/a.cram", "GCF_1.1", "66.89", "2.640", "中文,\nnotes"],
        ["s1", "/missing/b.cram", "GCF_1.1", "0", "0", ""],
        ["s1", "/missing/b.cram", "GCF_1.1", "0", "0", "duplicate key"],
        ["s3", "/missing/c.cram", "GCF_2.1", "99", "20", ""],
    ]
    write_input(source, rows)
    before = source.read_bytes()
    assert build_index(source, output)["records"] == 4
    with sqlite3.connect(output) as connection:
        assert connection.execute("PRAGMA integrity_check").fetchone() == ("ok",)
        assert connection.execute("SELECT * FROM reference_records ORDER BY id").fetchall() == [
            (i, *row) for i, row in enumerate(rows, 1)]
        selected = connection.execute(
            "SELECT id FROM reference_records WHERE GCFID=? ORDER BY Sample,id", ("GCF_1.1",)).fetchall()
        assert selected == [(2,), (3,), (1,)]
        plan = connection.execute(
            "EXPLAIN QUERY PLAN SELECT * FROM reference_records WHERE GCFID=? ORDER BY Sample,id",
            ("GCF_1.1",)).fetchall()
        assert any("USING INDEX reference_gcf_sample" in row[3] for row in plan)
        metadata = dict(connection.execute("SELECT * FROM metadata"))
        assert json.loads(metadata["alignment_files_checked"]) is False
    assert source.read_bytes() == before
    original = output.read_bytes()
    with pytest.raises(FileExistsError):
        build_index(source, output)
    assert output.read_bytes() == original


@pytest.mark.parametrize("change", [
    (0, ""), (1, "relative.cram"), (2, " padded "), (3, "-1"), (4, "nan"), (4, "bad"),
])
def test_bad_record_leaves_no_database(tmp_path, change):
    source, output = tmp_path / "input.csv", tmp_path / "index.sqlite"
    row = ["s1", "/missing/a.cram", "GCF_1", "90", "2", ""]
    row[change[0]] = change[1]
    write_input(source, [row])
    with pytest.raises(ValueError):
        build_index(source, output)
    assert not output.exists()


@pytest.mark.parametrize("content", [
    "Sample,GCFID\ns1,GCF_1\n",
    "Sample,CramFilePath,GCFID,Coverage,Depth\n",
    "Sample,CramFilePath,GCFID,Coverage,Depth\ns1,/a,GCF_1,90\n",
    "Sample,CramFilePath,GCFID,Coverage,Depth\ns1,/a,GCF_1,90,2,extra\n",
    "Sample,CramFilePath,GCFID,Coverage,Depth,sample\n",
])
def test_bad_csv(tmp_path, content):
    source, output = tmp_path / "input.csv", tmp_path / "index.sqlite"
    source.write_text(content)
    with pytest.raises(ValueError):
        build_index(source, output)
    assert not output.exists()


def write_catalogs(tmp_path, fna_rows, fai_rows):
    paths = [tmp_path / 'fna.csv', tmp_path / 'fai.csv']
    for path, column, rows in zip(paths, ['GCF_path', 'FAI_path'], [fna_rows, fai_rows]):
        with path.open('w', newline='', encoding='utf-8-sig') as handle:
            writer = csv.writer(handle)
            writer.writerow(['GCF_ID', column])
            writer.writerows(rows)
    return paths


def test_catalog_full_join_preserves_unused_and_cram_duplicates(tmp_path):
    source, output = tmp_path / 'input.csv', tmp_path / 'index.sqlite'
    rows = [['s1', '/a.cram', 'GCF_1.1', '90', '2', ''],
            ['s1', '/a.cram', 'GCF_1.1', '90', '2', 'duplicate'],
            ['s2', '/b.cram', 'GCA_2.1', '1', '0.1', '']]
    write_input(source, rows)
    fna, fai = write_catalogs(tmp_path,
        [('GCF_1.1', '/refs/GCF_1.1.fna'), ('GCF_1.1', '/refs/GCF_1.1.fna'),
         ('GCA_2.1', '/refs/GCA_2.1.fna'), ('unused', '/refs/unused.fna')],
        [('GCF_1.1', '/refs/GCF_1.1.fna.fai'), ('GCA_2.1', '/refs/GCA_2.1.fai'),
         ('only_fai', '/refs/only_fai.fai')])
    meta = build_index(source, output, fna, fai)
    assert meta['all_crams_catalog_matched'] is True
    assert meta['catalog_gcfs'] == 4
    assert meta['fna_rows'] == 4
    with sqlite3.connect(output) as con:
        assert con.execute('PRAGMA integrity_check').fetchone()[0] == 'ok'
        assert con.execute('SELECT * FROM reference_records ORDER BY id').fetchall() == [
            (i, *r) for i, r in enumerate(rows, 1)]
        assert con.execute('SELECT id,FnaFilePath,FaiFilePath FROM reference_cram_catalog ORDER BY id').fetchall() == [
            (1, '/refs/GCF_1.1.fna', '/refs/GCF_1.1.fna.fai'),
            (2, '/refs/GCF_1.1.fna', '/refs/GCF_1.1.fna.fai'),
            (3, '/refs/GCA_2.1.fna', '/refs/GCA_2.1.fai')]
        assert con.execute('SELECT GCFID,cram_records,status FROM gcf_catalog_audit ORDER BY GCFID').fetchall() == [
            ('GCA_2.1', 1, 'MATCHED'), ('GCF_1.1', 2, 'MATCHED'),
            ('only_fai', 0, 'MISSING_FNA'), ('unused', 0, 'MISSING_FAI')]
        # Existing query loader's exact schema guard still accepts this additive extension.
        assert json.loads(dict(con.execute('SELECT * FROM metadata'))['schema']) == 'coverage-raw-cram-index-1'


@pytest.mark.parametrize('case', ['missing_fna', 'missing_fai', 'version', 'conflict',
                                  'basename', 'relative', 'blank', 'padded', 'empty'])
def test_catalog_failure_no_output(tmp_path, case):
    source, output = tmp_path / 'input.csv', tmp_path / 'index.sqlite'
    write_input(source, [['s', '/a.cram', 'GCF_1.1', '90', '2', '']])
    fna_rows = [('GCF_1.1', '/refs/GCF_1.1.fna')]
    fai_rows = [('GCF_1.1', '/refs/GCF_1.1.fai')]
    if case == 'missing_fna':
        fna_rows = [('GCF_2.1', '/refs/GCF_2.1.fna')]
    elif case == 'missing_fai':
        fai_rows = [('GCF_2.1', '/refs/GCF_2.1.fai')]
    elif case == 'version':
        fna_rows = [('GCF_1.2', '/refs/GCF_1.2.fna')]
    elif case == 'conflict':
        fna_rows.append(('GCF_1.1', '/different/GCF_1.1.fna'))
    elif case == 'basename':
        fna_rows = [('GCF_1.1', '/refs/GCF_2.1.fna')]
    elif case == 'relative':
        fna_rows = [('GCF_1.1', 'GCF_1.1.fna')]
    elif case == 'blank':
        fna_rows = [('GCF_1.1', '')]
    elif case == 'padded':
        fna_rows = [('GCF_1.1 ', '/refs/GCF_1.1.fna')]
    else:
        fna_rows = []
    fna, fai = write_catalogs(tmp_path, fna_rows, fai_rows)
    with pytest.raises(ValueError):
        build_index(source, output, fna, fai)
    assert not output.exists()


@pytest.mark.parametrize('content', ['GCF_ID,wrong\nGCF_1.1,/a\n',
                                    'GCF_ID,GCF_path\nGCF_1.1\n',
                                    'GCF_ID,GCF_path\nGCF_1.1,/a,extra\n'])
def test_malformed_catalog(tmp_path, content):
    source, output = tmp_path / 'input.csv', tmp_path / 'index.sqlite'
    write_input(source, [['s', '/a.cram', 'GCF_1.1', '90', '2', '']])
    fna, fai = write_catalogs(tmp_path, [], [('GCF_1.1', '/refs/GCF_1.1.fai')])
    fna.write_text(content)
    with pytest.raises(ValueError):
        build_index(source, output, fna, fai)
    assert not output.exists()


def test_catalog_requires_both_and_protects_inputs(tmp_path):
    source = tmp_path / 'input.csv'
    write_input(source, [['s', '/a.cram', 'GCF_1.1', '90', '2', '']])
    fna, fai = write_catalogs(tmp_path, [('GCF_1.1', '/refs/GCF_1.1.fna')],
                             [('GCF_1.1', '/refs/GCF_1.1.fai')])
    before = fna.read_bytes()
    with pytest.raises(ValueError, match='both'):
        build_index(source, tmp_path / 'index.sqlite', fna)
    with pytest.raises(FileExistsError):
        build_index(source, fna, fna, fai)
    assert fna.read_bytes() == before
