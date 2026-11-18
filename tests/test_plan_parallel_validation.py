"""Planning correctness, bounded concurrency and content identity without decoding reads."""
import argparse
from collections import Counter
import csv
import hashlib
import json
import os
from pathlib import Path
import sys
import threading

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from bin import ref_build_pipeline as pipeline
from bin import ref_depth_catalog as catalog


COLUMNS = dict(sample='Sample', cram='CramFilePath', gcf='GCFID', ref_gcf='GCF_ID', fai='GCF_path')


def fixture(root, count=12):
    references = []
    for number in (1, 2):
        path = root / f'{number}.fai'
        path.write_text('chr2\t40\nchr1\t60\n')
        references.append(dict(GCF_ID=f'GCF_{number:09d}.1', GCF_path=str(path)))
    rows = []
    for number in range(count):
        extension = '.dat' if number == 1 else '.bam' if number == 2 else '.cram'
        path = root / f's{number:04d}{extension}'
        path.write_bytes((b'BAM\x01' if number == 2 else b'CRAM') + bytes([number % 250]) * 128)
        if number != 3:
            index = Path(str(path) + '.crai')
            if number == 4:
                target = root / 'index-target'
                target.write_bytes(b'index4')
                index.symlink_to(target)
            else:
                index.write_bytes(f'index{number}'.encode())
        rows.append(dict(Sample=f's{number:04d}', CramFilePath=str(path),
                         GCFID=f'GCF_{number % 2 + 1:09d}.1'))
    records, refs = root / 'records.csv', root / 'refs.csv'
    for path, fields, values in [(records, list(rows[0]), list(reversed(rows))),
                                  (refs, list(references[0]), references)]:
        with path.open('w', newline='') as handle:
            writer = csv.DictWriter(handle, fields)
            writer.writeheader()
            writer.writerows(values)
    return records, refs, rows


def plan(root, records, refs, *, cpus=1, smoke=0):
    args = argparse.Namespace(records=records, references=refs, outdir=root,
        cpus=cpus, chunk_size=3, smoke_records=smoke, window_size=20, bin_count=None)
    pipeline.plan_stage(args, eager_sources=True)
    return json.loads((root / 'database.json').read_text())


def chunks(root):
    return {p.name: p.read_bytes() for p in (root / 'chunks').glob('chunk-*.json')}


def test_full_plan_and_chunks_identical_across_cpu_counts(tmp_path):
    records, refs, rows = fixture(tmp_path)
    serial = plan(tmp_path / 'serial', records, refs)
    parallel = plan(tmp_path / 'parallel', records, refs, cpus=4)
    assert serial == parallel
    assert chunks(tmp_path / 'serial') == chunks(tmp_path / 'parallel')
    assert (tmp_path / 'serial/PLAN.sha256').read_bytes() == (tmp_path / 'parallel/PLAN.sha256').read_bytes()
    assert parallel['counts']['sample_gcf_pairs'] == len(rows)
    for group in parallel['groups'].values():
        assert [s['sample'] for s in group['samples']] == sorted(s['sample'] for s in group['samples'])
        for sample in group['samples']:
            path = Path(sample['cram'])
            assert sample['source_sha256'] == hashlib.sha256(path.read_bytes()).hexdigest()
            if sample['sample'] == 's0001':  # CRAM magic under an unknown extension.
                assert sample['cram_index']['sha256'] == hashlib.sha256(Path(str(path)+'.crai').read_bytes()).hexdigest()
            if sample['sample'] in ('s0002', 's0003'):
                assert sample['cram_index'] is None


def test_plan_cpus_reaches_actual_parallel_file_checks(tmp_path, monkeypatch):
    records, refs, _ = fixture(tmp_path)
    original = catalog.input_file_info
    barrier, lock = threading.Barrier(4), threading.Lock()
    state = dict(active=0, peak=0, entered=0)

    def inspect(value, manifest):
        if str(value).endswith('.fai'):
            return original(value, manifest)
        with lock:
            state['active'] += 1
            state['entered'] += 1
            state['peak'] = max(state['peak'], state['active'])
            first = state['entered'] <= 4
        try:
            if first:
                barrier.wait(timeout=5)
            return original(value, manifest)
        finally:
            with lock:
                state['active'] -= 1

    monkeypatch.setattr(catalog, 'input_file_info', inspect)
    plan(tmp_path / 'parallel', records, refs, cpus=4)
    assert state['peak'] == 4


def test_large_first_file_does_not_idle_the_other_worker():
    release, started = threading.Event(), threading.Event()
    state = dict(produced=0)

    def jobs():
        for number in range(100):
            state['produced'] += 1
            yield number

    def inspect(number):
        if number == 0:
            started.set()
            assert release.wait(timeout=5)
        return number

    results = catalog.bounded_plan_map(inspect, jobs(), 2)
    observed = []
    try:
        for _ in range(6):
            observed.append(next(results))
            assert state['produced'] <= len(observed) + 4
        assert started.is_set() and 0 not in observed
    finally:
        release.set()
    observed.extend(results)
    assert sorted(observed) == list(range(100))


def test_failure_does_not_submit_the_entire_manifest():
    state = dict(produced=0)

    def jobs():
        for number in range(10000):
            state['produced'] += 1
            yield number

    def inspect(number):
        if number == 0:
            raise ValueError('bad source')
        return number

    with pytest.raises(ValueError, match='bad source'):
        list(catalog.bounded_plan_map(inspect, jobs(), 2))
    assert state['produced'] <= 4


def test_same_stat_content_change_only_invalidates_its_chunk(tmp_path):
    records, refs, rows = fixture(tmp_path)
    plan(tmp_path / 'before', records, refs, cpus=4)
    before = chunks(tmp_path / 'before')
    source = Path(rows[0]['CramFilePath'])
    stat = source.stat()
    source.write_bytes(b'CRAM' + b'changed' + b'x' * (stat.st_size - 11))
    os.utime(source, ns=(stat.st_atime_ns, stat.st_mtime_ns))
    plan(tmp_path / 'after', records, refs, cpus=2)
    after = chunks(tmp_path / 'after')
    assert [name for name in before if before[name] != after[name]] == ['chunk-000000.json']


@pytest.mark.parametrize('target', ['cram', 'crai', 'fai'])
def test_source_mutation_during_validation_does_not_publish_plan(tmp_path, monkeypatch, target):
    records, refs, rows = fixture(tmp_path)
    source = Path(rows[0]['CramFilePath'])
    damaged = source if target == 'cram' else Path(str(source)+'.crai') if target == 'crai' else tmp_path/'1.fai'
    module = pipeline if target == 'cram' else catalog
    original = module.digest
    changed = False

    def digest(path):
        nonlocal changed
        value = original(path)
        if Path(path) == damaged and not changed:
            changed = True
            damaged.write_bytes(damaged.read_bytes() + b'changed')
        return value

    monkeypatch.setattr(module, 'digest', digest)
    with pytest.raises(ValueError, match='changed'):
        plan(tmp_path / 'failed', records, refs, cpus=4)
    assert not (tmp_path / 'failed').exists()


def test_smoke_still_validates_all_rows_but_hashes_only_selected(tmp_path, monkeypatch):
    records, refs, rows = fixture(tmp_path)
    original = pipeline.digest
    hashed = []

    def digest(path):
        if str(path).endswith(('.cram', '.bam', '.dat')):
            hashed.append(Path(path).name)
        return original(path)

    monkeypatch.setattr(pipeline, 'digest', digest)
    selected = plan(tmp_path / 'smoke', records, refs, cpus=4, smoke=2)
    assert sorted(hashed) == ['s0000.cram', 's0002.bam']
    assert selected['scope']['full_records'] == len(rows)
    Path(rows[-1]['CramFilePath']).unlink()
    with pytest.raises(FileNotFoundError):
        plan(tmp_path / 'bad-smoke', records, refs, cpus=4, smoke=2)
    assert not (tmp_path / 'bad-smoke').exists()


def test_known_cram_is_opened_once_for_full_content_hash(tmp_path, monkeypatch):
    records, refs, rows = fixture(tmp_path)
    original, lock, opened = Path.open, threading.Lock(), Counter()

    def open_file(path, *args, **kwargs):
        mode = args[0] if args else kwargs.get('mode', 'r')
        if mode == 'rb' and path.suffix.lower() == '.cram':
            with lock:
                opened[str(path)] += 1
        return original(path, *args, **kwargs)

    monkeypatch.setattr(Path, 'open', open_file)
    plan(tmp_path / 'full', records, refs, cpus=4)
    assert opened == Counter({row['CramFilePath']: 1 for row in rows if row['CramFilePath'].endswith('.cram')})


def test_index_changed_with_preserved_stat_during_cram_hash_is_rejected(tmp_path, monkeypatch):
    records, refs, rows = fixture(tmp_path)
    source = Path(rows[0]['CramFilePath'])
    index = Path(str(source)+'.crai')
    original = pipeline.digest

    def digest(path):
        result = original(path)
        if Path(path) == source:
            info = index.stat()
            index.write_bytes(b'X' * info.st_size)
            os.utime(index, ns=(info.st_atime_ns, info.st_mtime_ns))
        return result

    monkeypatch.setattr(pipeline, 'digest', digest)
    with pytest.raises(ValueError, match='CRAM index changed'):
        plan(tmp_path / 'failed', records, refs, cpus=4)
    assert not (tmp_path / 'failed').exists()


def test_atomic_json_streams_same_bytes_and_preserves_target_on_failure(tmp_path, monkeypatch):
    path = tmp_path / 'plan.json'
    value = dict(samples=[dict(sample=str(i), label='\u6d4b\u8bd5') for i in range(100)])
    expected = (json.dumps(value, sort_keys=True, indent=2, allow_nan=False)+'\n').encode()
    monkeypatch.setattr(catalog.json, 'dumps', lambda *args, **kwargs: pytest.fail('Whole-plan string allocation'))
    catalog.atomic_json(path, value)
    assert path.read_bytes() == expected
    with pytest.raises(ValueError):
        catalog.atomic_json(path, {'invalid': float('nan')})
    assert path.read_bytes() == expected


def test_dangling_optional_index_is_not_treated_as_absent(tmp_path):
    records, refs, rows = fixture(tmp_path)
    index = Path(rows[0]['CramFilePath'] + '.crai')
    index.unlink()
    index.symlink_to(tmp_path / 'absent')
    with pytest.raises(FileNotFoundError, match='unreadable') as error:
        plan(tmp_path / 'failed', records, refs, cpus=4)
    assert any('sample=s0000' in note for note in error.value.__notes__)
    assert not (tmp_path / 'failed').exists()


def test_duplicate_keys_fail_before_per_file_io(tmp_path, monkeypatch):
    records, refs, rows = fixture(tmp_path)
    with records.open('a') as handle:
        handle.write(f"{rows[0]['Sample']},{rows[0]['CramFilePath']},{rows[0]['GCFID']}\n")
    monkeypatch.setattr(catalog, 'input_file_info', lambda *args: pytest.fail('Unexpected per-file I/O'))
    with pytest.raises(ValueError, match='Duplicate sample/GCF'):
        plan(tmp_path / 'failed', records, refs, cpus=4)


def test_legacy_default_remains_metadata_only(tmp_path):
    records, refs, _ = fixture(tmp_path)
    value = catalog.create_plan(records, refs, tmp_path / 'legacy',
        dict(storage_mode='native', execution_model='nextflow-chunks', window_size=20), COLUMNS)
    assert all('source_sha256' not in sample for group in value['groups'].values() for sample in group['samples'])
    assert (tmp_path / 'legacy/PLAN.sha256').is_file()
