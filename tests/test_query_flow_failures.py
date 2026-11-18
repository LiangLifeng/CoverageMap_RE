"""Execution errors must not become successful cached scientific results."""
import csv
import json
from pathlib import Path
from types import SimpleNamespace

import pytest

from bin import coverage_query_flow as flow
from bin.query_manifest import DISTRIBUTION_COLUMNS, table


@pytest.mark.parametrize('child_exit,expected', [(0, 0), (1, 1), (2, 2), (75, 75),
                                               (137, 137), (-9, 137), (-15, 143)])
def test_score_retains_receipt_and_propagates_child_exit(tmp_path, monkeypatch, child_exit, expected):
    monkeypatch.chdir(tmp_path)
    source = tmp_path / 'input.bam'
    source.write_bytes(b'synthetic')
    row = dict(query_id='q00000001', Sample='sample', GCFID='gcf', Depth='1', CramPath=str(source))
    manifest = tmp_path / 'queries.tsv'
    table(manifest, ('Sample', 'GCFID', 'Depth', 'CramPath'),
          [{k: row[k] for k in ('Sample', 'GCFID', 'Depth', 'CramPath')}])
    result = tmp_path / 'query_result'

    def child(command, **kwargs):
        assert command[command.index('--reference-draws')+1] == '137'
        result.mkdir()
        summary = dict.fromkeys(DISTRIBUTION_COLUMNS, '')
        summary.update(query_id=row['query_id'], Sample='sample', GCFID='gcf',
                       status='COMPLETE', detection_status='DETECTED',
                       requested_draws='137', valid_draws='137')
        table(result / 'summary.tsv', DISTRIBUTION_COLUMNS, [summary])
        kwargs['stderr'].write(b'injected scorer diagnostic')
        return SimpleNamespace(returncode=child_exit)

    monkeypatch.setattr(flow.subprocess, 'run', child)
    args = flow.parser().parse_args([
        'score', '--manifest', str(manifest), '--index', 'unused.sqlite',
        '--query-id', row['query_id'], '--input-stamp', flow.input_stamp(row),
        '--query-outdir', str(result), '--receipt', str(tmp_path / 'receipt.json'),
        '--alpha', '.05', '--threads', '2', '--seed', '401', '--min-samples', '9',
        '--reference-draws', '137',
        '--depth-threshold', '5', '--cache-mib', '1024', '--work-mib', '4096', '--window-size', '1'])
    assert flow.score(args) == expected
    receipt = json.loads(args.receipt.read_text())
    assert receipt['query_exit_code'] == expected
    if expected:
        assert receipt['status'] == 'ERROR'
        assert receipt['detection_status'] == 'NOT_EVALUATED'
        assert 'injected scorer diagnostic' in receipt['error']
    else:
        assert receipt['detection_status'] == 'DETECTED'


@pytest.mark.parametrize('reported,valid,status,raises', [
    ('100', '100', 'COMPLETE', True),
    ('', '', 'COMPLETE', True),
    ('1000', '100', 'COMPLETE', True),
    ('1000', '1000', 'COMPLETE', False),
    ('1000', '999', 'COMPLETE_WITH_GAPS', False),
    ('1000', '', 'SKIPPED', False),
    ('', '', 'ERROR', False),
])
def test_reference_draw_contract(reported, valid, status, raises):
    result = dict(requested_draws=reported, valid_draws=valid, status=status)
    if raises:
        with pytest.raises(ValueError, match='mismatch'):
            flow.validate_reference_draws(result, 1000)
    else:
        flow.validate_reference_draws(result, 1000)


def test_changed_input_receipt_returns_nonzero(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    source = tmp_path / 'input.bam'
    source.write_bytes(b'synthetic')
    manifest = tmp_path / 'queries.tsv'
    table(manifest, ('Sample', 'GCFID', 'Depth', 'CramPath'),
          [dict(Sample='sample', GCFID='gcf', Depth='1', CramPath=str(source))])
    args = flow.parser().parse_args([
        'score', '--manifest', str(manifest), '--index', 'unused.sqlite',
        '--query-id', 'q00000001', '--input-stamp', 'wrong-stamp',
        '--query-outdir', str(tmp_path / 'result'), '--receipt', str(tmp_path / 'receipt.json'),
        '--alpha', '.05', '--threads', '2', '--seed', '401', '--min-samples', '9',
        '--depth-threshold', '5', '--cache-mib', '1024', '--work-mib', '4096', '--window-size', '1'])
    assert flow.score(args) == 2
    assert json.loads(args.receipt.read_text())['status'] == 'ERROR'


def test_all_failed_cohort_still_publishes_explicit_partial_ledger(tmp_path, capsys):
    roster = tmp_path / 'roster.json'
    roster.write_text(json.dumps([dict(query_id=f'q{i:08d}', Sample=f's{i}',
                                     GCFID='g', manifest_depth='1') for i in (1, 2)]))
    args = flow.parser().parse_args([
        'collect', '--roster', str(roster), '--receipts', str(tmp_path / 'no-receipts'),
        '--outdir', str(tmp_path / 'collected'), '--jaccard-methods', 'normal,empirical,minimum,bootstrap'])
    flow.collect(args)
    state = json.loads((args.outdir / 'completion.json').read_text())
    assert state['completion_status'] == 'PARTIAL_FAILURE'
    assert state['qc_status'] == 'EXECUTION_ERRORS'
    assert state['ERROR'] == state['missing_receipts'] == 2
    with (args.outdir / 'decision_methods.tsv').open() as handle:
        rows = list(csv.DictReader(handle, delimiter='\t'))
    assert len(rows) == 8
    assert all(r['method_status'] == 'query_task_error' and
               r['detection_status'] == 'NOT_EVALUATED' for r in rows)
    assert '2/2 queries failed' in capsys.readouterr().err
