"""Pure decision-method tests: same draws, no hidden CRAM reprocessing."""
import csv
import gzip
import json
from pathlib import Path
import subprocess
import sys

import numpy as np
import pytest

from bin.query_decision_methods import (
    COLUMNS, empirical_lower, evaluate, evaluate_saved, parse_methods, validate_options,
)


def scenario(values, *, score=.5, jaccard_p=.9, pearson_p=.8, groups=10, status='COMPLETE'):
    n = len(values)
    summary = dict(query_id='q00000001', Sample='q', GCFID='GCF_000000001.1',
                   requested_draws=n, valid_draws=n, status=status, arm='low',
                   jaccard=score, jaccard_p=jaccard_p, jaccard_status='estimable',
                   pearson_p=pearson_p, jaccard_null_mean=float(np.mean(values)),
                   jaccard_null_sd=float(np.std(values)))
    draws = [dict(draw_id=str(i), status='READY', reference_id=str(i % groups))
             for i in range(n)]
    nulls = [dict(metric='jaccard', left_draw=str(i), value=str(v))
             for i, v in enumerate(values)]
    return summary, draws, nulls


def test_empirical_is_fifth_order_statistic_not_interpolated_or_pairwise_minimum():
    values = [0.1, 0.2, 0.3, 0.4, 0.5] + [0.8]*95
    summary, draws, nulls = scenario(values, score=.25)
    rows = evaluate(summary, draws, nulls, methods='normal,empirical,minimum', alpha=.05)
    assert [row['method'] for row in rows] == ['normal', 'empirical', 'minimum']
    assert rows[0]['jaccard_decision'] == 'PASS'  # frozen P-value baseline
    assert rows[1]['reference_lower_limit'] == .5
    assert rows[1]['detection_status'] == 'NOT_DETECTED'
    assert rows[2]['reference_lower_limit'] == .1
    assert rows[2]['detection_status'] == 'DETECTED'
    assert empirical_lower(np.asarray(values), .01) == .1


def test_zero_variance_normal_unscored_but_minimum_and_empirical_remain_estimable():
    summary, draws, nulls = scenario([.7]*100, score=.7, jaccard_p=None)
    summary['jaccard_status'] = 'zero_variance_null'
    rows = evaluate(summary, draws, nulls, methods='normal,empirical,minimum')
    assert rows[0]['jaccard_decision'] == 'NOT_EVALUATED'
    assert rows[0]['detection_status'] == 'NOT_EVALUATED'
    assert all(row['jaccard_decision'] == 'PASS' for row in rows[1:])
    assert all(row['detection_status'] == 'DETECTED' for row in rows[1:])
    assert all(row['reference_lower_limit'] == .7 for row in rows[1:])


def test_nonestimable_null_never_drops_rows_or_makes_minimum():
    summary, draws, nulls = scenario([.8]*100, score=.8)
    nulls[7]['value'] = ''
    for row in evaluate(summary, draws, nulls, methods='normal,empirical,minimum,bootstrap'):
        assert row['jaccard_decision'] == 'NOT_EVALUATED'
        assert row['detection_status'] == 'NOT_EVALUATED'
        assert row['method_status'] == 'nonestimable_reference_null'
    summary, draws, nulls = scenario([.8]*100, score=.8)
    draws[5]['status'] = 'ERROR'
    assert all(row['method_status'] == 'incomplete_or_invalid_draws'
               for row in evaluate(summary, draws, nulls, methods='empirical,minimum'))


def test_bootstrap_is_clustered_by_distinct_source_and_is_seed_reproducible():
    values = [.5 + (i % 20)*.01 for i in range(100)]
    summary, draws, nulls = scenario(values, score=.57, groups=20)
    args = dict(methods='bootstrap', alpha=.05, bootstrap_replicates=100,
                bootstrap_confidence=.95, seed=101)
    a = evaluate(summary, draws, nulls, **args)[0]
    b = evaluate(summary, draws, nulls, **args)[0]
    assert a == b
    assert a['unique_references'] == 20
    assert a['bootstrap_lower_ci'] <= a['bootstrap_upper_ci']
    assert 0 <= a['bootstrap_pass_fraction'] <= 1
    assert a['jaccard_decision'] in ('PASS','FAIL','UNCERTAIN')
    summary, draws, nulls = scenario(values, score=.6, groups=1)
    assert evaluate(summary, draws, nulls, methods='bootstrap')[0]['method_status'] == (
        'insufficient_unique_references_for_bootstrap')
    assert evaluate(summary, draws, nulls, methods='minimum')[0]['jaccard_decision'] == 'PASS'


def test_pearson_remains_unchanged_and_unscored_controls_joint_result():
    values = [.5]*100
    summary, draws, nulls = scenario(values, score=.6, pearson_p=None)
    rows = evaluate(summary, draws, nulls, methods='empirical,minimum')
    assert all(row['jaccard_decision'] == 'PASS' for row in rows)
    assert all(row['pearson_decision'] == 'NOT_EVALUATED' for row in rows)
    assert all(row['detection_status'] == 'NOT_EVALUATED' for row in rows)


def test_invalid_cli_options_and_draw_identity_fail_closed():
    for bad in ('', 'normal,normal', 'normal,bogus', ' minimum,'):
        with pytest.raises(ValueError):
            parse_methods(bad)
    for alpha, reps, confidence in ((0,100,.95),(.05,0,.95),(.05,100,1)):
        with pytest.raises(ValueError):
            validate_options(alpha,reps,confidence)
    summary,draws,nulls=scenario([.7]*100,score=.7)
    nulls[0]['left_draw']='1'
    assert evaluate(summary, draws, nulls, methods='minimum')[0]['method_status'] == (
        'duplicate_or_missing_null_rows')


def test_saved_results_offline_cli_does_not_modify_inputs(tmp_path):
    summary, draws, nulls = scenario([.65+.001*(i%10) for i in range(100)], score=.65)
    query = tmp_path/'run/queries/q00000001'
    query.mkdir(parents=True)
    (query/'summary.json').write_text(json.dumps(summary))
    for name, rows, columns in (
        ('draws.tsv.gz', draws, ['draw_id','status','reference_id']),
        ('null.tsv.gz', nulls, ['metric','left_draw','value']),
    ):
        with gzip.open(query/name, 'wt', newline='') as handle:
            writer = csv.DictWriter(handle, fieldnames=columns, delimiter='\t')
            writer.writeheader()
            writer.writerows(rows)
    before = {p.name:p.read_bytes() for p in query.iterdir()}
    output = tmp_path/'decisions.tsv'
    run = subprocess.run(
        [sys.executable,'-B','-m','bin.query_decision_methods',
         '--run-dir',str(tmp_path/'run'),'--out',str(output),
         '--methods','normal,empirical,minimum,bootstrap','--bootstrap-replicates','20'],
        capture_output=True,text=True)
    assert run.returncode == 0, run.stderr
    with output.open() as handle:
        reader = csv.DictReader(handle,delimiter='\t')
        assert tuple(reader.fieldnames) == COLUMNS
        rows = list(reader)
    assert len(rows)==4
    assert {r['method'] for r in rows} == {'normal','empirical','minimum','bootstrap'}
    assert before == {p.name:p.read_bytes() for p in query.iterdir()}


def test_tp_fp_report_uses_total_denominators_and_checks_pairing(tmp_path):
    from bin.query_method_report import load, summarize
    fields=['query_id','Sample','GCFID','method','detection_status']
    def write(path, data):
        with path.open('w',newline='') as handle:
            writer=csv.DictWriter(handle,fields,delimiter='\t')
            writer.writeheader()
            writer.writerows([
                dict(query_id=f'q{i}', Sample=f's{i}', GCFID='GCF_000000001.1',
                     method=method,detection_status=status)
                for method, outcomes in data.items()
                for i,status in enumerate(outcomes)
            ])
    tp=tmp_path/'tp.tsv'; fp=tmp_path/'fp.tsv'
    write(tp,dict(normal=['DETECTED','DETECTED','NOT_EVALUATED'],
                  minimum=['DETECTED','DETECTED','DETECTED']))
    write(fp,dict(normal=['NOT_DETECTED','NOT_DETECTED','NOT_DETECTED'],
                  minimum=['NOT_DETECTED','DETECTED','NOT_DETECTED']))
    rows={r['method']:r for r in summarize(load(tp),load(fp))}
    assert rows['normal']['tp_total']==3 and rows['normal']['tpr_total']==pytest.approx(2/3)
    assert rows['minimum']['tpr_total']==1
    assert rows['minimum']['delta_tpr_vs_normal']==pytest.approx(1/3)
    assert rows['minimum']['fp_detected']==1
    assert rows['minimum']['fpr_total']==pytest.approx(1/3)
    assert rows['minimum']['delta_fpr_vs_normal']==pytest.approx(1/3)
    write(fp,dict(normal=['NOT_DETECTED','NOT_DETECTED','NOT_DETECTED'],
                  minimum=['NOT_DETECTED']))
    with pytest.raises(ValueError,match='rosters'):
        load(fp)


def test_nextflow_collector_rejects_wrong_method_identity_and_overrides_errors(tmp_path):
    from bin.coverage_query_flow import collect, parser
    roster = [dict(query_id='q00000001', Sample='sample-a',
                   GCFID='GCF_000000001.1', manifest_depth='1')]
    roster_path=tmp_path/'roster.json';roster_path.write_text(json.dumps(roster))
    receipts=tmp_path/'receipts';receipts.mkdir()
    entry=dict(query_id='q00000001', Sample='WRONG', GCFID='GCF_000000001.1',
               method='minimum', jaccard_decision='PASS', pearson_decision='PASS',
               detection_status='DETECTED')
    receipt=dict(roster[0], status='ERROR', error='failed task',
                 detection_status='NOT_EVALUATED', result_source='TASK_RECEIPT',
                 threshold_results=[], method_results=[entry], query_exit_code=2)
    receipt_file=receipts/'q00000001.json'
    receipt_file.write_text(json.dumps(receipt))
    args=parser().parse_args([
        'collect','--roster',str(roster_path),'--receipts',str(receipts),
        '--outdir',str(tmp_path/'collected'), '--jaccard-methods','minimum'])
    with pytest.raises(ValueError,match='method receipt identity'):
        collect(args)
    entry['Sample']='sample-a'
    receipt_file.write_text(json.dumps(receipt))
    collect(args)
    with (tmp_path/'collected/decision_methods.tsv').open() as handle:
        rows=list(csv.DictReader(handle,delimiter='\t'))
    assert len(rows)==1
    assert rows[0]['detection_status']=='NOT_EVALUATED'
    assert rows[0]['jaccard_decision']==rows[0]['pearson_decision']=='NOT_EVALUATED'
    assert rows[0]['method_status']=='query_task_error'


@pytest.mark.parametrize('n', [2, 17, 137])
def test_methods_accept_requested_reference_draws(n):
    values = np.linspace(.1, .9, n).tolist()
    summary, draws, nulls = scenario(values, score=.6)
    rows = evaluate(summary, draws, nulls, methods='normal,empirical,minimum,bootstrap')
    assert all(row['null_count'] == n for row in rows)
    assert rows[1]['reference_lower_limit'] == empirical_lower(values, .05)
    assert rows[2]['reference_lower_limit'] == .1
    summary['requested_draws'] = n+1
    rows = evaluate(summary, draws, nulls, methods='minimum')
    assert rows[0]['detection_status'] == 'NOT_EVALUATED'
