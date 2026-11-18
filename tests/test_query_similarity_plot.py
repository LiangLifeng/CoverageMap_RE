"""Verify saved statistics, output PNG metadata, and missing/error cases."""

import argparse
import csv
import gzip
import json
from pathlib import Path
from statistics import fmean, pstdev

import pytest
import numpy as np

from bin.diagnostics import query_similarity_plot as plot


def table(path, rows, compressed=False):
    opener = gzip.open if compressed else open
    with opener(path, 'wt', newline='') as handle:
        writer = csv.DictWriter(handle, list(rows[0]), delimiter='\t', lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)


def metric(row, name, values, *, observed, p, state='estimable'):
    finite = [x for x in values if x is not None]
    fields = dict(mean=observed, p=p, null_mean=fmean(finite) if finite else None,
                  null_sd=pstdev(finite) if finite else None, null_count=len(finite),
                  status=state, decision='NOT_EVALUATED' if p is None else 'PASS' if p>=.05 else 'FAIL')
    row.update({f'{name}_{key}': val for key, val in fields.items()})


def fixture(root):
    from PIL import Image
    rows, inputs = [], {}
    for i in range(1, 7):
        qid = f'q{i:08d}'
        row = dict(query_id=qid, Sample=f'Synthetic {i}', GCFID='GCF_000000001.1',
                   status='COMPLETE', detection_status='DETECTED', error='',
                   arm='low', alpha=.05, windows=60, requested_bin_count=60,
                   nominal_window_size=None, actual_query_depth=.12, manifest_depth=.13,
                   window_policy='gcf-concat-v1')
        jaccard, pearson = [.7,.8]*50, [-.3,.3]*50
        jp, pp, jm, pm, js, ps = .9, .4, .82, -.1, 'estimable', 'estimable'
        if i == 2:
            row.update(arm='high', requested_bin_count=None, nominal_window_size=100,
                       windows=28, detection_status='NOT_DETECTED')
            jp, jm = 1.24e-120, .17
        if i == 3:
            jaccard, jp, js = [1.]*100, None, 'zero_variance_null'
        if i == 4:
            pearson, pp, pm, ps = [None]*100, None, None, 'incomplete_or_nonestimable_draws'
        if i in (3,4):
            row.update(status='COMPLETE_WITH_GAPS', detection_status='NOT_EVALUATED')
        if i == 6:
            row.update(status='ERROR', detection_status='NOT_EVALUATED', error='Unreadable input',
                       manifest_depth='invalid')
            jaccard, pearson = [], []
            jp = pp = jm = pm = None
            js = ps = 'not_evaluated'
        row.update(actual_query_coverage=.25, gcf_reference_length=row['windows']*10)
        metric(row, 'jaccard', jaccard, observed=jm, p=jp, state=js)
        metric(row, 'pearson', pearson, observed=pm, p=pp, state=ps)
        rows.append(row)
        if i < 5:
            folder = root/'nulls'/qid
            folder.mkdir(parents=True)
            file = folder/'null.tsv.gz'
            table(file, [dict(metric=m, left_draw=j, right_draw='', value=v)
                         for m, values in [('jaccard',jaccard),('pearson',pearson)]
                         for j,v in enumerate(values)], compressed=True)
            np.save(folder/'window_lengths.npy', np.full(row['windows'], 10, dtype='int64'))
            np.save(folder/'reference_covered_bases.npy',
                    np.tile(np.array([2,4]*50)[:,None], (1,row['windows'])))
            table(folder/'draws.tsv.gz', [dict(query_id=qid, Sample=row['Sample'], GCFID=row['GCFID'],
                  draw_id=j, status='READY') for j in range(100)], compressed=True)
            inputs[qid] = file
    table(root/'summary.tsv', rows)
    return rows, inputs, Image


def args(tmp, **changes):
    return argparse.Namespace(**(dict(summary=tmp/'summary.tsv', outdir=tmp/'figures',
                 null_root=tmp/'nulls', work_dir=None, trace=None, staged_nulls=None,
                 staged_ids=None, query_id=None, jobs=2, stub=False) | changes))


def test_png_values_status_and_deterministic_overwrite(tmp_path):
    pytest.importorskip('matplotlib')
    rows, inputs, Image = fixture(tmp_path)
    config = args(tmp_path)
    assert plot.run(config) == 0
    index = plot.read_table(config.outdir/'plot_index.tsv', plot.PLOT_FIELDS)
    assert [x['plot_status'] for x in index] == ['PLOTTED']*2 + ['PLOTTED_WITH_GAPS']*4
    assert 'zero_variance_null' in index[2]['reason']
    assert index[3]['pearson_reference_n']=='0'
    assert 'null.tsv.gz missing' in index[4]['reason']
    assert 'Unreadable input' in index[5]['reason']
    for row, audit in zip(rows,index,strict=True):
        path = config.outdir/audit['plot_path']
        assert audit['plot_sha256']==plot.sha(path)
        with Image.open(path) as image:
            assert image.size == (1400,672)
            metadata = json.loads(image.info['Description'])
        assert metadata['summary']['query_id'] == row['query_id']
        for name,panel in metadata['panels'].items():
            assert panel['observed'] == row[name+'_mean'] and panel['p']==row[name+'_p']
        assert metadata['coverage']['query'] == .25
        if row['query_id'] in inputs:
            assert metadata['coverage']['reference_mean'] == pytest.approx(.3)
            assert metadata['coverage']['reference_min'] == .2
            assert metadata['coverage']['reference_max'] == .4
            assert metadata['coverage']['reference_n'] == 100
            assert metadata['panels']['jaccard']['values'] == plot.load_nulls(inputs[row['query_id']])['jaccard']
    assert plot.run(config) == 0
    assert plot.read_table(config.outdir/'plot_index.tsv',plot.PLOT_FIELDS) == index
    assert plot.p_label(1.24e-120) == 'p = 1.240e-120'
    assert plot.p_label(0) == 'p below float range'


def test_trace_uses_the_selected_work_task_and_pairwise_input_errors(tmp_path):
    pytest.importorskip('matplotlib')
    rows, inputs, _ = fixture(tmp_path)
    row = rows[0]
    qid = row['query_id']
    task = tmp_path/'work'/'af'/'abcded000000000000000000000000'
    qdir = task/'query_result'/'queries'/qid
    qdir.mkdir(parents=True)
    table(task/'query_result'/'summary.tsv', [row])
    for saved in inputs[qid].parent.iterdir():
        (qdir/saved.name).write_bytes(saved.read_bytes())
    trace = tmp_path/'trace.txt'
    table(trace, [dict(name=f'COVERAGE_QUERY:QUERY_SCORE ({qid})', hash='af/abcded', status='COMPLETED')])
    config = args(tmp_path, null_root=None, work_dir=tmp_path/'work', trace=trace, query_id=[qid])
    assert plot.run(config) == 0
    indexed, = plot.read_table(config.outdir/'plot_index.tsv',plot.PLOT_FIELDS)
    assert indexed['plot_status']=='PLOTTED'
    table(qdir/'null.tsv.gz', [dict(metric='jaccard', left_draw=0, right_draw=1, value=.6)], compressed=True)
    assert plot.run(config) == 1
    indexed, = plot.read_table(config.outdir/'plot_index.tsv',plot.PLOT_FIELDS)
    assert indexed['plot_status']=='ERROR' and 'pairwise' in indexed['reason'] and not indexed['plot_path']


def test_coverage_uses_base_lengths_and_ready_draws(tmp_path):
    rows, inputs, _ = fixture(tmp_path)
    row = rows[0]
    folder = inputs[row['query_id']].parent
    draws = plot.read_table(folder/'draws.tsv.gz', plot.IDENTITY, compressed=True)
    draws[0]['status'] = 'ERROR'
    table(folder/'draws.tsv.gz', draws, compressed=True)
    coverage, reasons = plot.coverage_data(row, (folder/'null.tsv.gz', None))
    assert coverage['reference_mean'] == pytest.approx((49*.2+50*.4)/99)
    assert coverage['reference_n'] == 99 and '99/100' in reasons[0]
    assert 'Query 25.00%' in plot.coverage_label(coverage)
    # Unequal window lengths: covered bases / full length, never average window fractions.
    row.update(gcf_reference_length=10, windows=2)
    np.save(folder/'window_lengths.npy', np.array([1,9]))
    np.save(folder/'reference_covered_bases.npy', np.tile([1,0], (100,1)))
    coverage, _ = plot.coverage_data(row, (folder/'null.tsv.gz', None))
    assert coverage['reference_mean'] == .1
    np.save(folder/'reference_covered_bases.npy', np.tile([2,0], (100,1)))
    with pytest.raises(ValueError, match='outside window lengths'):
        plot.coverage_data(row, (folder/'null.tsv.gz', None))
    row['actual_query_coverage'] = 25
    with pytest.raises(ValueError, match='outside its expected range'):
        plot.coverage_data(row, (folder/'null.tsv.gz', None))


def test_staged_query_directories_and_legacy_missing_coverage(tmp_path):
    rows, inputs, _ = fixture(tmp_path)
    qid = rows[0]['query_id']
    stage = tmp_path/'stage'/'input000001'
    stage.mkdir(parents=True)
    (stage/qid).symlink_to(inputs[qid].parent, target_is_directory=True)
    ids = tmp_path/'ids.txt'
    ids.write_text(qid+'\n')
    config = args(tmp_path, null_root=None, staged_nulls=tmp_path/'stage', staged_ids=ids)
    source = plot.null_sources(config, rows)[qid]
    assert source[0].resolve() == inputs[qid].resolve()
    coverage, reasons = plot.coverage_data(rows[0], source)
    assert coverage['reference_n'] == 100 and not reasons
    (inputs[qid].parent/'reference_covered_bases.npy').unlink()
    coverage, reasons = plot.coverage_data({}, source)
    assert coverage['query'] is None and coverage['reference_mean'] is None
    assert len(reasons) == 2 and 'NA' in plot.coverage_label(coverage)


@pytest.mark.parametrize('n', [2, 17, 137])
def test_plot_variable_reference_draw_count(tmp_path, n):
    rows, inputs, _ = fixture(tmp_path)
    row = rows[0]
    row['requested_draws'] = n
    folder = inputs[row['query_id']].parent
    values = np.linspace(.2, .8, n).tolist()
    for m in plot.METRICS:
        metric(row, m, values, observed=.5, p=.7)
    table(folder/'null.tsv.gz', [dict(metric=m,left_draw=i,right_draw='',value=v)
          for m in plot.METRICS for i,v in enumerate(values)], compressed=True)
    table(folder/'draws.tsv.gz', [dict(query_id=row['query_id'], Sample=row['Sample'],
          GCFID=row['GCFID'], draw_id=i,status='READY') for i in range(n)], compressed=True)
    np.save(folder/'reference_covered_bases.npy', np.full((n,row['windows']),3))
    panels, gaps, _ = plot.plot_data(row, (inputs[row['query_id']],None))
    coverage, coverage_gaps = plot.coverage_data(row, (inputs[row['query_id']],None))
    assert not gaps and not coverage_gaps
    assert len(panels['jaccard']['values']) == coverage['reference_n'] == n
    assert f'n={n}/{n}' in plot.coverage_label(coverage)
    target = tmp_path/'variable.png'
    plot.render(row,panels,[],target,coverage)
    from PIL import Image
    with Image.open(target) as image:
        metadata = json.loads(image.info['Description'])
    assert metadata['coverage']['reference_expected'] == n


def test_minimum_labels_follow_selected_method(tmp_path):
    pytest.importorskip('matplotlib')
    rows, _, image_reader = fixture(tmp_path)
    q = rows[1]
    methods = tmp_path/'decision_methods.tsv'
    common = {key: q[key] for key in plot.IDENTITY}
    common.update(jaccard_score=.17, bootstrap_lower_ci=None, bootstrap_upper_ci=None,
                  pearson_p=.4, pearson_decision='PASS', method_status='estimable')
    table(methods, [dict(common, method='normal', reference_lower_limit=.65,
                         jaccard_decision='FAIL', detection_status='NOT_DETECTED'),
                    dict(common, method='minimum', reference_lower_limit=.16,
                         jaccard_decision='PASS', detection_status='DETECTED')])
    config = args(tmp_path, methods='normal,minimum', decision_methods=methods,
                  query_id=[q['query_id']], jobs=1)
    assert plot.run(config) == 0
    normal, minimum = plot.read_table(config.outdir/'plot_index.tsv', plot.PLOT_FIELDS)
    assert (normal['method'], normal['detection_status']) == ('normal', 'NOT_DETECTED')
    assert (minimum['method'], minimum['detection_status']) == ('minimum', 'DETECTED')
    with image_reader.open(config.outdir/minimum['plot_path']) as image:
        payload = json.loads(image.info['Description'])
    assert payload['summary']['method'] == 'minimum'
    assert payload['panels']['jaccard']['decision'] == 'PASS'
    assert payload['panels']['jaccard']['p'] is None
    assert payload['panels']['jaccard']['limit'] == .16
    assert payload['panels']['pearson']['decision'] == 'PASS'
