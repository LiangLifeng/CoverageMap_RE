"""Opt-in native Nextflow wrapper/cache tests, with a synthetic scorer.

Set COVERAGEMAP_TEST_NEXTFLOW to the Nextflow executable. These tests exercise
the real process wrappers, not BAM/CRAM decoding or statistical calibration.
"""
import csv
import json
import os
from pathlib import Path
import subprocess
import sys

import pytest

from bin.coverage_query_flow import input_stamp
from bin.query_manifest import table


NEXTFLOW = os.environ.get('COVERAGEMAP_TEST_NEXTFLOW')
pytestmark = pytest.mark.skipif(not NEXTFLOW, reason='Set COVERAGEMAP_TEST_NEXTFLOW for native tests')
ROOT = Path(__file__).resolve().parents[1]
METHODS = 'normal,empirical,minimum,bootstrap'


BACKEND = '''#!/usr/bin/env python3
import argparse, json, os, signal, sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).absolute().parents[1]))
from bin.query_manifest import DISTRIBUTION_COLUMNS,table
from bin.query_decision_methods import COLUMNS
parser=argparse.ArgumentParser()
parser.add_argument('--version',action='version',version='synthetic-test-scorer')
for name in ('query-id','outdir','jaccard-methods'):
    parser.add_argument('--'+name,required=True)
parser.add_argument('--reference-draws',type=int,required=True)
args,_=parser.parse_known_args()
fault=Path(os.environ['COVERAGEMAP_TEST_FAULT'])
mode=fault.read_text() if fault.exists() else ''
if mode and mode!='wrong-draws' and (args.query_id=='q00000002' or mode.startswith('all:')):
    print('INJECTED_QUERY_FAILURE',file=sys.stderr)
    if mode=='missing-summary':
        sys.exit(0)
    code=int(mode.removeprefix('all:'))
    if code<0:
        os.kill(os.getpid(),-code)
    sys.exit(code)
root=Path(args.outdir);root.mkdir()
row=dict.fromkeys(DISTRIBUTION_COLUMNS,'')
row.update(query_id=args.query_id,Sample=args.query_id,GCFID='g',manifest_depth='1',
           status='COMPLETE',detection_status='DETECTED',
           requested_draws=100 if mode=='wrong-draws' else args.reference_draws,
           valid_draws=100 if mode=='wrong-draws' else args.reference_draws)
table(root/'summary.tsv',DISTRIBUTION_COLUMNS,[row])
table(root/'sensitivity.tsv',DISTRIBUTION_COLUMNS,[row])
table(root/'decision_methods.tsv',COLUMNS,[dict.fromkeys(COLUMNS,'')|dict(
    query_id=args.query_id,Sample=args.query_id,GCFID='g',method=method,
    status='COMPLETE',detection_status='DETECTED') for method in args.jaccard_methods.split(',')])
'''


def fixture(root, *, include_backend=True):
    source = root / 'input.bam'
    source.write_bytes(b'synthetic')
    manifest = root / 'manifest.tsv'
    rows = [dict(Sample=f'q{i:08d}', GCFID='g', Depth='1', CramPath=str(source)) for i in (1, 2)]
    table(manifest, ('Sample', 'GCFID', 'Depth', 'CramPath'), rows)
    roster = root / 'roster.json'
    roster.write_text(json.dumps([dict(query_id=row['Sample'], Sample=row['Sample'],
                                      GCFID='g', manifest_depth='1') for row in rows]))
    index = root / 'index.sqlite'
    index.write_bytes(b'synthetic index, never opened by the fixture scorer')
    backend = root / 'backend'
    backend.mkdir()
    script = backend / 'query_manifest_v2.py'
    script.write_text(BACKEND)
    script.chmod(0o755)
    code = sorted([p for p in (ROOT / 'bin').glob('*.py') if p.name != script.name]
                  + ([script] if include_backend else []))
    metas = [dict(id=row['Sample'], input_stamp=input_stamp(row)) for row in rows]
    main = root / 'main.nf'
    main.write_text(f"""include {{ QUERY_SCORE }} from '{ROOT}/modules/local/query_score/main'
include {{ QUERY_COLLECT }} from '{ROOT}/modules/local/query_collect/main'
workflow {{
    code = channel.value([{','.join('file('+json.dumps(str(p))+')' for p in code)}])
    tools = channel.value(file('{ROOT}/Tools.py'))
    queries = channel.of({','.join('[id:'+json.dumps(m['id'])+',input_stamp:'+json.dumps(m['input_stamp'])+']' for m in metas)})
    QUERY_SCORE(queries,channel.value(file('{manifest}')),channel.value(file('{index}')),code,tools)
    receipts = QUERY_SCORE.out.receipts.toSortedList {{ a,b -> a[0] <=> b[0] }}
        .map {{ rows -> rows.collect {{ row -> row[1] }} }}
    QUERY_COLLECT(channel.value(file('{roster}')),receipts,code,tools)
}}
""")
    config = root / 'local.config'
    config.write_text(f"""conda.enabled=false
process.memory=null
process.time=null
env.PATH='{Path(sys.executable).parent}:'+System.getenv('PATH')
env.COVERAGEMAP_TEST_FAULT='{root}/fault'
""")
    return main, config


def run(root, number, *, resume=False):
    command = [NEXTFLOW, '-log', str(root / f'nextflow-{number}.log'),
               '-C', f'{ROOT / "nextflow.config"},{root / "local.config"}',
               'run', str(root / 'main.nf'), '-profile', 'local,nhb',
               '-work-dir', str(root / 'work'), '-with-trace', str(root / f'trace-{number}.tsv'),
               '--outdir', str(root / 'results'), '--window_size', '1', '--alpha', '.05',
               '--query_cpus', '2', '--diagnostics', 'false', '--jaccard_methods', METHODS]
    if resume:
        command.append('-resume')
    env = dict(os.environ, NXF_SYNTAX_PARSER='v2', NXF_OFFLINE='true', PYTHONDONTWRITEBYTECODE='1')
    result = subprocess.run(command, cwd=root, env=env, capture_output=True, text=True, timeout=120)
    (root / f'command-{number}.json').write_text(json.dumps(command, indent=2))
    (root / f'console-{number}.log').write_text(result.stdout + result.stderr)
    with (root / f'trace-{number}.tsv').open() as handle:
        trace = list(csv.DictReader(handle, delimiter='\t'))
    return result, trace


def scores(trace, query_id):
    return [row for row in trace if row['name'].endswith(f'({query_id})')]


def test_native_1000_draws_and_stale_backend_rejection(tmp_path):
    fixture(tmp_path)
    with (tmp_path / 'local.config').open('a') as handle:
        handle.write('\nparams.reference_draws=1000\n')
    (tmp_path / 'fault').write_text('wrong-draws')
    failed, trace = run(tmp_path, 1)
    assert failed.returncode != 0
    assert all(row['exit']=='2' for row in trace if row['name'].startswith('QUERY_SCORE'))
    receipts = list((tmp_path / 'work').rglob('q00000001.json'))
    assert len(receipts)==1
    receipt=json.loads(receipts[0].read_text())
    assert receipt['status']=='ERROR' and receipt['detection_status']=='NOT_EVALUATED'
    assert 'requested 1000, backend reported' in receipt['error']
    (tmp_path / 'fault').unlink()
    passed, trace = run(tmp_path, 2, resume=True)
    assert passed.returncode == 0, passed.stdout+passed.stderr
    assert all(row['status']=='COMPLETED' for row in trace if row['name'].startswith('QUERY_SCORE'))
    summary=list(csv.DictReader((tmp_path/'results/collected/summary.tsv').open(),delimiter='\t'))
    assert len(summary)==2 and all(row['requested_draws']==row['valid_draws']=='1000' for row in summary)


def test_native_failure_resume_and_stable_collector(tmp_path):
    fixture(tmp_path)
    fault = tmp_path / 'fault'
    fault.write_text('1')
    first, trace = run(tmp_path, 1)
    assert first.returncode != 0
    assert [r['status'] for r in scores(trace, 'q00000001')] == ['COMPLETED']
    assert len(scores(trace, 'q00000002')) == 1  # deterministic error is not retried
    failed = scores(trace, 'q00000002')[0]
    assert failed['status'] == 'FAILED' and failed['exit'] == '1'
    state = json.loads((tmp_path / 'results/collected/completion.json').read_text())
    assert state['completion_status'] == 'PARTIAL_FAILURE'
    assert state['DETECTED'] == state['ERROR'] == 1
    task = next(p.parent for p in (tmp_path / 'work').rglob('q00000002.json'))
    assert (task / '.exitcode').read_text().strip() == '1'
    assert 'INJECTED_QUERY_FAILURE' in (task / '.command.err').read_text()
    assert '--jaccard-methods "normal,empirical,minimum,bootstrap"' in (task / '.command.sh').read_text()

    fault.unlink()  # same source, CLI, declared inputs, config, launch and work dirs
    second, trace = run(tmp_path, 2, resume=True)
    assert second.returncode == 0, second.stdout + second.stderr
    assert [r['status'] for r in scores(trace, 'q00000001')] == ['CACHED']
    assert [r['status'] for r in scores(trace, 'q00000002')] == ['COMPLETED']
    before = (tmp_path / 'results/collected/decision_methods.tsv').read_bytes()
    third, trace = run(tmp_path, 3, resume=True)
    assert third.returncode == 0, third.stdout + third.stderr
    assert all(r['status'] == 'CACHED' for r in trace)
    assert before == (tmp_path / 'results/collected/decision_methods.tsv').read_bytes()


@pytest.mark.parametrize('signal_code,exit_code', [(-9, 137), (-15, 143)])
def test_native_resource_failure_has_one_retry_and_keeps_healthy_query(tmp_path, signal_code, exit_code):
    fixture(tmp_path)
    (tmp_path / 'fault').write_text(str(signal_code))
    result, trace = run(tmp_path, 1)
    assert result.returncode != 0
    assert [r['status'] for r in scores(trace, 'q00000001')] == ['COMPLETED']
    failed = scores(trace, 'q00000002')
    assert len(failed) == 2 and all(r['exit'] == str(exit_code) for r in failed)
    state = json.loads((tmp_path / 'results/collected/completion.json').read_text())
    assert state['ERROR'] == 1 and state['DETECTED'] == 1


@pytest.mark.parametrize('mode', ['all:1', 'missing-summary'])
def test_native_empty_receipts_and_missing_output_are_explicit_errors(tmp_path, mode):
    fixture(tmp_path)
    (tmp_path / 'fault').write_text(mode)
    result, trace = run(tmp_path, 1)
    assert result.returncode != 0
    state = json.loads((tmp_path / 'results/collected/completion.json').read_text())
    assert state['ERROR'] == (2 if mode == 'all:1' else 1)
    assert state['completion_status'] == 'PARTIAL_FAILURE'
    assert any(row['name'].startswith('QUERY_COLLECT') and row['status'] == 'COMPLETED' for row in trace)


def test_native_missing_scorer_executable_is_not_cached_as_success(tmp_path):
    fixture(tmp_path, include_backend=False)
    result, trace = run(tmp_path, 1)
    assert result.returncode != 0
    assert all(row['status'] == 'FAILED' and row['exit'] == '2'
               for row in trace if row['name'].startswith('QUERY_SCORE'))
    errors = [p.read_text() for p in (tmp_path / 'work').rglob('.command.err')]
    assert sum('No such file or directory' in text and 'query_manifest_v2.py' in text for text in errors) == 2
    state = json.loads((tmp_path / 'results/collected/completion.json').read_text())
    assert state['ERROR'] == state['missing_receipts'] == 2
    assert (tmp_path / 'backend/query_manifest_v2.py').is_file()


@pytest.mark.parametrize('failure', ['exit', 'signal', 'missing-executable'])
def test_native_diagnostics_reports_original_failure(tmp_path, failure):
    fixture(tmp_path)
    collected = tmp_path / 'collected'
    collected.mkdir()
    (collected / 'summary.tsv').write_text('unused synthetic summary\n')
    nulls = tmp_path / 'null.tsv.gz'
    nulls.write_bytes(b'unused synthetic nulls')
    script = tmp_path / ('missing_plotter.py' if failure == 'missing-executable' else 'query_similarity_plot.py')
    script.write_text('''#!/usr/bin/env python3
import os,signal,sys
if '--version' in sys.argv:
    print('synthetic-test-plotter');sys.exit(0)
print('INJECTED_PLOT_FAILURE',file=sys.stderr)
''' + ('os.kill(os.getpid(),signal.SIGTERM)\n' if failure == 'signal' else 'sys.exit(7)\n'))
    script.chmod(0o755)
    (tmp_path / 'main.nf').write_text(f"""include {{ QUERY_DIAGNOSTICS }} from '{ROOT}/modules/local/query_diagnostics/main'
workflow {{
    QUERY_DIAGNOSTICS(channel.value(file('{collected}')),
        channel.value(tuple(['q00000001'],[file('{nulls}')])),channel.value(file('{script}')))
}}
""")
    result, trace = run(tmp_path, 1)
    expected = '143' if failure == 'signal' else '127' if failure == 'missing-executable' else '7'
    assert result.returncode != 0
    assert len(trace) == 1 and trace[0]['status'] == 'FAILED' and trace[0]['exit'] == expected
    task = next((tmp_path / 'work').rglob('plot_exit_code.txt'))
    assert task.read_text().strip() == expected
    assert ('command not found' if failure == 'missing-executable' else 'INJECTED_PLOT_FAILURE') in (
        task.parent / 'plot_stderr.log').read_text()


@pytest.mark.parametrize('option,value,message', [
    ('jaccard_methods', 'normal,normal', 'distinct names'),
    ('bootstrap_replicates', '9', '10..10000'),
    ('bootstrap_confidence', '1', 'between 0 and 1'),
    ('diagnostics', 'maybe', 'must be true or false'),
    ('reference_draws', '1', 'integer at least 2'),
    ('reference_draws', '2.5', 'integer at least 2')])
def test_native_entry_rejects_invalid_options_before_planning(tmp_path, option, value, message):
    fixture(tmp_path)
    command = [NEXTFLOW, '-log', str(tmp_path / 'preflight.log'),
               '-C', f'{ROOT / "nextflow.config"},{tmp_path / "local.config"}',
               'run', str(ROOT / 'main.nf'), '-profile', 'local,nhb',
               '--manifest', str(tmp_path / 'manifest.tsv'), '--index', str(tmp_path / 'index.sqlite'),
               '--outdir', str(tmp_path / 'results'), '--window_size', '1', '--alpha', '.05',
               '--' + option, value]
    result = subprocess.run(command, cwd=tmp_path,
                            env=dict(os.environ, NXF_SYNTAX_PARSER='v2', NXF_OFFLINE='true'),
                            capture_output=True, text=True, timeout=60)
    (tmp_path / 'preflight.console.log').write_text(result.stdout + result.stderr)
    assert result.returncode != 0
    assert message in result.stdout + result.stderr
    assert 'Submitted process' not in result.stdout + result.stderr
