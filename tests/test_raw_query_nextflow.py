"""User-run synthetic Raw Query DAG/Diagnostics/cache/resume acceptance.

Set COVERAGEMAP_TEST_NEXTFLOW to the official Nextflow executable to opt in.
No real database or clinical manifest is used.
"""
import csv
import json
import os
from pathlib import Path
import subprocess
import sys

import pytest

from test_raw_query_vnext import _fixture

NEXTFLOW = os.environ.get('COVERAGEMAP_TEST_NEXTFLOW')
pytestmark = pytest.mark.skipif(not NEXTFLOW, reason='Set COVERAGEMAP_TEST_NEXTFLOW for native tests')
ROOT = Path(__file__).resolve().parents[1]


def test_raw_native_methods_diagnostics_resume(tmp_path):
    spec, _ = _fixture(tmp_path)
    override = tmp_path / 'local.config'
    override.write_text('conda.enabled=false\n' +
                        f"env.PATH='{Path(sys.executable).parent}:'+System.getenv('PATH')\n")
    base = [NEXTFLOW, '-C', f'{ROOT}/raw_query/nextflow.config,{override}',
            'run', str(ROOT / 'raw_query/main.nf'), '-profile', 'local',
            '-work-dir', str(tmp_path / 'work'), '--outdir', str(tmp_path / 'results'),
            '--manifest', str(tmp_path / 'queries.tsv'), '--index', str(tmp_path / 'reference.sqlite'),
            '--alpha', '.05', '--window_size', '1', '--reference_draws', '11',
            '--query_cpus', '1', '--diagnostics_cpus', '1', '--min_samples', '0',
            '--diagnostics_bin_dir', str(Path(sys.executable).parent),
            '--geometry_cache_dir', str(tmp_path / 'shared-cache'),
            '--jaccard_methods', 'normal,empirical,minimum,bootstrap',
            '--bootstrap_replicates', '100', '--jaccard_coverage_thresholds', 'default',
            '--diagnostics', 'true']
    env = dict(os.environ, NXF_SYNTAX_PARSER='v2', PYTHONDONTWRITEBYTECODE='1')
    for number, extra in enumerate([[], ['-resume'], ['--sampling_unit', 'fragment']], 1):
        trace = tmp_path / f'trace-{number}.tsv'
        command = base + ['-with-trace', str(trace)] + extra
        result = subprocess.run(command, cwd=tmp_path, env=env, text=True,
                                capture_output=True, timeout=180)
        (tmp_path / f'command-{number}.json').write_text(json.dumps(command, indent=2))
        (tmp_path / f'console-{number}.log').write_text(result.stdout + result.stderr)
        assert result.returncode == 0, result.stdout + result.stderr
        with trace.open() as handle:
            scores = [r for r in csv.DictReader(handle, delimiter='\t') if 'RAW_QUERY_SCORE' in r['name']]
        assert len(scores) == 1
        assert scores[0]['status'] == ('CACHED' if number == 2 else 'COMPLETED')
        summary = json.loads((tmp_path / 'results/queries/q00000001/summary.json').read_text())
        assert summary['requested_draws'] == summary['valid_draws'] == 11
        assert list((tmp_path / 'results/diagnostics').rglob('*.png'))
        with (tmp_path / 'results/collected/decision_methods.tsv').open() as handle:
            methods = list(csv.DictReader(handle, delimiter='\t'))
        assert len(methods) == 4


def test_failed_score_audit_survives_scratch_cleanup(tmp_path):
    import base64
    spec, _ = _fixture(tmp_path)
    spec['error'] = 'Synthetic failure for scratch audit validation'
    payload = base64.b64encode(json.dumps(spec).encode()).decode()
    scratch = tmp_path/'scratch'
    scratch.mkdir()
    main = tmp_path/'failure.nf'
    main.write_text(f"""include {{ RAW_QUERY_SCORE }} from '{ROOT}/modules/local/raw_query_score/main'
workflow {{
    queries = channel.of(tuple([id: 'q00000001'], '{payload}', []))
    RAW_QUERY_SCORE(queries, channel.value(files('{ROOT}/bin/*.py')))
}}
""")
    config = tmp_path/'failure.config'
    config.write_text(f"""process.executor='local'
conda.enabled=false
process.scratch='{scratch}'
process.errorStrategy='ignore'
workflow.failOnIgnore=true
process.shell=['/bin/bash','-euo','pipefail']
env.PATH='{Path(sys.executable).parent}:'+System.getenv('PATH')
""")
    command = [NEXTFLOW, '-C', str(config), 'run', str(main), '-work-dir', str(tmp_path/'work')]
    result = subprocess.run(command, cwd=tmp_path, text=True, capture_output=True,
                            env=dict(os.environ, NXF_SYNTAX_PARSER='v2'), timeout=120)
    (tmp_path/'failure-console.log').write_text(result.stdout+result.stderr)
    assert result.returncode != 0
    audits = list((tmp_path/'work').rglob('failure_audit/summary.json'))
    assert len(audits)==1, result.stdout+result.stderr
    assert json.loads(audits[0].read_text())['error']==spec['error']
    assert not list(scratch.rglob('summary.json'))
