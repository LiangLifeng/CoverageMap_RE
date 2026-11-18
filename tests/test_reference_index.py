"""SQLite metadata integrity, bounded reuse, exact batch/legacy concordance."""
import csv
import json
from pathlib import Path
import shutil
import sqlite3
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from bin import concordance_batch as batch
from bin import coverage_concordance as cli
from bin.reference_index import ReferenceIndex, build_index, depth_key
from bin.ref_depth_catalog import collect, digest, create_plan
from bin.ref_depth_worker import run_gcf, load_depth
from bin.concordance_core import task_rng, evaluate_windows
from bin.coverage_core import select_units, window_sums
import test_ref_depth_database as fixtures


class IndexTests(unittest.TestCase):
    def setUp(self):
        self.fixture = fixtures.DepthDatabaseTests(methodName='runTest')
        self.fixture.setUp()
        self.root = self.fixture.root
        self.database, self.plan = self.fixture.plan()
        self.assertEqual(run_gcf(self.database, fixtures.GCF1, self.fixture.extractor), 0)
        collect(self.database)
        self.path = self.root / 'references.sqlite'
        build_index(self.database, self.path, verify_payloads=True)
        self.query = self.root / 'query.cram'; self.query.write_bytes(b'distinct query')

    def tearDown(self):
        if batch.CACHE is not None:
            batch.CACHE.close(); batch.INDEX.close(); batch.CACHE = batch.INDEX = None
        self.fixture.tearDown()

    def args(self, name='query-output', **changes):
        arguments = ['coverage_concordance.py', '--reference-sqlite', str(self.path),
            '--test-alignment', str(self.query), '--query-sample', 'query', '--gcf', fixtures.GCF1,
            '--depths', '.1,1', '--outdir', str(self.root / name), '--seed', '123']
        with patch.object(sys, 'argv', arguments):
            args = cli.parse_args()
        for key, value in changes.items():
            setattr(args, key, value)
        return args

    def execute(self, args):
        with (patch.object(batch, 'extract_pool', return_value=(self.fixture.big, self.fixture.layout)) as extraction,
              patch('bin.reference_store.load_reference_bundle', side_effect=AssertionError('Must not decode reference geometry'))):
            code = batch.run_indexed(args, cli.depth_values(args))
        return code, extraction.call_count

    def summaries(self, args):
        with (args.outdir / 'summary.tsv').open() as handle:
            return list(csv.DictReader(handle, delimiter='\t'))

    def test_index_counts_read_only_exact_keys_and_exclusive_publish(self):
        index = ReferenceIndex(self.path)
        try:
            self.assertEqual(index.metadata['counts'], dict(gcfs=1, references=2, depths=6, ready=5))
            self.assertEqual(depth_key('.10'), depth_key(.1))
            self.assertEqual(len(list(index.references(fixtures.GCF1, .1))), 2)
            rows = list(index.references(fixtures.GCF1, 10))
            self.assertEqual({r['depth_status'] for r in rows}, {'READY', 'UNSUPPORTED_DEPTH'})
            plan = index.connection.execute('''EXPLAIN QUERY PLAN SELECT r.id,d.status FROM
                references_table r INDEXED BY references_lookup LEFT JOIN depths d
                ON d.reference_id=r.id AND d.depth_key=? WHERE r.gcf=?''', (depth_key(.1), fixtures.GCF1)).fetchall()
            self.assertIn('SEARCH r USING', ' '.join(r[3] for r in plan))
            with self.assertRaises(sqlite3.OperationalError):
                index.connection.execute('DELETE FROM depths')
        finally:
            index.close()
        before = digest(self.path)
        with self.assertRaises(FileExistsError):
            build_index(self.database, self.path)
        self.assertEqual(digest(self.path), before)
        summary = json.loads((self.database / 'SUMMARY.json').read_text()); summary['complete'] = False
        (self.database / 'SUMMARY.json').write_text(json.dumps(summary))
        with self.assertRaisesRegex(ValueError, 'completed'):
            build_index(self.database, self.root / 'incomplete.sqlite')

    def test_one_extraction_multiple_references_matches_core_exactly(self):
        args = self.args()
        code, calls = self.execute(args)
        self.assertEqual((code, calls), (0, 1))
        rows = self.summaries(args)
        self.assertEqual(len(rows), 4)
        for row in rows:
            observed = json.loads(Path(row['result_path']).read_text())
            depth = float(row['depth'])
            ids = select_units(self.fixture.big, task_rng(123, depth, 1, 0), target_bases=depth * 31)
            matrix, edges, _ = load_depth(row['sample_manifest'], depth)
            expected = evaluate_windows(window_sums(self.fixture.big, ids, edges), matrix)
            self.assertEqual(observed['jaccard'], expected['jaccard'])
            self.assertEqual(observed['pearson'], expected['pearson'])
            self.assertFalse(observed['reference_geometry_loaded'])
        outcome = json.loads((args.outdir / 'provenance.json').read_text())['outcomes'][0]
        report = json.loads(Path(outcome['report']).read_text())
        self.assertTrue((args.outdir/'stage_timing.tsv').is_file())
        self.assertEqual(report['test_draws'], 2)
        stages = report['stage_timing']['stages']
        self.assertEqual(stages['query_extract']['calls'],1)
        self.assertEqual(stages['query_sample_count']['calls'],2)
        self.assertEqual(stages['reference_load_decompress']['calls'],4)
        self.assertEqual(stages['reference_payload_sha256']['calls'],4)
        self.assertLessEqual(sum(v['wall_seconds'] for v in stages.values()),report['elapsed_seconds'])
        for row in stages.values():
            self.assertTrue(all(np.isfinite(row[k]) and row[k]>=0 for k in ('wall_seconds','user_cpu_seconds','system_cpu_seconds')))
        for sample in self.plan['groups'][fixtures.GCF1]['samples']:
            value = next((self.database/fixtures.GCF1/'samples'/sample['folder']).glob('attempt-*/timing.json'))
            self.assertIn('final_payload_sha256',json.loads(value.read_text())['stages'])
        self.assertEqual(report['reference_geometry_loads'], 0)

    def test_two_queries_reuse_cache_and_zero_budget_is_bounded(self):
        manifest = self.root / 'queries.csv'
        fixtures.write_csv(manifest, ['Sample', 'CramFilePath', 'GCFID'],
            [dict(Sample=name, CramFilePath=str(self.query), GCFID=fixtures.GCF1) for name in ('query-one', 'query-two')])
        for cache_mb in (10, 0):
            args = self.args(name=f'cache-{cache_mb}', query_manifest=manifest, test_alignment=None,
                             query_sample=None, gcf=None, reference_cache_mb=cache_mb)
            code, calls = self.execute(args)
            self.assertEqual((code, calls), (0, 2))
            outcomes = json.loads((args.outdir / 'provenance.json').read_text())['outcomes']
            second = json.loads(Path(outcomes[1]['report']).read_text())
            self.assertEqual(second['cache']['hits'], 4 if cache_mb else 0)
            self.assertLessEqual(second['cache_after']['accounted_bytes'], cache_mb * 1024 * 1024)

    def test_cache_rejects_mutated_payload_and_registry(self):
        index = ReferenceIndex(self.path)
        cache = batch.MatrixCache(index, 10 * 1024 * 1024)
        try:
            row = next(index.references(fixtures.GCF1, .1))
            cache.load(row, .1); cache.load(row, .1)
            self.assertEqual(cache.hits, 1)
            payload = self.database / row['matrix_path']; payload.write_bytes(b'corrupt')
            with self.assertRaisesRegex(ValueError, 'checksum'):
                cache.load(row, .1)
        finally:
            cache.close(); index.close()

    def test_relocation_does_not_need_raw_reference_or_geometry(self):
        relocated = self.root / 'relocated'
        shutil.copytree(self.database, relocated)
        shutil.rmtree(self.database)
        self.fixture.a.unlink(); self.fixture.b.unlink(); self.fixture.fai.unlink()
        for directory in relocated.glob('*/samples/*/attempt-*/geometry'):
            shutil.rmtree(directory)
        args = self.args(database_root=relocated)
        self.assertEqual(self.execute(args), (0, 1))

    def test_self_exclusion_include_self_depth_gaps_and_no_usable_reference(self):
        args = self.args(test_alignment=self.fixture.a, query_sample='a')
        self.assertEqual(self.execute(args), (0, 1))
        self.assertEqual({r['reference_sample'] for r in self.summaries(args)}, {'b'})
        args = self.args('include', test_alignment=self.fixture.a, query_sample='a', include_self=True)
        self.assertEqual(self.execute(args), (0, 1))
        self.assertIn('same-source-self-consistency', {r['comparison_scope'] for r in self.summaries(args)})
        args = self.args('gaps', depths='10')
        self.assertEqual(self.execute(args), (0, 1))
        self.assertEqual({r['status'] for r in self.summaries(args)}, {'READY', 'UNSUPPORTED_REFERENCE_DEPTH'})
        args = self.args('none', reference_samples='not-in-library')
        self.assertEqual(self.execute(args), (0, 1))
        self.assertEqual({r['status'] for r in self.summaries(args)}, {'NO_ELIGIBLE_REFERENCE'})

    def test_missing_gcf_mixed_and_all_skipped_are_explicit(self):
        unknown = 'GCF_999999999.1'
        args = self.args('absent', gcf=unknown)
        self.assertEqual(self.execute(args),(0,0))
        self.assertEqual({r['status'] for r in self.summaries(args)}, {'NO_REFERENCE_GCF'})
        report=json.loads((args.outdir/'provenance.json').read_text())
        self.assertEqual((report['state'],report['ready_query_count'],report['skipped_query_count']),('COMPLETE_WITH_GAPS',0,1))
        item=report['outcomes'][0]
        self.assertTrue((Path(item['report']).parent/'status/SKIPPED').is_file())
        self.assertFalse((Path(item['report']).parent/'status/SUCCESS').exists())
        manifest=self.root/'mixed.csv'
        fixtures.write_csv(manifest,['Sample','CramFilePath','GCFID'],
            [dict(Sample='novel-query',CramFilePath=str(self.query),GCFID=g) for g in (fixtures.GCF1,unknown)])
        args=self.args('mixed',query_manifest=manifest,test_alignment=None,query_sample=None,gcf=None)
        self.assertEqual(self.execute(args),(0,1))
        report=json.loads((args.outdir/'provenance.json').read_text())
        self.assertEqual((report['ready_query_count'],report['skipped_query_count']),(1,1))
        self.assertEqual(len(self.summaries(args)),6)

    def test_only_self_and_no_supported_depth_are_skips(self):
        args=self.args('only-self',test_alignment=self.fixture.a,query_sample='renamed-query',reference_samples='a')
        self.assertEqual(self.execute(args),(0,1))
        report=json.loads(Path(json.loads((args.outdir/'provenance.json').read_text())['outcomes'][0]['report']).read_text())
        self.assertEqual((report['state'],report['reason']),('SKIPPED','NO_ELIGIBLE_REFERENCE'))
        args=self.args('unsupported',depths='1000')
        self.assertEqual(self.execute(args),(0,1))
        report=json.loads(Path(json.loads((args.outdir/'provenance.json').read_text())['outcomes'][0]['report']).read_text())
        self.assertEqual((report['state'],report['reason']),('SKIPPED','NO_COMPARABLE_DEPTH'))
        self.assertEqual({r['status'] for r in self.summaries(args)},{'UNSUPPORTED_QUERY_DEPTH'})

    def test_duplicate_manifest_and_conflicting_policies_fail_before_output(self):
        manifest = self.root / 'duplicate.csv'
        fixtures.write_csv(manifest, ['Sample', 'CramFilePath', 'GCFID'],
            [dict(Sample='same', CramFilePath=str(self.query), GCFID=fixtures.GCF1)] * 2)
        args = self.args(query_manifest=manifest, test_alignment=None, gcf=None, query_sample=None)
        with self.assertRaisesRegex(ValueError, 'Duplicate'):
            self.execute(args)
        self.assertFalse(args.outdir.exists())
        args = self.args(min_mapq=40)
        with self.assertRaisesRegex(ValueError, 'Policy'):
            self.execute(args)
        self.assertFalse(args.outdir.exists())


try:
    import pysam
except ImportError:
    pysam = None


@unittest.skipIf(pysam is None, 'Synthetic BAM/process integration requires pysam on NH')
class BatchIntegrationTests(unittest.TestCase):
    def test_serial_process_and_legacy_exact_and_failed_query_continues(self):
        with tempfile.TemporaryDirectory(prefix='indexed-query-') as temporary:
            root = Path(temporary).resolve(); fai = root / 'layout.fai'; fai.write_text('a\t20\nb\t30\n')
            files = {}
            for number, name in enumerate(('ref-a', 'ref-b', 'ref-c', 'query-one', 'query-two')):
                bam = root / (name + '.bam'); files[name] = bam
                sq = [dict(SN='a', LN=20), dict(SN='b', LN=30)]
                if name == 'ref-b': sq = sq[:1]
                with pysam.AlignmentFile(str(bam), 'wb', header=dict(HD=dict(VN='1.6'), SQ=sq)) as handle:
                    for i in range(80):
                        read = pysam.AlignedSegment(handle.header)
                        read.query_name = f'{name}-{i}'; read.query_sequence = 'A' * 5
                        read.flag = 0; read.reference_id = i % len(sq)
                        read.reference_start = (i * (number + 1)) % (sq[read.reference_id]['LN'] - 5)
                        read.mapping_quality = 40; read.cigar = [(0, 5)]; handle.write(read)
            gcf = fixtures.GCF1
            records, refs = root / 'records.csv', root / 'refs.csv'
            fixtures.write_csv(records, ['Sample', 'CramFilePath', 'GCFID'],
                [dict(Sample=n, CramFilePath=str(files[n]), GCFID=gcf) for gcf in (fixtures.GCF1, fixtures.GCF2) for n in ('ref-a', 'ref-b', 'ref-c')])
            fixtures.write_csv(refs, ['GCF_ID', 'GCF_path'], [dict(GCF_ID=g, GCF_path=str(fai)) for g in (fixtures.GCF1, fixtures.GCF2)])
            config = dict(depths=[.1, 1., 10.], cpus=2, memory='1G', python=sys.executable,
                replicates=5, bins=10, seed=401, batch_size=2, sampler='permutation',
                geometry_format='packed', matrix_format='npz', insufficient_depth='record',
                policies=dict(mate_policy='both', sampling_unit='read', min_mapq=0, exclude_duplicates=False, exclude_qcfail=False))
            database = root / 'database'; create_plan(records, refs, database, config, fixtures.COLUMNS)
            self.assertEqual(run_gcf(database, fixtures.GCF1), 0)
            self.assertEqual(run_gcf(database, fixtures.GCF2), 0); collect(database)
            index = root / 'references.sqlite'; build_index(database, index, True)
            queries = root / 'queries.csv'
            fixtures.write_csv(queries, ['Sample', 'CramFilePath', 'GCFID'],
                [dict(Sample=n, CramFilePath=str(files[n]), GCFID=gcf) for gcf in (fixtures.GCF1, fixtures.GCF2) for n in ('query-one', 'query-two')])
            entry = Path(__file__).resolve().parents[1] / 'bin/coverage_concordance.py'
            tables = []
            for jobs, scheduling in ((1,'gcf'), (2,'gcf'), (2,'query')):
                output = root / f'jobs-{jobs}-{scheduling}'
                subprocess.run([sys.executable, '-B', str(entry), '--reference-sqlite', str(index),
                    '--query-manifest', str(queries), '--depths', '.1,1', '--seed', '123',
                    '--query-jobs', str(jobs), '--query-scheduling', scheduling, '--outdir', str(output)], check=True)
                with (output / 'summary.tsv').open() as handle: rows = list(csv.DictReader(handle, delimiter='\t'))
                self.assertEqual(len(rows), 24)
                tables.append({(r['gcf'], r['query_sample'], r['reference_sample'], r['depth']):json.loads(Path(r['result_path']).read_text()) for r in rows})
                outcomes = json.loads((output / 'provenance.json').read_text())['outcomes']
                reports = [json.loads(Path(o['report']).read_text()) for o in outcomes]
                for report in reports:
                    self.assertEqual(report['extraction_calls'],1)
                if scheduling=='gcf':
                    for g in (fixtures.GCF1, fixtures.GCF2):
                        selected = [r for r in reports if r['query']['gcf']==g]
                        self.assertEqual(selected[0]['worker_pid'],selected[1]['worker_pid'])
                        self.assertEqual(selected[1]['cache']['hits'],6)
                    self.assertEqual(json.loads((output/'provenance.json').read_text())['active_workers'], min(jobs,2))
                else:
                    self.assertEqual(json.loads((output/'provenance.json').read_text())['active_workers'],2)
            for key in tables[0]:
                for field in ('jaccard', 'pearson', 'test_selected_ids_sha256'):
                    self.assertEqual(tables[0][key][field], tables[1][key][field])
                    self.assertEqual(tables[0][key][field], tables[2][key][field])
            with (database / 'SAMPLE_INDEX.tsv').open() as handle: saved = list(csv.DictReader(handle, delimiter='\t'))
            for record in saved:
                old = root / ('legacy-' + record['gcf'] + '-' + record['sample'])
                subprocess.run([sys.executable, '-B', str(entry), '--reference-depth-manifest', record['sample_manifest'],
                    '--test-alignment', str(files['query-one']), '--depths', '.1,1', '--seed', '123', '--outdir', str(old)], check=True)
                for number, depth in enumerate((.1, 1.)):
                    expected = json.loads((old / f'depth-{number:05d}/concordance.json').read_text())
                    observed = tables[0][(record['gcf'], 'query-one', record['sample'], str(depth))]
                    for field in ('jaccard', 'pearson', 'test_selected_ids_sha256'):
                        self.assertEqual(observed[field], expected[field])
            # Unknown GCF must skip before parsing the alignment, in process mode too.
            unknown='GCF_999999999.1'
            fixtures.write_csv(queries, ['Sample','CramFilePath','GCFID'],
                [dict(Sample='new-query',CramFilePath=str(files['query-one']),GCFID=g)
                 for g in (fixtures.GCF1,unknown)])
            output=root/'missing-process'
            subprocess.run([sys.executable,'-B',str(entry),'--reference-sqlite',str(index),
                '--query-manifest',str(queries),'--depth','1','--query-jobs','2','--outdir',str(output)],check=True)
            report=json.loads((output/'provenance.json').read_text())
            self.assertEqual((report['ready_query_count'],report['skipped_query_count']),(1,1))
            skipped=json.loads(Path(report['outcomes'][1]['report']).read_text())
            self.assertEqual((skipped['reason'],skipped['extraction_calls']),('NO_REFERENCE_GCF',0))
            bad = root / 'bad.bam'; bad.write_bytes(b'not alignment')
            fixtures.write_csv(queries, ['Sample', 'CramFilePath', 'GCFID'],
                [dict(Sample='bad', CramFilePath=str(bad), GCFID=gcf),
                 dict(Sample='good', CramFilePath=str(files['query-one']), GCFID=gcf)])
            output = root / 'partial'
            result = subprocess.run([sys.executable, '-B', str(entry), '--reference-sqlite', str(index),
                '--query-manifest', str(queries), '--depth', '1', '--query-jobs', '2', '--outdir', str(output)])
            self.assertEqual(result.returncode, 1)
            outcomes = json.loads((output / 'provenance.json').read_text())['outcomes']
            self.assertEqual([r['state'] for r in outcomes], ['FAILED', 'COMPLETE'])


if __name__ == '__main__':
    unittest.main()
