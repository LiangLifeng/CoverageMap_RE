"""Independent numerical baseline and real multi-GCF CLI serial/process replay."""
import csv
import json
import math
from pathlib import Path
from statistics import NormalDist
import subprocess
import sys
import tempfile
from types import SimpleNamespace
import typing
import unittest
from unittest.mock import patch

import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from bin.coverage_core import Reference, BlockPool, dense_coverage, select_units, window_sums
from bin.query_manifest import (windows, draw_count, source_draws, rng, table, read_manifest,
                                high_null, pair_population_size, parser, run)
from bin.query_manifest_stats import similarity, low_null, aggregate, normal_test
from bin.ref_depth_catalog import digest, create_plan, collect
from bin.ref_depth_worker import run_gcf
from bin.reference_index import build_index
from bin.native_concordance import NativeCache
import test_ref_depth_database as fixtures


def baseline():
    namespace = dict(np=np, Optional=typing.Optional, Union=typing.Union, Tuple=typing.Tuple,
                     norm=SimpleNamespace(cdf=NormalDist().cdf))
    exec((Path(__file__).parent/'data/tools_statistics_original.py').read_text(), namespace)
    return namespace


class NumericalTests(unittest.TestCase):
    def test_unlimited_cache_retains_large_entries_and_still_invalidates_stamps(self):
        cache = NativeCache(None, None)
        for number in range(3):
            cache.put(number, ('payload', number), 1024**3, ['original'])
        self.assertEqual(cache.bytes, 3*1024**3)  # Accounted sizes, no large allocation.
        self.assertEqual(cache.evictions, 0)
        self.assertIsNone(cache.stats()['budget_bytes'])
        self.assertEqual(cache.get(0, ['original']), ('payload', 0))
        self.assertIsNone(cache.get(0, ['changed']))
        self.assertEqual(cache.bytes, 2*1024**3)
        cache.close()
        self.assertEqual(cache.bytes, 0)
        self.assertFalse(cache.entries)
        bounded = NativeCache(None, 10)
        bounded.put('a', 'a', 6, [])
        bounded.put('b', 'b', 6, [])
        bounded.put('oversized', 'large', 11, [])
        self.assertEqual(list(bounded.entries), ['b'])
        self.assertEqual(bounded.evictions, 1)
        disabled = NativeCache(None, 0)  # Preserve internal zero-budget callers.
        disabled.put('a', 'a', 1, [])
        self.assertFalse(disabled.entries)

    def test_negative_memory_budgets_fail_before_reading_inputs(self):
        base = ['--manifest','missing','--index','missing','--outdir','missing','--window-size','100']
        defaults = parser().parse_args(base)
        self.assertEqual((defaults.cache_mib, defaults.work_mib), (131072, 131072))
        for flag in ('--cache-mib', '--work-mib'):
            with self.subTest(flag=flag), self.assertRaisesRegex(ValueError, 'Invalid numeric configuration'):
                run(parser().parse_args(base+[flag, '-1']))

    def test_global_gcf_width_matches_original_size_mode_and_independent_dense(self):
        original = baseline()
        layout = Reference(('one', 'two'), (20, 11), ('', ''))
        pool = BlockPool.from_intervals([[0, 7], [3, 19], [20, 23], [25, 31]], 31)
        dense = dense_coverage(pool, np.arange(len(pool)))
        for size, count in ((7, None), (None, 3), (100, None), (None, 100)):
            edges, rows = windows(layout, size, count)
            values = window_sums(pool, np.arange(len(pool)), edges)
            width = size if size is not None else math.ceil(layout.length/count)
            expected = np.concatenate([original['bin_coverage'](v, bin_count=None, bin_size=width)
                                       for v in (dense[:20], dense[20:])])
            np.testing.assert_array_equal(values, expected)
            np.testing.assert_array_equal(values, [dense[a:b].sum() for a,b in zip(edges[:-1], edges[1:])])
            self.assertEqual(sum(r['width'] for r in rows), 31)
            self.assertTrue(all(r['end'] <= (20 if r['contig'] == 'one' else 11) for r in rows))
        self.assertEqual(windows(layout, count=6)[1][0]['nominal_width'], 6)
        self.assertEqual(len(windows(layout, count=6)[1]), 6)
        with self.assertRaises(ValueError): windows(layout, 7, 3)

    def test_global_divisor_uses_full_reference_and_retains_whole_short_contigs(self):
        layout=Reference(('a','b','c'),(1800,900,300),('','',''))
        edges,rows=windows(layout,count=3)
        np.testing.assert_array_equal(edges,[0,1000,1800,2700,3000])
        self.assertEqual([r['width'] for r in rows],[1000,800,900,300])
        self.assertEqual({r['nominal_width'] for r in rows},{1000})
        self.assertEqual([r['start'] for r in rows],[0,1000,0,0])
        self.assertEqual(len(rows),4)  # B sets width, not the actual number after contig splits.
        # Missing/zero-covered contigs still contribute their full reference lengths.
        missing=Reference(layout.names,layout.lengths,layout.md5,missing_header_contigs=('b','c'))
        np.testing.assert_array_equal(windows(missing,count=3)[0],edges)
        self.assertEqual([r['width'] for r in windows(layout,count=1)[1]],[1800,900,300])

    def test_jaccard_signed_pearson_zeros_constants_and_order(self):
        x, y = np.array([0., 1, 4, 0, 2]), np.array([3., 0, 0, 0, 1])
        observed, original = similarity(x, y), baseline()
        self.assertEqual(observed['intersection'], 1)
        self.assertEqual(observed['union'], 4)
        self.assertEqual(observed['jaccard'], original['jaccard_similarity'](x, y))
        self.assertAlmostEqual(observed['pearson'], original['pearson_corr'](x, y), places=14)
        self.assertLess(observed['pearson'], 0)
        self.assertNotEqual(observed['pearson'], similarity(np.sort(x), np.sort(y))['pearson'])
        zero = similarity(np.zeros(3), np.zeros(3))
        self.assertIsNone(zero['jaccard'])
        self.assertEqual(zero['jaccard_status'], 'zero_union')
        self.assertEqual(zero['union'], 0)
        self.assertIsNone(zero['pearson'])
        self.assertIsNone(normal_test(.5, [1., 1.])['p'])
        with self.assertRaises(ValueError): similarity([1, -1], [1, 0])

    def test_original_low_formula_and_no_failed_draw_denominator_change(self):
        matrix = np.array([[0,2,4,3,0], [1,0,3,2,0], [0,2,1,0,2], [3,0,0,1,1]])
        query = np.array([0,1,3,1,1])
        original = baseline()
        nulls = low_null(matrix)
        draws = [dict(status='READY', **similarity(query, row)) for row in matrix]
        result = aggregate(draws, nulls, 'low', 4)
        for metric, function in (('jaccard', 'jaccard_similarity'), ('pearson', 'pearson_corr')):
            full = original[function](matrix)
            expected = np.array([np.delete(row, i).mean() for i, row in enumerate(full)])
            np.testing.assert_allclose(nulls[metric], expected, atol=1e-14)
            observed = original[function](query, matrix).mean()
            self.assertAlmostEqual(result[metric]['p'], original['z_test'](observed, expected), places=12)
        draws[1]['status'] = 'INSUFFICIENT_REFERENCE_DEPTH'
        self.assertIsNone(aggregate(draws, nulls, 'low', 4)['jaccard']['p'])

    def test_replacement_percent_count_and_independent_streams(self):
        self.assertEqual(draw_count('20%', 11), 3)
        self.assertEqual(draw_count('100', 11), 100)
        query = dict(query_id='q1', Sample='s', GCFID='g')
        first, token = source_draws([1,2,3], 100, 91, query)
        second, _ = source_draws([1,2,3], 100, 91, query)
        np.testing.assert_array_equal(first, second)
        self.assertEqual(len(set(first)), 3)
        self.assertEqual(len(first), 100)
        self.assertNotEqual(token, rng(91, 'read-units', 0)[1])
        for spec in ('0', '101%', 'nan%', '-1', 'garbage'):
            with self.assertRaises(ValueError): draw_count(spec, 11)

    def test_high_optional_null_excludes_same_sample_preserves_draw_weights(self):
        matrix = np.array([[0,1,4,0], [1,0,2,4], [0,1,4,0], [2,1,0,3]])
        draws = [dict(status='READY', reference_sample=s, **similarity([0,2,3,1], row))
                 for s,row in zip(('a','b','a','c'), matrix)]
        config = dict(high_null='draw-pairs', seed=9, null_limit=1000)
        nulls, records = high_null(matrix, draws, config, dict(query_id='q1'))
        self.assertEqual(len(nulls['jaccard']), 5)
        self.assertEqual(len(records), 10)
        self.assertEqual(pair_population_size(draws), 5)
        self.assertTrue(all(r['left_reference_sample'] != r['right_reference_sample'] for r in records))
        self.assertTrue(all(r['intersection'] <= r['union'] for r in records))
        self.assertEqual(sum(r['left_draw'] == 0 and r['right_draw'] == 2 for r in records), 0)
        result = aggregate(draws, nulls, 'high', 4)
        original = baseline()
        for metric, func in (('jaccard','jaccard_similarity'), ('pearson','pearson_corr')):
            value = original[func](np.array([0,2,3,1]), matrix).max()
            self.assertAlmostEqual(result[metric]['p'], original['z_test'](value, np.array(nulls[metric])))
        config['null_limit'] = 2
        self.assertEqual(high_null(matrix, draws, config, dict(query_id='q1')),
                         high_null(matrix, draws, config, dict(query_id='q1')))
        self.assertEqual(len(high_null(matrix, draws, config, dict(query_id='q1'))[0]['jaccard']), 2)

    def test_high_null_1000_pair_cap_default_and_degenerate_distributions(self):
        draws = [dict(status='READY', reference_sample=str(i % 2)) for i in range(64)]
        matrix = np.asarray([[0,1,4,0] if i % 2 else [1,0,2,4] for i in range(64)])
        config = dict(high_null='draw-pairs', seed=11, null_limit=1000)
        values, records = high_null(matrix, draws, config, dict(query_id='cap'))
        pairs = {(r['left_draw'], r['right_draw']) for r in records}
        self.assertEqual(pair_population_size(draws), 1024)
        self.assertEqual(len(pairs), 1000)
        self.assertEqual(len(values['jaccard']), 1000)
        self.assertEqual((values, records), high_null(matrix, draws, config, dict(query_id='cap')))
        self.assertTrue(all(draws[i]['reference_sample'] != draws[j]['reference_sample'] for i,j in pairs))
        self.assertIsNone(normal_test(.5, values['jaccard'])['p'])  # Zero variance is not manufactured p.
        constant = normal_test(.7, [.7]*1000)
        self.assertEqual(constant['null_sd'], 0)
        self.assertEqual(constant['status'], 'zero_variance_null')
        self.assertIsNone(constant['p'])
        for d in draws: d['reference_sample']='same'
        self.assertEqual(pair_population_size(draws), 0)
        self.assertEqual(high_null(matrix, draws, config, dict(query_id='cap')), ({'jaccard': [], 'pearson': []}, []))
        args=parser().parse_args(['--manifest','missing','--index','missing','--outdir','missing','--window-size','1000'])
        self.assertEqual(args.high_null, 'draw-pairs')
        args.null_limit=1001
        with self.assertRaises(ValueError): run(args)

    def test_complete_unit_sampling_matches_dense_no_truncation_or_reuse(self):
        pool = BlockPool.from_intervals([[0,7], [3,10], [12,17], [20,29]], 31)
        edges = windows(Reference(('a','b'), (20,11), ('','')), size=6)[0]
        for number in range(5):
            ids = select_units(pool, rng(50, number)[0], target_bases=10, depth_sampler='permutation')
            self.assertEqual(len(ids), len(set(ids)))
            dense = dense_coverage(pool, ids)
            np.testing.assert_array_equal(window_sums(pool, ids, edges), [dense[a:b].sum() for a,b in zip(edges[:-1], edges[1:])])
            self.assertGreaterEqual(int(dense.sum()), 10)
            self.assertEqual(int(dense.sum()), int(pool.unit_bases[ids].sum()))
        with self.assertRaises(ValueError): select_units(pool, rng(50)[0], target_bases=100)

    def test_manifest_validation_keeps_bad_rows_and_deterministic_compression(self):
        with tempfile.TemporaryDirectory() as temporary:
            root=Path(temporary)
            (root/'x.cram').write_bytes(b'fixture')
            rows = [dict(Sample='s', GCFID='g', Depth='1', CramPath=str(root/'x.cram')),
                    dict(Sample='bad', GCFID='g', Depth='nan', CramPath='relative')]
            table(root/'q.tsv', ('Sample','GCFID','Depth','CramPath'), rows)
            loaded = read_manifest(root/'q.tsv')
            self.assertFalse(loaded[0]['_error']); self.assertTrue(loaded[1]['_error'])
            for name in ('a.gz','b.gz'): table(root/name, ('x',), [dict(x=1)])
            self.assertEqual((root/'a.gz').read_bytes(), (root/'b.gz').read_bytes())
            table(root/'coverage.tsv', ('Sample','GCFID','Depth','CramPath','Coverage'), [dict(rows[0],Coverage='44.03')])
            self.assertTrue(read_manifest(root/'coverage.tsv')[0]['_error'])
            self.assertFalse(read_manifest(root/'coverage.tsv','percent')[0]['_error'])
            self.assertTrue(read_manifest(root/'coverage.tsv','fraction')[0]['_error'])
            (root/'bad.tsv').write_text('Sample,GCFID,Depth,CramPath\ns,g,1,x\n')
            with self.assertRaises(ValueError): read_manifest(root/'bad.tsv')
            (root/'fields.tsv').write_text('Sample\tGCFID\tDepth\tCramPath\ns\tg\t1\tx\textra\n')
            malformed = read_manifest(root/'fields.tsv')[0]
            self.assertEqual(malformed['_extra_fields'], ['extra'])
            self.assertTrue(malformed['_error'])
            json.dumps(malformed, sort_keys=True)


class ProcessTests(unittest.TestCase):
    def test_exact_parent_mass_does_not_round_target_above_available_bases(self):
        from bin.query_manifest import parser, run
        fx=fixtures.DepthDatabaseTests('runTest');fx.setUp()
        try:
            pool=BlockPool.from_intervals([[1,11]]*24+[[1,10]],31)
            self.assertEqual(pool.total_bases,249)
            self.assertGreater((249/31)*31,249)  # Concrete regression, not a tolerance choice.
            fx.small=fx.big=pool
            config=dict(fx.config,storage_mode='native',depths=[],window_size=7,bins=None)
            database=fx.root/'exact-parent'
            create_plan(fx.records,fx.refs,database,config,fixtures.COLUMNS)
            self.assertEqual(run_gcf(database,fixtures.GCF1,fx.extractor),0)
            self.assertTrue(collect(database)['complete'])
            index=fx.root/'exact.sqlite';build_index(database,index,True)
            manifest=fx.root/'q.tsv'
            table(manifest,('Sample','GCFID','Depth','CramPath'),
                  [dict(Sample='a',GCFID=fixtures.GCF1,Depth=249/31,CramPath=str(fx.a))])
            args=parser().parse_args(['--manifest',str(manifest),'--index',str(index),'--outdir',str(fx.root/'out'),
                                      '--window-size','7','--depth-threshold','10','--min-samples','0','--simulations','3'])
            with patch('bin.query_manifest.extract_pool',return_value=(pool,fx.layout)):
                self.assertEqual(run(args),0)
            row=json.loads((args.outdir/'queries/q00000001/summary.json').read_text())
            self.assertEqual(row['valid_draws'],3)
            self.assertEqual(row['status'],'COMPLETE')
            # The bounded guard rejects just below the required bytes; both exact
            # boundary and unlimited execute the same real synthetic calculation.
            from bin.query_manifest import calculate
            from bin.reference_index import ReferenceIndex
            cfg = json.loads((args.outdir/'run.json').read_text())['config']
            query = read_manifest(manifest)[0]
            required = 3*len(windows(fx.layout, size=7)[1])*8*5
            idx = ReferenceIndex(index)
            cache = NativeCache(idx, None)
            try:
                for label, budget in (('below',required-1),('boundary',required),('unlimited',None)):
                    output = fx.root/label
                    output.mkdir()
                    with patch('bin.query_manifest.extract_pool',return_value=(pool,fx.layout)):
                        if label == 'below':
                            with self.assertRaisesRegex(ValueError,'exceed work budget'):
                                calculate(query, dict(cfg,work_bytes=budget), idx, cache, output)
                        else:
                            actual, _ = calculate(query, dict(cfg,work_bytes=budget), idx, cache, output)
                            self.assertEqual(actual, row)
            finally:
                cache.close()
                idx.close()
        finally:fx.tearDown()

    def test_actual_bam_serial_process_boundary_resume_and_source_immutability(self):
        try: import pysam
        except ImportError: self.skipTest('pysam required')
        with tempfile.TemporaryDirectory() as temporary:
            root=Path(temporary).resolve()
            fai=root/'ref.fai'; fai.write_text('one\t300\ntwo\t100\n')
            def bam(name, depth, offset):
                path=root/(name+'.bam')
                with pysam.AlignmentFile(path, 'wb', header={'HD':{'VN':'1.6'}, 'SQ':[{'SN':'one','LN':300}]}) as handle:
                    for i in range(depth*20):
                        read=pysam.AlignedSegment();read.query_name=f'{name}-{i}';read.query_sequence='A'*20
                        read.flag=0;read.reference_id=0;read.reference_start=(i*offset)%280
                        read.mapping_quality=60;read.cigar=[(0,20)];handle.write(read)
                return path
            sources={d:bam('r'+str(d),d,13+d) for d in (2,6,8)}
            refs=root/'refs.csv'; records=root/'records.csv'
            gcfs=(fixtures.GCF1,fixtures.GCF2)
            fixtures.write_csv(refs, ['GCF_ID','GCF_path'], [dict(GCF_ID=g,GCF_path=str(fai)) for g in gcfs])
            fixtures.write_csv(records, ['Sample','CramFilePath','GCFID'],
                               [dict(Sample='r'+str(d),CramFilePath=str(p),GCFID=g) for g in gcfs for d,p in sources.items()])
            config=dict(storage_mode='native',window_size=70,cpus=1,python=sys.executable,memory=None,geometry_format='packed',
                        policies=dict(mate_policy='both',sampling_unit='read',min_mapq=0,exclude_duplicates=False,exclude_qcfail=False))
            database=root/'db';create_plan(records,refs,database,config,fixtures.COLUMNS)
            for g in gcfs: self.assertEqual(run_gcf(database,g),0)
            self.assertTrue(collect(database)['complete'])
            index=root/'index.sqlite';build_index(database,index,True)
            low=bam('q1',1,23); exact=bam('q5',5,17)
            rows=[dict(Sample='q1',GCFID=gcfs[0],Depth=1,CramPath=str(low)),
                  dict(Sample='q5',GCFID=gcfs[0],Depth=5,CramPath=str(exact)),
                  dict(Sample='r6',GCFID=gcfs[1],Depth=6,CramPath=str(sources[6])),
                  dict(Sample='missing-gcf',GCFID='missing',Depth=1,CramPath=str(low)),
                  dict(Sample='bad',GCFID=gcfs[0],Depth='nan',CramPath=str(low))]
            manifest=root/'q.tsv';table(manifest, ('Sample','GCFID','Depth','CramPath'), rows)
            before={str(p):digest(p) for p in root.rglob('*') if p.is_file()}
            entry=Path(__file__).resolve().parents[1]/'bin/query_manifest.py'
            summaries=[]
            budgets = ((1,512,512),(2,512,512),(1,0,0),(2,0,0),(2,0,512),(2,512,0))
            for case, (jobs, cache_mib, work_mib) in enumerate(budgets):
                out=root/f'out{case}'
                command=[sys.executable,'-B',str(entry),'--manifest',str(manifest),'--index',str(index),
                         '--outdir',str(out),'--window-size','70','--min-samples','2','--simulations','12',
                         '--reference-draws','20','--jobs',str(jobs),'--seed','91',
                         '--cache-mib',str(cache_mib),'--work-mib',str(work_mib)]
                result=subprocess.run(command, capture_output=True, text=True)
                self.assertEqual(result.returncode,1,result.stderr)  # retained bad input row
                with (out/'summary.tsv').open() as handle: summary=list(csv.DictReader(handle,delimiter='\t'))
                self.assertEqual(len(summary),5)
                self.assertEqual([r['arm'] for r in summary[:3]],['low','low','high'])
                self.assertEqual(summary[0]['status'],'COMPLETE')
                self.assertEqual(summary[1]['status'],'COMPLETE_WITH_GAPS')
                self.assertEqual(summary[3]['status'],'INSUFFICIENT_REFERENCE_SAMPLES')
                self.assertEqual(summary[4]['status'],'ERROR')
                self.assertGreater(int(summary[2]['duplicate_draws']),0)
                self.assertEqual(summary[2]['high_null_policy'],'draw-pairs')
                self.assertLessEqual(int(summary[2]['high_null_sampled_pairs']),1000)
                self.assertNotEqual(summary[2]['pearson_p'],'')
                import gzip
                with gzip.open(out/'queries/q00000003/draws.tsv.gz','rt') as handle:
                    draws=list(csv.DictReader(handle,delimiter='\t'))
                self.assertTrue(any(d['same_source']=='True' and d['same_sample']=='True' for d in draws))
                summaries.append((out,summary))
                recorded = json.loads((out/'run.json').read_text())['config']
                for name, amount in (('cache',cache_mib),('work',work_mib)):
                    self.assertEqual(recorded[name+'_mib'], amount)
                    self.assertEqual(recorded[name+'_bytes'], None if amount == 0 else amount*1024**2//jobs)
                audit = json.loads((out/'queries/q00000001/audit.json').read_text())
                self.assertEqual(audit['cache']['budget_bytes'], recorded['cache_bytes'])
                self.assertEqual(subprocess.run(command+['--resume'],capture_output=True).returncode,1)
                changed_budget = '1' if cache_mib == 0 else '0'
                changed = subprocess.run(command+['--cache-mib',changed_budget,'--resume'],capture_output=True,text=True)
                self.assertEqual(changed.returncode,2)
                self.assertIn('Resume configuration/source/input/database identity differs', changed.stderr)
                self.assertEqual(subprocess.run(command,capture_output=True).returncode,2)
            for _, summary in summaries[1:]:
                self.assertEqual(summaries[0][1], summary)
            for p in (summaries[0][0]/'queries').rglob('*'):
                if p.is_file() and p.name not in ('audit.json','checksums.json'):
                    for other_out, _ in summaries[1:]:
                        other=other_out/p.relative_to(summaries[0][0])
                        self.assertEqual(p.read_bytes(),other.read_bytes(),str(p))
            self.assertEqual(before,{p:digest(p) for p in before})
            # Additional N-mode real execution and explicit missing high null.
            minimal=root/'one.tsv';table(minimal, ('Sample','GCFID','Depth','CramPath'), [rows[2]])
            command=[sys.executable,'-B',str(entry),'--manifest',str(minimal),'--index',str(index),
                     '--outdir',str(root/'nmode'),'--bin-count','6','--min-samples','2','--reference-draws','100%',
                     '--high-null','none']
            result=subprocess.run(command,capture_output=True,text=True)
            self.assertEqual(result.returncode,0,result.stderr)
            with (root/'nmode/summary.tsv').open() as handle: nrow=next(csv.DictReader(handle,delimiter='\t'))
            self.assertEqual(nrow['requested_draws'],'3')
            self.assertEqual(nrow['gcf_reference_length'],'400')
            self.assertEqual(nrow['nominal_window_size'],'67')
            self.assertEqual(nrow['requested_bin_count'],'6')
            self.assertEqual(nrow['jaccard_status'],'HIGH_NULL_NOT_DEFINED')
            self.assertEqual(nrow['jaccard_p'],'')
            # Strict gate: exactly 3 available does not pass >3.
            command=[sys.executable,'-B',str(entry),'--manifest',str(manifest),'--index',str(index),
                     '--outdir',str(root/'gate'),'--windows-per-contig','6','--min-samples','3']
            self.assertEqual(subprocess.run(command,capture_output=True).returncode,1)
            with (root/'gate/summary.tsv').open() as handle: gated=list(csv.DictReader(handle,delimiter='\t'))
            self.assertEqual(gated[0]['status'],'INSUFFICIENT_REFERENCE_SAMPLES')
            self.assertIn('deprecated',subprocess.run(command+['--help'],capture_output=True,text=True).stdout)
            # Both spellings select the same corrected global divisor.
            parsed=parser().parse_args(['--manifest',str(manifest),'--index',str(index),'--outdir','unused',
                                       '--windows-per-contig','6'])
            self.assertEqual(parsed.bin_count,6)
            # Mutation blocks resume before any saved query is accepted.
            manifest.write_text(manifest.read_text()+'\n')
            self.assertEqual(subprocess.run(command+['--resume'],capture_output=True).returncode,2)


if __name__ == '__main__': unittest.main()
