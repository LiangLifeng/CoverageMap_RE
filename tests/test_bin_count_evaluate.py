"""Independent geometry/distance references and synthetic frozen-index CLI acceptance."""
from contextlib import closing, redirect_stdout
import csv
import hashlib
import io as textio
import json
from pathlib import Path
import sqlite3
import subprocess
import sys
import tempfile
import unittest

import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from bin import bin_count_eval_io as io
from bin import bin_count_eval_math as maths
from bin.coverage_core import BlockPool, Reference, window_sums
from bin.reference_store import build_bundle
from bin.bin_count_evaluate import main, evaluate


def fixture(root, samples_per_gcf=4):
    """Four GCFs, two independent source roots, three codecs, no alignment files."""
    roots={1:root/'source1',2:root/'source2'}
    for p in roots.values():p.mkdir()
    db=root/'references.sqlite'
    pre=root/'preflight';pre.mkdir()
    rows=[];gcf_rows=[]
    with closing(sqlite3.connect(db)) as c, c:
        c.executescript('''CREATE TABLE metadata(key TEXT,value TEXT);
        CREATE TABLE sources(id INTEGER,root TEXT);
        CREATE TABLE gcfs(gcf TEXT,layout_json TEXT);
        CREATE TABLE references_table(id INTEGER,gcf TEXT,sample TEXT,genome_len INTEGER,
          full_mean_depth REAL,covered_fraction REAL,aligned_bases INTEGER,
          manifest_path TEXT,manifest_sha TEXT,checksums_sha TEXT,geometry_path TEXT,source_id INTEGER);''')
        c.execute('INSERT INTO metadata VALUES (?,?)',('schema',json.dumps('coverage-reference-sqlite-3')))
        c.executemany('INSERT INTO sources VALUES (?,?)',[(i,str(p)) for i,p in roots.items()])
        number=0
        for gi,(gcf,length) in enumerate([('A',23),('B',27),('C',31),('D',37)]):
            split=length//2
            layout=dict(names=['contig1','contig2'],lengths=[split,length-split],md5=['',''])
            ref=Reference(*(tuple(layout[k]) for k in ('names','lengths','md5')))
            c.execute('INSERT INTO gcfs VALUES (?,?)',(gcf,json.dumps(layout)))
            gcf_rows.append(dict(gcf=gcf,genome_length_bp=length,contig_count=2,
                                 total_records=samples_per_gcf,retained_records=samples_per_gcf))
            for si in range(samples_per_gcf):
                number+=1;sample=f'sample{si}'
                sid=1+gi%2
                attempt=roots[sid]/gcf/sample
                attempt.mkdir(parents=True)
                intervals=[[0,split],[split,length]]*3
                if gi%2==0:
                    intervals += [[si%3,split]]*(si+1)
                else:
                    intervals += [[split,min(length,split+si+4)]]*(si+1)
                pool=BlockPool.from_intervals(np.asarray(intervals,dtype=np.int64),length)
                build_bundle(pool,ref,attempt/'geometry',sample,{},[],workers=1,chunk_units=3,
                             storage_format=['packed','npz','npy'][si%3])
                record=dict(gcf=gcf,sample=sample,genome_len=length,aligned_bases=pool.total_bases,
                            full_mean_depth=pool.total_bases/length,layout=layout)
                (attempt/'status').mkdir()
                (attempt/'status/SUCCESS').touch()
                (attempt/'status/exit_code.txt').write_text('0\n')
                io.write_json(attempt/'result.json',record)
                checks={str(p.relative_to(attempt)):io.sha(p) for p in (attempt/'geometry').rglob('*') if p.is_file()}
                io.write_json(attempt/'checksums.json',checks)
                mp=str((attempt/'result.json').relative_to(roots[sid]));gp=str((attempt/'geometry').relative_to(roots[sid]))
                c.execute('INSERT INTO references_table VALUES (?,?,?,?,?,?,?,?,?,?,?,?)',
                          (number,gcf,sample,length,record['full_mean_depth'],1.,pool.total_bases,mp,
                           io.sha(attempt/'result.json'),io.sha(attempt/'checksums.json'),gp,sid))
                rows.append(dict(reference_id=number,gcf=gcf,sample=sample,genome_length_bp=length,
                                 full_mean_depth=record['full_mean_depth'],covered_fraction=1.,source_id=sid,
                                 manifest_path=mp,geometry_path=gp))
    io.write_tsv(pre/'filtered_records.tsv',list(rows[0]),rows)
    io.write_tsv(pre/'filtered_gcf_lengths.tsv',list(gcf_rows[0]),gcf_rows)
    samples = sorted({r['sample'] for r in rows})
    (pre/'filtered_samples.txt').write_text(''.join(f'{sample}\n' for sample in samples))
    io.write_json(pre/'summary.json',dict(counts=dict(retained_records=4*samples_per_gcf,
        retained_gcfs=4,retained_unique_samples=samples_per_gcf), input_sha256=io.sha(db),
        filters=dict(depth_gt=0,cover_fraction_gt=0,genome_length_bp_gt=0,min_samples_per_gcf=samples_per_gcf)))
    (pre/'SHA256SUMS').write_text(''.join(f'{io.sha(p)}  {p.name}\n' for p in sorted(pre.iterdir())))
    (pre/'SUCCESS').write_text('COMPLETE\n')
    return pre,db


class MathTests(unittest.TestCase):
    def test_mapped_matrix_closes_on_success_and_exception_with_live_references(self):
        with tempfile.TemporaryDirectory() as tmp:
            path=Path(tmp)/'mapped.npy'
            with io.mapped_matrix(path,mode='w+',shape=(2,3)) as writer:
                writer[:]=1
                kept_writer=writer
            self.assertTrue(kept_writer._mmap.closed)
            with self.assertRaisesRegex(RuntimeError,'injected'):
                with io.mapped_matrix(path) as reader:
                    kept_reader=reader
                    np.testing.assert_array_equal(reader,np.ones((2,3)))
                    raise RuntimeError('injected')
            self.assertTrue(kept_reader._mmap.closed)
            path.unlink()
            self.assertFalse(path.exists())

    def test_remainder_and_cross_contig_dense_reference(self):
        edges=maths.exact_edges(23,5)
        np.testing.assert_array_equal(edges,[0,4,8,12,16,23])
        # Contig boundary at10 is crossed by bin[8,12); no block spans boundary.
        blocks=np.array([[0,10],[10,23],[7,10],[10,13],[8,10],[10,12]],dtype=np.int64)
        pool=BlockPool.from_intervals(blocks,23)
        dense=np.zeros(23,dtype=np.int64)
        for begin,end in blocks:dense[begin:end]+=1
        observed=window_sums(pool,np.arange(len(pool)),edges)
        expected=np.array([dense[a:b].sum() for a,b in zip(edges[:-1],edges[1:])])
        np.testing.assert_array_equal(observed,expected)
        self.assertEqual(observed.sum(),pool.total_bases)
        np.testing.assert_allclose(maths.normalized(observed/np.diff(edges)).sum(),1e6)
        with self.assertRaises(ValueError):maths.exact_edges(4,5)
        with self.assertRaises(ValueError):maths.normalized([0,0])

    def test_strata_boundaries_aliases_and_range(self):
        layouts={g:dict(lengths=[n]) for g,n in [('a',10),('b',29),('c',30),('d',100)]}
        rows=[dict(gcf=g) for g in layouts]
        strata=io.strata_plan(rows,layouts,[10,30,100,1000],[1,3,4])
        self.assertEqual(strata[0]['gcfs'],['a','b'])
        self.assertEqual(strata[1]['gcfs'],['c'])
        self.assertEqual(strata[2]['gcfs'],['d'])
        self.assertEqual(strata[3]['status'],'SKIPPED_EMPTY')
        self.assertEqual(io.sizes(start=1,stop=100,step=50),[1,51])
        self.assertEqual(io.sizes(),[1,50,100])
        alias=io.strata_plan([dict(gcf='x')],{'x':dict(lengths=[10])},[0],[4,5])
        self.assertEqual(alias[0]['candidates'][0]['bin_sizes_bp'],[4,5])
        with self.assertRaises(ValueError):io.sizes([1],1,10,1)
        with self.assertRaises(ValueError):io.strata_plan(rows,layouts,[11],[1])

    def test_metric_references_and_macro_weighting(self):
        from scipy.spatial.distance import cosine,jensenshannon
        self.assertEqual(maths.distances([0,0],[3,4],'euclidean')[0,0],5)
        x=np.array([1.,3.]);y=np.array([3.,1.])
        self.assertAlmostEqual(maths.distances(x,y,'cosine')[0,0],cosine(x,y))
        self.assertAlmostEqual(maths.distances(x,y,'jensenshannon')[0,0],jensenshannon(x,y))
        matrix=np.array([[1.,0.],[3.,0.],[10.,0.],[10.,0.],[10.,0.],[10.,0.]])
        s,_,_,within,_,pairs=maths.score(matrix,['a','a','b','b','b','b'],'euclidean')
        self.assertEqual(s['within_mean'],.5) # macro mean, NOT pooled2/6
        self.assertEqual(s['between_mean'],8)
        self.assertEqual(s['ratio'],16)
        self.assertEqual(maths.score(matrix[:2],['a','a'],'euclidean')[0]['status'],'SINGLE_GCF')
        self.assertEqual(maths.score(np.ones((4,2)),['a','a','b','b'],'euclidean')[0]['status'],'ZERO_WITHIN_DISTANCE')

    def test_pcoa_reconstructs_distances_and_records_negative_spectrum(self):
        from scipy.spatial.distance import pdist,squareform
        d=squareform(pdist(np.array([[0.,0.],[1.,0.],[0.,2.],[1.,2.]])))
        points,meta=maths.pcoa(d)
        np.testing.assert_allclose(squareform(pdist(points)),d,atol=1e-12)
        _,negative=maths.pcoa(np.array([[0.,1.,1.],[1.,0.,3.],[1.,3.,0.]]))
        self.assertEqual(negative['warning'],'NEGATIVE_EIGENVALUES')

    def test_parallel_pairwise_and_scores_match_serial_exactly(self):
        matrix = np.random.default_rng(31).uniform(size=(151, 137))
        matrix[matrix < .3] = 0
        matrix[6] = matrix[5]
        labels = [f'g{i % 13}' for i in range(len(matrix))]
        for metric in maths.METRICS:
            with self.subTest(metric=metric):
                for values in (matrix, matrix.astype(np.float32)):
                    serial = maths.pairwise(values, metric)
                    np.testing.assert_array_equal(maths.pairwise(values, metric, jobs=7), serial)
                expected = maths.score(matrix, labels, metric)
                observed = maths.score(matrix, labels, metric, jobs=7)
                for index in (0, 1, 3, 5):
                    self.assertEqual(expected[index], observed[index])
                for index in (2, 4):
                    np.testing.assert_array_equal(expected[index], observed[index])

    def test_tsne_deterministic_and_perplexity_rejection(self):
        d=maths.pairwise(np.random.default_rng(5).uniform(size=(16,4)),'euclidean')
        a,_=maths.tsne(d,19,3,300,'auto',1)
        b,_=maths.tsne(d,19,3,300,'auto',1)
        np.testing.assert_array_equal(a,b)
        with self.assertRaisesRegex(ValueError,'perplexity'):maths.tsne(d,19,16,300,'auto',1)

    def test_blocked_tsne_preserves_exact_neighbors_and_coordinates(self):
        from scipy.spatial.distance import pdist, squareform
        from sklearn.neighbors import NearestNeighbors
        from threadpoolctl import threadpool_limits

        rng = np.random.default_rng(24)
        values = rng.integers(1, 4, (275, 7)).astype(np.float64)
        values[:, 0] = 4
        values[:40] = values[40]
        perplexity = 5
        size = maths.tsne_neighbor_count(len(values), perplexity) + 1
        with threadpool_limits(limits=1):
            for metric in maths.METRICS:
                with self.subTest(metric=metric):
                    dense = squareform(pdist(values, metric=metric))
                    if metric == 'correlation':
                        np.testing.assert_array_equal(maths.pairwise(values, metric, jobs=3), dense)
                    oracle = NearestNeighbors(metric='precomputed').fit(dense)
                    expected = oracle.kneighbors_graph(dense, n_neighbors=size, mode='distance')
                    graph = maths.tsne_neighbor_graph(values, metric, perplexity, jobs=3)
                    np.testing.assert_array_equal(graph.indptr, expected.indptr)
                    np.testing.assert_array_equal(graph.indices, expected.indices)
                    np.testing.assert_array_equal(graph.data.view(np.uint64), expected.data.view(np.uint64))
                    dense_xy, dense_meta = maths.tsne(dense, 19, perplexity, 300, 'auto', 1)
                    graph_xy, graph_meta = maths.tsne(graph, 19, perplexity, 300, 'auto', 1)
                    np.testing.assert_array_equal(dense_xy.view(np.uint32), graph_xy.view(np.uint32))
                    self.assertEqual(dense_meta, graph_meta)

    def test_streamed_edges_keep_legacy_bytes_across_chunks(self):
        lengths = {'B': 40003, 'A': 21887, '非ASCII': 18010}
        bins = 9001
        with tempfile.TemporaryDirectory() as tmp:
            streamed = Path(tmp) / 'streamed.json'
            legacy = Path(tmp) / 'legacy.json'
            io.write_edges(streamed, lengths, bins)
            io.write_json(legacy, {g: maths.exact_edges(length, bins).tolist()
                                  for g, length in lengths.items()})
            self.assertEqual(streamed.read_bytes(), legacy.read_bytes())


class PipelineTests(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory();self.addCleanup(self.tmp.cleanup)
        self.root=Path(self.tmp.name);self.pre,self.db=fixture(self.root)

    def plan(self,name='plan'):
        return io.create_plan(self.pre,self.db,self.root/name,[10,30,100],[1,5,10])

    def test_grid_budget_check_and_tsne_only(self):
        config = self.root / 'strata.json'
        io.write_json(config, {'strata': [
            dict(min_genome_bp=10, max_genome_bp_exclusive=None,
                 nominal_bin_widths_bp=[10, 23]),
        ]})
        common = ['--ordination', 'tsne', '--distance', 'euclidean']
        check = ['check-grid', '--preflight-dir', str(self.pre), '--strata-config', str(config),
                 *common, '--jobs', '1']
        with redirect_stdout(textio.StringIO()) as report:
            main([*check, '--index', str(self.db), '--max-matrix-gib', '1', '--max-working-gib', '1'])
        self.assertIn('Validated frozen SQLite', report.getvalue())
        self.assertIn('peak evaluation per job', report.getvalue())
        wrong_index = self.root / 'wrong.sqlite'
        wrong_index.write_bytes(self.db.read_bytes() + b'changed')
        with self.assertRaises(SystemExit) as error:
            main([*check, '--index', str(wrong_index), '--max-matrix-gib', '1', '--max-working-gib', '1'])
        self.assertEqual(error.exception.code, 1)
        for matrix_limit, work_limit in [('0.00000001', '1'), ('1', '0.00000001')]:
            with self.assertRaises(SystemExit) as error:
                main([*check, '--max-matrix-gib', matrix_limit, '--max-working-gib', work_limit])
            self.assertEqual(error.exception.code, 1)

        plan = self.root / 'plan'
        main(['plan', '--preflight-dir', str(self.pre), '--index', str(self.db),
              '--strata-config', str(config), '--outdir', str(plan)])
        matrices = self.root / 'matrices'
        io.build(plan, matrices)
        output = self.root / 'plots'
        main(['evaluate', '--matrices', str(matrices), '--outdir', str(output), *common,
              '--plots-only', '--perplexity', '3', '--iterations', '300', '--seeds', '7'])
        _, result = io.load(output, 'evaluation')
        self.assertTrue(result['plots_only'])
        self.assertEqual([(e['bin_count'], e['ordination'], e['status']) for e in result['embeddings']],
                         [(2, 'tsne', 'COMPLETE')])
        view = output / result['embeddings'][0]['folder']
        self.assertEqual(len(io.read_tsv(view / 'coordinates.tsv')), 16)
        self.assertEqual(len(list(view.glob('*.pdf'))), 1)
        for unwanted in ('scores*', 'centers.npy', 'pairwise_distances.npy', '*curves*'):
            self.assertFalse(list(output.rglob(unwanted)))

    def test_end_to_end_geometry_build_and_reevaluation(self):
        before=io.sha(self.db)
        self.plan();io.build(self.root/'plan',self.root/'matrices')
        root,m=io.load(self.root/'matrices','matrices')
        self.assertEqual(len(m['candidates']),6)
        plan=json.loads((self.root/'plan/result.json').read_text())
        for c in m['candidates']:
            folder=root/c['folder'];raw=np.load(folder/'raw_mean_depth.npy');norm=np.load(folder/'normalized.npy')
            order=io.read_tsv(folder/'rows.tsv');edges=json.loads((folder/'edges.json').read_text())
            for i,r in enumerate(order):
                row=next(x for x in plan['rows'] if x['id']==int(r['reference_id']))
                pool=io.load_pool(row,plan['layouts'][r['gcf']])
                dense=np.zeros(pool.genome_len)
                for a,b in zip(pool.starts,pool.ends):dense[a:b]+=1
                e=edges[r['gcf']]
                expected=[dense[a:b].mean() for a,b in zip(e[:-1],e[1:])]
                np.testing.assert_array_equal(raw[i],expected)
            np.testing.assert_allclose(norm.sum(axis=1),1e6)
        evaluate(root,self.root/'eval',list(maths.METRICS),['pcoa'],[7],3,300,'auto',1,1)
        evaluate(root,self.root/'eval2',list(maths.METRICS),['none'],[99],3,300,'auto',1,1)
        self.assertEqual(io.sha(self.root/'eval/scores.json'),io.sha(self.root/'eval2/scores.json'))
        self.assertEqual(before,io.sha(self.db))
        with self.assertRaises(FileExistsError):self.plan()
        with self.assertRaisesRegex(ValueError,'matrix disk'):io.build(self.root/'plan',self.root/'small',1e-12)
        with self.assertRaisesRegex(ValueError,'working memory'):
            evaluate(root,self.root/'lowmem',['euclidean'],['pcoa'],[7],3,300,'auto',1,1e-12)
        with (root/m['candidates'][0]['folder']/'normalized.npy').open('ab') as h:h.write(b'changed')
        with self.assertRaisesRegex(ValueError,'Changed artifact'):io.load(root,'matrices')

    def test_cli_run_tsne_and_all_help(self):
        script=Path(__file__).parents[1]/'bin/bin_count_evaluate.py'
        for sub in ([],['plan'],['build'],['evaluate'],['run']):
            r=subprocess.run([sys.executable,'-B',str(script),*sub,'--help'],capture_output=True,text=True)
            self.assertEqual(r.returncode,0,r.stderr)
        r=subprocess.run([sys.executable,'-B',str(script),'run','--preflight-dir',str(self.pre),
            '--index',str(self.db),'--outdir',str(self.root/'run'),'--genome-size-breaks-bp','10','100',
            '--bin-sizes-bp','10','--ordination','tsne','pcoa','--distance','cosine',
            '--perplexity','3','--iterations','300','--seeds','7','8'],capture_output=True,text=True)
        self.assertEqual(r.returncode,0,r.stdout+r.stderr)
        _,result=io.load(self.root/'run/evaluation','evaluation')
        self.assertEqual(len(result['embeddings']),3)
        self.assertEqual({e['status'] for e in result['embeddings']},{'COMPLETE'})
        self.assertTrue((self.root/'run/SUCCESS').exists())

    def test_strata_json_run_and_equivalent_legacy_payloads(self):
        script = Path(__file__).parents[1] / 'bin/bin_count_evaluate.py'
        config = self.root / 'strata.json'
        strata = [
            dict(min_genome_bp=10, max_genome_bp_exclusive=30, nominal_bin_widths_bp=[5]),
            dict(min_genome_bp=30, max_genome_bp_exclusive=None, nominal_bin_widths_bp=[7]),
        ]
        io.write_json(config, {'strata': strata})
        output = self.root / 'configured'
        command = [sys.executable, '-B', str(script), 'run', '--preflight-dir', str(self.pre),
                   '--index', str(self.db), '--outdir', str(output), '--strata-config', str(config),
                   '--ordination', 'none', '--jobs', '4', '--distance', *maths.METRICS]
        run = subprocess.run(command, capture_output=True, text=True)
        self.assertEqual(run.returncode, 0, run.stdout + run.stderr)
        self.assertTrue((output/'SUCCESS').is_file())
        _, planned = io.load(output/'plan', 'plan')
        self.assertEqual(planned['strata_config'], strata)
        _, configured = io.load(output/'matrices', 'matrices')
        self.assertEqual(configured['build_jobs'], 4)
        self.assertEqual([(c['stratum'], c['bin_count']) for c in configured['candidates']],
                         [('stratum_01', 4), ('stratum_02', 4)])
        io.create_plan(self.pre, self.db, self.root/'legacy_plan', [10, 30], [5, 7])
        io.build(self.root/'legacy_plan', self.root/'legacy_matrices')
        for candidate in configured['candidates']:
            folder = candidate['folder']
            for artifact in ('raw_mean_depth.npy', 'normalized.npy', 'edges.json', 'rows.tsv'):
                self.assertEqual(io.sha(output/'matrices'/folder/artifact),
                                 io.sha(self.root/'legacy_matrices'/folder/artifact))
        evaluate(self.root/'legacy_matrices', self.root/'legacy_eval', list(maths.METRICS),
                 ['none'], [7], 3, 300, 'auto', 1, 1)
        result = json.loads((output/'evaluation/scores.json').read_text())
        previous = json.loads((self.root/'legacy_eval/scores.json').read_text())
        self.assertEqual(result, [r for r in previous if r['bin_count'] == 4])
        self.assertEqual(len(result), 8)

        conflict = subprocess.run(command + ['--bin-sizes-bp', '5'], capture_output=True, text=True)
        self.assertNotEqual(conflict.returncode, 0)
        self.assertIn('cannot be combined', conflict.stderr)

        strata[1]['min_genome_bp'] = 31
        io.write_json(config, {'strata': strata})
        with self.assertRaisesRegex(ValueError, 'previous exclusive upper bound'):
            io.read_strata_config(config)

    def test_input_corruption_and_missing_geometry(self):
        self.plan()
        plan=json.loads((self.root/'plan/result.json').read_text());row=plan['rows'][0]
        geometry=Path(row['source_root'])/row['geometry_path']
        (geometry/'manifest.json').write_text('{}')
        with self.assertRaisesRegex(ValueError,'Geometry metadata'):
            io.build(self.root/'plan',self.root/'bad',jobs=4)
        self.assertFalse((self.root/'bad/SUCCESS').exists())
        (self.pre/'filtered_records.tsv').write_text('changed')
        with self.assertRaisesRegex(ValueError,'checksum'):self.plan('badplan')

    def test_root_map_is_explicit_and_index_mismatch_fatal(self):
        mapping={'1':str(self.root/'relocated')}
        io.create_plan(self.pre,self.db,self.root/'mapped',[10],[10],mapping)
        with self.assertRaisesRegex(ValueError,'Source-root-map'):
            io.create_plan(self.pre,self.db,self.root/'unknown',[10],[10],{'99':'/unknown'})
        with closing(sqlite3.connect(self.db)) as c,c:c.execute("UPDATE references_table SET sample='different' WHERE id=1")
        with self.assertRaisesRegex(ValueError,'hash differs'):self.plan('badidentity')


if __name__=='__main__':unittest.main()
