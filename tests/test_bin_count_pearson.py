"""Independent Pearson references and immutable saved-matrix CLI acceptance."""
from contextlib import redirect_stderr
import io as textio
import json
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest

import numpy as np
from scipy.stats import pearsonr

sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from bin import bin_count_eval_io as io
from bin import bin_count_eval_math as maths
from bin.bin_count_evaluate import evaluate
from tests.test_bin_count_evaluate import fixture


class PearsonMathTests(unittest.TestCase):
    def test_independent_pearson_reference_and_affine_invariance(self):
        x=np.array([1.,2.,4.,7.])
        vectors=np.array([x, -x, [1.,5.,2.,4.], 3*x+9])
        expected=[1-float(pearsonr(x,y).statistic) for y in vectors]
        observed=maths.distances(x,vectors,'correlation')[0]
        np.testing.assert_allclose(observed,expected,atol=1e-14)
        self.assertAlmostEqual(observed[0],0)
        self.assertAlmostEqual(observed[1],2)
        self.assertAlmostEqual(observed[3],0)
        zero=maths.distances([-1.,0.,1.],[1.,-2.,1.],'correlation')[0,0]
        self.assertAlmostEqual(zero,1)
        np.testing.assert_allclose(maths.distances(maths.normalized(x),
            maths.normalized(vectors[2]),'correlation')[0,0],expected[2],atol=1e-14)

    def test_center_distances_and_equal_gcf_weighting(self):
        matrix=np.array([[1.,2.,4.,8.],[2.,1.,5.,7.],[8.,4.,2.,1.],
                         [7.,5.,1.,2.],[6.,4.,3.,2.]])
        matrix=np.array([maths.normalized(row) for row in matrix])
        labels=['a','a','b','b','b']
        centers=[matrix[:2].mean(axis=0),matrix[2:].mean(axis=0)]
        within=np.array([1-pearsonr(x,centers[0 if i<2 else 1]).statistic
                         for i,x in enumerate(matrix)])
        macro=(within[:2].mean()+within[2:].mean())/2
        between=1-pearsonr(*centers).statistic
        result,groups,actual,_,per_row,pairs=maths.score(matrix,labels,'correlation')
        np.testing.assert_allclose(actual,centers)
        np.testing.assert_allclose(per_row,within,atol=1e-14)
        self.assertEqual(groups,['a','b'])
        self.assertAlmostEqual(result['within_mean'],macro)
        self.assertNotAlmostEqual(macro,float(within.mean()),places=7)
        self.assertAlmostEqual(result['between_mean'],between)
        self.assertAlmostEqual(result['ratio'],between/macro)
        self.assertAlmostEqual(result['within_pearson_r_mean'],1-macro)
        self.assertAlmostEqual(result['between_pearson_r_mean'],1-between)
        self.assertEqual(len(pairs),1)

    def test_pairwise_bounds_symmetry_and_independent_reference(self):
        matrix=np.random.default_rng(391).uniform(.1,9,size=(7,11))
        expected=np.array([[1-pearsonr(a,b).statistic for b in matrix] for a in matrix])
        actual=maths.pairwise(matrix,'correlation')
        np.testing.assert_allclose(actual,expected,atol=1e-14)
        np.testing.assert_array_equal(actual,actual.T)
        np.testing.assert_array_equal(actual.diagonal(),np.zeros(7))
        self.assertTrue(np.all((actual>=0)&(actual<=2)))

    def test_undefined_rows_centers_and_one_bin_are_never_imputed(self):
        for row in ([.1,.1,.1],[1.],[0.,0.]):
            with self.assertRaises(maths.UndefinedCorrelationError):
                maths.distances(row,row,'correlation')
            with self.assertRaises(maths.UndefinedCorrelationError):
                maths.pairwise(np.array([row,row]),'correlation')
        # Both rows vary, but their mean-profile center is exactly flat.
        matrix=np.array([[1.,2.,3.],[3.,2.,1.],[2.,1.,4.]])
        with self.assertRaises(maths.UndefinedCorrelationError) as caught:
            maths.score(matrix,['a','a','b'],'correlation')
        self.assertEqual(caught.exception.issues,
                         [dict(kind='center',index=0,reason='CONSTANT_PROFILE')])

    def test_near_constants_warn_without_changing_values_or_population(self):
        row=np.array([1000.,1000.+1e-12,1000.-1e-12])
        before=row.copy()
        issues=maths.correlation_diagnostics(row)
        self.assertEqual(issues,[dict(kind='row',index=0,reason='NEAR_CONSTANT_PROFILE')])
        distance=maths.distances(row,[1.,2.,3.],'correlation')[0,0]
        self.assertTrue(np.isfinite(distance))
        np.testing.assert_array_equal(before,row)

    def test_single_gcf_and_zero_within_keep_explicit_null_ratio(self):
        rows=np.array([[1.,2.,3.,7.],[1.,2.,3.,7.],[7.,3.,2.,1.],[7.,3.,2.,1.]])
        result=maths.score(rows,['a','a','b','b'],'correlation')[0]
        self.assertEqual(result['status'],'ZERO_WITHIN_DISTANCE')
        self.assertIsNone(result['ratio'])
        result=maths.score(rows[:2],['a','a'],'correlation')[0]
        self.assertEqual(result['status'],'SINGLE_GCF')
        self.assertIsNone(result['between_mean'])
        self.assertIsNone(result['ratio'])


class PearsonSavedMatrixTests(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory();self.addCleanup(self.tmp.cleanup)
        self.root=Path(self.tmp.name)
        pre,db=fixture(self.root)
        io.create_plan(pre,db,self.root/'plan',[10],[5])
        io.build(self.root/'plan',self.root/'matrices')
        # Matrix evaluation must not reopen any original biological input.
        self.input_hashes={str(p.relative_to(self.root/'matrices')):io.sha(p)
                           for p in (self.root/'matrices').rglob('*') if p.is_file()}
        db.unlink()
        shutil.rmtree(self.root/'source1');shutil.rmtree(self.root/'source2')
        self.script=Path(__file__).resolve().parents[1]/'bin/bin_count_evaluate.py'

    def test_cli_all_four_metrics_and_embeddings_reuse_sealed_matrices(self):
        command=[sys.executable,'-B',str(self.script),'evaluate','--matrices',str(self.root/'matrices'),
                 '--outdir',str(self.root/'evaluation'),'--distance',*maths.METRICS,
                 '--ordination','tsne','pcoa','--perplexity','3','--iterations','300','--seeds','7']
        run=subprocess.run(command,capture_output=True,text=True)
        self.assertEqual(run.returncode,0,run.stdout+run.stderr)
        _,manifest=io.load(self.root/'evaluation','evaluation')
        self.assertEqual(manifest['distance'],list(maths.METRICS))
        self.assertEqual(len(manifest['embeddings']),8)
        self.assertEqual({r['status'] for r in manifest['embeddings']},{'COMPLETE'})
        self.assertEqual(len(io.read_tsv(self.root/'evaluation/scores.tsv')),4)
        candidate='stratum_01/B4/correlation'
        folder=self.root/'evaluation'/candidate
        for table in ('within_rows.tsv','between_centers.tsv'):
            records=io.read_tsv(folder/table)
            self.assertTrue(records)
            for r in records:
                self.assertAlmostEqual(float(r['pearson_r']),1-float(r['distance']))
        self.assertEqual(io.read_tsv(folder/'pearson_diagnostics.tsv'),[])
        for relative,sha in self.input_hashes.items():
            self.assertEqual(io.sha(self.root/'matrices'/relative),sha)
        evaluate(self.root/'matrices',self.root/'repeat',['correlation'],['none'],[9],3,300,'auto',1,8)
        scores=json.loads((self.root/'evaluation/scores.json').read_text())
        pearson=[r for r in scores if r['distance']=='correlation']
        self.assertEqual(pearson,json.loads((self.root/'repeat/scores.json').read_text()))
        second=subprocess.run(command,capture_output=True,text=True)
        self.assertNotEqual(second.returncode,0)
        self.assertIn('File exists',second.stderr)

    def test_constant_profiles_preserve_identifiers_and_fail_without_success(self):
        root,manifest=io.load(self.root/'matrices','matrices')
        folder=root/manifest['candidates'][0]['folder']
        with io.mapped_matrix(folder/'normalized.npy',mode='r+') as matrix:
            matrix[0]=1e6/matrix.shape[1]
        # Reseal a synthetic deliberately constant input, not a production file.
        io.seal(root,manifest)
        stderr=textio.StringIO()
        with redirect_stderr(stderr),self.assertRaises(maths.UndefinedCorrelationError):
            evaluate(root,self.root/'bad',['correlation'],['none'],[7],3,300,'auto',1,8)
        dest=self.root/'bad'/manifest['candidates'][0]['folder']/'correlation'
        issues=io.read_tsv(dest/'pearson_diagnostics.tsv')
        first=io.read_tsv(folder/'rows.tsv')[0]
        matching=[r for r in issues if r['kind']=='row' and r['index']=='0']
        self.assertEqual(len(matching),1)
        for name in ('row_id','row_key','sample','gcf'):
            self.assertEqual(matching[0][name],first[name])
        failure=json.loads((dest/'FAILURE.json').read_text())
        self.assertEqual(failure['records_removed'],0)
        self.assertEqual(failure['status'],'UNDEFINED_PEARSON')
        self.assertFalse((self.root/'bad/SUCCESS').exists())
        command=[sys.executable,'-B',str(self.script),'evaluate','--matrices',str(root),
                 '--outdir',str(self.root/'cli_bad'),'--distance','correlation','--ordination','none']
        run=subprocess.run(command,capture_output=True,text=True)
        self.assertNotEqual(run.returncode,0)
        self.assertIn('Pearson correlation undefined',run.stderr)
        self.assertIn('pearson_diagnostics.tsv',run.stderr)
        self.assertFalse((self.root/'cli_bad/SUCCESS').exists())

    def test_near_constant_rows_are_reported_and_retained(self):
        root,manifest=io.load(self.root/'matrices','matrices')
        folder=root/manifest['candidates'][0]['folder']
        with io.mapped_matrix(folder/'normalized.npy',mode='r+') as matrix:
            matrix[0]=[250000.,250000.+1e-9,250000.-1e-9,250000.]
        io.seal(root,manifest)
        stderr=textio.StringIO()
        with redirect_stderr(stderr):
            evaluate(root,self.root/'near',['correlation'],['none'],[7],3,300,'auto',1,8)
        dest=self.root/'near'/manifest['candidates'][0]['folder']/'correlation'
        issues=io.read_tsv(dest/'pearson_diagnostics.tsv')
        self.assertTrue(any(r['kind']=='row' and r['index']=='0' and
                            r['reason']=='NEAR_CONSTANT_PROFILE' for r in issues))
        self.assertIn('WARNING Pearson',stderr.getvalue())
        self.assertEqual(len(io.read_tsv(dest/'within_rows.tsv')),16)
        self.assertTrue((self.root/'near/SUCCESS').exists())


if __name__=='__main__':
    unittest.main()
