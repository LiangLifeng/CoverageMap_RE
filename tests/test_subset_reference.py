"""Full-FAI geometry with subset CRAM headers, persistence and downstream use."""
import csv
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from bin.coverage_core import Reference, merge_references, extract_pool, dense_coverage
from bin.ref_depth_worker import merge_layout, load_depth
from bin.reference_store import build_bundle, load_reference_bundle
try:
    import pysam
except ImportError:
    pysam = None


class CompatibilityTests(unittest.TestCase):
    def test_known_identities_accumulate_and_conflicts_remain_errors(self):
        first = Reference(('a', 'b'), (10, 20), ('one', ''))
        second = Reference(('a', 'b'), (10, 20), ('', 'two'))
        merged = merge_references(first, second)
        self.assertEqual(merged.md5, ('one', 'two'))
        with self.assertRaisesRegex(ValueError, 'M5'):
            merge_references(merged, Reference(('a', 'b'), (10, 20), ('', 'wrong')))
        with self.assertRaisesRegex(ValueError, 'layouts differ'):
            merge_references(first, Reference(('b', 'a'), (20, 10), ('', '')))
        layouts = [dict(names=list(r.names), lengths=list(r.lengths), md5=list(r.md5)) for r in (first, second)]
        combined = merge_layout(*layouts)
        self.assertEqual(combined['md5'], ['one', 'two'])
        with self.assertRaisesRegex(ValueError, 'M5'):
            merge_layout(combined, dict(names=['a', 'b'], lengths=[10, 20], md5=['', 'wrong']))


@unittest.skipIf(pysam is None, 'Synthetic CRAM tests require pysam on NH')
class SubsetIntegrationTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory(prefix='subset-reference-')
        self.root = Path(self.temporary.name).resolve()
        self.fasta = self.root / 'full.fna'
        self.fasta.write_text('>a\n' + 'A' * 20 + '\n>b\n' + 'C' * 30 + '\n>c\n' + 'G' * 40 + '\n')
        pysam.faidx(str(self.fasta))
        self.fai = Path(str(self.fasta) + '.fai')
        self.b = self.make_cram('b', 30, 'C')
        self.c = self.make_cram('c', 40, 'G')

    def tearDown(self):
        self.temporary.cleanup()

    def make_cram(self, name, length, base):
        path = self.root / (name + '.cram')
        with pysam.AlignmentFile(str(path), 'wc', header=dict(HD=dict(VN='1.6'), SQ=[dict(SN=name, LN=length)]), reference_filename=str(self.fasta)) as handle:
            for i in range(80):
                read = pysam.AlignedSegment(handle.header)
                read.query_name = f'{name}-{i}'; read.query_sequence = base * 5
                read.flag = 0; read.reference_id = 0; read.reference_start = i % (length - 5)
                read.mapping_quality = 40; read.cigar = [(0, 5)]
                handle.write(read)
        return path

    def test_full_coordinates_zero_regions_and_persisted_absence(self):
        pool, ref = extract_pool(self.b, reference_index=self.fai)
        self.assertEqual(ref.names, ('a', 'b', 'c'))
        self.assertEqual(ref.length, 90)
        self.assertEqual(pool.total_bases, 400)
        self.assertEqual(ref.missing_header_contigs, ('a', 'c'))
        self.assertEqual(ref.md5[0], '')
        expected = np.zeros(90, dtype=np.int64)
        for i in range(80):
            start = 20 + i % 25
            expected[start:start + 5] += 1
        np.testing.assert_array_equal(dense_coverage(pool, np.arange(len(pool))), expected)
        policies = dict(mate_policy='both', sampling_unit='read', min_mapq=0, exclude_duplicates=False, exclude_qcfail=False)
        bundle = self.root / 'bundle'
        build_bundle(pool, ref, bundle, 'subset', policies, [], workers=1, storage_format='packed')
        saved_pool, saved_ref, manifest = load_reference_bundle(bundle, 1)
        self.assertEqual(saved_ref.missing_header_contigs, ('a', 'c'))
        self.assertEqual(manifest['missing_header_contigs'], ['a', 'c'])
        np.testing.assert_array_equal(saved_pool.starts, pool.starts)
        from bin.pool_bridge import extract_with_python
        _, bridged, _ = extract_with_python(sys.executable, self.b, None, self.fai, policies, self.root)
        self.assertEqual(bridged.missing_header_contigs, ('a', 'c'))
        with self.assertRaisesRegex(ValueError, 'contig set'):
            extract_pool(self.b, self.fasta)  # FASTA mode remains strict.

    def test_extra_contigs_and_length_conflicts_rejected(self):
        for name, length in (('unknown', 30), ('b', 29)):
            path = self.root / (name + '.bam')
            with pysam.AlignmentFile(str(path), 'wb', header=dict(SQ=[dict(SN=name, LN=length)])):
                pass
            with self.assertRaisesRegex(ValueError, 'contig set|header mismatch'):
                extract_pool(path, reference_index=self.fai)

    def test_database_two_different_subsets_and_cached_query(self):
        # Preserve the optional-index fixture while verifying full-FAI native coordinates.
        for path in (self.b,self.c):
            sorted_path=path.with_name(path.stem+'-sorted.cram')
            pysam.sort('-O','CRAM','--reference',str(self.fasta),'-o',str(sorted_path),str(path))
            sorted_path.replace(path)
            pysam.index(str(path))
        records, refs = self.root / 'records.csv', self.root / 'refs.csv'
        gcf = 'GCF_000000001.1'
        with records.open('w', newline='') as handle:
            writer = csv.writer(handle); writer.writerow(['Sample', 'CramFilePath', 'GCFID'])
            writer.writerows([['b', str(self.b), gcf], ['c', str(self.c), gcf]])
        with refs.open('w', newline='') as handle:
            writer = csv.writer(handle); writer.writerow(['GCF_ID', 'GCF_path']); writer.writerow([gcf, str(self.fai)])
        database = self.root / 'database'
        scripts = Path(__file__).resolve().parents[1] / 'bin'
        subprocess.run([sys.executable, '-B', str(scripts / 'ref_depth_database.py'), 'run', '--mode', 'local',
            '--records', str(records), '--references', str(refs), '--outdir', str(database), '--window-size', '10', '--cpus', '2'], check=True, stdout=subprocess.DEVNULL)
        summary = json.loads((database / 'SUMMARY.json').read_text())
        self.assertTrue(summary['complete']); self.assertEqual(summary['ready_samples'], 2)
        with (database / 'SAMPLE_INDEX.tsv').open() as handle:
            rows = {r['sample']: r for r in csv.DictReader(handle, delimiter='\t')}
        for sample,missing in (('b',['a','c']),('c',['a','b'])):
            manifest=Path(rows[sample]['sample_manifest']);result=json.loads(manifest.read_text())
            self.assertEqual(result['missing_header_contigs'],missing)
            self.assertEqual(result['missing_header_contig_count'],2)
            self.assertEqual(result['full_mean_depth'],400/90)
            edges=np.load(manifest.parent/'edges.npy');self.assertEqual(edges[-1],90)
            with np.load(manifest.parent/'native_coverage.npz') as data:
                pointers = data['segment_ptr']
                for i, mass in enumerate(data['aligned_bases']):
                    contigs = data['segment_contig_index'][pointers[i]:pointers[i+1]]
                    if all(result['layout']['names'][number] in missing for number in contigs):
                        self.assertEqual(mass, 0)
        from bin.reference_index import build_index
        index=self.root/'reference.sqlite';build_index(database,index,True)
        query = self.root / 'query'
        subprocess.run([sys.executable, '-B', str(scripts / 'coverage_concordance.py'),
            '--comparison-policy','nearest','--reference-sqlite', str(index), '--test-alignment', str(self.c),
            '--query-sample','query','--gcf',gcf,'--reference-samples','b','--depth','1','--replicates','5','--outdir', str(query)], check=True, stdout=subprocess.DEVNULL)
        evidence = json.loads(next(query.glob('queries/*/*/comparisons/*/*/concordance.json')).read_text())
        self.assertEqual(evidence['reference_missing_header_contigs'], ['a', 'c'])
        self.assertEqual(evidence['test_missing_header_contigs'], ['a', 'b'])
        self.assertEqual(evidence['identity_unverified_contigs'], ['a', 'b', 'c'])


if __name__ == '__main__':
    unittest.main()
