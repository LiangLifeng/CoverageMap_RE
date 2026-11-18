import unittest
from unittest.mock import patch
from bin.stage_timing import StageTiming
from bin.concordance_batch import query_tasks


class TimingTests(unittest.TestCase):
    def test_nested_exclusive_cpu_wall_and_exception(self):
        with patch.object(StageTiming, 'clock', side_effect=[(i, i*2, i*3) for i in range(6)]):
            timer = StageTiming()
            with timer.measure('outer'):
                with self.assertRaises(ValueError):
                    with timer.measure('inner'):
                        raise ValueError('failed stage')
            result = timer.snapshot()['stages']
        self.assertEqual({k:v['wall_seconds'] for k,v in result.items()},dict(other=2., outer=2., inner=1.))
        for row in result.values():
            self.assertEqual(row['user_cpu_seconds'], 2*row['wall_seconds'])
            self.assertEqual(row['system_cpu_seconds'], 3*row['wall_seconds'])

    def test_gcf_tasks_cover_queries_once_and_keep_affinity(self):
        queries = [dict(gcf=gcf, sample=str(i)) for i,gcf in enumerate(['b','a','b','c','b','a'])]
        tasks = query_tasks(queries, 'gcf')
        self.assertEqual([len(t) for t in tasks], [3,2,1])
        self.assertEqual(sorted(n for t in tasks for n,_ in t), list(range(6)))
        self.assertTrue(all(len({q['gcf'] for _,q in t})==1 for t in tasks))
        self.assertEqual(len(query_tasks(queries, 'query')),6)
