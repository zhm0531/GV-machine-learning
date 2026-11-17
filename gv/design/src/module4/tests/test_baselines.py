# -*- coding: utf-8 -*-
"""模块4快速回归测试；从 gv/design/src 运行：python -m unittest module4.tests.test_baselines"""
import unittest

from module4.baselines import (
    REF_SEQ, MASK, check, random_mutation_search, greedy_hill_climb
)


class BaselineTests(unittest.TestCase):
    def _assert_result(self, result):
        self.assertEqual(len(result.best_sequence), len(REF_SEQ))
        self.assertTrue(all(result.best_sequence[i] == REF_SEQ[i] for i in MASK))
        self.assertEqual(check(result.best_sequence), [])
        self.assertGreaterEqual(result.best_score, result.reference_score)
        self.assertGreaterEqual(result.eval_count, 1)
        self.assertLessEqual(result.eval_count, result.budget)
        bests = [float(row["best_so_far"]) for row in result.history]
        self.assertTrue(all(bests[i] <= bests[i + 1] for i in range(len(bests) - 1)))

    def test_random_mutation(self):
        result = random_mutation_search(budget=20, seed=42)
        self._assert_result(result)

    def test_greedy_hill_climb(self):
        result = greedy_hill_climb(budget=20, seed=42)
        self._assert_result(result)


if __name__ == "__main__":
    unittest.main()
