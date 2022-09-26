

from celest.schedule.alns import ALNS
from unittest import TestCase
import numpy as np


class TestALNS(TestCase):

    def setUp(self):
        self.alns = ALNS([1, 1, 1])

        self.alns.destroy_weights = [1, 2, 3]
        self.alns.destroy_functions = [0, 0, 0]

        self.alns.repair_weights = [1, 2, 3, 4]
        self.alns.repair_functions = [0, 0, 0, 0]

    def test_add_cos_func(self):
        self.alns.add_cost_function(True)
        self.assertTrue(self.alns.cost_function)

    def test_add_destroy_functions(self):
        self.alns.add_destroy_functions([True])
        self.assertTrue(np.all(self.alns.destroy_functions))

    def test_add_repair_functions(self):
        self.alns.add_repair_functions([True])
        self.assertTrue(np.all(self.alns.repair_functions))

    def test_add_completeness_function(self):
        self.alns.add_completeness_function(True)
        self.assertTrue(self.alns.completeness_function)

    def test_solve(self):

        def cost_function(x):
            return -len(x)

        def destroy_function(x, q):
            while q > 0:
                if len(x) > 0:
                    x.pop()
                q -= 1

        def repair_function(x, q):
            while q > 0:
                x.append(1)
                q -= 1

        def completeness_function(x):
            return len(x) == 5

        self.alns.add_cost_function(cost_function)
        self.alns.add_destroy_functions([destroy_function])
        self.alns.add_repair_functions([repair_function])
        self.alns.add_completeness_function(completeness_function)

        self.assertTrue(self.alns.solve(10, 0.5, 0.5, 1, 2) == [1, 1, 1, 1, 1])

    def test_destroy_index_is_int(self):
        self.assertIsInstance(self.alns._get_destroy_index(), int)

    def test_destroy_index_in_range(self):
        self.assertTrue(0 <= self.alns._get_destroy_index() <= 2)

    def test_get_probabilities(self):
        weights = [1, 2, 3]
        true_probabilities = [1 / 6, 2 / 6, 3 / 6]
        test_probabilities = self.alns._get_probabilities(weights)

        self.assertTrue(np.array_equal(true_probabilities, test_probabilities))

    def test_repair_index_is_int(self):
        self.assertIsInstance(self.alns._get_repair_index(), int)

    def test_repair_index_in_range(self):
        self.assertTrue(0 <= self.alns._get_repair_index() <= 3)

    def test_get_destroy_function(self):
        self.alns.destroy_functions = [True, False, False]
        self.assertTrue(self.alns._get_destroy_function(0))

    def test_get_repair_function(self):
        self.alns.repair_functions = [True, False, False]
        self.assertTrue(self.alns._get_repair_function(0))

    def test_accept_if_improved_solution(self):
        self.alns.cost_function = lambda x: len(x)
        self.assertTrue(self.alns._accept([1, 2, 3], [1, 2, 3, 4], 1))

    def test_accept_if_non_improved_solution(self):
        self.alns.cost_function = lambda x: len(x)
        self.assertIsInstance(self.alns._accept([1, 2, 3, 4], [1, 2, 3], 1), bool)

    def test_update_destroy_weights(self):
        self.alns._update_destroy_weights(0.5, 2, 0)
        self.assertTrue(self.alns.destroy_weights[0] == 0.5 + 1 / 6)

    def test_update_repair_weights(self):
        self.alns._update_repair_weights(0.5, 2, 0)
        self.assertTrue(self.alns.repair_weights[0] == 0.5 + 1 / 6)
