"""A nominal signal must not be erased by a fixed 'neither direction' verdict."""
import contextlib
import importlib.util
import io
import itertools
import numpy as np
from pathlib import Path
import shutil
import tempfile
import unittest
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location('permutation_reading', ROOT / 'scripts/permutation_test.py')
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)


class AssociationReadingTests(unittest.TestCase):
    def status(self, score, rho, p):
        reader = getattr(module, 'nominal_association', None)
        self.assertTrue(callable(reader), 'The report needs a data-driven, direction-aware verdict')
        return reader(score, rho, p)

    def test_negative_affinity_signal_has_the_expected_direction(self):
        self.assertEqual(self.status('affinity_pred', -0.537, 0.0171), 'expected_direction')

    def test_positive_affinity_signal_has_the_opposite_direction(self):
        self.assertEqual(self.status('affinity_pred', 0.537, 0.0171), 'opposite_direction')

    def test_nonsignificant_published_subset_is_not_a_signal(self):
        self.assertEqual(self.status('affinity_pred', -0.543, 0.0864), 'not_detected')

    def test_vina_and_probability_scores_use_the_opposite_convention(self):
        for score in ('vina_neg', 'binary_bind_prob'):
            with self.subTest(score=score):
                self.assertEqual(self.status(score, 0.7, 0.01), 'expected_direction')
                self.assertEqual(self.status(score, -0.7, 0.01), 'opposite_direction')

    def test_report_preserves_the_full_set_signal_without_a_utility_claim(self):
        with tempfile.TemporaryDirectory() as td:
            results = Path(td)
            shutil.copyfile(ROOT / 'results/ext_ranking.csv', results / 'ext_ranking.csv')
            with patch.object(module, 'RESULTS', results), contextlib.redirect_stdout(io.StringIO()):
                module.main()
            report = (results / 'permutation.md').read_text(encoding='utf-8')
        self.assertIn('未校正、方向正确的关联信号', report)
        self.assertNotIn('两个方向都没有证据', report)
        self.assertIn('不能据此证明预测可用性', report)

    def test_zero_and_boundary_p_are_not_positive_evidence(self):
        self.assertEqual(self.status('affinity_pred', 0, 0.001), 'not_detected')
        self.assertEqual(self.status('affinity_pred', -0.7, 0.05), 'not_detected')


class PermutationBoundaryTests(unittest.TestCase):
    def test_all_720_permutations_include_equal_tail_values(self):
        # For these six distinct ranks, rho = (35 - sum squared rank differences)/35.
        # Count that integer numerator independently of the production correlation code.
        x = np.array([1, 3, 5, 4, 2, 6])
        y = np.array([6, 5, 3, 4, 2, 1])
        permutations = np.array(list(itertools.permutations(range(6))))
        observed = abs(35 - sum((int(a) - int(b)) ** 2 for a, b in zip(x, y)))
        extreme = sum(abs(35 - sum((int(a) - int(b)) ** 2 for a, b in zip(x, y[idx])))
                      >= observed for idx in permutations)
        self.assertEqual(extreme, 126)

        class ExhaustiveKeys:
            def random(self, shape):
                if shape != (720, 6):
                    raise AssertionError(shape)
                # argsort of these keys returns each index permutation exactly once.
                return np.argsort(permutations, axis=1)

        with patch.object(module, 'N_PERM', 720):
            result = module.perm_test(x, y, ExhaustiveKeys())
        # The function retains its Monte Carlo +1 rule even for this exhaustive fixture.
        self.assertEqual(result['p_two_sided'], 127 / 721)

    def test_published_active_vina_keeps_all_3502_extreme_shuffles(self):
        frame = module.pd.read_csv(ROOT / 'results/ext_ranking.csv')
        active = frame[frame['true_pkd'] > 5]
        result = module.perm_test(-active['vina_best'].to_numpy(float),
                                  active['true_pkd'].to_numpy(float),
                                  np.random.default_rng(42))
        self.assertEqual(result['p_two_sided'], 3503 / 20001)


if __name__ == '__main__':
    unittest.main()
