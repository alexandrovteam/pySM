from unittest import TestCase

import numpy as np
import pandas as pd
from numpy.testing import assert_array_equal, assert_array_almost_equal

from fdr_measures import find_crossing, calc_fdr_df, calc_fdr_arr, decoy_adducts


class CalcFdrDfTest(TestCase):
    def test_calc_fdr_df(self):
        test_cases = (
            {
                'target': [0.9, 0.92, 0.97],
                'decoy': [0.1, 0.2, 0.5],
                'ascending': False,
                'fdr_curve': [0, 0, 0, 1. / 3, 2. / 3, 1],
                'target_hits': [1, 2, 3, 3, 3, 3],
                'score_vect': [0.97, 0.92, 0.9, 0.5, 0.2, 0.1]
            }, {
                'target': [0.9, 0.92, 0.97, 0.99, 0.99, 0.99, 0.99],
                'decoy': [0.1, 0.2, 0.5, 0.3, 0.3, 0.3, 0.3],
                'ascending': False,
                'fdr_curve': [0, 0, 0, 0, 0, 0, 0, 1. / 7, 2. / 7, 3. / 7, 4. / 7, 5. / 7, 6. / 7, 1],
                'target_hits': [1, 2, 3, 4, 5, 6, 7, 7, 7, 7, 7, 7, 7, 7],
                'score_vect': [0.99, 0.99, 0.99, 0.99, 0.97, 0.92, 0.9, 0.5, 0.3, 0.3, 0.3, 0.3, 0.2, 0.1]
            }, {
                'target': [0.9, 0.92, 0.97, 0, 0, 0, 0],
                'decoy': [0.1, 0.2, 0.5, 0, 0, 0, 0],
                'ascending': False,
                'fdr_curve': [0, 0, 0, 1./3, 2./3, 1, 1, 1, 1, 1, 1, 1, 1, 1],
                'target_hits': [1, 2, 3, 3, 3, 3, 3.5, 4, 4.5, 5, 5.5, 6, 6.5, 7],
                'score_vect': [0.97, 0.92, 0.9, 0.5, 0.2, 0.1, 0, 0, 0, 0, 0, 0, 0, 0]
            },)
        for case in test_cases:
            target_arr, decoy_arr = case['target'], case['decoy']
            ascending = case['ascending']
            target_df = pd.DataFrame({'msm': target_arr})
            decoy_df = pd.DataFrame({'msm': decoy_arr})

            np.random.seed(0)
            arr_result = calc_fdr_arr(target_arr, decoy_arr, ascending=ascending)
            np.random.seed(0)
            df_result = calc_fdr_df(target_df, decoy_df, col='msm', ascending=ascending)

            assert_array_almost_equal(arr_result[0], case['fdr_curve'])
            assert_array_almost_equal(df_result[0], case['fdr_curve'])
            assert_array_equal(arr_result[1], case['target_hits'])
            assert_array_equal(df_result[1], case['target_hits'])
            assert_array_equal(arr_result[2], case['score_vect'])
            assert_array_equal(df_result[2], case['score_vect'])


class FindCrossingTest(TestCase):
    def test_valid_inputs(self):
        test_cases = (
            (([0, 1], 0.5), 0),
            (([0, 0, 0.8, 1], 0.85), 2),
            (([0, 0, 0.8, 1], 0.8), 1),
            (([0, 0.7, 0.8, 0.8, 0.8, 0.8, 0.8, 1], 0.8), 1),
            (([0, 0.2, 0.8, 1], 0.1), 0),
            (([0, 0, 0.8, 1], 0), -1),
            (([0, 0, 0.8, 1], 1), -1),
            (([0, 0, 0.8, 1], 1.1), -1),
            (([0, 0, 0.8, 1], -0.1), -1),
            (([0, 1, 0.5, 0.7, 0.2, 1], 0.6), 4),
        )
        for (fdr_curve, fdr_target), expected in test_cases:
            fdr_curve = np.asarray(fdr_curve)
            self.assertEqual(find_crossing(fdr_curve, fdr_target), expected)

    def test_invalid_inputs(self):
        inputs = (
            ([], 0.5),
            ([0], 0.5),
            ([1], 0.5),
            ([0, 0.99], 0.5),
            ([0, 0.1, 1.01], 0.5),
            ([0.1, 0.2, 1], 0.5),
            ([-0.1, 0.2, 1], 0.5),
        )
        for args in inputs:
            self.assertRaises(ValueError, find_crossing, *args)
