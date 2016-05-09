import collections
from collections import defaultdict
import numpy as np
import pandas as pd

__author__ = 'palmer'

def get_adducts(json_filename):
    import json
    config = json.loads(open(json_filename).read())
    target_adducts = [a['adduct'] for a in config['fdr']['pl_adducts']]
    decoy_adducts = list(set([a['adduct'] for a in config['isotope_generation']['adducts']])-set(target_adducts))
    return target_adducts, decoy_adducts


def calc_fdr_df(target_df, decoy_df, col='mult', ascending=False):
    # possible backward-incompatibility: calc_fdr_df previously returned the last argument as a series, not as an array
    return calc_fdr_arr(target_df[col], decoy_df[col], ascending=ascending)


def calc_fdr_arr(target_arr, decoy_arr, ascending=False):
    """
    Calculate the FDR curve for arrays of target scores and decoy scores.

    :param target_arr: sequence of scores from sum formulae with target adducts
    :param decoy_arr: sequence of scores from sum formulae with decoy adducts
    :param ascending: True if smaller scores are better, False if higher scores are better
    :raise TypeError: if the two sequences do not have the same length
    :return:
     fdr_curve: 1d array, twice as long as each of the input sequences, where the i-th element is the number of
     decoy hits divided by the number of target hits among the i highest scores (combined).

     target_hits: 1d array, twice as long as each of the input sequences, where the i-th element is the number of
     target hits among the i highest scores (combined)

     sorted_values: 1d array, twice as long as each of the input sequences, containing the sorted values from both
     input sequences
    """
    n, m = len(target_arr), len(decoy_arr)
    if n != m:
        raise TypeError('target should be same length as decoy {} {}'.format(n, m))
    ordering = 1 if ascending else -1  # reversed sorting if score is not ascending
    combined = np.concatenate((target_arr, decoy_arr))
    combined.sort()
    # count how often each value occurs
    target_bag, decoy_bag = _count(target_arr), _count(decoy_arr)
    unique_sorted = np.unique(combined)[::ordering]
    target_hits, decoy_hits = zip(*_iter_hits(target_bag, decoy_bag, unique_sorted))
    target_hits, decoy_hits = np.cumsum(target_hits), np.cumsum(decoy_hits)
    fdr_curve = decoy_hits / target_hits
    fdr_curve[target_hits == 0] = 0
    return fdr_curve, target_hits, combined[::ordering]


def _count(it):
    bag = defaultdict(int)
    for x in it:
        bag[x] += 1
    return bag


def _iter_hits(bag1, bag2, unique_sorted):
    # For each unique value, determine how often the value occurs in the target array and decoy array respectively
    # and yield the fractions (into the hits arrays).
    # Example: if there are two target scores and one decoy score with the same value, append 2/3, 2/3, 2/3 to the
    # target_hits array and 1/3, 1/3, 1/3 to the decoy_hits array.
    # This is to avoid the previously randomized approach of randomly assigning the value to one of the arrays.
    for x in unique_sorted:
        tc, dc = bag1[x], bag2[x]
        tc_part, dc_part = float(tc) / (tc + dc), float(dc) / (tc + dc)
        for _ in range(int(tc + dc)):
            yield (tc_part, dc_part)


def get_decoy_df(decoy_adducts, results_data, sf):
    import numpy as np
    import bisect
    import pandas as pd
    # decoy selection
    nDecoy = len(sf)
    decoy_adducts = np.random.choice(decoy_adducts, size=nDecoy, replace=True)
    dc_df = []
    for s, da in zip(sf, decoy_adducts):
        dc_df.append(results_data[da].iloc[bisect.bisect_left(sf, s)])
        dc_df = pd.DataFrame.from_records(dc_df)
    return dc_df


def select_passes(score_data, t_holds):
    import numpy as np
    moc_t, spec_t, spat_t = t_holds
    pass_tf = np.all([score_data['moc'] > moc_t, score_data['spat'] > spat_t, score_data['spec'] > spec_t], axis=0)
    pass_data = score_data[pass_tf]
    return pass_data


def select_passes_linear(score_data, t_hold):
    import numpy as np
    pass_tf = np.asarray(score_data['moc'] + score_data['spat'] + score_data['spec']) > t_hold
    pass_data = score_data[pass_tf]
    return pass_data


def select_passes_mult(score_data, t_hold):
    import numpy as np
    pass_tf = np.asarray(score_data['moc'] * score_data['spat'] * score_data['spec']) > t_hold
    pass_data = score_data[pass_tf]
    return pass_data


def select_passes_l2(score_data, t_hold):
    import numpy as np
    pass_tf = np.sqrt(np.asarray(score_data['moc'] ** 2 + score_data['spat'] ** 2 + score_data['spec'])) > t_hold
    pass_data = score_data[pass_tf]
    return pass_data


def count_adducts(pass_data, adducts):
    ## number of false hits by adduct
    fd_count = pass_data['adduct'].value_counts()  # count adduct occurances
    set_diff = set(adducts) - set(fd_count.keys())
    for s in set_diff:  # add zero occurance back in
        fd_count[s] = 0
    return fd_count.loc[[f in adducts for f in fd_count.keys()]]


def calc_fdr_adducts(fd_count, adducts, plausible_adducts, average='mean'):
    import numpy as np
    implausible_adducts = [a for a in adducts if a not in plausible_adducts]
    implausible_list = fd_count[implausible_adducts]
    implausible_hits = sum(implausible_list)
    if average == "mean":
        mean_implausible = float(implausible_hits) / len(implausible_list)
    elif average == 'median':
        mean_implausible = np.percentile(implausible_list, 50)
    else:
        raise ValueError('average type not recognised: {}'.format(average))
    plausible_hits = np.sum(fd_count[list(plausible_adducts)])
    if plausible_hits == 0:
        return 1, plausible_hits, implausible_hits
    fdr = (len(plausible_adducts) * mean_implausible) / plausible_hits
    return fdr, plausible_hits, implausible_hits


def is_fdr_curve(fdr_curve):
    """
    A FDR curve has these properties:

    - It starts at 0 and ends at 1
    - All points are greater than or equal to zero

    :param fdr_curve: the sequence to be tested
    :return: True if fdr_curve satisfies the conditions above, False otherwise
    """
    if len(fdr_curve) < 2:
        raise ValueError("FDR length < 2")
    if fdr_curve[0] != 0:
        raise ValueError("FDR start not 0")
    if not np.allclose(fdr_curve[-1], 1.0):
        raise ValueError("FDR end not 1")
    if any(fdr_curve[1:-1] < 0):
        raise ValueError("FDR values not all > 0")
    return True


def find_crossing(fdr_curve, fdr_target):
    """
    Find the index of the point before the rightmost crossing point between an FDR curve and a FDR target value.

    Formally speaking, given an array fdr_curve and a number fdr_target, find the smallest index i such that
    fdr_curve[j] >= fdr_target for all j > i

    :param fdr_curve: the FDR curve as a 1d array
    :param fdr_target: the target FDR
    :type fdr_target: float between 0 and 1
    :return: the index i or -1 if no crossing was found
    :raise ValueError: if fdr_curve is not a valid FDR curve
    """
    #if not is_fdr_curve(fdr_curve):
    #    raise ValueError("Not a valid FDR curve") #ADP - need to review is_fdr_curve criteria +noise means can start above 0
    if not 0 < fdr_target < 1:
        return -1

    less_zero_indices = np.where(fdr_curve <= fdr_target)[0]
    if len(less_zero_indices) == 0:
        return len(fdr_curve)-1
    i = less_zero_indices[-1]
    return i


def score_msm(score_data_df):
    return score_data_df['moc'] * score_data_df['spat'] * score_data_df['spec']


class decoy_adducts():
    #### Some class for dealing with FDR results ####
    def __init__(self, fname, target_adducts, decoy_adducts, shuf=None):
        self.target_adducts = target_adducts
        self.decoy_adducts = decoy_adducts
        self.shuf = shuf or _shuffle_rows
        # read in raw score file and calculate metabolite signal match
        with open(fname) as f_in:
            self.score_data_df = pd.read_csv(f_in, quotechar='"').fillna(0)
        self.score_data_df["msm"] = score_msm(self.score_data_df)
        self.score_data_df.sort_values(by="sf", inplace=True)
        self._count_formulas_per_adduct()

    def _count_formulas_per_adduct(self):
        # store some data info
        self.sf_l = {}
        self.n_sf = {}
        for a in self.target_adducts:
            self.sf_l[a] = np.unique(self.score_data_df.ix[self.score_data_df['adduct'] == a]["sf"])
            self.n_sf[a] = len(self.sf_l[a])

    def decoy_adducts_get_pass_list(self, fdr_target, n_reps, col='msm', return_decoy=False):
        # Get MSM threshold @ target fdr
        # Return molecules with higher MSM value
        msm_vals = self.get_msm_threshold_per_adduct(fdr_target, n_reps)
        pass_list = {}
        for a in self.target_adducts:
            target_df = self.score_data_df.ix[self.score_data_df["adduct"] == a]
            pass_list[a] = target_df.ix[target_df[col] > msm_vals[a]]['sf'].values
        if not return_decoy:
            return pass_list
        else:
            pass_list_decoy = {}
            for a in self.target_adducts:
                decoy_df = self.score_data_df.ix[~self.score_data_df["adduct"].isin(self.target_adducts)]
                pass_list_decoy[a] = decoy_df.ix[decoy_df[col] > msm_vals[a]][['sf', 'adduct']].values
            return pass_list, pass_list_decoy

    def get_msm_thresholds(self, adduct, fdr_target, n_reps=10, col='msm'):
        """
            Calculate the MSM crossing point at a given target fdr
        """
        fdr_curves, _, score_vects = self.get_fdr_curve(adduct, n_reps, col)
        msm_vals = []
        for fdr_curve, score_vect in zip(fdr_curves, score_vects):
            crossing_idx = find_crossing(fdr_curve, fdr_target)
            if crossing_idx > -1:
                msm_vals.append(score_vect[crossing_idx])
            else:
                msm_vals.append(0)
        return msm_vals

    def get_msm_threshold_per_adduct(self, fdr_target, n_reps=10, col='msm'):
        # Repeatedly calcualte FDR curves
        #   Find target crossing point -> find correspdoning msm score
        # return average score per adduct
        msm_vals = {}
        for a in self.target_adducts:
            msm_vals[a] = self.get_msm_thresholds(a, fdr_target, n_reps=n_reps, col=col)
            # calculate average
            msm_vals[a] = np.median(msm_vals[a])
        return msm_vals

    def get_fdr_curve(self, adduct, n_reps=10, col='msm'):
        # for a particular adduct, calcualte n_reps fdr curves
        target_df = self.score_data_df.ix[self.score_data_df["adduct"] == adduct]
        col_vector_decoy = self.score_data_df.ix[self.score_data_df['adduct'].isin(self.decoy_adducts) &
                                                 self.score_data_df['sf'].isin(self.sf_l[adduct])][col].values
        data_reps = len(col_vector_decoy) / len(self.sf_l[adduct])
        try:
            col_vector_decoy = col_vector_decoy.reshape((self.n_sf[adduct], data_reps))
        except ValueError as e:
            print len(col_vector_decoy) , len(self.sf_l[adduct]), np.shape(self.score_data_df), adduct, self.n_sf[adduct], data_reps
            raise
        self.shuf(col_vector_decoy)
        fdr_curves = []
        target_hits = []
        score_vects = []
        for col_vector in col_vector_decoy.T[:n_reps]:
            fdr_curve, target_hit, score_vect = calc_fdr_arr(target_df[col].values, col_vector, ascending=False)
            fdr_curves.append(fdr_curve)
            target_hits.append(target_hit)
            score_vects.append(score_vect)
        return fdr_curves, target_hits, score_vects


def _shuffle_rows(scores):
    # consume the iterator without using the values by "inserting" into an empty immutable queue;
    # see https://docs.python.org/2/library/itertools.html#recipes
    collections.deque((np.random.shuffle(row) for row in scores), maxlen=0)
