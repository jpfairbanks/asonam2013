import pandas as pd
import numpy as np
import scipy as sp
import scipy.stats as stats
#frame = bcf.join(oldbcf[[4]], how='outer', lsuffix='new', rsuffix='old')

#TODO: add an analysis for clustering vertices based on there kernel features
#TODO: try and find journalists based on this
#TODO: PCA the deg, cc, bc data for some batch

def crosstabs(value_series, rank_series, epsilon):
    """ uses a +1 to avoid dividing by zero later

    Arguments:
    - `value_series`:
    - `rank_series`:
    - `epsilon`: the threshold for significance
    """
    mask = value_series.abs().dropna() > epsilon
    sigvals = mask*value_series
    ct = pd.crosstab(np.sign(sigvals), np.sign(rank_series), margins=False) + 1
    return ct

def exec_crosstabs(df,timestamps,epsilon):
    """

    Arguments:
    - `df`:
    - `timestamps`:
    - `epsilon`:
    """
    frame = df[timestamps]
    last_col = timestamps[-1]
    valsframe = frame
    vals = valsframe[last_col] - valsframe[timestamps[0]]
    ranks = frame.rank(method='min',
                       na_option='top').T.diff().T[last_col]
    ranks.name = 'delta_rank'
    vals.name = 'delta_val'
    ct = crosstabs(vals, ranks, epsilon)
    marg_vals = ((ct.T+0.0)/ct.sum(axis=1)).T
    print('P(deltaRank|deltaValue):')
    print(marg_vals)
    marg_ranks = (ct+0.0)/ct.sum(axis=0)
    return ct, marg_vals, marg_ranks
