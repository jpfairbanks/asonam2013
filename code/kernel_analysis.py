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

    Returns:
    - `ct`: the crosstabulation dataframe

    """
    mask = value_series.abs().dropna() > epsilon
    sigvals = mask*value_series
    ct = pd.crosstab(np.sign(sigvals), np.sign(rank_series), margins=False) + 1
    return ct

def exec_crosstabs(df,timestamps,epsilon):
    """ use the crosstabs function to make the tables corresponding to the joint
    probability distribution of changes in value and changes in rank and the
    conditional distribution of change in rank given change in value. This is
    helps one understand the relationship between these changes. The conditional
    probability tables tells us what the probability of the change in rank
    occuring once we know what the change in value was.

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

def rhotk(df,sample_times, display_starts, method='pearson'):
    """ compute the correlation going forward. This corresponds to the
    rho(t,t+k) for t+k in sample_times. We only return the frame containing
    columns for each t in display_starts.

    Arguments:
    - `df`: the data from must contain sample_times as columns
    - `sample_times`: the points in time at which to compute the correlations
    - `display_starts`: must be a subset of sample_times
    - `method`: pearson, spearman, kendall. pearson is the fastest and default

    Returns:
    -`rhos`: a dataframe containing a column for each starting point and the
             rows are all t+k in sample_times. You can call plot on this to
             visualize it quickly

    Note:
    pearson correlation is the fastest but only captures linear correlation.
    Spearman's correlation is a rank correlation and so is for any monotonic
    relationship.
    """
    rhos = df[sample_times].corr(method=method)[display_starts]
    return rhos
