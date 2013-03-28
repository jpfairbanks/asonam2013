import pandas as pd
import numpy as np
import scipy as sp
import scipy.stats as stats
import matplotlib.mlab as mlab
#frame = bcf.join(oldbcf[[4]], how='outer', lsuffix='new', rsuffix='old')

#TODO: add an analysis for clustering vertices based on there kernel features
#TODO: try and find journalists based on this
#TODO: PCA the deg, cc, bc data for some batch
def summarize_vertices_many(df):
    """ Try and extract some features from the data by summarizing each vertex 
    according to some reduction.

    Arguments:
    - `df`: Vertices by timesteps data frame
    
    Returns:
    - `mframe`: Vertices by features dataframe

    """
    dmag, dcnt  = np.square(df).sum(axis=1), df.count(axis=1)
    dlcnt = np.log1p(df.count(axis=1))
    dsqvar, dvar = np.square(df).var(axis=1), df.var(axis=1)
    dsqmean = np.square(df).mean(axis=1)
    dmean = df.mean(axis=1)
    mframe = pd.DataFrame({'sqsum':dmag,'count':dcnt,'logcount':dlcnt,'sqvar':dsqvar,
                           'var':dvar,'sqmean':dsqmean,'mean':dmean})
    return mframe

def summarize_vertices(df, pos_filter=False, whiten=False):
    """ Compute the mean and std deviation for each vertex accross time and optionally whiten it.
    That is to subtract off the mean and divide by the standard deviation
    Arguments:
    - `df`: the DataFrame columns are time rows are vertices
    - `pos_filter`: only use vertices with positive mean
    - `whiten`: boolean defaults to false
    Returns:
    - `whitf`: columns are 'mu', 'sigma' rows are vertices

    Note: pos_filter is useful for things graph statistics that return only nonnegative values.
    We don't want to include vertices that have a 0 for every time step. They will only skew our 
    statistics.
    """
    sumf = pd.DataFrame({'mu':df.mean(axis=1), 'sigma':df.std(axis=1),
	                 #'kurt':lf.kurt(axis=1),
			 }).dropna()
    if pos_filter:
        posf  = sumf[sumf['mu']>0]
    else:
        posf  = sumf
    if whiten:
        centf = posf-posf.mean()
        whitf = centf/centf.std()
    else:
        whitf = posf
    return whitf

def deriv(frame):
    """Return the centered derivative of the frame.

    Arguments:
    - `df`: the frame you care about columns are time rows are vertices

    """
    ds = dict()
    indices = range(1,len(frame.columns)-1)
    #print(indices)
    #print(frame.columns)
    for i in indices:
        ds[i] = (frame[frame.columns[i+1]] - frame[frame.columns[i-1]])/2.0
    return pd.DataFrame(ds)

def flatten(frame):
    """Flattens a frame into a big sequence.
    does not make a hierarchical index or anything.
    You can call hist on the result in order to see what all of the data in a frame look like.
    This came about because I was trying to histogram a distribution and wanted to consider
    the data from all time steps together.

    Arguments:
    - `frame`: A dataframe that has columns of time steps and rows of entities
    - `returns`: The flattened series of all data in one pool

    """
    ser = pd.concat(diffframe[i] for i in diffframe.columns)
    return ser


def count_change_directions(df, eps=None):
    """Count the number of vertices that increase, stay constant or decrease
    df is expected to be items in the rows and time steps in the columns

    Arguments:
    - `df`: the derivative of the frame you care about columns are time rows are vertices
    - `eps`: threshold to use, not implemented yet
    """
    if eps is not None:
        print('\nthresholding not implemented yet. We are using threshold of 0')
    return np.sign(df).apply(pd.value_counts).T

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
    rf = df[sample_times].corr(method=method)
    rhos = pd.DataFrame(np.tril(rf),
                        index=rf.index,columns=rf.columns)[display_starts]
    rhos = rhos.replace(0,np.nan)
    return rhos

# # Hypothesis Testing
# ====================

# TODO: use this
def rank_sums_test(treatment1, treatment2):
    """ See if the distribution of treatmen1 is different than
    the distribution treatment2

    Arguments:
    - `treatment1`:
    - `treatment2`:
    """
    z_stat, p_val = stats.ranksums(treatment1, treatment2)
    print "Mann-Whitney-Wilcoxon RankSum P for treatments 1 and 2 =", p_val
    return z_stat, p_val

# # Dimmensionality Reduction
# =========================

def pca(df):
    """

    Arguments:
    - `df`:
    """
    pika = mlab.PCA(df.dropna())
    pj = np.array([pika.project(df.T[r]) for r in df.T])
    pjf = pd.DataFrame(pj)
    return pjf

# # Comparing Rank and Value
#   ============================
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
