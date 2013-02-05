"""Processes the temporal streams of histograms into visualizations for the paper
"""
from __future__ import print_function
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scipy.stats as stats
from time import time as time
def load_data_histogram(pathfmt, batch):
    """Gets a batch into a series to be used in a the analysis

    Arguments:
    - `pathfmt`:
    - `batch`:
    """
    path = pathfmt.format(batch)
    dframe = pd.read_csv(path, columns=frozenset('bin','count'))
    return dframe

def load_sparse_vec(pathfmt, batch):
    """reads in a sparse vector that is represented as key value pairs.

    Arguments:
    - `pathfmt`:
    - `batch`:
    - `names`:
    - `'val']`:
    """
    path = pathfmt % (batch)
    dframe = pd.read_csv(path, index_col=0,names=[batch])
    return dframe

def load_batches(pathfmt, batches):
    """Load a set of batches into a big data frame

    Arguments:
    - `pathfmt`:
    - `batches`:
    """
    series = [load_sparse_vec(pathfmt, b) for b in batches]
    return pd.DataFrame.join(series[0], series[1:],how='outer')

def lognormal_estimates_histogram(dframe):
    """Calculates the log-mean and log-standard deviation assuming that the data
    is log normal.

    Arguments:
    - `dframe`:
    """
    dframe = dframe[1:] #drop the 0 item because log(0) = -inf
    logs = np.log(dframe['bin'])
    weights = dframe['count']
    mu = np.inner(logs, weights)/len(weights)
    centered = logs-mu
    sigma_sq = np.inner(centered, centered)/len(weights)
    return (mu, sigma_sq)

def run_lognormal_analysis(directory, nsamples, kernel_name):
    """load the data analyse it and make the plots.

    Arguments:
    - `directory`:
    - `nsamples`:
    - `kernel_name`:
    """
    muestimates = np.zeros((nsamples, 1), dtype=float)
    sigestimates = np.zeros((nsamples, 1), dtype=float)

    for i in range(0, nsamples):
        DF = load_data_histogram(directory + kernel_name + ".%d.csv", i)
        print(str(i)+":"+ str(len(DF)))
        est = lognormal_estimates_histogram(DF)
        muestimates[i] = est[0]
        sigestimates[i] = est[1]

    plt.plot(range(nsamples), muestimates, label=kernel_name+'mu')
    #plt.plot(range(nsamples), np.exp(muestimates), label=kernel_name+'exp(mu)')
    #plt.plot(range(nsamples), muestimates * muestimates, label=kernel_name+'mu^2')
    plt.plot(range(nsamples), np.sqrt(sigestimates), label=kernel_name+'sigma')
    plt.legend()
    plt.show()

def normalize(seq):
    m = np.max(seq)
    return seq/m

def percentiles(sequence, queries):
    """returns the percentile of each element in the sequence

    Arguments:
    - `sequence`:
    - `queries`:
    """
    sequence = np.sort(sequence)
    return [stats.percentileofscore(sequence, sequence[d]) for d in queries]

def plot_kernel(dframe, targets, kernel_name,
                figure_path, save=False, **kwargs):
    """ Makes a plot  of the kernel values for the targets over time

    Arguments:
    - `dframe`: the  DataFrame of kernel values at time steps in the columns and vertices as rows
    - `targets`: the set of vertices to display
    - `kernel_name`: the text that will be used in the title and filename
    - `figure_path`: where to save the figure ending in a slash
    - `save`: boolean indicating whether to save or not default to False
    - `**kwargs`: matplotlib keyword args passed to plot

    Returns:
    - `percs`: the kernel values restricted to the targets transposed for plotting
    """
    # TODO: normalize
    percs =  dframe.ix[targets].T
    ax = percs.plot(kwargs)
    ax.set_title(kernel_name+" @mentions by vertex")
    ax.set_xlabel("batch number")
    ax.set_ylabel("percentile of "+kernel_name)
    if save:
        plt.savefig(figure_path+kernel_name+str(time())+'.png')
    return percs

def random_targets(dframe, nplots, nseries, **kwargs):
    """ Make nplots figures that each display nseries samples from the data
        in dframe
    """
    targetcollection  = [np.random.permutation(dframe.index)[:nseries]
                         for i in range(nplots)]
    data = [plot_kernel(dframe, t, KERNEL_NAME, FIGUREPATH, kwargs)
            for t in targetcollection]
    return data

DATA_DIR = u'/shared/users/jfairbanks/sandyfull/'
FIGUREPATH = u'/shared/users/jfairbanks/smisc.sandystudy/output/'
NSAMPLES = 1000
STRIDE = 10
NTARGETS = 8
TARGETSV = [328,457,5056,340807,12859,12860]
KERNEL_NAME = "bc"
KERNEL_NAMES = ['bc', 'communities', 'components']
#run_lognormal_analysis(DATA_DIR, NSAMPLES, KERNEL_NAME)
df = load_batches(DATA_DIR+ KERNEL_NAME+".%d.vec", range(1,NSAMPLES,STRIDE))
#count how many entries are not defined,
#vertices who have not  appeared yet or have fallen below epsilon
num_nans = df.apply(np.isnan).sum()
print(num_nans)

#replace nans  with 0 to fix sorting order
rf = df.rank(ascending=False)
print (rf)
signal = rf.ix[TARGETSV]
#print(signal)
#signal.T.apply(normalize).plot()

#show that the path to the end is jagged and jumpy
randframes = random_targets(rf, 6, 8)
# TODO: quantiles

#demonstrate that the leaders  at the  beginning are the leaders at the end.
#the rich get richer
# targets_top_end = rf.sort(columns=[rf.columns[-1]], ascending=False)[:10].index
# targets_top_start = rf.sort(columns=[rf.columns[5]], ascending=False)[:10].index
# print(targets_top_end)
# toppe = plot_kernel(rf, list(targets_top_end), KERNEL_NAME, FIGUREPATH)
# topps = plot_kernel(rf, list(targets_top_start), KERNEL_NAME, FIGUREPATH)

# TODO: sample from the set of vertices that have ever been above median
#rankf = df.rank().T.ix[targetsv]
plt.show()
