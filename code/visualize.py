"""Processes the temporal streams of histograms into visualizations for the paper
"""
from __future__ import print_function
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scipy.stats as stats
from time import time as time
import timedict as td

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
    frame = pd.DataFrame.join(series[0], series[1:],how='outer')
    frame.save
    return frame

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
    ax = percs.plot(**kwargs)
    ax.set_title(kernel_name+" @mentions by vertex")
    ax.set_xlabel("batch number")
    ax.set_ylabel("percentile of "+kernel_name)
    if save:
        plt.savefig(figure_path+kernel_name+str(time())+'.png')
    return percs

def random_targets(dframe, nplots, nseries, pool=None, **kwargs):
    """ Make nplots figures that each display nseries samples from the data
        in dframe

        draws from pool or the index of the dframe if no  pool is given
    """
    if pool==None:
        pool = dframe.index
    targetcollection  = [np.random.permutation(pool)[:nseries]
                         for i in range(nplots)]
    [t.sort() for t in targetcollection]
    data = [plot_kernel(dframe, t, KERNEL_NAME, FIGUREPATH, **kwargs)
            for t in targetcollection]
    return data

def describe_vertex(df,lf,rf,v,plot=False):
    """
    """
    si = pd.DataFrame({'val': df.ix[v].dropna(),
                       'diff':df.ix[v].dropna().diff(),
                       'ratio':lf.ix[v].dropna().diff(),
                       'rank':rf.ix[v].dropna()})
    if plot:
        si.plot(subplots=True)
    return si

DATA_DIR = u'/shared/users/jfairbanks/sandyfull/'
FIGUREPATH = u'/shared/users/jfairbanks/smisc.sandystudy/output/'
NSAMPLES = 1000 #number of batches
STRIDE = 10 #resolution of time in batches
NTARGETS = 8 #number of series for the plots
BINCOUNT = 50 #how many bins for histograms
TARGETSV = [328,457,5056,340807,12859,12860] #vertices that I picked by hand
KERNEL_NAME = "bc"
KERNEL_NAMES = ['bc', 'communities', 'components']
timer = td.timedict()
#run_lognormal_analysis(DATA_DIR, NSAMPLES, KERNEL_NAME)
timer.tic('loading data')
df = load_batches(DATA_DIR+ KERNEL_NAME+".%d.vec", range(1,NSAMPLES,STRIDE))
timer.toc('loading data')
#count how many entries are not defined,
#vertices who have not appeared yet or have fallen below epsilon
num_nans = df.apply(np.isnan).sum()
print("The number of NANs by timestep:\n", num_nans)
#initial vertex set
TARGETSV = df[1].dropna().index
#replace nans  with 0 to fix s#orting order
timer.tic('replacing nans')
#df = df.replace(np.nan, 0)
timer.toc('replacing nans')

#compute the ranks for each vertex at each timethod='minm, na_option='top'e step
timer.tic('ranking')
rf = df.rank(method='min', na_option='top')
timer.toc('ranking')
signal = rf.ix[TARGETSV]
#print(signal)
#signal.T.apply(normalize).plot()

# TODO: show that the path to the end is jagged and jumpy
#randframes = random_targets(rf, 1, 6)

# Quantile analysis
timer.tic('quants')
quants = [.1,.25,.5,.75,.9,.95,.99,1]
tf = df.ix[TARGETSV]
#use both the initial vertices and the entire vertex set for quantiles
#qf = pd.DataFrame({q:tf.quantile(q) for q in quants})
# TODO: come back to the quantile analysis
#qf = pd.DataFrame({q:df.quantile(q) for q in quants})
#qf.plot(logy=True)
timer.toc('quants')
# TODO: demonstrate that the leaders  at the  beginning are the leaders at the end.
#the rich get richer
# targets_top_end = rf.sort(columns=[rf.columns[-1]], ascending=False)[:10].index
# targets_top_start = rf.sort(columns=[rf.columns[5]], ascending=False)[:10].index
# print(targets_top_end)
# toppe = plot_kernel(rf, list(targets_top_end), KERNEL_NAME, FIGUREPATH)
# topps = plot_kernel(rf, list(targets_top_start), KERNEL_NAME, FIGUREPATH)
# TODO: what happens to the people who start at the bottom.
smallfish = df[df[21] <.1][21].index
# TODO: sample from the set of vertices that have ever been above a high quantile
q=.500
qmedians = df.quantile(q)
topq = df[df > qmedians]
targets = random_targets(rf, 1, 6, pool=topq.index)[0].columns

numzeros = (rf[rf<=1].count())
numzeros.plot(label='rank($\epsilon)$')
nnz = rf[rf>1].count()
timer.tic('diff')
# show the derivative using logs to account for exponential dist of BC
lf = df.apply(np.log)
diffs = lf.T.diff().T
difftarg = diffs.ix[targets].T
lftarg = lf.ix[targets].T
difftarg.plot(title='diffs in log '+ KERNEL_NAME)
plt.figure()
#TODO: what do the diffs look like en masse
timer.tic('density')
diffsums = diffs.sum().dropna()
ax = diffsums.hist(bins=BINCOUNT, normed=True)
diffsums.plot(kind='kde',
              title='distribution of total change in log bc',)
#diffs.unstack().dropna().plot(kind='kde', title='kde of diffs flattened')
timer.toc('density')
# percent change is kind of useless
#lftarg.T.pct_change().T.plot(title='percent change of diffs in log'+ KERNEL_NAME)
# TODO: find the peaks of the derivative, these are the jumps in BC
w = lf.ix[targets].T
pos = w[w.diff()>0]
timer.toc('diff')
v = targets[0]
summary = describe_vertex(df,lf,rf,v,plot=True)
print(summary.head(40))
#TODO: how many jumps does each vertex have
# TODO: find the distribution of the jumps over time.
#peaks in the differences
plt.figure()
peakslocs = diffs.replace(np.nan,-1).T.apply(np.argmax)*STRIDE
peakslocs.hist(bins=NSAMPLES/STRIDE)
# TODO: measure slopes as a function distance from rank(e)
print(timer)
plt.show()
# TODO: can we segment the lines into jumps and linear segments?
# the difference between the slope of val and rank lines is the attrition rate
# TODO: autocorrelation
correlogram = pd.tools.plotting.autocorrelation_plot
plt.axes()
correlogram(summary['diff'].dropna())
plt.show()
