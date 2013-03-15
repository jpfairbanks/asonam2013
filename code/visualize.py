"""Processes the temporal streams of histograms and dense vectors
into visualizations for the paper
"""
from __future__ import print_function
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scipy.stats as stats
from time import time as time
import timedict as td
import util
import kernelio as kio
import kernel_analysis as ka
import plotting as plg
import paper_figures as pf
from kernelio import *
from plotting import *
from paper_figures import *

def lognormal_estimates_histogram(dframe):
    """Calculates the log-mean and log-standard deviation assuming that the data
    is log normal.

    Arguments:
    - `dframe`:
    """
    dframe = dframe[1:] #drop the 0 item because log(0) = -inf
    logs = np.log(dframe['bin'])
    weights = dframe['count']
    mu = np.inner(logs, weights)/np.sum(weights)
    centered = logs-mu
    sigma_sq = np.inner(centered, centered)/np.sum(weights)
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
    plt.plot(range(nsamples), np.sqrt(sigestimates), label=kernel_name+'sigma')
    plt.legend()
    plt.show()


def compare_kings_seq(rf, seq, start_col=0, end_col=-1,
                  ascending=False, plot=False, **kwargs):
    ''' Compare the entries at the top at the start to
    the entries at the top at the end.

    Uses start_col and end_col to define start and end.
    Arguments:
    - seq: use the top s vertices for s in seq
    - start_col: avoid burn in time
    - end_col: trim the slack at the end
    - ascending: reverse the sort
    - plot: display a figure of this ranging over seq
    - **kwargs: matplolib args

    '''
    #the rich stay richer
    ser = pd.Series(seq)
    targets_top_end = rf.sort(columns=[rf.columns[end_col]],
                              ascending=ascending)
    targets_top_start = rf.sort(columns=[rf.columns[start_col]],
                                ascending=ascending)
    similarity = lambda s: util.jaccard(targets_top_end[:s].index,
                                        targets_top_start[:s].index)
    ser = ser.apply(similarity)
    if plot:
        ser.plot(**kwargs)
    return ser

def plot_kernel(dframe, rows, kernel_name,
                figure_path, save=False, **kwargs):
    """ Makes a plot  of the kernel values for the rows over time

    Arguments:
    - `dframe`: the  DataFrame of kernel values at time steps in the columns and vertices as rows
    - `rows`: the set of vertices to display
    - `kernel_name`: the text that will be used in the title and filename
    - `figure_path`: where to save the figure ending in a slash
    - `save`: boolean indicating whether to save or not default to False
    - `**kwargs`: matplotlib keyword args passed to plot

    Returns:
    - `percs`: the kernel values restricted to the rows transposed for plotting
    """
    # TODO: normalize
    percs =  dframe.ix[rows].T
    ax = percs.plot(**kwargs)
    ax.set_title("Traces in %s" % kernel_name)
    ax.set_xlabel("batch number")
    ax.set_ylabel(kernel_name)
    if save:
        plt.savefig(figure_path+kernel_name+str(time())+'.'+FIGURE_EXTENSION)
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

def describe_vertex(df, lf, rf, v, plot=False, **kwargs):
    """ Gives a dataframe describing the vertex over time for quick inspection
    and the option to plot.
    """
    si = pd.DataFrame({'logval':lf.ix[v].dropna(), 'diff':lf.ix[v].dropna().diff(),
                       'rank':rf.ix[v].dropna()},
                      columns=['rank', 'logval', 'diff'])
    if plot:
        si.plot(subplots=True, **kwargs)
    return si

def distribution_describe(df, colindex=None, plot=False,
                          transform=None, **kwargs):
    """Makes a simple description using pandas.describe

    Arguments:
    - `df`:
    - `colindex`:
    - `plot`:

    """
    if transform is None:
        lf = df
    else:
        lf = transform(df)
    if colindex == None:
        colindex = df.columns
    fr = lf[colindex]
    desc = fr.describe().ix[[1, 2, 4, 5, 6]]
    if plot:
        desc.T.plot(**kwargs)
    return desc

def quantile_analysis(df, vertices, quants, timer=None, plot=False, **kwargs):
    """show how the quantiles of the distribution change over time

    Arguments:
    - `df`:
    - `vertices`:
    - `quants`:
    - `timer`:
    - `plot`
    - `**kwargs`: matplotlib kwargs
    """
    # Quantile analysis
    if timer:
        timer.tic('quants')
    quants = [.1, .25, .5, .75, .9, .95, .99, 1]
    tf = df.ix[vertices]
    qf = pd.DataFrame({q:tf.quantile(q) for q in quants})
    if plot:
        qf.plot(**kwargs)
    if timer:
        timer.toc('quants')
    return qf

def auto_correlate(df, columns='diff'):
    """makes a autocorrelation plot of the columns

    Arguments:
    - `df`:
    - `columns`:
    """
    # TODO: autocorrelation
    correlogram = pd.tools.plotting.autocorrelation_plot
    plt.figure()
    ax = correlogram(df[columns].dropna())
    return ax

def run_bc_analysis(df, timer):
    """

    Arguments:
    - `df`:
    """
    TARGETS_INITIAL = df[1].dropna().index #initial vertex set
    #count how many entries are not defined,
    #vertices who have not appeared yet or have fallen below epsilon
    num_nans = df.apply(np.isnan).sum()
    #print("The number of NANs by timestep:\n", num_nans)
    #compute the ranks for each vertex at each timethod='minm, na_option='top'e step
    timer.tic('ranking')
    rf = df.rank(method='min', na_option='top')
    timer.toc('ranking')

    # TODO: demonstrate that the leaders at the beginning are the leaders at the end
    compare_kings_seq(rf, [10, 30, 300], start_col=5)
    #the similarity of the top at begin and end is large and then falls very quickly
    # TODO: what happens to the people who start at the bottom.
    compare_kings_seq(rf, [100, 300, 3000], start_col=5,
                      ascending=True,)#the ylabel should be similarity
    #smallfish = df[df[21] < 0.1][21].index
    # TODO: can we find anyone making a smooth climb to the top

    timer.tic('density')
    distribution_describe(df, df.columns[::10], transform=np.log,
                          plot=True, title='changing distribution of log(bc)')
    timer.toc('density')
    # Tracing some vertices though time either random or seeded
    lf = np.log(df)
    if TARGETSV:
        targets = random_targets(lf, 1, NTARGETS, pool=TARGETSV, save=True)[0].columns
    else:
        # sample from the set of vertices that have ever been above a high quantile
        q=.500
        qmedians = df.quantile(q)
        topq = df[df > qmedians]
        targets = random_targets(lf, 1, NTARGETS, pool=topq.index)[0].columns

    numzeros = (rf[rf<=1].count())
    #numzeros.plot(label='rank($\epsilon)$')
    v = targets[0]
    lf = df.apply(np.log)
    summary = describe_vertex(df, lf, rf, v, plot=True,
                              title='profile of vertex %d'%v,)
    difftarg = lf.ix[targets].T.diff()
    diffax = difftarg.plot(title='diffs in log '+ KERNEL_NAME)
    diffax.set_xlabel('batch number')
    diffax.set_ylabel('change')
    return targets, summary

def derivative_analysis(lf, vertices, timer=None):
    """

    Arguments:
    - `lf`:
    - `vertices`:
    """

    ## DERIVATIVES of log BC
    if timer:
        timer.tic('diff')
    # show the derivative using logs to account for exponential dist of BC
    diffs = lf.T.diff().T

    # TODO: find the peaks of the derivative, these are the jumps in BC
    if timer:
        timer.toc('diff')

    # TODO: find significant jumps and count them
    # TODO: find the distribution of the jumps over time.
    #peaks in the differences
    peakslocs = diffs.replace(np.nan, -100).T.apply(np.argmax)*STRIDE
    ax = peakslocs.hist(bins=NSAMPLES/STRIDE,
                   normed=False,)
    ax.set_title('location of maximum change in bc')
    ax.set_ylabel('count')
    ax.set_xlabel('batch')
    # TODO: why do some batches have so many vertices with peaks there
    #david thinks these are connected component merges
    # TODO: measure slopes as a function distance from rank(e)
    # TODO: can we segment the lines into jumps and linear segments?
    # the difference between the slope of val and rank lines is the attrition rate
    #TODO: how many jumps does each vertex have
    return peakslocs


def correlation_analysis_from_top(sorted_frame, seq, corrmethod, ):
    """ Always compare the vertices that are highest at the end to the vertices
    that are highest at the current batch

    Arguments:
    - `sorted_frame`:
    - `seq`: the number of vertices to include in 'top'
    - `corrmethod`: ['spearman', 'pearson', 'kendall']
    """
    fun = lambda i: sorted_frame[:i].corr(method=corrmethod)
    return seq.map(fun)

def scatter_matrix_topp(sorted_frame, selected_axes, percentile=1):
    """

    Arguments:
    - `sorted_frame`:
    - `selected_axes`: the axes to include in .the scatterplot matrix
    - `percentile`:
    """
    pd.scatter_matrix(
        np.log(sorted_frame[selected_axes]+1)[:(percentile*len(sorted_frame)/100)]
        )

def exec_correlation_analysis(frame, selected_axes,  seq=None,
                              corrmethod=None,plot=False, **kwargs):
    """

    Arguments:
    - `frame`:
    - `seq`:
    - `plot`:
    - `corrmethod`:
    - `**kwargs`:
    """
    methods = ['spearman', 'pearson', 'kendall']
    if not (corrmethod in methods):
        corrmethod = methods[0]
    if seq==None:
        seq = pd.Index(range(10,2000,100))
    sorted_frame = frame.sort(columns=frame.columns[-1], ascending=False).dropna()
    print(sorted_frame)
    cseq = correlation_analysis_from_top(sorted_frame,
                                         seq, corrmethod)
    if plot:
        scatter_matrix_topp(sorted_frame, selected_axes, 10)
    return (cseq)


def load_mean(count):
    df = load_batches('/home/users/jfairbanks/Projects/morerad/'+ KERNEL_NAME+".%d.csv",
                      range(1,count), column=-1)
    return df, df.mean()

def main(df, t, timer=None):
    """

    Arguments:
    - `df`:
    - `t`: a reference column for all single time studies
    - `timer`:
    """
    lf = np.log(df)
    run_bc_analysis(df, timer)
    #peakslocs = derivative_analysis(df.apply(np.log), None, timer)
    #peakloc_counts = peakslocs.value_counts()
    #peakpeaks = peakslocs.value_counts()[:3]
    #cf = load_batches(DATA_DIR+'components.%d.csv', TIMEINDEX,)
    #cf[[421,431,441]].dropna(how='all')
    #trif = load_batches(DATA_DIR+ 'triangles'+".%d.csv", TIMEINDEX, column=-1)
    return True

def get_me_data():
    """

    Arguments:
    - `store_name`:
    - `frame_name`:
    - `filename`:
    """
    FILENAME = 'kernelframe.%s.%d.%d.%d.csv' %(KERNEL_NAME,
                                               INIT_SAMPLE, END_SAMPLE, STRIDE)
    store_name, frame_name = format_hdf_names(DATA_DIR,
                                              KERNEL_NAME,
                                              INIT_SAMPLE, END_SAMPLE, STRIDE)
    #get the data no matter what it takes, save it to hdf5 for future speed
    try:
        df = load_hdf_table(store_name, frame_name)
    except KeyError as e:
        print(e)
        print('we could not find the data in an hdfstore.')
        df = load_csv_frame(DATA_DIR, KERNEL_NAME, FILENAME, TIMEINDEX)
        try:
            write_hdf_table(store_name, frame_name, df)
        except Exception as ex:
            print(ex)
            print('we failed to write to hdf, is the hdf library installed?')
    return df
if __name__ == '__main__':
    FIGUREPATH = u'./figures/'
    FIGURE_EXTENSION = u'png'
    DATA_DIR = u'../data/kernels/' #symlink this for portability
    NSAMPLES = 1000 #number of batches
    STRIDE = 10 #resolution of time in batches
    INIT_SAMPLE = 1
    END_SAMPLE = NSAMPLES #TODO subtract one from this
    TIMEINDEX = pd.Index(np.arange(INIT_SAMPLE,END_SAMPLE,STRIDE))
    NTARGETS = 8 #number of series for the plots
    BINCOUNT = 50 #how many bins for histograms
    TARGETSV = []
    '''TARGETSV = [688773,
                756680,
                984640,
                1067645,
                3030528,
                3035516,
                ]'''
    #[3784858, 2357671, 2975930, 359724, 2124973, 3732925,]
    #vertices that I picked by hand
    KERNEL_NAME = "betweenness_centrality"
    KERNEL_NAMES = ['bc', 'communities', 'components', 'triangles']
    timer = td.timedict()

    #run_lognormal_analysis(DATA_DIR, NSAMPLES, KERNEL_NAME)

    timer.tic('loading data')
    # figure out where data should be found
    df = get_me_data()
    print('data loaded')
    timer.toc('loading data')
    t = 981
    lf = np.log(df)
    main_out = main(df, t,timer)
    #show that we should use the median instead of the mean
    print('comparing')
    #====Parametric Analysis================#
    ##=======show density at fixed time t=============
    #show_histogram_logbc(lf[t], t, median=True, fitter=stats.expon)
    #show_histogram_logbc(lf[t], t, median=True, fitter=stats.beta)
    print('filtering')
    ##===========Show the CDF for the statistic at a fixed time========#
    print('starting CDF static')
    filt = lf[lf>(lf.median())]
    pf.cdf_plot_save(filt[[t-(10*STRIDE),t]])
    print('ending CDF static')

    #======Analysis of Derivatives==============
    #filt[t].hist(bins=BINCOUNT,normed=True)
    print('starting derivative analysis')
    frame = df#[df>df.median()]
    diffs = (frame[t+STRIDE] - frame[t-STRIDE])/2 #use a centered derivative
    seq = np.log(diffs[(diffs)>0].dropna())
    seq_neg = np.log(diffs[(diffs)<0].abs().dropna())
    diffr = pd.DataFrame({'pos':seq, 'neg':seq_neg})
    ##======Show the CDF for the diffs separating positive and negative====
    pf.cdf_plot_save_diffs(diffr,)
    ka.rank_sums_test(seq, seq_neg)
    ksp = stats.kstest(seq, 'norm')
    print(ksp)
    ##=====Show the density estimates for pos and neg separately
    show_histogram_diffs(seq,t, fitter=stats.beta, name='pos-beta')
    #show_histogram_diffs(pd.Series(stats.trimboth(seq.order(),.025,)), t,
    #                     fitter=stats.norm, name='pos-trimmed-norm')
    show_histogram_diffs(seq_neg,t, fitter=stats.beta, name='neg-beta')
    #seq_combined = np.log(diffs[diffs!=0].abs()).dropna()
    #not useful because most of the vertices lose BC each round
    #show_histogram_diffs(seq_combined,t, fitter=stats.beta, name='neg-beta')
    print('ending derivative analysis')

    #show how the mean of the distribution changes over time


    #======Correlation Analysis===========#
    print('starting correlation analysis')
    corr_plot(df[[801,811,901]])
    tmpframe = df[df.columns[10::2]]
    #pearson is linear and spearman in any monotonic correlation
    ax, rhoframe = corr_model(tmpframe, degree=1, method='pearson')
    #ax, rhoframe = corr_model(tmpframe, degree=2, method='spearman')
    print('ending correlation analysis')

    #=====Crosstabs for conditional probability
    pf.save_crosstabs(df, t, STRIDE=STRIDE, eps=1)
