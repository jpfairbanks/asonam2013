"""Processes the temporal streams of histograms and dense vectors
into visualizations for the paper
"""
from __future__ import print_function
from multiprocessing import Process
#from threading import Thread as Process
import argparse
import numpy as np
import pandas as pd
import matplotlib
#matplotlib.use('Agg')
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
    """ Uses global parameters KERNEL_NAME, INIT_SAMPLE, END_SAMPLE, STRIDE
    DATA_DIR

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
    names_file = DATA_DIR+('betweenness_centrality.%d.csv'%END_SAMPLE)
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
    try:
        namesf = pd.read_csv(names_file, header=None)
        names = namesf[[0,1]].set_index(0)
    except:
        print("we failed to read names: does file %s exist?" % namesfile)
        names = None
    return df, names

def get_args():
    """
    This functions parses the options and returns them as an object
    use args.name to access the value of argument <name>

    add new options here
    """
    parser = argparse.ArgumentParser()

    parser.add_argument('-u', '--summary',
                        help='summarize each vertex accross time',
                        action='store_true')
    parser.add_argument('-s', '--static',
                        help='show static density estimates',
                        action='store_true')

    parser.add_argument('-t', '--temporal',
                        help='show how statistics change over time',
                        action='store_true')

    parser.add_argument( '--traces',
                         help='show traces for random vertices over time+derivative',
                        action='store_true')

    parser.add_argument('-d', '--derivative',
                        help='show an analysis of the derivatives of the data',
                        action='store_true')

    parser.add_argument('-r', '--correlation',
                        help='plot the function rho(t,k)',
                        action='store_true')

    parser.add_argument('--scatter',
                        help='show scatter plot when doing correlation analysis.',
                        action='store_true')

    parser.add_argument('-x','--crosstabs',
                        help='Create a latex table of the crosstabs',
                        action='store_true')

    args = parser.parse_args()
    return args

if __name__ == '__main__':
    args = get_args()
    FIGUREPATH = u'./figures/'
    FIGURE_EXTENSION = u'png'
    DATA_DIR         = u'../data/kernels/' #symlink this for portability
    POST_PROCESS_DIR = u'../data/post_process/'#symlink this for portability
    NSAMPLES = 100 #number of batches
    STRIDE = 1 #resolution of time in batches
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
    df, names = get_me_data()
    print('data loaded')
    timer.toc('loading data')
    t = 98
    lf = np.log1p(df)
    if args.summary:
        timer.tic('summary')
        print("summarizing vertices")
        whitf = ka.summarize_vertices(lf, pos_filter=True, whiten=True)
        whitf = whitf.join(names)
        whitf.to_csv(POST_PROCESS_DIR+'whitened_BC_stats.csv')
        covm  = whitf.cov()
        print(covm)
        print('a sample of the vertices as shown by the mean and std of their time series')
        fig, ax = plg.scatter_vertices(whitf, alpha=.3)
        timer.toc('summary')

    timer.tic('differentiate')
    diffframe = pd.DataFrame({t:(lf[t+STRIDE]-lf[t-STRIDE])/2
                              for t in lf.columns[1:-1]})
    timer.toc('differentiate')
    print('Beginning the Analysis')
    #main_out = main(df, t,timer)
    #show that we should use the median instead of the mean
    #====Parametric Analysis================#

    ##===========Show the CDF for the statistic at a fixed time========#
    if args.static:
        ##=======show density at fixed time t=============
        filt = lf[lf!=0]
        show_histogram_logbc(filt[t], t, median=True, fitter=stats.expon)
        #show_histogram_logbc(filt[t], t, median=False, fitter=stats.norm)
        print('starting CDF static')
        timer.tic('static')
        #staticf = lambda t: pf.cdf_plot_save(filt[[t,]],fitter=stats.norm)
        #staticf(t)
        filtmed = lf[lf>(filt.median())]
        staticf = lambda t: pf.cdf_plot_save(filtmed[[t,]],fitter=stats.expon)
        staticf(t)
        print('ending CDF static')
        timer.toc('static')

    if args.traces:
        timer.tic('traces')
        q=.500
        qmedians = df.quantile(q)
        topq = df[df > qmedians].dropna()
        tracef = lambda : pf.bc_traces(lf,diffframe, topq)
        tracef()
        timer.toc('traces')

    #======Correlation Analysis===========#
    if args.correlation:
        print('starting correlation analysis')
        timer.tic('correlation')
        tmpframe = df[df.columns[10::2]]
        #pearson is linear and spearman in any monotonic correlation
        corr_model(tmpframe,1,'pearson')
        #ax, rhoframe = corr_model(tmpframe, degree=1, method='pearson')
        #ax, rhoframe = corr_model(tmpframe, degree=2, method='spearman')
        timer.toc('correlation')
        print('ending correlation analysis')
    if args.scatter:
        timer.tic('scatter')
        print('starting scatter plot')
        scatter_frame = df[[t-2*STRIDE, t-STRIDE, t]]
        corr_plot(scatter_frame)
        #corr_plot(scatter_frame)
        timer.toc('scatter')
        print('ending scatter plot')

    if args.temporal:
        #show how the distribution changes over time
        print('beginning temporal analysis')
        print('we make two plots here one for all vertices as they show up')
        print('and one for only the initial vertex set')

        timer.tic('temporal')
        temporal_kwargs = {"title":"Distribution for all vertices"}
        tfilt = lf[lf>0]
        plg.distribution_describe(tfilt, colindex=None, plot=True,
                              transform=None, **temporal_kwargs)
        temporal_kwargs['title'] = "Distribution for initial vertex set"
        plg.distribution_describe(tfilt.dropna(), colindex=None, plot=True,
                              transform=None, **temporal_kwargs)
        timer.toc('temporal')
        print('ending temporal analysis')

    #======Analysis of Derivatives==============
    if args.derivative:
    #filt[t].hist(bins=BINCOUNT,normed=True)
        print('starting derivative analysis')
        timer.tic('derivative')
        frame = df#[df>df.median()]
        diffs = diffframe[t]
        seq = (diffs[(diffs)>0].dropna())
        seq_neg = (diffs[(diffs)<0].abs().dropna())
        diffr = pd.DataFrame({'pos':seq, 'neg':seq_neg})
        diffr = np.log1p(diffr)
        ##======Show the CDF for the diffs separating positive and negative====
        pf.cdf_plot_save_diffs(diffr,)
        ka.rank_sums_test(seq, seq_neg)
        ksp = stats.kstest(seq, 'norm')
        print(ksp)
        ##=====Show the density estimates for pos and neg separately
        show_histogram_diffs(diffr.pos,t, fitter=stats.beta, name='pos-beta')
        #show_histogram_diffs(pd.Series(stats.trimboth(seq.order(),.025,)), t,
        #                     fitter=stats.norm, name='pos-trimmed-norm')
        show_histogram_diffs(diffr.neg,t, fitter=stats.beta, name='neg-beta')
        timer.toc('derivative')
        print('ending derivative analysis')


    #=====Crosstabs for conditional probability
    if args.crosstabs:
        timer.tic('crosstabs')
        pf.save_crosstabs(df, t, STRIDE=STRIDE, eps=1)
        timer.toc('crosstabs')
    #TODO: Implement Isotonic regresion on the values
    print(timer.ends)
    print('total time: %f' % sum(timer.ends.values()))
    print('\n\tDONE')
