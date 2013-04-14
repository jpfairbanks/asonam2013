""" This module uses the plotting module (written by James Fairbanks) in order
to make the publication ready figures for the paper about this research. The
goal is to run a script that will generate all of the figures with titles,
captions, and axis labels. This will allow the maximum agility in publication.

These functions should just call functions from plotting.py with the propoer
arguments for the paper. Then they annotate the resulting plots and save them to disk
with the filenames that are in the latex document for includegraphics.

Any function that makes a plot should be in the plotting module written to be
reused and the specialized in this module.

"""

from time import time as time
import pandas as pd
import matplotlib
#if running without a display use this for generating to files
matplotlib.use('svg')
import matplotlib.pyplot as plt
import scipy.stats as stats
import numpy as np
import plotting as plg
from plotting import *
import kernel_analysis as ka
FIGUREPATH = u'./figures/'
FIGURE_EXTENSION = u'png'

def bc_traces(lf, diffframe, vpool):
    """
    Make and save a figure that shows the trace of 6 vertices over time.
    This is like a seismograph for BC.

    Writes out a plot using the current time to prevent over writing.
    Arguments:
    - `lf`: the values of the statistic
    - `diffframe`: the derivative frame
    - `topq`: the vertices to sample uniformly at random from

    Returns:
    - `targets`: the vertices that were selected.
    """

    kernel_name = 'betweenness centrality'
    traces_kwargs = {'title':"Traces in $\log(BC)$",}
    targets = plg.random_targets_trace(lf, 1, 6,
                                 pool=vpool.index, **traces_kwargs)[0].columns
    fig = plt.gcf()
    ax = fig.axes[0]
    ax.set_xlabel("batch number")
    ax.set_ylabel(kernel_name)
    plt.savefig(FIGUREPATH+'bc'+str(time())+'.'+FIGURE_EXTENSION)
    print(targets)
    traces_kwargs['title'] = "Traces in $d\log(BC)/dt$"
    plg.random_targets_trace(diffframe, 1, 6,
                       pool=targets, **traces_kwargs)
    fig = plt.gcf()
    ax = fig.axes[0]
    ax.set_title(traces_kwargs['title'])
    ax.set_xlabel("batch number")
    ax.set_ylabel(kernel_name)
    plt.savefig(FIGUREPATH+'bc'+str(time())+'.'+FIGURE_EXTENSION)
    return targets

def show_histogram_logbc(ls, t, median=False, fitter=stats.expon):
    """ Histogram some data and show a best fit distribution on top.

    Arguments:
    - `ls`:
    - `t`:
    - `median`:
    - `fitter`: the scipy.stats class to try and fit to the top half
    """
    plt.figure()
    if median:
        q = 0.5
    else:
        q = False
    show_histogram_parameteric_fit(ls, t, q, fitter)
    plt.title('Density of log(BC) at batch %d'%t)
    plt.xlabel('log(BC)')
    plt.ylabel('frequency')
    plt.legend()
    plt.savefig('%slogbc-histogram.%s'% (FIGUREPATH, FIGURE_EXTENSION))
    plt.show()

def cdf_plot_save(df, fitter=stats.norm):
    """
    Customize this function to chose colors, cdf vs survival func, kernel_name
    """
    kernel_name='logBC'
    fig, ax = cdf_plot(df, fitter=fitter,)
    ax.legend(loc='best')
    ax.set_title('CDF of %s(v)'%kernel_name)
    ax.set_xlabel('%s(v)'%kernel_name)
    ax.set_ylabel('CDF(%s(v))'%kernel_name)
    fig.savefig(FIGUREPATH+'cdf-logbc.%s'%FIGURE_EXTENSION,)
    return fig, ax

def cdf_plot_save_diffs(df):
    """
    Customize this function to chose colors, cdf vs survival func, kernel_name
    """
    kernel_name='|dlog(BC)/dt|'
    fig, ax = cdf_plot(df, fitter=stats.beta,)
    ax.legend(loc='best')
    ax.set_title('CDF of %s(v)'%kernel_name)
    ax.set_xlabel('log(%s(v))'%kernel_name)
    ax.set_ylabel('CDF(%s(v))'%kernel_name)
    fig.savefig(FIGUREPATH+'cdf-deriv-logbc.%s'%FIGURE_EXTENSION,)
    print('saved the right figure')
    return fig, ax


def show_histogram_diffs(seq,t, q=0, fitter=stats.norm, name=None):
    """ Histograms the diffs and show a best fit distribution on top.

    Arguments:
    - `seq`:
    - `t`:
    - `q`:
    - `fitter`:
    """
    plt.figure()
    show_histogram_parameteric_fit(seq, t, q, fitter)
    plt.title('Density of log(dBC/dt) at batch %d'%t)
    plt.xlabel('log(dBC/dt)')
    plt.ylabel('frequency')
    plt.legend()
    plt.savefig('%sdiff-histogram-%s.%s'% (FIGUREPATH, name, FIGURE_EXTENSION))
    plt.show()

def save_crosstabs(df, t, STRIDE=10, eps=1):
    """ Makes tables showing the probabilities of
    changes in ranks against changes in values.

    Arguments:
    - `df`:
    - `t`:

    """
    ct, margv, margr = ka.exec_crosstabs(df, [t, t+STRIDE], eps)
    fix_str = lambda s: s.replace('delta\_','$\Delta$')
    ctstr = fix_str(ct.to_latex())
    margvstr = fix_str(margv.to_latex())
    print(ctstr, margvstr)
    f = open('crosstab.tex', 'w')
    f.write(ctstr)
    f.close()
    f = open('margtab.tex', 'w')
    f.write(margvstr)
    f.close()
    return True

def corr_plot(df):
    """ Show a scatter plot between the first column and all of the other
    columns. This allows an analyst to see if there is a linear relationship
    between the values and how far that trend goes.

    Arguments:
    - `df`:

    """
    times = df.columns
    frame = df
    FILENAME = 'figures/correlation-scatter.%s' % FIGURE_EXTENSION
    titles = ['Correlation of Betweenness Centrality over time', '','']
    xlabels = ['logBC after batch %d']*2
    ylabels = ['logBC %d']*2
    print(xlabels)
    fig, axes = correlation_changes_over_time(df, times, color1='b', color2='b')
    for i,axis in enumerate(axes):
        axis.set_title(titles[i])
        axis.set_xlabel(xlabels[i] % times[0])
        axis.set_ylabel(ylabels[i] % times[i+1])
        axis.axvline(x=np.log(frame[times[0]].quantile(0.5)), color='k',
                    linestyle='-', label='median %d' % times[0])
        axis.axhline(y=np.log(frame[times[i+1]].quantile(0.5)),
                    color='k', linestyle='--',
                    label='median %d' % times[i+1])
        axis.legend()
    fig.savefig(FILENAME)
    return fig, axes

def corr_model(df, degree=1, method='pearson'):
    """ Make and save a figure showing how the correlation changes over time.
    Makes a fixed number of lines and only uses predictions into the future.

    The figure name uses the correlation method and the degree of the model.
    Arguments:
    - `df`:
    -`degree`: of the model to fit with
    -`method`: correlation method

    """
    FILENAME = 'figures/%s-correlation-deg%d-model%s' % (method,
                                                         degree,
                                                         FIGURE_EXTENSION)
    portion = len(df.columns)/2
    rhoframe = ka.rhotk(df, df.columns, df.columns[:portion:portion/4],
                               method=method)
    ax  = polyfit_plot(rhoframe, degree=degree, residuals=False,)
    ax.legend(ncol=2)
    ax.set_xlabel("$t+k$ batches")
    ax.set_ylabel("correlation")
    ax.set_title('correlation decays as gapsize increases')
    ax.figure.savefig(FILENAME)
    return ax, rhoframe

def rank_vals(df,tend=911, tstart=891, ):
    """

    Arguments:
    - `df`:
    """
    drf = (df[tend].rank()-df[tstart].rank())/(tend-tstart)
    ctdiff = (df[tend]-df[tstart])/(tend-tstart)
    table = pd.DataFrame({'dval':ctdiff,'drank':drf})
    ptab = table[table.dval>0]
    ntab = table[table.dval<0]
    fig, ax = plt.subplots(1,1,1)
    ax.scatter(np.log(-1 * ntab.dval), ntab.drank, color='r', alpha=.1)
    ax.scatter(np.log( 1 * ptab.dval), ptab.drank, color='b', alpha=.1)
    return fig, ax
