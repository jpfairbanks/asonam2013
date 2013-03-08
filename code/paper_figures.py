""" This module uses the plotting module (written by James Fairbanks) in order
to make the publication ready figures for the paper about this research. The
goal is to run a script that will generate all of the figures with titles,
captions, and axis labels. This will allow the maximum agility in publication.

These functions should just call functions from plotting.py with the propoer
arguments for the paper. Then annotate the resulting plots and save them to disk
with the filenames that are in the latex document for includegraphics.

Any function that makes a plot should be in the plotting module written to be
reused and the specialized in this module.


"""
import pandas as pd
import matplotlib
#if running without a display use this for generating to files
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import scipy.stats as stats
import numpy as np
from plotting import *
import kernel_analysis as ka
FIGUREPATH = u'./figures/'
FIGURE_EXTENSION = u'png'

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
        q = .5
    else:
        q = False
    show_histogram_parameteric_fit(ls, t, q, fitter)
    plt.title('Density of log(BC) at batch %d'%t)
    plt.xlabel('log(BC)')
    plt.ylabel('frequency')
    plt.legend()
    plt.savefig('%slogbc-histogram.%s'% (FIGUREPATH, FIGURE_EXTENSION))
    plt.show()

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


def corr_plot(df):
    """

    Arguments:
    - `df`:

    """
    FILENAME = 'figures/correlation-scatter.png'
    titles = ['Correlation of Betweenness Centrality over time', '','']
    xlabels = ['logBC after batch %d']*2
    ylabels = ['logBC %d']*2
    print(xlabels)
    times = [801,811,901]
    frame = df[times]
    fig, axes = correlation_changes_over_time(frame, times, color1='b', color2='b')
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
    FILENAME = 'figures/%s-correlation-deg%d-model.png' % (method, degree)
    portion = len(df.columns)/2
    rhoframe = ka.rhotk(df, df.columns, df.columns[:portion:portion/4],
                               method=method)
    ax  = polyfit_plot(rhoframe, degree=degree, residuals=False,)
    ax.legend(ncol=2)
    ax.set_xlabel("$t+k$ batches")
    ax.set_ylabel("correlation")
    ax.set_title('correlation decays quadratically in gapsize')
    ax.figure.savefig(FILENAME)
    return ax, rhoframe
