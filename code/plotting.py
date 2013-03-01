"""This module is for functions that create plots for analysing the data from
frunning a kernel on a graph. The module paper_figures uses these functions for
making figures with titles labels and captions for the paper about this research

"""

import pandas as pd
import matplotlib.pyplot as plt
import scipy.stats as stats
import numpy as np

BINCOUNT = 50

def show_histogram_parameteric_fit(seq, t, quantile=0, fitter=stats.norm):
    """ Show the histogram of a sequence along with a parametric fit. Allows for
    filtering using a quantile in case the fit only applies to the tail of the
    distibution.

    Arguments:
    - `seq`:
    - `t`:
    - `fitter`:
    """
    seq.hist(bins=BINCOUNT, normed=True)
    scaling = 1
    if quantile:
        plt.axvline(x=seq.quantile(quantile), color='k',label='median')
        filtered = seq[seq>seq.quantile(quantile)]
        scaling = 1-quantile
    else:
        filtered = seq
    if fitter is not None:
        params = fitter.fit(filtered)
        rv = fitter(*params)
        r = (filtered.min(), filtered.max())
        domain = np.arange(*r,step=(r[1]-r[0])/1000)
        pdf = pd.Series(rv.pdf(domain)*scaling, index=domain, name='pdf')
        pdf.plot(color='r')

def correlation_changes_over_time(df, times, log=True,
                                  condition=np.median, color1='b', color2='r'
                                  ):
    """
    Show how we have a strong correlation between two adjacent batches but a
    weaker correlation as the time goes further away. The condition and colors
    are to show that the correlation is better when condition=True and worse
    when it is False. greater than condition is color1. The colors will blend
    for records that change their conditional

    Arguments:
    - `df`:
    - `times`: the batches to show
    - `log`:default=True
    - `condition`:=np.median
    - `color1`:='b'
    - `color2`:='r'

    """

    frame = df[times]
    if log:
        frame = np.log(frame)
    func = lambda s: condition(s.dropna())
    color_mask = np.where(frame > frame.apply(func), color1, color2)
    fig, axes = plt.subplots(len(times)-1,1)
    for i,t in enumerate(frame.columns[1:]):
        axis = axes[i]
        axis.scatter(x=frame[frame.columns[0]], y=frame[t],
                    s=10, alpha=.25, c=color_mask[:,0])
        if color2:
            axis.scatter(x=frame[frame.columns[0]], y=frame[t],
                         s=10, alpha=.25, c=color_mask[:,i])
    return fig, axes
