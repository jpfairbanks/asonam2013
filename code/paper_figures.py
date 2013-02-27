import pandas as pd
import matplotlib.pyplot as plt
import scipy.stats as stats
import numpy as np
from plotting import *
FIGUREPATH = u'./figures/'
FIGURE_EXTENSION = u'png'
def show_histogram_logbc(ls, t, median=False, fitter=stats.expon):
    """

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
    """

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
