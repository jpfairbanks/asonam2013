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
