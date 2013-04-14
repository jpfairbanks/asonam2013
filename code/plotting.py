"""This module is for functions that create plots for analysing the data from
frunning a kernel on a graph. The module paper_figures uses these functions for
making figures with titles labels and captions for the paper about this research

"""

import pandas as pd
import matplotlib.pyplot as plt
import scipy.stats as stats
import numpy as np

BINCOUNT = 50
# #Traces
# ============================================
def plot_kernel_traces(dframe, rows, **kwargs):
    """ Makes a plot  of the kernel values for the rows over time

    Arguments:
    - `dframe`: the  DataFrame of kernel values at time steps in the columns and vertices as rows
    - `rows`: the set of vertices to display
    - `**kwargs`: matplotlib keyword args passed to plot

    Returns:
    - `percs`: the kernel values restricted to the rows transposed for plotting
    """
    # TODO: normalize
    percs =  dframe.ix[rows].T
    ax = percs.plot(**kwargs)
    return percs


def random_targets_trace(dframe, nplots, nseries, pool=None, **kwargs):
    """ Make nplots figures that each display nseries samples from the data
        in dframe

        draws from pool or the index of the dframe if no  pool is given
    """
    if pool==None:
        pool = dframe.index
    targetcollection  = [np.random.permutation(pool)[:nseries]
                         for i in range(nplots)]
    [t.sort() for t in targetcollection]
    data = [plot_kernel_traces(dframe, t, **kwargs)
            for t in targetcollection]
    return data


def show_PCA(df):
    pca = mlab.PCA(df)
    plt.plot(pca.Y)
    return pca

# # Density Estimation
# ================================================================

# ## Nonparametric
# ----------------
def distribution_describe(df, colindex=None, plot=False,
                          transform=None, **kwargs):
    """Makes a simple description using pandas.describe
    If df is DataFrame with vertices in rows and time steps in columns
    then this will give a rough picture of how the distribution is changing over
    time.

    Arguments:
    - `df`:
    - `colindex`:
    - `plot`:

    Returns:
    - `desc`: the description frame

    Note:
    You can select elements from the description frame that is return if you
    would rather manipulate the data by hand.
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

# ## Parametric
# -------------
def cdf_plot(df, fitter=stats.norm,
             cdf=True, colors = ['b','g','r','c','m','y','k','w'],):
    """

    Arguments:
    - `df`:
    - `fitter`: a statistical distribution with a fit method
    - `cdf`: defaults to True, showing the CDF version.
             If false use the survival function.
    """
    ords = map(lambda s: df[s].dropna().order(ascending=cdf), df.columns)
    names = [s.name for s in ords]
    oframes = map(lambda s: pd.DataFrame({'x':s,
                                            'CDF(x)':s.rank(ascending=cdf)/s.count()
                                            }).set_index('x'),
                  ords)
    fig, ax = plt.subplots(1,1,1)
    figs = [ax.plot(seq.index, seq, color=col, label='%s empirical'%name)
            for seq, col, name in zip(oframes, colors, names)]
    #using a model if one is suggested to us.
    if fitter is not None:
        models = map(lambda s: fitter(*fitter.fit(s)), ords)
        domains = map(lambda s: np.linspace(s.min(),s.max(), 1000), ords)
        if cdf:
            yvals = [model.cdf(domain) for model, domain in zip(models,domains)]
        else:
            yvals = [model.sf(domain) for model, domain in zip(models,domains)]
        for x, y, col, name  in zip(domains, yvals, colors, names):
            ax.plot(x, y, label='%s model' %name,
                    color=col, linestyle='--')
        for model in models:
            print(model.args)
    return fig, ax

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
        plt.axvline(x=seq.quantile(quantile), color='k',
                    label='quantile %.2f' % quantile)
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


# # Correlation
#======================================================
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
                         s=10, alpha=.7, c=color_mask[:,i])
    return fig, axes



def polyfit_plot(frame, degree=1, residuals=True,):
    """ Fit a dataframe and return the polynomials then show them on a plot.

    Arguments:
    - `frame`: the data you want to fit
    - `degree`: of the polynomial
    - `residuals`: do you want to residuals defaults to True
    """

    paramframe = pd.DataFrame({s:np.polyfit(y=frame[s].ix[s::],
                                            x=frame[s].ix[s::].index,deg=degree)
                               for s in frame})
    modelframe = pd.DataFrame({s:pd.Series(np.polyval(p=paramframe[s],
                                                      x=frame[s].ix[s::].index),
                                           index=frame[s].ix[s::].index)
                               for s in frame})
    ax = modelframe.plot(style='--', legend=False)
    frame.plot(ax=ax, style='+-', legend=False)
    ax.legend(ncol=2)
    if residuals:
        resids = (modelframe-frame)
        rax = resids.plot(kind='bar')
    return ax

# # Putting vertices into Feature Space
def scatter_vertices(df, alpha=.3):
    """ make a scatter plot whose data elements are vertices and whos axes
    are summary statistics of kernel values over time.
    See kernel_analysis.summararize_vertices for a description of arguments.
    """
    plt.scatter(df[df.columns[0]], df[df.columns[1]], alpha=alpha)
    fig = plt.gcf()
    ax = fig.axes[0]
    ax.set_xlabel(df.columns[0])
    ax.set_ylabel(df.columns[1])
    return fig, ax

# # To see if we can apply a model
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
