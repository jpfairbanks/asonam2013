"""Clustering coefficients analyzes both the number of triangles and the 
clustering coefficient which is the number of triangles divided by degree squared 
for each vertex.
"""

import kernelio as kio
import kernel_analysis as ka
import plotting as plg
import paper_figures as pf
import numpy as np
import scipy as scipy
import scipy.stats as stats
import pandas as pd
import matplotlib.pyplot as plt
import argparse
# Input Output
DATA_DIR      = u'/scratch/jfairbanks/tri_bc_sandy/'
POST_PROC_DIR = u'/scratch/jfairbanks/tri_bc_sandy/post_processing/'
FIGUREPATH = u'./figures/'
TIMEINDEX=pd.Index(np.arange(1,100,1))
t = 99 #the time point of interest
def load_data(data_dir, post_proc_dir, timeindex, timepoint):
    """
    Arguments:
    =========
    - `data_dir` : where the data lives
    - `post_proc_dir` : where to put my output/cache
    - `timeindex`: the indices to load
    - `timepoint`: the time stamp of interest
    """
    ccf    = kio.load_csv_frame(DATA_DIR, 'local_clustering_coefficients',
                             POST_PROC_DIR+'localcc.1.100.10.csv', TIMEINDEX)
    trif   = kio.load_csv_frame(DATA_DIR, 'triangles',
                             POST_PROC_DIR+'triangles.1.100.10.csv', TIMEINDEX)
    namesf = pd.read_csv(DATA_DIR+'triangles.100.csv', header=None)
    names = namesf[[0,1]].set_index(0)
    return trif, ccf, namesf
# ========
def model_func():
    """ summarize each vertex into a NV x Properties frame
    Calling scatter plot on the return value will show some patterns in the data.
    """
    sumframe = ka.summarize_vertices_many(ccf)[['mean','var',]]#'logcount']]
    nzd = diffframe[diffframe!=0]
    mframe = ka.summarize_vertices_many(nzd)[['count', 'mean','var']]
    #mframe['count'] = (mframe['count']+0.0).replace(0.0,np.nan).dropna()
    frame = sumframe.join(mframe, lsuffix='value',rsuffix='deriv',how='outer')
    print(frame.head())
    return frame

# Plotting
# ========
def count_changes_plot(change_frame, ):
    """ Plots the number of vertices changing in each direction

    Arguments:
    - `change_frame`: use count_change_directions to compute this

    Returns:
    - `axes`: the axes that contains the plot
    """
    vals = change_frame[[1,-1]]
    ax = plt.axes()
    ax.scatter(vals[-1], vals[1])
    ax.set_title('batches by counting positive and negative derivatives')
    ax.set_xlabel('number of decreasing vertices')
    ax.set_ylabel('number of increasing vertices')
    #ratios = ((vals[1]+0.0)/vals[-1])
    #ratios.plot()
    axes = vals.plot(title='change in Local Clustering Coefficient over time')
    axes.set_xlabel('batch number')
    axes.set_ylabel('number of vertices')
    fig = plt.gcf()
    fig.savefig(FIGUREPATH+'tri-deriv-sign.png')
    return axes

def save_globalcc_over_time(ccf, title=""):
    """

    Arguments:
    - `ccf`:
    """
    fig, axes = plt.subplots(1,1,1)
    globalcc = ccf.mean()
    axes.plot(globalcc.index, globalcc)
    axes.set_ylabel("mean Clustering Coefficient")
    axes.set_xlabel("time (batches)")
    axes.set_title(title)
    return fig, axes

def save_cdf(cdf, title=""):
    """

    Arguments:
    - `cdf`:
    """
    fig, axes = plt.subplots(1,1,1)
    axes.plot(cdf, np.arange(cdf.count()))
    axes.set_ylabel("number of vertices")
    axes.set_xlabel("log(triangles)")
    aces.set_title(title)
    return fig, axes

def get_args():
    """
    This functions parses the options and returns them as an object
    use args.name to access the value of argument <name>

    add new options here
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('-s', '--static', help='show static density estimates',
                        action='store_true')
    parser.add_argument('-t', '--temporal', help='show how statistics change over time',
                        action='store_true')
    parser.add_argument('-u', '--summary',
                        help='summarize each vertex accross time, write into the post_process dir',
                        action='store_true')
    parser.add_argument('-d', '--derivative', help='show an analysis of the derivatives of the data',
                        action='store_true')
    parser.add_argument('-m', '--model', help='extract some features from the derivative of CC and plot them',
                        action='store_true')
    parser.add_argument('-r', '--correlation', help='plot the function rho(t,k)',
                        action='store_true')
    parser.add_argument('--scatter',
                        help='show scatter matrix of the data',
                        action='store_true')
    parser.add_argument('--pca',
                        help='show scatter matrix using  principal components of the data',
                        action='store_true')
    args = parser.parse_args()
    return args

if __name__ == '__main__':
    args = get_args()
    t = 99
    #get me some data
    trif, ccf, namesf, = load_data(DATA_DIR,
                                   POST_PROC_DIR, TIMEINDEX, t)

    logtrif = np.log1p(trif)
    #degf = kio.load_csv_frame(DATA_DIRm 'degree',
    #'degree.1.100.10', TIMEINDEX)

    # Static analyses
    if args.static:
        #density estimation static of number of triangles
        plg.cdf_plot(logtrif[[t]], fitter=stats.expon)
        #density estimation of CC
        plg.cdf_plot(ccf[[t]], fitter=stats.gamma)


    # Temporal analysis
    if args.temporal:
        filtitle = 'For vertices with at least one triangle'
        save_globalcc_over_time(ccf[ccf>0], title=filtitle)
        save_globalcc_over_time(ccf, title='All vertices')

    #summarize vertices over time
    if args.summary:
        #timer.tic('summary')
        print("summarizing vertices")
        whitf = ka.summarize_vertices(logtrif, pos_filter=True, whiten=True)
        whitf = whitf.join(namesf)
        whitf.to_csv(POST_PROC_DIR+'whitened_tri_stats.csv')
        covm  = whitf.cov()
        print(covm)
        print('a sample of the vertices as shown by the mean and std of their time series')
        #fig, ax = plg.scatter_vertices(whitf, alpha=.3)
        #timer.toc('summary')

    if args.correlation:
    # Correlation analysis
        sample_times = ccf.columns
        display_starts = sample_times[:len(sample_times)/2:10]
        #rhoframe = ka.rhotk(ccf, sample_times, display_starts, method='pearson')
        ax, rhoframe = pf.corr_model(ccf,degree=1, method='pearson')
        print(rhoframe)
        rhoframe.plot(title="Rho_cc(t,k)")
    # Differences Analysis
    diffframe = ka.deriv(ccf)
    if args.derivative:
        #how many vertices change their CC in each direction over time
        change_frame = ka.count_change_directions(ka.deriv(ccf))
        count_changes_plot(change_frame)
        diffst = (ccf[t]-ccf[t-2])/2
        #plt.figure()
        #diffframe[diffframe>0].mean(axis=0).plot()
        #diffframe[diffframe<0].mean(axis=0).plot()
        #ser = ka.flatten(diffframe)
        #ser.hist()
        print("we only model vertices that have a change bigger than 0.01")
        nzdiffst = diffst[diffst.abs() > 0.01].dropna()
        nzdiffst.plot(kind='kde')
        fig, axes = plt.subplots(1,1,1)
        plg.show_histogram_parameteric_fit(nzdiffst, t-1,)
        plg.cdf_plot(pd.DataFrame(nzdiffst))

    if args.model:

        mframe = model_func()
        maxes = np.log(trif[trif>0].max(axis=1)).dropna()
        maxes.name = 'maxtri'
        jf = mframe.join(maxes)
        maxcc = ccf.max(axis=1).replace(0,np.nan).dropna()
        maxcc.name = 'maxcc'
        jf = jf.join(maxcc)
        jf.save("ccfeatures.data")
        #pd.scatter_matrix(jf, alpha=.3,color='g')

    inframe = pd.DataFrame({'tri': trif[t],
                            'cc' : ccf[t]})
    if args.scatter:
        pd.scatter_matrix(inframe)
    if args.pca:
        pjf = ka.pca(inframe)
        pd.scatter_matrix(pjf)
