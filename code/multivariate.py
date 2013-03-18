import kernelio as kio
import kernel_analysis as ka
import plotting as plg
import numpy as np
import scipy as scipy
import scipy.stats as stats
import pandas as pd
import matplotlib.pyplot as plt
import argparse
# Input Output
DATA_DIR      = u'/scratch/jfairbanks/tri_bc_sandy/'
POST_PROC_DIR = u'/scratch/jfairbanks/tri_bc_sandy/post_processing/'
TIMEINDEX=pd.Index(np.arange(1,100,1))
t = 99 #the time point of interest
def multivariate_load(data_dir, post_proc_dir, timeindex, timepoint):
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
    bcf    = kio.load_csv_frame(DATA_DIR, 'betweenness_centrality',
                             POST_PROC_DIR+'bc.1.100.10.csv', TIMEINDEX)
    namesf = pd.read_csv(DATA_DIR+'triangles.100.csv', header=None)
    names = namesf[[0,1]].set_index(0)
    names[1].name = 'name'
    triseq      = trif[t]
    ccseq       = ccf[t]
    bcseq       =  bcf[t]
    triseq.name = 'tri'
    ccseq.name  = 'cc'
    bcseq.name  = 'bc'
    bcseq.dtype = np.float
    tf = names.join(triseq)
    tf = tf.join(ccseq)
    tf = tf.join(bcseq)
    return trif, ccf, bcf, namesf, tf

# Analysis
# ========

# Plotting
# ========


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
    parser.add_argument('-d', '--derivative', help='show an analysis of the derivatives of the data',
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
    trif, ccf, bcf, namesf, tf = multivariate_load(DATA_DIR,
                                                   POST_PROC_DIR, TIMEINDEX, t)
    #
    #tf = tf[tf!=0].dropna()
    print(tf.head())
    bctf = tf[['bc','tri']].ix[::10]
    lf = np.log1p(bctf)
    inframe = lf.join(tf['cc'].ix[::10])
    inframe = inframe[inframe>0]
    if args.scatter:
        pd.scatter_matrix(inframe)
    if args.pca:
        pjf = ka.pca(inframe)
        pd.scatter_matrix(pjf)
