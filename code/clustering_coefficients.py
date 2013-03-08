#import kernel_analysis as ka
import kernelio as kio
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
def count_change_directions(df, eps=None):
    """Count the number of vertices that increase, stay constant or decrease

    Arguments:
    - `df`: the frame you care about columns are time rows are vertices
    - `eps`: threshold to use, not implemented yet
    """
    if eps is not None:
        print('\nthresholding not implemented yet. We are using threshold of 0')
    return np.sign(df.T.diff().T).apply(pd.value_counts).T

def change_counting_plot(change_frame, ):
    """ Plots the number of vertices changing in each direction

    Arguments:
    - `change_frame`: use count_change_directions to compute this

    Returns:
    - `axes`: the axes that contains the plot
    """
    vals = change_frame[[1,-1]]
    axes = vals.plot(title='change in Local Clustering Coefficient over time')
    axes.set_xlabel('batch number')
    axes.set_ylabel('number of vertices')
    return axes
DATA_DIR =u'/scratch/jfairbanks/sandy_triangles10x/'
TIMEINDEX=pd.Index(np.arange(1,100,1))
ccf = kio.load_csv_frame(DATA_DIR, 'local_clustering_coefficients',
                         'localcc.1.100.10.csv', TIMEINDEX)
trif = kio.load_csv_frame(DATA_DIR, 'triangles',
                         'triangles.1.100.10.csv', TIMEINDEX)
#how many vertices change their CC in each direction over time
change_counting_plot(count_change_directions(df,))
fig = plt.figure()
np.log1p(trif[99]).order().plot(fig=fig)
fig = plt.figure()
globalcc = ccf.mean()
globalcc.plot(fig=fig)

"""bcf = kio.load_hdf_table(*kio.format_hdf_names(DATA_DIR,
                                               'betweenness_centrality',
                                               1,1000, 10))
df = ccf.join(bcf, how='inner', lsuffix='cc', rsuffix='bc')
"""
