#import kernel_analysis as ka
import kernelio as kio
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

DATA_DIR =u'../data/kernels/'
TIMEINDEX=pd.Index(np.arange(1,1000,10))
ccf = kio.load_csv_frame(DATA_DIR, 'local_clustering_coefficients',
                         'localcc.1.1000.10.csv', TIMEINDEX)
#how many vertices change their CC in each direction over time
vals = np.sign(ccf.T.diff().T).apply(pd.value_counts).T
axes = vals[[-1,1]].plot(title='change in Local Clustering Coefficient over time')
axes.set_xlabel('batch')
axes.set_ylabel('number of vertices')
