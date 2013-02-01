"""Processes the temporal streams of histograms into visualizations for the paper
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def load_data(pathfmt, batch):
    """Gets a batch into a series to be used in a the analysis

    Arguments:
    - `pathfmt`:
    - `batch`:
    """
    path = pathfmt % batch
    dframe = pd.read_csv(path, names=['bin', 'count'])
    return dframe

def lognormal_estimates_histogram(dframe):
    """Calculates the log-mean and log-standard deviation assuming that the data
    is log normal.

    Arguments:
    - `dframe`:
    """
    dframe = dframe[1:] #drop the 0 item because log(0) = -inf
    logs = np.log(dframe['bin'])
    weights = dframe['count']
    mu = np.inner(logs, weights)/len(weights)
    centered = logs-mu
    sigma_sq = np.inner(centered, centered)/len(weights)
    return (mu, sigma_sq)

DATA_DIR = u'/shared/users/jfairbanks/smisc.sandystudy/output/temporal/sandy.results/'
nsamples = 1000
muestimates = np.zeros((nsamples, 1), dtype=float)
sigestimates = np.zeros((nsamples, 1), dtype=float)
kernel_name = "bc"
kernel_names = ['bc', 'communities', 'components']
for i in range(0, nsamples):
    DF = load_data(DATA_DIR + kernel_name + ".%d.csv", i)
    print(str(i)+":"+ str(len(DF)))
    est = lognormal_estimates_histogram(DF)
    muestimates[i] = est[0]
    sigestimates[i] = est[1]

plt.plot(range(nsamples), muestimates, label=kernel_name+'mu')
#plt.plot(range(nsamples), np.exp(muestimates), label=kernel_name+'exp(mu)')
#plt.plot(range(nsamples), muestimates * muestimates, label=kernel_name+'mu^2')
plt.plot(range(nsamples), np.sqrt(sigestimates), label=kernel_name+'sigma')
plt.legend()
plt.show()
