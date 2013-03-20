import numpy as np
import scipy as sp
import matplotlib.pyplot as plt

def mvnrand(covariance,sample_size=10000, mean=None):
    """	
    Draw from a multivariate gaussian preserving the axes labels on the 
    covariance matrix.
    Arguments:
    - `covariance`:
    - `sample_size`:
    - `mean`:
    """
    if mean is None:
        mean = np.zeros(covariance.shape[0])
        data = rand.multivariate_normal(mean, covariance, sample_size)
    return pd.DataFrame(data,columns=covariance.columns)


def jaccard(A, B):
    num = (len(A.intersection(B))+0.0)
    denom = (len(A.union(B))+0.0)
    return num/denom

def normalize(seq, prob=False):
    ''' Takes a sequence and normalizes it

    Arguments:
    - seq
    - prob: if true then it will normalize to sum to 1
    '''
    if prob:
        return seq/np.sum(seq)
    else:
        m = np.max(seq)
        return seq/m
