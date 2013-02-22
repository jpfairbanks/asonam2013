import numpy as np
import scipy as sp
import matplotlib.pyplot as plt

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
