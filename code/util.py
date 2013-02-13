import numpy as np
import scipy as sp
import matplotlib.pyplot as plt

def jaccard(A, B):
    num = (len(A.intersection(B))+0.0)
    denom = (len(A.union(B))+0.0)
    return num/denom
