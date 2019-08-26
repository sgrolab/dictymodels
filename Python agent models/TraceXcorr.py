# -*- coding: utf-8 -*-
"""
Created on Wed Jul 24 13:45:46 2019

@author: Chuqiao
"""
import operator as op
from functools import reduce
import numpy as np

def ncr(n, r):
    r = min(r, n-r)
    numer = reduce(op.mul, range(n, n-r, -1), 1)
    denom = reduce(op.mul, range(1, r+1), 1)
    return numer / denom

def TracesXcorr(SingleTraces):
    TraceLen,NumTraces = SingleTraces.shape
    r=[] #np.zeros(int(ncr(NumTraces,2)))
    for i in range(NumTraces):
        CurrTrace=SingleTraces[:,i]
        k=i+1
        while k<= (NumTraces-1):
            FollowTrace=SingleTraces[:,k]
            R = np.corrcoef(CurrTrace,FollowTrace)
            r.append(R[1,0])
            k=k+1
    
    MeanR = np.mean(r)
    return MeanR


