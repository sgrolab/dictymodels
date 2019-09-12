#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep  6 09:10:43 2019

@author: bgregor
"""

import os
import platform

def all_indices(a_array, b_array):
    ''' Take two incoming 1-D arrays. Generate
        a list of tuples containing indices to
        access all combinations of the values in 
        the two arrays.  For use with a grid search
        over the values spanned by the two arrays.
        
        Returns:  ((0,0),(0,1),...,(N,M)) '''
    indices=[None] * len(a_array) * len(b_array)
    i = 0
    # simple way with a pair of for loops.
    for a in range(len(a_array)):
        for b in range(len(b_array)):
            indices[i] = (a,b)
            i += 1
    return indices


def get_n_workers(verbose=False):
    ''' Get the assigned cores on the SCC and subtract 1 so one 
        core is reserved to manage the workers. If NSLOTS is not found
        return 1. '''
    # If on Linux check for NSLOTS
    if platform.system() == 'Linux':
        nslots = os.getenv('NSLOTS')
        if nslots is not None:
          nslots = int(nslots) 
          if nslots > 1:
              return nslots - 1
        # oops, job run with 1 core.
        if verbose:
            print('using 1 worker as NSLOTS=1 or is not defined.')
        return 1
    # Otherwise get the number of CPUs and subtract 1.
    if verbose:
        print('On Windows or Mac. Using all available cores.')
    return os.cpu_count() - 1

 