# -*- coding: utf-8 -*-
"""
Created on Wed Apr 22 15:28:00 2020

@author: cqhuyan
"""
import numpy as np

def import_npz(npz_file, namespace):
    Data = np.load(npz_file)
    for varName in Data:
        namespace[varName] = Data[varName]