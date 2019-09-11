# -*- coding: utf-8 -*-
"""
Created on Sat Dec 08 15:51:50 2018

@author: Alex
"""

from __future__ import division
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import copy
from multiprocessing import Pool
from functools import partial
from activeGelFramework import Grid

T=10
dt=.01
size=1

grid=Grid()