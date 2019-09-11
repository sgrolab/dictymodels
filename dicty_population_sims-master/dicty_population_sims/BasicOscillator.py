#import stuff
from __future__ import division
import numpy as np

def makeActivatorDynamics(nonlinearity='Cubic', externalInput='Log'):
    #Define activator dynamics
    Adot = lambda A, R, cAMPe, params: (A-(A**3)/3 - R + params['a']*np.log(1+cAMPe/params['Kd'])/params['tA'])/params['tA']
    return Adot

def makeRepressorDynamics():
    Rdot = lambda A, R, params: ((A-params['gamma']*R)+params['c0'])/params['tR']
    return Rdot

def addNoise():
    noise = lambda dt, params, NDisc: params['sigma']*np.random.normal(loc = 0.0, scale = np.sqrt(dt), size = NDisc)
    return noise