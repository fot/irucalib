# irudefs.py

import numpy as np 
from pylab import *
from math import *

def irusigns(Umat, vectors):
    """function to sign codes for an array of vector rotations
       input  Umat     : array (numaxes,3) of iru rotation axes
              vectors  : array (3,num) of rotation vectors
       output signcode : array (num,) of sign codes for each rotation vector
    """
    signs = np.dot(Umat, vectors)
    signs = (signs >= 0)
    signcode = np.dot(np.array([1, 2, 4, 8]), signs)
    return signcode
# Maneuver signs and codes
# ---- 0    ++++ 15
# ---+ 1    +++- 14
# --+- 2    ++-+ 13
# --++ 3    ++-- 12
# -+-- 4    +-++ 11
# -+-+ 5    +-+- 10
# -++- 6    +--+  9
# -+++ 7    +---  8

def irurates(Dmat, Mmat, Gmat, SFac, Bias4, Bias3, ratecnts):
    """function to compute iru 3-vector rate
       input  Dmat : 3x3 array of coefficients to adjust on-board Mmat
              Mmat : 3x3 array of on-board misalignment/scale-adjust matrix
              Gmat : 3x4 array to combine IRU channels to obtain 3-vector
              SFac : 2x4 array of scale factors (radians/count) for each sign of eah channel
              Bias4 : 4-vector bias (counts/sec)
              Bias3 : 3-vector bias (radians/sec)
              ratecnts : 5xnum array of time & count rate (counts/sec) for each channel at each time
       output
              chanrate :
              angrate  :
    """
    chanrate = ratecnts.copy()
    Bias4 = Bias4.copy()
    Bias3 = Bias3.copy()
    if Bias4.ndim == 1:
        Bias4 = reshape(Bias4,(4,1))
    if Bias3.ndim == 1:
        Bias3 = reshape(Bias3,(3,1))
    
    I3 = np.eye(3)
    ProdMat = dot(dot(I3 + Dmat, I3 + Mmat), Gmat)
    
    chanrate[1:, :] = chanrate[1:,:] - Bias4
    chanrate[1, :] = chanrate[1, :] * ((chanrate[1, :] > 0.0) * SFac[0, 1] + (chanrate[1, :] < 0.0) * SFac[1, 1])
    chanrate[2, :] = chanrate[2, :] * ((chanrate[2, :] > 0.0) * SFac[0, 2] + (chanrate[2, :] < 0.0) * SFac[1, 2])
    chanrate[3, :] = chanrate[3, :] * ((chanrate[3, :] > 0.0) * SFac[0, 3] + (chanrate[3, :] < 0.0) * SFac[1, 3])
    chanrate[4, :] = chanrate[4, :] * ((chanrate[4, :] > 0.0) * SFac[0, 4] + (chanrate[4, :] < 0.0) * SFac[1, 4])
    
    angrate = zeros((4,chanrate.shape[1]))
    angrate[0, :] = chanrate[0, :]
    angrate[1:,:] = dot(ProdMat, chanrate[1:,:]) - Bias3
    
    return (chanrate, angrate)


def irucounts(accumcnts):
    """function to compute delta iru channel counts, count rate, and
       adjust accumulated angle for roll over
       input  accumcnts : array(0 to numchan, 0 to numtime)
                          row 0, time in sec for numtime number of time steps (numtime+1 values)
                          rows 1 to numchan, accumulated angle counts for each channel, 1 to numchan
       output accumcnts : array(0 to numchan, 0 to numtime)
                          row 0, time in sec
                          rows 1 to numchan, accumulated angle counts for each channel, 1 to numchan
                          adjusted for rollover, first value zero
              deltacnts : row 0, delta time in sec for numtime number of time steps, numtime+1 values
                          rows 1-numchan, change of counts for each channel, 1 to numchan
              ratecnts  : row 0, time in sec for numtime-1 number of time steps (numtime values)
                          rows 1 to numchan, count rate (cnts/sec) for each channel, 1 to numchan
    """
    # change in iru counts per channel for each time step with first values zero
    accumcnts = accumcnts.copy()
    deltacnts = np.zeros(accumcnts.shape)
    deltacnts[0, 1:] = accumcnts[0, 1:] - accumcnts[0, 0:-1]
    deltacnts[1:, 1:] = int16(accumcnts[1:, 1:] - accumcnts[1:, :-1])
    # accumulated angle counts per channel, adjusted for rollover, with first value zero
    accumcnts[1:, :] = deltacnts[1:, :].cumsum(axis=1)
    # count rate per channel for each time step, number of time values one less than accumcnts
    ratecnts = zeros((accumcnts.shape[0], accumcnts.shape[1] - 1))
    ratecnts[0, :] = accumcnts[0, 1:]
    ratecnts[1:, :] = deltacnts[1:, 1:] / deltacnts[0, 1:]
    # 
    return (accumcnts, deltacnts, ratecnts)
    