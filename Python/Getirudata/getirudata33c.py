#!/usr/bin/env python
# getirudata.py
# python script to find time spans of NPNT iru data with other constraints
# and compute quantities for each maneuver for subsequent calibration 
# Version 27: add selection of interval with a number
#             use average bias
# Version 28: put quaternion functions in quatfunc.py
#             do matrix computations without constant maneuver axis approximation
# Version 29: change default M-matrix to all zeros
#             move D-matrix to each interval
#             output pcad bias at start of maneuver
# Version 30: option to bypass intermediate batch computations
#             auto select NMAN based on before and after NPNT bias ave diff & stdev
#             add option to apply M-matrix
#             edited and added plots
#             get AOUNLOAD and reject NPNT & NMAN with AOUNLOAD = 'GRND' (remove mups unloads)
#             get AORWBIAS and reject NPNT & NMAN with AORWBIAS = 'DISA' (remove SCS-107 unloads)
# Version 31: style and syntax changes
# Version 32: add option for aberration adjustment
#             add option for rejecting bad times from file
# Version 33: a) 
#             b) Determine & display rotation signs for each maneuver
#                Option to output signs to file
#                Replace IRU rate calc with defs from irudefs.py
#                Fix SFave, add "/ 2", adjust parameters to match
#                Add option to bypass filter for bias stsev & diff
#                Add more config messages
#             c) Add separate options for filter for bias stsev & diff
#             

import Ska.engarchive.fetch as fetch
from Ska.Matplotlib import plot_cxctime
import Ska.Numpy
from pylab import *
import Chandra.Time
from arraydata import getstrstartstop
#import sys
import os
from math import *
from quatdefs import *
from scipy.interpolate import lagrange
import irudefs as iru

#######################################################################
# Initialization

version = 'v33c'
interval = 29 # default time interval
plot_man_flag = False # True for single maneuver mode, False for multi-maneuver mode
use_ave_bias = False # True for use of computed average bias, False for use PCAD bias (use False)
use_zero_Mmat = True # True for computing calib param, False for evaluating current ones
auto_save_plot = True
compute_batch = True # compute and output intermediate batch quantities
adj_aber = False # Adjust attitudes for aberration is True
filter_bad_times = True # Remove AOATTER bad times
filter_stdev_bias_limits = False # for True maneuvers filtered for bias limits
filter_diff_bias_limits = True # for True maneuvers filtered for bias limits
write_signs = False # Option to write signs to output file
two15 = 2**15
two16 = 2**16
rad2deg = 180.0 / pi # radians to degrees
deg2rad = pi / 180.0 # degrees to radians
rad2asec = 180.0 / pi * 3600.0 # radians to arcsec
rps2dph  = 180.0 / pi * 3600.0 # radians/sec to deg/hr
dph2rps  = 1.0 / rps2dph # deg/hr to radians/sec

# print "sys.argv=", sys.argv
if (size(sys.argv) > 1):
    interval = float(sys.argv[1])

if (interval == 1):
    tstart = '2003:274:14:00:00.000'  # start time for interval 01
    tstop  = '2004:093:00:00:00.000'  # stop time for interval 01
    interval = 'i01'
    Dmat = array([[ 0.0,  0.0,  0.0], # correction to Mmat
                  [ 0.0,  0.0,  0.0],
                  [ 0.0,  0.0,  0.0]])
elif (interval == 2):
    tstart = '2004:093:00:00:00.000'  # start time for interval 02
    tstop  = '2004:276:00:00:00.000'  # stop time for interval 02
    interval = 'i02'
    Dmat = array([[ 0.0,  0.0,  0.0], # correction to Mmat
                  [ 0.0,  0.0,  0.0],
                  [ 0.0,  0.0,  0.0]])
elif (interval == 3):
    tstart = '2004:276:00:00:00.000'  # start time for interval 03
    tstop  = '2005:092:00:00:00.000'  # stop time for interval 03
    interval = 'i03'
    Dmat = array([[ 0.0,  0.0,  0.0], # correction to Mmat
                  [ 0.0,  0.0,  0.0],
                  [ 0.0,  0.0,  0.0]])
elif (interval == 4):
    tstart = '2005:092:00:00:00.000'  # start time for interval 04
    tstop  = '2005:275:00:00:00.000'  # stop time for interval 04
    interval = 'i04'
    Dmat = array([[ 0.0,  0.0,  0.0], # correction to Mmat
                  [ 0.0,  0.0,  0.0],
                  [ 0.0,  0.0,  0.0]])
elif (interval == 5):
    tstart = '2005:275:00:00:00.000'  # start time for interval 05
    tstop  = '2006:090:00:00:00.000'  # stop time for interval 05
    interval = 'i05'
    Dmat = array([[ 0.0,  0.0,  0.0], # correction to Mmat
                  [ 0.0,  0.0,  0.0],
                  [ 0.0,  0.0,  0.0]])
elif (interval == 6):
    tstart = '2006:090:00:00:00.000'  # start time for interval 06
    tstop  = '2006:271:00:00:00.000'  # stop time for interval 06
    interval = 'i06'
    Dmat = array([[ 0.0,  0.0,  0.0], # correction to Mmat
                  [ 0.0,  0.0,  0.0],
                  [ 0.0,  0.0,  0.0]])
elif (interval == 7):
    tstart = '2006:271:00:00:00.000'  # start time for interval 07
    tstop  = '2006:352:00:00:00.000'  # stop time for interval 07
    interval = 'i07'
    Dmat = array([[ 0.0,  0.0,  0.0], # correction to Mmat
                  [ 0.0,  0.0,  0.0],
                  [ 0.0,  0.0,  0.0]])
elif (interval == 8):
    tstart = '2006:352:16:00:00.000'  # start time for interval 08
    tstop  = '2007:148:00:00:00.000'  # stop time for interval 08
    interval = 'i08'
    Dmat = array([[ 0.0,  0.0,  0.0], # correction to Mmat
                  [ 0.0,  0.0,  0.0],
                  [ 0.0,  0.0,  0.0]])
elif (interval == 9):
    tstart = '2007:148:00:00:00.000'  # start time for interval 09
    tstop  = '2007:306:00:00:00.000'  # stop time for interval 09
    interval = 'i09'
    Dmat = array([[ 0.0,  0.0,  0.0], # correction to Mmat
                  [ 0.0,  0.0,  0.0],
                  [ 0.0,  0.0,  0.0]])
elif (interval == 10):
    tstart = '2007:306:00:00:00.000'  # start time for interval 10
    tstop  = '2008:060:00:00:00.000'  # stop time for interval 10
    interval = 'i10'
    Dmat = array([[ 0.0,  0.0,  0.0], # correction to Mmat
                  [ 0.0,  0.0,  0.0],
                  [ 0.0,  0.0,  0.0]])
elif (interval == 11):
    tstart = '2008:060:00:00:00.000'  # start time for interval 11
    tstop  = '2008:186:00:00:00.000'  # stop time for interval 11
    interval = 'i11'
    Dmat = array([[ 0.0,  0.0,  0.0], # correction to Mmat
                  [ 0.0,  0.0,  0.0],
                  [ 0.0,  0.0,  0.0]])
elif (interval == 12):
    tstart = '2008:186:00:00:00.000'  # start time for interval 12
    tstop  = '2008:319:00:00:00.000'  # stop time for interval 12
    interval = 'i12'
    Dmat = array([[ 0.0,  0.0,  0.0], # correction to Mmat
                  [ 0.0,  0.0,  0.0],
                  [ 0.0,  0.0,  0.0]])
elif (interval == 13):
    tstart = '2008:319:00:00:00.000'  # start time for interval 13
    tstop  = '2009:052:00:00:00.000'  # stop time for interval 13
    interval = 'i13'
    Dmat = array([[ 0.0,  0.0,  0.0], # correction to Mmat
                  [ 0.0,  0.0,  0.0],
                  [ 0.0,  0.0,  0.0]])
elif (interval == 14):
    tstart = '2009:052:00:00:00.000'  # start time for interval 14
    tstop  = '2009:156:00:00:00.000'  # stop time for interval 14
    interval = 'i14'
    Dmat = array([[ 0.0,  0.0,  0.0], # correction to Mmat
                  [ 0.0,  0.0,  0.0],
                  [ 0.0,  0.0,  0.0]])
elif (interval == 15):
    tstart = '2009:156:00:00:00.000'  # start time for interval 15
    tstop  = '2009:275:00:00:00.000'  # stop time for interval 15
    interval = 'i15'
    Dmat = array([[ 0.0,  0.0,  0.0], # correction to Mmat
                  [ 0.0,  0.0,  0.0],
                  [ 0.0,  0.0,  0.0]])
elif (interval == 16):
    tstart = '2009:275:00:00:00.000'  # start time for interval 16
    tstop  = '2010:001:00:00:00.000'  # stop time for interval 16
    interval = 'i16'
    Dmat = array([[ 0.0,  0.0,  0.0], # correction to Mmat
                  [ 0.0,  0.0,  0.0],
                  [ 0.0,  0.0,  0.0]])
elif (interval == 17):
    tstart = '2010:001:00:00:00.000'  # start time for interval 17
    tstop  = '2010:106:00:00:00.000'  # stop time for interval 17
    interval = 'i17'
    Dmat = array([[ 0.0,  0.0,  0.0], # correction to Mmat
                  [ 0.0,  0.0,  0.0],
                  [ 0.0,  0.0,  0.0]])
elif (interval == 18):
    tstart = '2010:106:00:00:00.000'  # start time for interval 18
    tstop  = '2010:204:00:00:00.000'  # stop time for interval 18
    interval = 'i18'
    Dmat = array([[ 0.0,  0.0,  0.0], # correction to Mmat
                  [ 0.0,  0.0,  0.0],
                  [ 0.0,  0.0,  0.0]])
elif (interval == 19):
    tstart = '2010:204:00:00:00.000'  # start time for interval 19
    tstop  = '2010:302:00:00:00.000'  # stop time for interval 19
    interval = 'i19'
    Dmat = array([[ 0.0,  0.0,  0.0], # correction to Mmat
                  [ 0.0,  0.0,  0.0],
                  [ 0.0,  0.0,  0.0]])
elif (interval == 20):
    tstart = '2010:302:00:00:00.000'  # start time for interval 20
    tstop  = '2010:350:20:00:00.000'  # stop time for interval 20
    interval = 'i20'
    Dmat = array([[ 0.0,  0.0,  0.0], # correction to Mmat
                  [ 0.0,  0.0,  0.0],
                  [ 0.0,  0.0,  0.0]])
elif (interval == 21):
    tstart = '2010:351:00:00:00.000'  # start time for interval 21
    tstop  = '2011:105:21:20:00.000'  # stop time for interval 21
    interval = 'i21'
    Dmat = array([[ 0.0,  0.0,  0.0], # correction to Mmat
                  [ 0.0,  0.0,  0.0],
                  [ 0.0,  0.0,  0.0]])
elif (interval == 22):
    tstart = '2011:105:21:20:00.000'  # start time for interval 22
    tstop  = '2011:187:08:00:00.000'  # stop time for interval 22
    interval = 'i22'                  # stops before Safe Mode 4
    Dmat = array([[ 0.0,  0.0,  0.0], # correction to Mmat
                  [ 0.0,  0.0,  0.0],
                  [ 0.0,  0.0,  0.0]])
elif (interval == 23):                # starts after Safe Mode 4
    tstart = '2011:192:03:00:00.000'  # start time for interval 23
    tstop  = '2011:257:00:00:00.000'  # stop time for interval 23
    interval = 'i23'
    Dmat = array([[ 0.0,  0.0,  0.0], # correction to Mmat
                  [ 0.0,  0.0,  0.0],
                  [ 0.0,  0.0,  0.0]])
elif (interval == 24):
    tstart = '2011:257:00:00:00.000'  # start time for interval 24
    tstop  = '2011:319:00:00:00.000'  # stop time for interval 24
    interval = 'i24'
    Dmat = array([[ 0.0,  0.0,  0.0], # correction to Mmat
                  [ 0.0,  0.0,  0.0],
                  [ 0.0,  0.0,  0.0]])
elif (interval == 25):
    tstart = '2011:319:00:00:00.000'  # start time for interval 25
    tstop  = '2012:022:00:00:00.000'  # stop time for interval 25
    interval = 'i25'
    Dmat = array([[ 0.0,  0.0,  0.0], # correction to Mmat
                  [ 0.0,  0.0,  0.0],
                  [ 0.0,  0.0,  0.0]])
elif (interval == 26):
    tstart = '2012:022:00:00:00.000'  # start time for interval 26
    tstop  = '2012:062:15:00:00.000'  # stop time for interval 26
    interval = 'i26'                  # stops before M-mat uplink 5
    Dmat = array([[ 0.0,  0.0,  0.0], # correction to Mmat
                  [ 0.0,  0.0,  0.0],
                  [ 0.0,  0.0,  0.0]])
elif (interval == 27):                # starts after M-mat uplink 5
    tstart = '2012:062:16:00:00.000'  # start time for interval 27
    tstop  = '2012:150:03:33:00.000'  # stop time for interval 27
    interval = 'i27'                  # stops before Safe Mode 5
    Dmat = array([[ 0.0,  0.0,  0.0], # correction to Mmat
                  [ 0.0,  0.0,  0.0],
                  [ 0.0,  0.0,  0.0]])
elif (interval == 28):                # starts after Safe Mode 5
    tstart = '2012:152:00:00:00.000'  # start time for interval 28
    tstop  = '2012:215:00:00:00.000'  # stop time for interval 28
    interval = 'i28'
    Dmat = array([[ 0.0,  0.0,  0.0], # correction to Mmat
                  [ 0.0,  0.0,  0.0],
                  [ 0.0,  0.0,  0.0]])
elif (interval == 29):
    tstart = '2012:336:00:00:00.000'  # start time for interval 29
    tstop  = '2013:021:00:00:00.000'  # stop time for interval 29
    interval = 'i29c'
    Dmat = array([[ 0.0,  0.0,  0.0], # correction to Mmat
                  [ 0.0,  0.0,  0.0],
                  [ 0.0,  0.0,  0.0]])
elif (interval == 30):
    tstart = '2012:062:16:00:00.000'  # start time for interval 30
    tstop  = '2012:364:00:00:00.000'  # stop time for interval 30
    interval = 'i30a'
    Dmat = array([[ 0.0,  0.0,  0.0], # correction to Mmat
                  [ 0.0,  0.0,  0.0],
                  [ 0.0,  0.0,  0.0]])
elif (interval == 99):
    tstart = '2003:207:00:00:00.000'  # start time for interval 99 (all)
    tstop  = '2099:365:23:59:59.999'  # stop time for interval 99 (all)
    interval = 'i99'
    Dmat = array([[ 0.0,  0.0,  0.0], # correction to Mmat
                  [ 0.0,  0.0,  0.0],
                  [ 0.0,  0.0,  0.0]])
else:
    tstart = '2012:338:00:00:00.000'  # start time for interval 00 (custom)
    tstop  = '2012:339:00:00:00.000'  # stop time for interval 00 (custom)
    interval = 'i00'
    Dmat = array([[ 0.0,  0.0,  0.0], # correction to Mmat
                  [ 0.0,  0.0,  0.0],
                  [ 0.0,  0.0,  0.0]])

# correction to Mmat
#   Dmat = array([[  8.792448e-05,  1.409469e-04,  3.321078e-05], 
#                 [ -1.200405e-04,  9.856613e-05,  2.051482e-05],
#                 [ -3.014042e-05,  2.017329e-05,  1.311966e-04]])
npnt_min_dur = 20.0 * 60.0 # minimum duration of NPNT before and after maneuver (sec)
conv_time = 360.0 # Kalman filter converge time for small errors (sec)
acisfreqydeg = 0.3600 # ACIS Y-dither frequency
acisfreqzdeg = 0.5091 # ACIS Z-dither frequency
hrcfreqydeg = 0.3312 # HRC Y-dither frequency
hrcfreqzdeg = 0.4684 # HRC Z-dither frequency
man_ang_min = 5.0 * deg2rad # minimum maneuver angle in radians
outputfile = 'getirudata_' + interval + '_' + version + '.out'
bias_diff_lim = 0.025 * dph2rps # maximum bias difference in rad/sec
bias_stdev_lim = 0.080 * dph2rps # maximum bias standard dev in rad/sec
dump_damp = 180.0 # damping time for momentum dump (sec)
eph_pad = 600.0 # Extend ephem start and stop by ephpad (sec)
summaryfile = 'getirudata_' + interval + '_' + version + '.sum'
print 'output file = %s' % outputfile
print 'summary file = %s' % summaryfile

# On-board calibration
# pseudo-inverse G-matrix
Gmat = array([[-0.499539493,  0.500015266,  0.500455729, -0.500504173],
              [-0.254059137,  0.609733116, -0.253191860,  0.610258254],
              [-0.557983976, -0.053139506, -0.556465488, -0.053843139]])

# U-matrix, rows axes of each channel, Gmat is pseudo-inverse of Umat
Umat = array([[-0.498768681599350, -0.076240169039052, -0.863375491725958],
              [ 0.500265748681156,  0.788096859372320, -0.358660734576184],
              [ 0.500711245008005, -0.075089767304433, -0.862351305901781],
              [-0.499738137314404,  0.788337746634606, -0.358866815375450]])

# Low-rate positive scale factors (rad/cnt)
SFpos = array([ 1.0, 0.1555267e-05, 0.1564191e-05, 
                     0.1552225e-05, 0.1569751e-05]) * 0.25625 / 4.0

# Low-rate negative scale factors (rad/cnt)
SFneg = array([-1.0, 0.1555571e-05, 0.1564383e-05, 
                     0.1552371e-05, 0.1570025e-05]) * 0.25625 / 4.0

# Average low-rate scale factor
SFave = (SFpos + SFneg) / 2.0

# Combine SFpos and SFneg for function call
SFact = vstack((SFpos,SFneg))

# Misalignment/scale-factor adjustment matrix, M-matrix
if (use_zero_Mmat):
    Mmat = zeros((3, 3))
else:
    MmatTimes = np.array([Chandra.Time.DateTime('2003:203:00:00:00').secs,  # initial M-matrix after IRU swap
                          Chandra.Time.DateTime('2003:274:13:19:00').secs,  # uplinked to IRU-2, CAP 891
                          Chandra.Time.DateTime('2006:352:14:25:00').secs,  # uplinked to IRU-2, CAP 1021
                          Chandra.Time.DateTime('2010:350:22:10:00').secs,  # uplinked to IRU-2, CAP 1170, PR-283
                          Chandra.Time.DateTime('2011:105:21:20:00').secs,  # uplinked to IRU-2, CAP 1179, PR-289
                          Chandra.Time.DateTime('2012:062:15:26:00').secs]) # uplinked to IRU-2, CAP 1227, PR-309

    MmatArrays = zeros((3, 3, 6))

    MmatArrays[:, :, 1] = array([[ 3.3203451e-06,  9.3199606e-05,  1.3764573e-05],
                                 [-1.3894030e-04, -1.3754384e-05,  9.8274797e-06],
                                 [-3.8983014e-05,  4.5790156e-06,  3.3727364e-06]])

    MmatArrays[:, :, 2] = array([[ 0.392656E-4,  0.920426E-4,  0.047589E-4],
                                 [-1.509408E-4,  0.429630E-4,  0.155405E-4],
                                 [-0.477453E-4,  0.168878E-4,  0.679124E-4]])

    MmatArrays[:, :, 3] = array([[ 8.792448e-05,  1.409469e-04,  3.321078e-05],
                                 [-1.200405e-04,  9.856613e-05,  2.051482e-05],
                                 [-3.014042e-05,  2.017329e-05,  1.311966e-04]])

    MmatArrays[:, :, 4] = array([[ 1.651433e-04,  1.888956e-04,  6.763121e-05],
                                 [-5.143320e-05,  1.669320e-04,  2.689127e-05],
                                 [-2.455693e-06,  1.191769e-05,  2.150693e-04]])

    MmatArrays[:, :, 5] = array([[ 2.484911e-04,  1.933052e-04,  5.790450e-05],
                                 [ 2.882515e-05,  2.505313e-04,  4.593649e-05],
                                 [ 4.250966e-05,  8.943458e-06,  2.930867e-04]])

print 'Get IRU Calibration Data'
print 'Minimum NPNT duration = %0.3f sec' % npnt_min_dur
print 'Kalman filter converge time = %0.3f sec' %  conv_time
print 'Minimum maneuver angle = %0.3f deg' % (man_ang_min * rad2deg)
print 'Maximum bias difference limit = %0.3f deg/hr' % (bias_diff_lim * rps2dph)
print 'Maximum bias std dev limit = %0.3f deg/hr' % (bias_stdev_lim * rps2dph)
print 'Requested time interval is %s to %s' %(tstart, tstop)

if adj_aber:
    print 'Aberration adjustment is true.'
else:
    print 'Aberration adjustment is false.'

if use_zero_Mmat:
    print 'Use zero Mmat adjustment.'
else:
    print 'Use onboard Mmat adjustment.'

if filter_diff_bias_limits:
    print 'Filter with diff bias limits.'
else:
    print 'Filter without diff bias limits.'

if filter_stdev_bias_limits:
    print 'Filter with stdev bias limits.'
else:
    print 'Filter without stdev bias limits.'

# Get IRU data
# AOPCADMD, PCAD mode flag (...NMAN, NPNT...)
# AOAUTTXN, autonamous mode transistion enable disable flag (AMT) (DISA, ENAB)
# AOACASEQ, Aspect camera processing sequence indicator (BRIT, AQXN, GUID, KALM)
# AOUNLOAD, Angular momentum unloading state flag (MON ,GRND, AUTO)
# AORWBIAS, reaction wheel bias (DISA, ENAB)

# Fetch data

print 'Fetch AOPCADMD, AOAUTTXN, AOACASEQ, AOUNLOAD, and AORWBIAS'
data = fetch.MSIDset(['AOPCADMD', 'AOAUTTXN', 'AOACASEQ', 'AOUNLOAD', 'AORWBIAS'], 
                     tstart, tstop, filter_bad=True)

aopcadmd_vals = np.array(data['AOPCADMD'].vals)
aopcadmd_times = data['AOPCADMD'].times[0:] 
aoauttxn_vals = np.array(data['AOAUTTXN'].vals) 
aoauttxn_times = data['AOAUTTXN'].times[0:] 
aoacaseq_vals = np.array(data['AOACASEQ'].vals)
aoacaseq_times = data['AOACASEQ'].times[0:] 
aounload_vals = np.array(data['AOUNLOAD'].vals)
aounload_times = data['AOUNLOAD'].times[0:]
aorwbias_vals = np.array(data['AORWBIAS'].vals)
aorwbias_times = data['AORWBIAS'].times[0:]

print 'PCAD mode times %s to %s' %(Chandra.Time.DateTime(aopcadmd_times[0]).date,
                                    Chandra.Time.DateTime(aopcadmd_times[-1]).date)

## Find start and stop times and indices for NPNT, NMAN, DISA, KALM, GRND, AUTO
## Put data in arrays... np.array(2, num)   
## array(0, num) is start time, array(1, num) is stop time

# NPNT data (normal point mode)
(npnt_indices, npnt_times) = getstrstartstop(aopcadmd_vals, aopcadmd_times, 'NPNT')
print "npnt_indices.shape[1] = %d, npnt_times.shape[1] = %d" % (npnt_indices.shape[1],
                                                                npnt_times.shape[1])
if (npnt_indices.shape[1] == npnt_times.shape[1]):
    num_npnt = npnt_indices.shape[1]
else:
    sys.exit(0)

# NMAN data (normal maneuver mode)
(nman_indices, nman_times) = getstrstartstop(aopcadmd_vals, aopcadmd_times, 'NMAN')
print "nman_indices.shape[1] = %d, nman_times.shape[1] = %d" % (nman_indices.shape[1],
                                                                nman_times.shape[1])
if (nman_indices.shape[1] == nman_times.shape[1]):
    num_nman = nman_indices.shape[1]
    num_nman0 = num_nman # initial number of maneuvers
else:
    sys.exit(0)

# DISA data (segmented maneuver, no attitude update, no NPNT)
(disa_indices, disa_times) = getstrstartstop(aoauttxn_vals, aoauttxn_times, 'DISA')
if (disa_indices.size > 0):
    if (disa_indices.shape[1] == disa_times.shape[1]):
        num_disa = disa_indices.shape[1]
    else:
        sys.exit(0)
else:
    num_disa = 0

print "num_disa = %d" % (num_disa)

# KALM data (Kalman filter is running)
(kalm_indices, kalm_times) = getstrstartstop(aoacaseq_vals, aoacaseq_times, 'KALM')
print "kalm_indices.shape[1] = %d, kalm_times.shape[1] = %d" % (kalm_indices.shape[1],
                                                                kalm_times.shape[1])
if (kalm_indices.shape[1] == kalm_times.shape[1]):
    num_kalm = kalm_indices.shape[1]
else:
    sys.exit(0)

# GRND data (ground momentum unload in progress)
(grnd_indices, grnd_times) = getstrstartstop(aounload_vals, aounload_times, 'GRND')
if (grnd_indices.size > 0):
    print "grnd_indices.shape[1] = %d, grnd_times.shape[1] = %d" % (grnd_indices.shape[1],
                                                                    grnd_times.shape[1])
    if (grnd_indices.shape[1] == grnd_times.shape[1]):
        num_grnd = grnd_indices.shape[1]
    else:
        sys.exit(0)
else:
    num_grnd = 0
print "num_grnd = %d" % (num_grnd)

# add dump damping time
if num_grnd > 0:
    grnd_times[1,:] = grnd_times[1,:] + dump_damp

# rwbi data (for SCS-107)
(rwbi_indices, rwbi_times) = getstrstartstop(aorwbias_vals, aorwbias_times, 'DISA')
if (rwbi_indices.size > 0):
    print "rwbi_indices.shape[1] = %d, rwbi_times.shape[1] = %d" % (rwbi_indices.shape[1],
                                                                    rwbi_times.shape[1])
    if (rwbi_indices.shape[1] == rwbi_times.shape[1]):
        num_rwbi = rwbi_indices.shape[1]
    else:
        sys.exit(0)
else:
    num_rwbi = 0
print "num_rwbias = %d" % (num_rwbi)

# Get bad times from file from 'aoatter_bad_times.dat'
if filter_bad_times:
    textin = open('aoatter_bad_times.dat','rU').readlines()
    timelist = [x.split() for x in textin]
    startstr = array([x[0] for x in timelist])
    stopstr = array([x[1] for x in timelist])
    num_bad = startstr.shape[0]
    print 'Number of bad times = %d' % num_bad
    bad_times = zeros((2,num_bad))
    bad_times[0, :] = Chandra.Time.DateTime(startstr).secs
    bad_times[1, :] = Chandra.Time.DateTime(stopstr).secs
else:
    num_bad = 0

## Select maneuvers for calibration

# 1. maneuver not segmented: disa_times start or stop not between nman start and stop
disa_start_notin_nman = ones(num_nman, dtype=bool) # pre-allocate array
disa_stop_notin_nman = ones(num_nman, dtype=bool) # pre-allocate array
disa_notin_nman = ones(num_nman, dtype=bool) # pre-allocate array
print 'num_disa = %d, num_nman = %d' % (num_disa, num_nman)
if (num_disa > 0):
    print 'Remove NMAN with DISA start and stop'
    for n in range(num_nman):
        disa_start_notin_nman[n] = not(((nman_times[0, n] <= disa_times[0, :]) & (disa_times[0, :] <= nman_times[1, n])).any())
        disa_stop_notin_nman[n] = not(((nman_times[0, n] <= disa_times[1, :]) & (disa_times[1, :] <= nman_times[1, n])).any())

disa_notin_nman = disa_start_notin_nman & disa_stop_notin_nman
nman_indices = nman_indices[:, disa_notin_nman]
nman_times = nman_times[:, disa_notin_nman]
num_nman = nman_times.shape[1]
print 'num_nman = %d' % (num_nman)

# 2. Select NPNT intervals which have at least one KALM start, except first NPNT
num_kalm = kalm_times.shape[1]
num_npnt = npnt_times.shape[1]
kalm_in_npnt = ones(num_npnt, dtype=bool) # pre-allocate array
print 'num_kalm = %d, num_npnt = %d' % (num_kalm, num_npnt)
print 'Remove NPNT without KALM start'
for n in range(1, num_npnt):
    kalm_in_npnt[n] = (((npnt_times[0, n] <= kalm_times[0, :]) & (kalm_times[0, :] <= npnt_times[1, n])).any())

npnt_indices = npnt_indices[:, kalm_in_npnt]
npnt_times = npnt_times[:, kalm_in_npnt]
num_npnt = npnt_times.shape[1]
print 'num_npnt = %d' % (num_npnt)

# 3. Each NPNT has minimum duration of npnt_min_dur 
# For each NPNT interval selected so far, select those with a minimum duration.  
print 'Remove NPNT with durations lower than the minimum: %f sec' % npnt_min_dur
npnt_dur = npnt_times[1, :] - npnt_times[0, :]
npnt_with_min = (npnt_dur >=  npnt_min_dur)
npnt_indices = npnt_indices[:, npnt_with_min]
npnt_times = npnt_times[:, npnt_with_min]
npnt_dur = npnt_times[1, :] - npnt_times[0, :]
print 'All NPNT intervals have duration >= min duration: %s' % (npnt_dur >=  npnt_min_dur).all()
num_npnt = npnt_times.shape[1]
print 'num_npnt = %d' % (num_npnt)

# 4. maneuver has NPNT before and after
# For each NMAN interval, select those with an NPNT interval 
# before and after, out of those NPNT intervals selected above.  
print 'Remove NMAN without NPNT before and after'
nman_with_npnt = zeros(num_nman, dtype=bool) # pre-allocate array
print 'num_nman = %d' % (num_nman)
for n in range(num_nman):
    nman_with_npnt[n] = ((nman_indices[0, n] - 1) == npnt_indices[1, :]).any() & ((nman_indices[1, n] + 1) == npnt_indices[0, :]).any()

nman_indices = nman_indices[:, nman_with_npnt]
print 'All NMAN starts are after an NPNT stop: %s' % (aopcadmd_vals[nman_indices[0, :] - 1, ] == 'NPNT').all()
print 'All NMAN stops are before an NPNT start: %s' % (aopcadmd_vals[nman_indices[1, :] + 1, ] == 'NPNT').all()
nman_times = nman_times[:, nman_with_npnt]
num_nman = nman_times.shape[1]
print 'num_nman = %d' % num_nman 

# 5. find NPNT before and after NMAN
# for each NPNT interval selected so far, determine which NMAN interval it
# preceeds and/or succeeds (or not).  
print 'Find NPNT before and after NMAN'
npnt_before_nman = zeros(num_npnt, dtype=bool) # pre-allocate array
npnt_after_nman =  zeros(num_npnt, dtype=bool) # pre-allocate array
for n in range(num_npnt): # search NMAN indices in PCAD mode 
    npnt_before_nman[n] = ((nman_indices[0, :] - 1) == npnt_indices[1, n]).any()
    npnt_after_nman[n] =  ((nman_indices[1, :] + 1) == npnt_indices[0, n]).any()

npnt_before_nman_indices = npnt_indices[:, npnt_before_nman]
npnt_before_nman_times =   npnt_times[:, npnt_before_nman]
print 'All NPNT before NMAN values are before the start of an NMAN interval: %s' % (aopcadmd_vals[npnt_before_nman_indices[1, :] + 1, ] == 'NMAN').all()
npnt_after_nman_indices = npnt_indices[:, npnt_after_nman]
npnt_after_nman_times = npnt_times[:, npnt_after_nman]
print 'All NPNT after NMAN values are after the stop of an NMAN interval: %s' % (aopcadmd_vals[npnt_after_nman_indices[0, :] - 1, ] == 'NMAN').all()
print 'num of NPNT before values equals num of NPNT after values equals num of NMAN intervals: %s' % (npnt_after_nman_indices.shape[1] == npnt_after_nman_indices.shape[1] ==  nman_indices.shape[1])

# 6. adjust npnt_before_nman_times and npnt_after_nman_times to npnt_min_dur
npnt_before_nman_times[0, :] = npnt_before_nman_times[1, :] - npnt_min_dur
npnt_after_nman_times[1, :] = npnt_after_nman_times[0, :] + npnt_min_dur

# 7. remove maneuvers with MUPS (GRND) start or stop time within NPNT or NMAN intervals
grnd_start_notin_nmannpnt = ones(num_nman, dtype=bool) # pre-allocate array
grnd_stop_notin_nmannpnt = ones(num_nman, dtype=bool) # pre-allocate array
grnd_notin_nmannpnt = ones(num_nman, dtype=bool) # pre-allocate array
print 'num_grnd = %d, num_nman = %d' % (num_grnd, num_nman)
if (num_grnd > 0):
    print 'Remove NMAN with GRND start and stop in NMAN or NPNT'
    for n in range(num_nman):
        grnd_start_notin_nmannpnt[n] = not(((npnt_before_nman_times[0, n] <= grnd_times[0, :]) & (grnd_times[0, :] <= npnt_after_nman_times[1, n])).any())
        grnd_stop_notin_nmannpnt[n] = not(((npnt_before_nman_times[0, n] <= grnd_times[1, :]) & (grnd_times[1, :] <= npnt_after_nman_times[1, n])).any())

grnd_notin_nmannpnt = grnd_start_notin_nmannpnt & grnd_stop_notin_nmannpnt
nman_indices = nman_indices[:, grnd_notin_nmannpnt]
nman_times = nman_times[:, grnd_notin_nmannpnt]
num_nman = nman_times.shape[1]
print 'num_nman = %d' % (num_nman)
npnt_before_nman_indices = npnt_before_nman_indices[:, grnd_notin_nmannpnt]
npnt_before_nman_times =   npnt_before_nman_times[:, grnd_notin_nmannpnt]
print 'All NPNT before NMAN values are before the start of an NMAN interval: %s' % (aopcadmd_vals[npnt_before_nman_indices[1, :] + 1, ] == 'NMAN').all()
npnt_after_nman_indices =  npnt_after_nman_indices[:, grnd_notin_nmannpnt]
npnt_after_nman_times =    npnt_after_nman_times[:, grnd_notin_nmannpnt]
print 'All NPNT after NMAN values are after the stop of an NMAN interval: %s' % (aopcadmd_vals[npnt_after_nman_indices[0, :] - 1, ] == 'NMAN').all()
print 'num of NPNT before values equals num of NPNT after values equals num of NMAN intervals: %s' % (npnt_before_nman_indices.shape[1] == npnt_after_nman_indices.shape[1] ==  nman_indices.shape[1])

# 8. remove maneuvers with SCS-107 (AORWBIAS = 'DISA') start or stop time within NPNT or NMAN intervals
rwbi_start_notin_nmannpnt = ones(num_nman, dtype=bool) # pre-allocate array
rwbi_stop_notin_nmannpnt = ones(num_nman, dtype=bool) # pre-allocate array
rwbi_notin_nmannpnt = ones(num_nman, dtype=bool) # pre-allocate array
print 'num_rwbi = %d, num_nman = %d' % (num_rwbi, num_nman)
if (num_rwbi > 0):
    print 'Remove NMAN with AORWBIAS/DISA start and stop in NMAN or NPNT'
    for n in range(num_nman):
        rwbi_start_notin_nmannpnt[n] = not(((npnt_before_nman_times[0, n] <= rwbi_times[0, :]) & (rwbi_times[0, :] <= npnt_after_nman_times[1, n])).any())
        rwbi_stop_notin_nmannpnt[n] = not(((npnt_before_nman_times[0, n] <= rwbi_times[1, :]) & (rwbi_times[1, :] <= npnt_after_nman_times[1, n])).any())

rwbi_notin_nmannpnt = rwbi_start_notin_nmannpnt & rwbi_stop_notin_nmannpnt
nman_indices = nman_indices[:, rwbi_notin_nmannpnt]
nman_times = nman_times[:, rwbi_notin_nmannpnt]
num_nman = nman_times.shape[1]
print 'num_nman = %d' % (num_nman)
npnt_before_nman_indices = npnt_before_nman_indices[:, rwbi_notin_nmannpnt]
npnt_before_nman_times =   npnt_before_nman_times[:, rwbi_notin_nmannpnt]
print 'All NPNT before NMAN values are before the start of an NMAN interval: %s' % (aopcadmd_vals[npnt_before_nman_indices[1, :] + 1, ] == 'NMAN').all()
npnt_after_nman_indices =  npnt_after_nman_indices[:, rwbi_notin_nmannpnt]
npnt_after_nman_times =    npnt_after_nman_times[:, rwbi_notin_nmannpnt]
print 'All NPNT after NMAN values are after the stop of an NMAN interval: %s' % (aopcadmd_vals[npnt_after_nman_indices[0, :] - 1, ] == 'NMAN').all()
print 'num of NPNT before values equals num of NPNT after values equals num of NMAN intervals: %s' % (npnt_before_nman_indices.shape[1] == npnt_after_nman_indices.shape[1] ==  nman_indices.shape[1])

# 9. remove maneuvers with bad times from file, start or stop time within NPNT or NMAN intervals or encloses entire interval
if (filter_bad_times and (num_bad > 0)):
    bad_start_notin_nmannpnt = ones(num_nman, dtype=bool) # pre-allocate array
    bad_stop_notin_nmannpnt = ones(num_nman, dtype=bool) # pre-allocate array
    nmannpnt_notin_bad = ones(num_nman, dtype=bool) # pre-allocate array
    nmannpnt_not_bad = ones(num_nman, dtype=bool) # pre-allocate array
    print 'num_bad = %d, num_nman = %d' % (num_bad, num_nman)
    print 'Remove NMAN with bad start and stop in NMAN or NPNT or all'
    for n in range(num_nman):
        bad_start_notin_nmannpnt[n] = not(((npnt_before_nman_times[0, n] <= bad_times[0, :]) & (bad_times[0, :] <= npnt_after_nman_times[1, n])).any())
        bad_stop_notin_nmannpnt[n] = not(((npnt_before_nman_times[0, n] <= bad_times[1, :]) & (bad_times[1, :] <= npnt_after_nman_times[1, n])).any())
        nmannpnt_notin_bad[n] = not(((npnt_before_nman_times[0, n] >= bad_times[0, :]) & (bad_times[1, :] >= npnt_after_nman_times[1, n])).any())
    nmannpnt_not_bad = bad_start_notin_nmannpnt & bad_stop_notin_nmannpnt & nmannpnt_notin_bad
    nman_indices = nman_indices[:, nmannpnt_not_bad]
    nman_times = nman_times[:, nmannpnt_not_bad]
    num_nman = nman_times.shape[1]
    print 'num_nman = %d' % (num_nman)
    npnt_before_nman_indices = npnt_before_nman_indices[:, nmannpnt_not_bad]
    npnt_before_nman_times =   npnt_before_nman_times[:, nmannpnt_not_bad]
    print 'All NPNT before NMAN values are before the start of an NMAN interval: %s' % (aopcadmd_vals[npnt_before_nman_indices[1, :] + 1, ] == 'NMAN').all()
    npnt_after_nman_indices =  npnt_after_nman_indices[:, nmannpnt_not_bad]
    npnt_after_nman_times =    npnt_after_nman_times[:, nmannpnt_not_bad]
    print 'All NPNT after NMAN values are after the stop of an NMAN interval: %s' % (aopcadmd_vals[npnt_after_nman_indices[0, :] - 1, ] == 'NMAN').all()
    print 'num of NPNT before values equals num of NPNT after values equals num of NMAN intervals: %s' % (npnt_before_nman_indices.shape[1] == npnt_after_nman_indices.shape[1] ==  nman_indices.shape[1])

# sys.exit()

# 10. find last kalman start time in npnt after each maneuver
#    find KALM start in npnt_after_nman_times
print 'Find Kalman start times in NPNT after each maneuver'
kalm_start_in_npnt = zeros(num_kalm, dtype = bool) # pre-allocate array
print "num_kalm = %d, num_nman = %d" % (num_kalm, num_nman)
for n in range(num_kalm):
    kalm_start_in_npnt[n] = ((npnt_after_nman_times[0, :] <= kalm_times[0, n]) & (kalm_times[0, n] <= npnt_after_nman_times[1, :])).any()

kalm_indices = kalm_indices[:, kalm_start_in_npnt]
kalm_times = kalm_times[:, kalm_start_in_npnt]
num_kalm = kalm_times.shape[1]
print "num_kalm = %d, num_nman = %d" % (num_kalm, num_nman)

#11. remove all but last KALM start times from npnt_after_nman_times
print 'Remove all but last Kalman time in NPNT intervals after NMAN'
kalm_start_in_npnt = ones(num_kalm, dtype=bool) # pre-allocate array
for n in range(num_nman):
    kalm_index = find((npnt_after_nman_times[0, n] <= kalm_times[0, :]) & (kalm_times[0, :] <= npnt_after_nman_times[1, n]))
    if (kalm_index.shape[0] > 1):
        kalm_start_in_npnt[kalm_index[:-1]] = False

kalm_indices = kalm_indices[:, kalm_start_in_npnt]
kalm_times = kalm_times[:, kalm_start_in_npnt]
num_kalm = kalm_times.shape[1]
print "num_kalm = %d, num_nman = %d" % (num_kalm, num_nman)

## Computations for each maneuver

# Preallocate arrays for each maneuver
initquat  = zeros((5, num_nman)) # initial pcad quaternion for each maneuver, time, q1, q2, q3, q4
finalquat = zeros((5, num_nman)) # final pcad quaternion for each maneuver, time, q1, q2, q3, q4
manvrquat = zeros((5, num_nman)) # observed rotation for each maneuver, time-diff, q1, q2, q3, q4
manvrtime = zeros((2, num_nman)) # maneuver time interval start (0) and stop (1)
intratebody = zeros((4, num_nman)) # integrated body rates over each maneuver
diffchancnts = zeros((5, num_nman)) # difference in counts across maneuver
ave_bias_before_nman = zeros((5, num_nman)) # average of rate in time interval before maneuver
std_bias_before_nman = zeros((5, num_nman)) # standard deviation of rate in time interval before maneuver
ave_bias_after_nman = zeros((5, num_nman)) # average of rate in time interval after maneuver
std_bias_after_nman = zeros((5, num_nman)) # standard of rate in time interval after maneuver
ave_cnt_bias = zeros((5, num_nman)) # average value of bias for maneuver (cnts/sec)
dif_cnt_bias = zeros((5, num_nman)) # average value of bias for maneuver (cnts/sec)
ini2finquat = zeros((5, num_nman)) # rotation quaternion from pcad initial & final attitudes
ini2finvect = zeros((4, num_nman)) # 0 index is stop time; 1, 2, 3 indices are rotation vector 
ini2finang = zeros((2, num_nman)) # 0 index is stop time; 1 index is rotation vector 
finalpropquat = zeros((5, num_nman)) # time, q1, q2, q3, q4; propagated initial quaternion
deltaquat = zeros((5, num_nman)) # time, q1, q2, q3, q4; difference between propagation and solution, 
deltavect = zeros((4, num_nman)) # time, v1, v2, v3; difference between propagation and solution, 
deltaYZ = zeros((2, num_nman)) # time & YZ magnitude of delta vector
pcadbias_start = zeros((4, num_nman)) # PCAD bias at start of maneuver
if compute_batch:
    sumprop = zeros((3, 3, num_nman)) # sum of prop-mat for maneuver
    sumproprot = zeros((3, 9, num_nman)) # sum of prop-mat-func(rate) for maneuver

if plot_man_flag:
    mannum = input('input number of maneuver: ')
    rng = [mannum]
else:
    rng = range(num_nman)

print "Begin loop over maneuvers for n = 0 to %d" % (num_nman - 1)
# Computations for each maneuver or for single specified maneuver (plot_man_flag == True)
for n in rng:
#   fetch all quaternions in pre, during, and post maneuver interval
    os.write(0, '+') # indicates start of loop on console
    data = fetch.MSIDset(['AOATTQT1', 'AOATTQT2', 'AOATTQT3', 'AOATTQT4'], 
                         npnt_before_nman_times[0, n], 
                         npnt_after_nman_times[1, n], filter_bad = True)
    pcadquat = np.array([data['AOATTQT1'].times[0:], # time of quaternion
                         data['AOATTQT1'].vals,      # q1
                         data['AOATTQT2'].vals,      # q2
                         data['AOATTQT3'].vals,      # q3
                         data['AOATTQT4'].vals])     # q4
    num_quat = pcadquat.shape[1]
    
#   fetch all CXO and Earth velocities in pre, during, and post maneuver interval
    if adj_aber:
        data = fetch.MSIDset(['orbitephem1_vx', 'orbitephem1_vy','orbitephem1_vz'], 
                             npnt_before_nman_times[0, n] - eph_pad, 
                             npnt_after_nman_times[1, n] + eph_pad, 
                             filter_bad = True)
        cxovel = np.array([np.array(data['orbitephem1_vx'].times),
                           np.array(data['orbitephem1_vx'].vals)/1000.0,  # km/sec
                           np.array(data['orbitephem1_vy'].vals)/1000.0,  # km/sec
                           np.array(data['orbitephem1_vz'].vals)/1000.0]) # km/sec
        
#       Get Sun position wrt ECI
        data = fetch.MSIDset(['solarephem1_vx', 'solarephem1_vy', 'solarephem1_vz'], 
                             npnt_before_nman_times[0, n] - eph_pad, 
                             npnt_after_nman_times[1, n] + eph_pad, 
                             filter_bad=True)
        sunvel = np.array([np.array(data['solarephem1_vx'].times),
                           np.array(data['solarephem1_vx'].vals)/1000.0,  # km/sec
                           np.array(data['solarephem1_vy'].vals)/1000.0,  # km/sec
                           np.array(data['solarephem1_vz'].vals)/1000.0]) # km/sec
        
#       CXO velocity with respect to Sun
        cxovel[1:, ] = cxovel[1:, ] - sunvel[1:, ]
        num_vel = cxovel.shape[1]
        
#       Velocity at the time of each PCAD quaternion
        quatvel = zeros((4,num_quat))
        for m in range(2,num_vel - 1):
            idx = find((cxovel[0,m-1] < pcadquat[0,]) & (pcadquat[0,] <= cxovel[0,m]))
            quatvel[0, idx] = pcadquat[0, idx]
            P = lagrange(cxovel[0, m-2:m+2] - cxovel[0, 0], cxovel[1, m-2:m+2])
            quatvel[1, idx] = P(quatvel[0, idx] - cxovel[0, 0])
            P = lagrange(cxovel[0, m-2:m+2] - cxovel[0, 0], cxovel[2, m-2:m+2])
            quatvel[2, idx] = P(quatvel[0, idx] - cxovel[0, 0])
            P = lagrange(cxovel[0, m-2:m+2] - cxovel[0, 0], cxovel[3, m-2:m+2])
            quatvel[3, idx] = P(quatvel[0, idx] - cxovel[0, 0])
#       Adjust pcad quaternion for velocity averration
        pcadquat = quatxaber(pcadquat,quatvel)
    
#   find last NPNT attitude quaternion before NMAN
    idx = find(pcadquat[0, :] < nman_times[0, n])
    initquat[:, n] = pcadquat[:, idx.max()]
    manvrtime[0, n] = initquat[0, n]

#   find final NPNT attitude quaternion after Kalman converges
    idx = find(pcadquat[0, :] > (kalm_times[0, n] + conv_time))
    finalquat[:, n] = pcadquat[:, idx.min()]
    manvrtime[1, n] = finalquat[0, n]
    
#   compute rotation quaternion and maneuver eigen axis & angle
    ini2finquat[:, n] = quatnorm(quatmult(quatconj(initquat[:, n]), finalquat[:, n])).T
    ini2finquat[0, n] = manvrtime[1, n] - manvrtime[0, n]
    ini2finvect[:, n] = quat2vect(ini2finquat[:, n]).T
    ini2finang[0, n] = ini2finvect[0, n]
    ini2finang[1, n] = vectmag(ini2finvect[1:, n])

#   obtain IRU channel accum cnts data 30 min before maneuver entire time interval (NPNT, NMAN, NPNT)
    data = fetch.MSIDset(['AOGYRCT1', 'AOGYRCT2', 'AOGYRCT3', 'AOGYRCT4'], npnt_before_nman_times[0, n], npnt_after_nman_times[1, n], filter_bad=True)
    accumcnts = np.array([data['AOGYRCT1'].times[0:], # time of accumulated counts
                          data['AOGYRCT1'].vals,      # channel-1 accum-cnts with roll-over
                          data['AOGYRCT2'].vals,      # channel-2 accum-cnts with roll-over
                          data['AOGYRCT3'].vals,      # channel-3 accum-cnts with roll-over
                          data['AOGYRCT4'].vals])     # channel-4 accum-cnts with roll-over
    num_cnts = accumcnts.shape[1]

#   plot raw counts
    if plot_man_flag:
        figure(1)
        clf()
        title('Raw IRU Counts for each Channel')
        xlabel('Time')
        ylabel('Counts')
        plot_cxctime(accumcnts[0, :], accumcnts[1, :], '-r')
        plot_cxctime(accumcnts[0, :], accumcnts[2, :], '-g')
        plot_cxctime(accumcnts[0, :], accumcnts[3, :], '-b')
        plot_cxctime(accumcnts[0, :], accumcnts[4, :], '-m')
        grid('on')
        legend(('Channel-1', 'Channel-2', 'Channel-3', 'Channel-4'), loc = 'best')
        draw()
        show(block = False)
        if auto_save_plot:
            figfilename = 'Fig01_RawCounts_%s_m%03d_%s.png' % (interval, n, version)
            savefig(figfilename)

#   compute delta-time and delta-counts & adjust for roll-over
    (accumcnts, deltacnts, ratecnts) = iru.irucounts(accumcnts)

#   plot counts adjusted for roll-over
    if plot_man_flag:
        figure(2)
        clf()
        title('Counts Adjusted for Rollover')
        xlabel('Time')
        ylabel('Counts')
        plot_cxctime(accumcnts[0, :], accumcnts[1, :], '-r')
        plot_cxctime(accumcnts[0, :], accumcnts[2, :], '-g')
        plot_cxctime(accumcnts[0, :], accumcnts[3, :], '-b')
        plot_cxctime(accumcnts[0, :], accumcnts[4, :], '-m')
        grid('on')
        legend(('Channel-1', 'Channel-2', 'Channel-3', 'Channel-4'), loc = 'best')
        draw()
        show(block = False)
        if auto_save_plot:
            figfilename = 'Fig02_AdjCounts_%s_m%03d_%s.png' % (interval, n, version)
            savefig(figfilename)

#   Fetch PCAD bias for entire interval
    data = fetch.MSIDset(['AOGBIAS1', 'AOGBIAS2', 'AOGBIAS3'], 
                          npnt_before_nman_times[0, n], 
                          npnt_after_nman_times[1, n], filter_bad=True)
    pcadbias = np.array([data['AOGBIAS1'].times[0:], # time of accumulated counts
                         data['AOGBIAS1'].vals,      # X-axis bias
                         data['AOGBIAS2'].vals,      # Y-axis bias
                         data['AOGBIAS3'].vals])     # Z-axis bias
    num_bias = accumcnts.shape[1]

#   compute the difference in channel counts (and time) across maneuver
    start_index = max(find(accumcnts[0, :] <= (initquat[0, n] + 0.01)))
    stop_index = min(find(accumcnts[0, :] >= (finalquat[0, n] - 0.01)))
    diffchancnts[:, n] = accumcnts[:, stop_index] - accumcnts[:, start_index]

    if plot_man_flag:
        figure(3)
        clf()
        title('Count Rate for Each Channel')
        xlabel('Time')
        ylabel('Count Rate (counts/sec)')
        plot_cxctime(ratecnts[0, :], ratecnts[1, :], '-r')
        plot_cxctime(ratecnts[0, :], ratecnts[2, :], '-g')
        plot_cxctime(ratecnts[0, :], ratecnts[3, :], '-b')
        plot_cxctime(ratecnts[0, :], ratecnts[4, :], '-m')
        grid('on')
        legend(('Channel-1', 'Channel-2', 'Channel-3', 'Channel-4'), loc = 'best')
        draw()
        show(block = False)
        if auto_save_plot:
            figfilename = 'Fig03_CountRate_%s_m%03d_%s.png' % (interval, n, version)
            savefig(figfilename)

#   compute average rate per channel & std (bias in cnts/sec) before maneuver
    start_index = 0
    stop_index = max(find(accumcnts[0, :] < nman_times[0, n]))
    ave_bias_before_nman[0, n] = accumcnts[0, stop_index] # stop time of bias
    ave_bias_before_nman[1:, n] = ratecnts[1:, start_index:stop_index].mean(axis = 1)
    std_bias_before_nman[0, n] = accumcnts[0, stop_index] # stop time of bias
    std_bias_before_nman[1:, n] = ratecnts[1:, start_index:stop_index].std(axis = 1)
    
#   compute average rate per channel & std (bias in cnts/sec) after maneuver
    start_index = min(find(accumcnts[0, :] > finalquat[0, n]))
    stop_index = num_cnts - 1
    ave_bias_after_nman[0, n] = accumcnts[0, start_index] # start time of bias
    ave_bias_after_nman[1:, n] = ratecnts[1:, start_index:stop_index].mean(axis = 1)
    std_bias_after_nman[0, n] = accumcnts[0, start_index] # start time of bias
    std_bias_after_nman[1:, n] = ratecnts[1:, start_index:stop_index].std(axis = 1)
    
#   Compute average bias for maneuver (cnts/sec)
    ave_cnt_bias[0, n] = (ave_bias_before_nman[0, n] + ave_bias_after_nman[0, n]) / 2.0
    ave_cnt_bias[1:, n] = (ave_bias_before_nman[1:, n] + ave_bias_after_nman[1:, n]) / 2.0
#   Compute the difference in bias between before and after
    dif_cnt_bias[0, n] = (ave_bias_before_nman[0, n] + ave_bias_after_nman[0, n]) / 2.0
    dif_cnt_bias[1:, n] = ave_bias_before_nman[1:, n] - ave_bias_after_nman[1:, n]
    
#   get first NMAN bias
    idx = find(pcadbias[0, :] >= nman_times[0, n])[0]
    pcadbias_start[:, n] = pcadbias[:, idx]

#   compute 3-vector angular rate (rad/sec)
#   set bias
    if use_ave_bias:
        Bias4 = ave_cnt_bias[1:, n]
        Bias3 = zeros((3,1))
    else:
        Bias4 = zeros((4,1))
        Bias3 = pcadbias_start[1:, n]
    
#   select M-matrix for time span n
    if (use_zero_Mmat):
        Mmat = zeros((3, 3))
    else:
        if (ratecnts[0, -1] < MmatTimes[0]):
            Mmat = zeros((3, 3))
        elif (ratecnts[0, -1] < MmatTimes[1]):
            Mmat = MmatArrays[:, :, 0]
        elif (ratecnts[0, -1] < MmatTimes[2]):
            Mmat = MmatArrays[:, :, 1]
        elif (ratecnts[0, -1] < MmatTimes[3]):
            Mmat = MmatArrays[:, :, 2]
        elif (ratecnts[0, -1] < MmatTimes[4]):
            Mmat = MmatArrays[:, :, 3]
        elif (ratecnts[0, -1] < MmatTimes[5]):
            Mmat = MmatArrays[:, :, 4]
        else:
            Mmat = MmatArrays[:, :, 5]
    
#   compute angular rate
    (angratechan, angratebody) = iru.irurates(Dmat, Mmat, Gmat, SFact, Bias4, Bias3, ratecnts)
    
#   plot angular IRU channel rates (maybe adjusted for bias) (deg/hr)
    if plot_man_flag:
        figure(4) 
        clf()
        title('Angular Rate for Each Channel')
        xlabel('Time')
        ylabel('Ang Rate (deg/hr)')
        plot_cxctime(angratechan[0, :], angratechan[1, :] * rps2dph, '-r')
        plot_cxctime(angratechan[0, :], angratechan[2, :] * rps2dph, '-g')
        plot_cxctime(angratechan[0, :], angratechan[3, :] * rps2dph, '-b')
        plot_cxctime(angratechan[0, :], angratechan[4, :] * rps2dph, '-m')
        grid('on')
        legend(('Channel-1', 'Channel-2', 'Channel-3', 'Channel-4'), loc = 'best')
        draw()
        show(block = False)
        if auto_save_plot:
            figfilename = 'Fig04_ChanAngRate_%s_m%03d_%s.png' % (interval, n, version)
            savefig(figfilename)

#   adjust time of rate
    angratebody[1, :-1] = angratebody[1, :-1] * 0.75 + angratebody[1, 1:] * 0.25
    angratebody[2, :-1] = angratebody[2, :-1] * 0.75 + angratebody[2, 1:] * 0.25
    angratebody[3, :-1] = angratebody[3, :-1] * 0.75 + angratebody[3, 1:] * 0.25

#   plot angular S/C 3-vector rates adjusted for bias (deg/hr)
    if plot_man_flag:
        figure(5)
        clf()
        title('Angular Rates Adjusted for Bias')
        xlabel('Time')
        ylabel('Ang Rate (deg/hr)')
        plot_cxctime(angratebody[0, :], angratebody[1, :] * rps2dph, '-r')
        plot_cxctime(angratebody[0, :], angratebody[2, :] * rps2dph, '-g')
        plot_cxctime(angratebody[0, :], angratebody[3, :] * rps2dph, '-b')
        legend(('X-axis (Roll))', 'Y-axis (Pitch)', 'Z-axis (Yaw)'), loc = 'best')
        grid('on')
        draw()
        show(block = False)
        if auto_save_plot:
            figfilename = 'Fig05_SCAngRate_%s_m%03d_%s.png' % (interval, n, version)
            savefig(figfilename)

#   propagated maneuver rates 
    idx = find(initquat[0, n] < angratebody[0, :])
    idx_begin = idx.min() # index of first angratebody at time after initial quaternion
    idx = find(finalquat[0, n] >= (angratebody[0, :] - 0.01))
    idx_end = idx.max() # index of last angratebody upto time of final quaternion
    manvrquat[-1, n] = 1 # sets initial maneuver rotation quaternion to [0.; 0.; 0.; 1.]
    if plot_man_flag:
        propquat =  initquat[:, n]
        propdeltquat = zeros((5, (idx_end - idx_begin + 1)))
        ix = 0
    for idx in range(idx_begin, (idx_end + 1)):
        deltatime = angratebody[0, idx] - angratebody[0, idx - 1]
        intratebody[0, n] = intratebody[0, n] + deltatime # accumulated time
        intratebody[1:, n] = intratebody[1:, n] + angratebody[1:, idx] * deltatime # integrated body rates
        manvrquat[0, n] = manvrquat[0, n] + deltatime # accumulated time
        rotvec = angratebody[:, idx]
        rotvec[1:, ] = rotvec[1:, ] * deltatime # rotation vector for interval
        rotquat = vect2quat(rotvec).T # rotation quaternion for interval
        manvrquat[:, n] = quatmult(manvrquat[:, n], rotquat).T # propagate maneuver quat
        manvrquat[:, n] = quatnorm(manvrquat[:, n]).T
        matrot3x9 = array([[rotvec[1, ], rotvec[2, ], rotvec[3, ], 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                           [0.0, 0.0, 0.0, rotvec[1, ], rotvec[2, ], rotvec[3, ], 0.0, 0.0, 0.0],
                           [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, rotvec[1, ], rotvec[2, ], rotvec[3, ]]])
        rotmat = quat2mat(quatconj(manvrquat[:, n])) # rotmat from current to initial
        if compute_batch:
            sumprop[:, :, n] = sumprop[:, :, n] + rotmat * deltatime
            sumproprot[:, :, n] = sumproprot[:, :, n] + dot(rotmat, matrot3x9)
        
        if plot_man_flag:
            propquat = quatmult(propquat, vect2quat(rotvec)).T  # propagate attitude quaternion
            timidx = find(abs(pcadquat[0, :] - angratebody[0, idx]) < 0.1)
            if (timidx.size == 1):
                propdeltquat[:, ix] = quatmult(quatconj(propquat), pcadquat[:, timidx[0]]).T
                propdeltquat[:, ix] = quatnorm(propdeltquat[:, ix]).T
                ix = ix + 1
    if plot_man_flag:
        ix = (propdeltquat[0, :] > 0)
        propdeltquat = propdeltquat[:, ix]

#   Compute final quaternion by propagation with rates
    finalpropquat[:, n] = quatmult(initquat[:, n], manvrquat[:, n]).T
    finalpropquat[:, n] = quatnorm(finalpropquat[:, n]).T

#   Compute final sums of propagation...
    rotmat = quat2mat(manvrquat[:, n]) # rotation matrix from initial to final
    if compute_batch:
        sumprop[:, :, n] = dot(rotmat, sumprop[:, :, n])
        sumproprot[:, :, n] = dot(rotmat, sumproprot[:, :, n])

#   Compute delta quaternion between propagation and solution
    deltaquat[:, n] = quatmult(quatconj(finalpropquat[:, n]), finalquat[:, n]).T
    deltaquat[:, n] = quatnorm(deltaquat[:, n]).T
    deltavect[:, n] = quat2vect(deltaquat[:, n]).T
    deltaYZ[0, n] = manvrtime[1, n]
    deltaYZ[1, n] = np.sqrt((deltavect[2:, n] * deltavect[2:, n]).sum(axis = 0))

#   Plot delta quat through maneuver
    if plot_man_flag:
        figure(6)
        clf()
        title('Delta Vector')
        xlabel('Time')
        ylabel('Delta Vec (arcsec)')
        plot_cxctime(propdeltquat[0, :], propdeltquat[1, :] * rad2asec * 2.0, '-r')
        plot_cxctime(propdeltquat[0, :], propdeltquat[2, :] * rad2asec * 2.0, '-g')
        plot_cxctime(propdeltquat[0, :], propdeltquat[3, :] * rad2asec * 2.0, '-b')
        grid('on')
        legend(('X-axis (Roll))', 'Y-axis (Pitch)', 'Z-axis (Yaw)'), loc='best')
        draw()
        show(block = False)
        if auto_save_plot:
            figfilename = 'Fig06_DeltaVector_%s_m%03d_%s.png' % (interval, n, version)
            savefig(figfilename)

#   end of loop, print '-'
    os.write(0, '-') # indicates end of each loop on console

print '.' # indicates end of maneuver for-loop

if plot_man_flag: # exit if plot single maneuver
    sys.exit()
    
num_man = nman_times.shape[1]
print "Number of maneuvers before angle & bias constraints applied = %d" % (num_man)

# select maneuvers above specified angle, bias diff below limit, bias stdev below lim
constraints = (ini2finang[1, :] >= man_ang_min)
print "Number of maneuvers after maneuver angle constraint = %d" % (constraints.sum())
if filter_stdev_bias_limits:
    constraints = constraints & ((std_bias_before_nman[1, :] * SFave[1]) < bias_stdev_lim) & ((std_bias_after_nman[1, :] * SFave[1]) < bias_stdev_lim)
    constraints = constraints & ((std_bias_before_nman[2, :] * SFave[2]) < bias_stdev_lim) & ((std_bias_after_nman[2, :] * SFave[2]) < bias_stdev_lim)
    constraints = constraints & ((std_bias_before_nman[3, :] * SFave[3]) < bias_stdev_lim) & ((std_bias_after_nman[3, :] * SFave[3]) < bias_stdev_lim)
    constraints = constraints & ((std_bias_before_nman[4, :] * SFave[4]) < bias_stdev_lim) & ((std_bias_after_nman[4, :] * SFave[4]) < bias_stdev_lim)
    print "Number of maneuvers after bias stdev constraint = %d" % (constraints.sum())
if filter_diff_bias_limits:
    constraints = constraints & (abs(dif_cnt_bias[1, :] * SFave[1]) < bias_diff_lim)
    constraints = constraints & (abs(dif_cnt_bias[2, :] * SFave[2]) < bias_diff_lim)
    constraints = constraints & (abs(dif_cnt_bias[3, :] * SFave[3]) < bias_diff_lim)
    constraints = constraints & (abs(dif_cnt_bias[4, :] * SFave[4]) < bias_diff_lim)
    print "Number of maneuvers after bias difference constraint = %d" % (constraints.sum())
idx = find(constraints)
num_man = idx.shape[0]
print('Number of Maneuvers with angle >= %7.3f deg is %d') % (man_ang_min * rad2deg, num_man)

# select maneuvers with idx
initquat  = initquat[:, idx]
finalquat = finalquat[:, idx]
manvrquat = manvrquat[:, idx]
manvrtime = manvrtime[:, idx]
intratebody = intratebody[:, idx]
diffchancnts = diffchancnts[:, idx]
ave_bias_before_nman = ave_bias_before_nman[:, idx]
std_bias_before_nman = std_bias_before_nman[:, idx]
ave_bias_after_nman = ave_bias_after_nman[:, idx]
std_bias_after_nman = std_bias_after_nman[:, idx]
ave_cnt_bias = ave_cnt_bias[:, idx]
dif_cnt_bias = dif_cnt_bias[:, idx]
ini2finquat = ini2finquat[:, idx]
ini2finvect = ini2finvect[:, idx]
ini2finang = ini2finang[:, idx]
finalpropquat = finalpropquat[:, idx]
deltaquat = deltaquat[:, idx]
deltavect = deltavect[:, idx]
deltaYZ = deltaYZ[:, idx]
pcadbias_start = pcadbias_start[:, idx]
sumprop = sumprop[:, :, idx]
sumproprot = sumproprot[:, :, idx]

# compute sign codes for each maneuver
signcode = iru.irusigns(Umat, ini2finvect[1:, :])
numcodes = np.array([(signcode == 0).sum(), (signcode == 1).sum(), (signcode == 2).sum(),
                     (signcode == 3).sum(), (signcode == 4).sum(), (signcode == 5).sum(),
                     (signcode == 6).sum(), (signcode == 7).sum(), (signcode == 8).sum(),
                     (signcode == 9).sum(), (signcode == 10).sum(), (signcode == 11).sum(),
                     (signcode == 12).sum(), (signcode == 13).sum(), (signcode == 14).sum(),
                     (signcode == 15).sum()])
print(' code signs num  code signs num')
print('   0  ----  %2d    15  ++++  %2d') % (numcodes[0], numcodes[15])
print('   1  ---+  %2d    14  +++-  %2d') % (numcodes[1], numcodes[14])
print('   2  --+-  %2d    13  ++-+  %2d') % (numcodes[2], numcodes[13])
print('   3  --++  %2d    12  ++--  %2d') % (numcodes[3], numcodes[12])
print('   4  -+--  %2d    11  +-++  %2d') % (numcodes[4], numcodes[11])
print('   5  -+-+  %2d    10  +-+-  %2d') % (numcodes[5], numcodes[10])
print('   6  -++-  %2d     9  +--+  %2d') % (numcodes[6], numcodes[9])
print('   7  -+++  %2d     8  +---  %2d') % (numcodes[7], numcodes[8])

figure(7)
clf()
title('STDev of Bias before NMAN')
xlabel('Time')
ylabel('Bias StDev (deg/hr)')
plot_cxctime(std_bias_before_nman[0, :], std_bias_before_nman[1, :] * SFave[1] * rps2dph, '-r')
plot_cxctime(std_bias_before_nman[0, :], std_bias_before_nman[2, :] * SFave[2] * rps2dph, '-g')
plot_cxctime(std_bias_before_nman[0, :], std_bias_before_nman[3, :] * SFave[3] * rps2dph, '-b')
plot_cxctime(std_bias_before_nman[0, :], std_bias_before_nman[4, :] * SFave[4] * rps2dph, '-m')
axis_arr = array(axis())
axis_arr[2] = 0.0
axis(axis_arr)
grid('on')
legend(('Channel-1', 'Channel-2', 'Channel-3', 'Channel-4'), loc = 'best')
draw()
show(block = False)
if auto_save_plot:
    figfilename = 'Fig07_STDevBiasBefore_%s_%s.png' % (interval, version)
    savefig(figfilename)

figure(8)
clf()
title('STDev of Bias after NMAN')
xlabel('Time')
ylabel('Bias StDev (deg/hr)')
plot_cxctime(std_bias_after_nman[0, :], std_bias_after_nman[1, :] * SFave[1] * rps2dph, '-r')
plot_cxctime(std_bias_after_nman[0, :], std_bias_after_nman[2, :] * SFave[2] * rps2dph, '-g')
plot_cxctime(std_bias_after_nman[0, :], std_bias_after_nman[3, :] * SFave[3] * rps2dph, '-b')
plot_cxctime(std_bias_after_nman[0, :], std_bias_after_nman[4, :] * SFave[4] * rps2dph, '-m')
axis_arr = array(axis())
axis_arr[2] = 0.0
axis(axis_arr)
grid('on')
legend(('Channel-1', 'Channel-2', 'Channel-3', 'Channel-4'), loc = 'best')
draw()
show(block = False)
if auto_save_plot:
    figfilename = 'Fig08_STDevBiasAfter_%s_%s.png' % (interval, version)
    savefig(figfilename)

figure(9)
clf()
title('Average Bias before & after in NPNT')
xlabel('Time')
ylabel('Ave Bias (deg/hr)')
plot_cxctime(ave_cnt_bias[0, :], ave_cnt_bias[1, :] * SFave[1] * rps2dph, '-r')
plot_cxctime(ave_cnt_bias[0, :], ave_cnt_bias[2, :] * SFave[2] * rps2dph, '-g')
plot_cxctime(ave_cnt_bias[0, :], ave_cnt_bias[3, :] * SFave[3] * rps2dph, '-b')
plot_cxctime(ave_cnt_bias[0, :], ave_cnt_bias[4, :] * SFave[4] * rps2dph, '-m')
grid('on')
legend(('Channel-1', 'Channel-2', 'Channel-3', 'Channel-4'), loc = 'best')
draw()
show(block = False)
if auto_save_plot:
    figfilename = 'Fig09_AveChanBiasNPNT_%s_%s.png' % (interval, version)
    savefig(figfilename)

figure(10)
clf()
title('Difference of Bias before & after in NPNT')
xlabel('Time')
ylabel('Bias Diff (deg/hr)')
plot_cxctime(dif_cnt_bias[0, :], dif_cnt_bias[1, :] * SFave[1] * rps2dph, '-r')
plot_cxctime(dif_cnt_bias[0, :], dif_cnt_bias[2, :] * SFave[2] * rps2dph, '-g')
plot_cxctime(dif_cnt_bias[0, :], dif_cnt_bias[3, :] * SFave[3] * rps2dph, '-b')
plot_cxctime(dif_cnt_bias[0, :], dif_cnt_bias[4, :] * SFave[4] * rps2dph, '-m')
axis_arr = array(axis())
axis_arr[3] = axis_arr[3] + bias_diff_lim * rps2dph
axis(axis_arr)
legend(('Channel-1', 'Channel-2', 'Channel-3', 'Channel-4'), loc = 'best')
draw()
grid('on')
show(block = False)
if auto_save_plot:
    figfilename = 'Fig10_DiffChanBiasNPNT_%s_%s.png' % (interval, version)
    savefig(figfilename)

figure(11)
clf()
title('Histogram of Bias Difference')
xlabel('Bias Diff (deg/hr)')
ylabel('Number')
bias_diff = np.append(dif_cnt_bias[1, :] * SFave[1], dif_cnt_bias[2, :] * SFave[2])
bias_diff = np.append(bias_diff, dif_cnt_bias[3, :] * SFave[3])
bias_diff = np.append(bias_diff, dif_cnt_bias[4, :] * SFave[4]) * rps2dph
hist(bias_diff, 20);
draw()
grid('on')
show(block = False)
if auto_save_plot:
    figfilename = 'Fig11_HistBiasDiff_%s_%s.png' % (interval, version)
    savefig(figfilename)

figure(12)
clf()
title('YZ-Error vs Maneuver Angle')
xlabel('Angle (deg)')
ylabel('DeltaYZ (arcsec')
plot(ini2finang[1, :] * rad2deg, deltaYZ[1, :] * rad2asec, '.')
#axis([0.0, 180.0, 0.0, 120.0])
draw()
grid('on')
show(block = False)
if auto_save_plot:
    figfilename = 'Fig12_YZManErr_%s_%s.png' % (interval, version)
    savefig(figfilename)

figure(13)
clf()
title('Y-Error vs Y-Maneuver Angle')
xlabel('Y-Angle (deg)')
ylabel('Delta-Y (arcsec')
plot(ini2finvect[2, :] * rad2deg, deltavect[2, :] * rad2asec, '.')
#axis([-90.0, 90.0, -60.0, 60.0])
draw()
grid('on')
show(block = False)
if auto_save_plot:
    figfilename = 'Fig13_YManErrY_%s_%s.png' % (interval, version)
    savefig(figfilename)

figure(14)
clf()
title('Y-Error vs Z-Maneuver Angle')
xlabel('Z-Angle (deg)')
ylabel('Delta-Y (arcsec')
plot(ini2finvect[3, :] * rad2deg, deltavect[2, :] * rad2asec, '.')
#axis([-90.0, 90.0, -60.0, 60.0])
draw()
grid('on')
show(block = False)
if auto_save_plot:
    figfilename = 'Fig14_YManErrZ_%s_%s.png' % (interval, version)
    savefig(figfilename)

figure(15)
clf()
title('Z-Error vs Y-Maneuver Angle')
xlabel('Y-Angle (deg)')
ylabel('Delta-Z (arcsec')
plot(ini2finvect[2, :] * rad2deg, deltavect[3, :] * rad2asec, '.')
#axis([-90.0, 90.0, -60.0, 60.0])
draw()
grid('on')
show(block = False)
if auto_save_plot:
    figfilename = 'Fig15_ZManErrY_%s_%s.png' % (interval, version)
    savefig(figfilename)

figure(16)
clf()
title('Z-Error vs Z-Maneuver Angle')
xlabel('Z-Angle (deg)')
ylabel('Delta-Z (arcsec')
plot(ini2finvect[3, :] * rad2deg, deltavect[3, :] * rad2asec, '.')
#axis([-90.0, 90.0, -60.0, 60.0])
draw()
grid('on')
show(block = False)
if auto_save_plot:
    figfilename = 'Fig16_ZManErrZ_%s_%s.png' % (interval, version)
    savefig(figfilename)

# Write summary file
fsum = open(summaryfile, 'w')
fsum.write('getirudata.py version %s\n' % (version))
fsum.write('processing start time = %s\n' % (tstart))
fsum.write('processing stop  time = %s\n' % (tstop))
fsum.write('minimum NPNT duration = %f sec\n' % (npnt_min_dur))
fsum.write('Kalman filter converge time = %f sec\n' % (conv_time))
Tinitialfirst = Chandra.Time.DateTime(manvrtime[0, 0])
Tfinallast  = Chandra.Time.DateTime(manvrtime[1, -1])
fsum.write('Number of  input maneuvers = %d\n' % (num_nman0))
fsum.write('Number of output maneuvers = %d\n' % (num_man))
fsum.write('Initial time of first maneuver = %s\n' % (Tinitialfirst.date))
fsum.write('Final   time of last  maneuver = %s\n' % (Tfinallast.date))
fsum.write('output file = %s' % (outputfile))
fsum.close()

# write output to file 
#   finalpropquat:float(5, num_nman); propagated initial quaternion, time, q1, q2, q3, q4
#   deltaquat:    float(5, num_nman); difference between propagation and solution, time, q1, q2, q3, q4
#   initquat:     float(5, num_nman); initial pcad quaternion for each maneuver, time, q1, q2, q3, q4
#   finalquat:    float(5, num_nman); final pcad quaternion for each maneuver, time, q1, q2, q3, q4
#   manvrquat:    float(5, num_nman); observed rotation for each maneuver, time-diff, q1, q2, q3, q4
#   manvrtime:    float(2, num_nman); maneuver time interval start (0) and stop (1)
#   intratebody:  float(4, num_nman); integrated body rates over each maneuver
#   diffchancnts: float(5, num_nman); difference in counts across maneuver
#   ini2finquat:  float(5, num_nman); rotation quaternion from pcad initial & final attitudes, time, q1, q2, q3, q4
#   ini2finvect:  float(4, num_nman); rotation vector from initial to final quaternion, time, v1, v2, v3  
#   ini2finang:   float(2, num_nman); rotation angle; time, angle
#   finalpropquat:float(5, num_nman); propagated initial quaternion, time, q1, q2, q3, q4
#   deltaquat:    float(5, num_nman); difference between propagation and solution, time, q1, q2, q3, q4
#   ini2fintim:   float(1, num_nman); duration of maneuver interval (sec)
#   ini2finang:   float(1, num_nman); maneuver angle (deg)
#   sumprop:      float(3, 3, num_nman); sum of propagation matrices over maneuver
#   sumproprot: float(3, 9, num_nman); sum of prop-mat-func(rate) for maneuver
#   pcadbias_start:  float(4, num_nman); PCAD bias at start of maneuver (rad/sec)
#   ave_bias_before_nman: float(4, num_nman), average 4-vector bias before nman (cnts/sec)
#   ave_bias_after_nman:  float(4, num_nman), average 4-vector bias after nman (cnts/sec)
#   signcode : int(num_nman), integer coding signs of rotation around each channel axis

# Write column names
fout = open(outputfile, 'w')
fout.write(' num            start_time             stop_time ')
fout.write('   initquat1    initquat2    initquat3    initquat4 ')
fout.write('  finalquat1   finalquat2   finalquat3   finalquat4 ')
fout.write('  manvrquat1   manvrquat2   manvrquat3   manvrquat4 ')
fout.write('intratebody1 intratebody2 intratebody3 ')
fout.write('diffchancnt1 diffchancnt2 diffchancnt3 diffchancnt4 ')
fout.write('ini2finquat1 ini2finquat2 ini2finquat3 ini2finquat4 ')
fout.write('ini2finvect1 ini2finvect2 ini2finvect3 ')
fout.write('finpropquat1 finpropquat2 finpropquat3 finpropquat4 ')
fout.write('  deltaquat1   deltaquat2   deltaquat3   deltaquat4 ')
fout.write('  ini2fintim ')
fout.write('    ini2finang ')
if compute_batch:
    fout.write('   sumprop[1,1]    sumprop[1,2]    sumprop[1,3] ')
    fout.write('   sumprop[2,1]    sumprop[2,2]    sumprop[2,3] ')
    fout.write('   sumprop[3,1]    sumprop[3,2]    sumprop[3,3] ')
    fout.write('sumproprot[1,1] sumproprot[1,2] sumproprot[1,3] ')
    fout.write('sumproprot[1,4] sumproprot[1,5] sumproprot[1,6] ')
    fout.write('sumproprot[1,7] sumproprot[1,8] sumproprot[1,9] ')
    fout.write('sumproprot[2,1] sumproprot[2,2] sumproprot[2,3] ')
    fout.write('sumproprot[2,4] sumproprot[2,5] sumproprot[2,6] ')
    fout.write('sumproprot[2,7] sumproprot[2,8] sumproprot[2,9] ')
    fout.write('sumproprot[3,1] sumproprot[3,2] sumproprot[3,3] ')
    fout.write('sumproprot[3,4] sumproprot[3,5] sumproprot[3,6] ')
    fout.write('sumproprot[3,7] sumproprot[3,8] sumproprot[3,9] ')

fout.write('pcadbias_start1 pcadbias_start2 pcadbias_start3 ')
fout.write('cntratebiasbef1 cntratebiasbef2 cntratebiasbef3 cntratebiasbef4 ')
if write_signs:
    fout.write('cntratebiasaft1 cntratebiasaft2 cntratebiasaft3 cntratebiasaft4 signcode\n')
else:
    fout.write('cntratebiasaft1 cntratebiasaft2 cntratebiasaft3 cntratebiasaft4\n')

# Write numeric data for each maneuver
for n in range(num_man):
    start = Chandra.Time.DateTime(manvrtime[0, n])
    stop  = Chandra.Time.DateTime(manvrtime[1, n])
    fout.write('%4d %21s %21s ' % (n, start.date, stop.date))
    fout.write('%12.9f %12.9f %12.9f %12.9f ' % (initquat[1, n], initquat[2, n], initquat[3, n], initquat[4, n]))
    fout.write('%12.9f %12.9f %12.9f %12.9f ' % (finalquat[1, n], finalquat[2, n], finalquat[3, n], finalquat[4, n]))
    fout.write('%12.9f %12.9f %12.9f %12.9f ' % (manvrquat[1, n], manvrquat[2, n], manvrquat[3, n], manvrquat[4, n]))
    fout.write('%12.8f %12.8f %12.8f ' % (intratebody[1, n], intratebody[2, n], intratebody[3, n]))
    fout.write('%12d %12d %12d %12d ' % (diffchancnts[1, n], diffchancnts[2, n], diffchancnts[3, n], diffchancnts[4, n]))
    fout.write('%12.9f %12.9f %12.9f %12.9f ' % (ini2finquat[1, n], ini2finquat[2, n], ini2finquat[3, n], ini2finquat[4, n]))
    fout.write('%12.8f %12.8f %12.8f ' % (ini2finvect[1, n], ini2finvect[2, n], ini2finvect[3, n]))
    fout.write('%12.9f %12.9f %12.9f %12.9f ' % (finalpropquat[1, n], finalpropquat[2, n], finalpropquat[3, n], finalpropquat[4, n]))
    fout.write('%12.9f %12.9f %12.9f %12.9f ' % (deltaquat[1, n], deltaquat[2, n], deltaquat[3, n], deltaquat[4, n]))
    fout.write('%12.6f ' % (ini2finang[0, n]))
    fout.write('%15.6f ' % (ini2finang[1, n] * rad2deg))
    if compute_batch:
        fout.write('%15.8e %15.8e %15.8e ' % (sumprop[0, 0, n], sumprop[0, 1, n], sumprop[0, 2, n]))
        fout.write('%15.8e %15.8e %15.8e ' % (sumprop[1, 0, n], sumprop[1, 1, n], sumprop[1, 2, n]))
        fout.write('%15.8e %15.8e %15.8e ' % (sumprop[2, 0, n], sumprop[2, 1, n], sumprop[2, 2, n]))
        fout.write('%15.8e %15.8e %15.8e ' % (sumproprot[0, 0, n], sumproprot[0, 1, n], sumproprot[0, 2, n]))
        fout.write('%15.8e %15.8e %15.8e ' % (sumproprot[0, 3, n], sumproprot[0, 4, n], sumproprot[0, 5, n]))
        fout.write('%15.8e %15.8e %15.8e ' % (sumproprot[0, 6, n], sumproprot[0, 7, n], sumproprot[0, 8, n]))
        fout.write('%15.8e %15.8e %15.8e ' % (sumproprot[1, 0, n], sumproprot[1, 1, n], sumproprot[1, 2, n]))
        fout.write('%15.8e %15.8e %15.8e ' % (sumproprot[1, 3, n], sumproprot[1, 4, n], sumproprot[1, 5, n]))
        fout.write('%15.8e %15.8e %15.8e ' % (sumproprot[1, 6, n], sumproprot[1, 7, n], sumproprot[1, 8, n]))
        fout.write('%15.8e %15.8e %15.8e ' % (sumproprot[2, 0, n], sumproprot[2, 1, n], sumproprot[2, 2, n]))
        fout.write('%15.8e %15.8e %15.8e ' % (sumproprot[2, 3, n], sumproprot[2, 4, n], sumproprot[2, 5, n]))
        fout.write('%15.8e %15.8e %15.8e ' % (sumproprot[2, 6, n], sumproprot[2, 7, n], sumproprot[2, 8, n]))
    
    fout.write('%15.8e %15.8e %15.8e ' % (pcadbias_start[1, n], pcadbias_start[2, n], pcadbias_start[3, n]))
    fout.write('%15.8e %15.8e %15.8e %15.8e ' % (ave_bias_before_nman[1, n], ave_bias_before_nman[2, n], ave_bias_before_nman[3, n], ave_bias_before_nman[4, n]))
    if write_signs:
        fout.write('%15.8e %15.8e %15.8e %15.8e %2d\n' % (ave_bias_after_nman[1, n], ave_bias_after_nman[2, n], ave_bias_after_nman[3, n], ave_bias_after_nman[4, n], signcode[n]))
    else:
        fout.write('%15.8e %15.8e %15.8e %15.8e\n' % (ave_bias_after_nman[1, n], ave_bias_after_nman[2, n], ave_bias_after_nman[3, n], ave_bias_after_nman[4, n]))

fout.close()


