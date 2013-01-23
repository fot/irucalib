#!/usr/bin/env python
# aoatter.py

print 'Running aoatter.py'

import Ska.engarchive.fetch as fetch
from Ska.Matplotlib import plot_cxctime
import Ska.Numpy
from pylab import *
import Chandra.Time
from arraydata import getstrstartstop
#import sys
import os
from math import *
from quatfunc import *

asec2rad = pi / 180.0 / 3600.0 # arcsec to radians multiplier
bad_aoatter_limit = 1.0 * asec2rad # aoatter bad if > bad_aoatter_limit, after other effects removed
settle_time = 300.0 # Kalman settling time in sec
radpersec2degperhr = 180.0 / pi * 3600 # rate units multiplier
dump_damp = 180.0 # dump damping time (sec)

#tstart = '2012:062:16:00:00.000' # Just after Mmat(5) update
#tstop  = '2012:150:03:33:00.000' # Just before safe mode 5
#tstart = '2012:152:00:00:00.000' # Just after safe mode 5
#tstop  = '2012:363:00:00:00.000' # 
tstart = '2012:336:00:00:00.000' # 1 Dec 2012
tstop  = '2013:021:00:00:00.000'

adj_kalm = True # Adjust start time of NPNT to KALM time + settle time, if present
filter_aoatter_bad_times = True # Remove AOATTER bad times
adj_mups_dump = True # Adjust  MUPS firing end time to remove damping motion
filter_mups_dump = True # Remove MUPS dumps
filter_rwbi_disable = True # Remove RW bias disable (SCS-107)
save_plots = False # Write plots to png file

# Read PCAD mode flag and get NPNT & NMAN start and stop times
print 'Requested time interval is %s to %s' %(tstart, tstop)
print 'Fetch AOPCADMD and AOACASEQ'
data = fetch.MSIDset(['AOPCADMD', 'AOACASEQ'], tstart, tstop, filter_bad=True)
aopcadmd_vals = np.array(data['AOPCADMD'].vals)
aopcadmd_times = data['AOPCADMD'].times[0:]
print 'PCAD mode times %s to %s' %(Chandra.Time.DateTime(aopcadmd_times[0]).date, Chandra.Time.DateTime(aopcadmd_times[-1]).date)
aoacaseq_vals = np.array(data['AOACASEQ'].vals)
aoacaseq_times = data['AOACASEQ'].times[0:]
(npnt_indices, npnt_times) = getstrstartstop(aopcadmd_vals, aopcadmd_times, 'NPNT')
npnt_num = npnt_times.shape[1]
print 'NPNT number = %d' % npnt_num
(nman_indices, nman_times) = getstrstartstop(aopcadmd_vals, aopcadmd_times, 'NMAN')
nman_num = nman_times.shape[1]
print 'NMAN number = %d' % nman_num
(kalm_indices, kalm_times) = getstrstartstop(aoacaseq_vals, aoacaseq_times, 'KALM')
kalm_num = kalm_indices.shape[1]
print 'KALMAN number = %d' % kalm_num

# Adjust start time of NPNT to KALM time, if present
kalm_time_in_npnt = zeros(kalm_num, dtype=bool)
if adj_kalm:
    for n in range(npnt_num):
        kalm_time_in_npnt = ((npnt_times[0,n] < kalm_times[0,:]) & (kalm_times[0,:] < npnt_times[1,n]))
        if kalm_time_in_npnt.any():
            npnt_times[0,n] = kalm_times[0,kalm_time_in_npnt].max() + settle_time

# Read PCAD attitude errors and times, compute YZ error
print 'Fetch AOATTER1, AOATTER2, and AOATTER3'
data = fetch.MSIDset(['AOATTER1', 'AOATTER2', 'AOATTER3'], tstart, tstop, filter_bad=True)
aoatter = np.array([data['AOATTER1'].times[0:], # time of attitude error
                    data['AOATTER1'].vals,      # roll attitude error
                    data['AOATTER2'].vals,      # pitch attitude error
                    data['AOATTER3'].vals])     # yaw attitude error
print 'Attitude error times %s to %s' %(Chandra.Time.DateTime(aoatter[0, 0]).date, Chandra.Time.DateTime(aoatter[0, -1]).date)
aoatter_num = aoatter.shape[1]
print 'AOATTER number = %d' % aoatter_num
aoatteryz = np.sqrt(aoatter[2, :]**2 + aoatter[3, :]**2)

figure(1)
clf()
title('YZ Attitude Error, NMAN & NPNT')
xlabel('Time')
ylabel('Attitude Err (arcsec)')
plot_cxctime(aoatter[0, :], aoatteryz * 180.0 / pi * 3600, '.r')
grid('on')
draw()
show(block=False)
if save_plots:
    savefig('YZ_Attitude_Error_'+tstart[0:4]+'-'+tstart[5:8]+'-'+tstart[9:11]+'.png')

# Find aoatter in NPNT
aoatter_after_npnt_start = zeros(aoatter_num, dtype=bool) # pre-allocate array
aoatter_before_npnt_stop = zeros(aoatter_num, dtype=bool) # pre-allocate array
aoatter_in_npnt = zeros(aoatter_num, dtype=bool) # pre-allocate array
for n in range(npnt_num):
    aoatter_after_npnt_start = (aoatter[0,:] >= npnt_times[0,n])
    aoatter_before_npnt_stop = (aoatter[0,:] <= npnt_times[1,n])
    aoatter_in_npnt = aoatter_in_npnt | (aoatter_after_npnt_start & aoatter_before_npnt_stop)
print 'Number of aoatter in NPNT = %d' % aoatter_in_npnt.sum()

figure(2)
clf()
title('YZ Attitude Error during NPNT')
xlabel('Time')
ylabel('Attitude Err (arcsec)')
plot_cxctime(aoatter[0, aoatter_in_npnt], aoatteryz[aoatter_in_npnt] * 180.0 / pi * 3600, '.r')
grid('on')
draw()
show(block=False)
if save_plots:
    savefig('YZ_Attitude_Error_NPNT_'+tstart[0:4]+'-'+tstart[5:8]+'-'+tstart[9:11]+'.png')

# Get data for mups dump and RW bias disable
print 'Fetch AOUNLOAD, and AORWBIAS'
data = fetch.MSIDset(['AOUNLOAD', 'AORWBIAS'], tstart, tstop, filter_bad=True)
aounload_vals = np.array(data['AOUNLOAD'].vals)
aounload_times = data['AOUNLOAD'].times[0:]
aorwbias_vals = np.array(data['AORWBIAS'].vals)
aorwbias_times = data['AORWBIAS'].times[0:]
(dump_indices, dump_times) = getstrstartstop(aounload_vals, aounload_times, 'GRND')
dump_num = dump_indices.shape[1]
print 'Number of momentum dumps = %d' % dump_num
# Adjust  MUPS firing end time to remove damping motion
if adj_mups_dump:
    dump_times[1,:] = dump_times[1,:] + dump_damp
(rwbi_indices, rwbi_times) = getstrstartstop(aorwbias_vals, aorwbias_times, 'DISA')
rwbi_num = rwbi_indices.shape[1]
print 'Number of RW bias disables = %d' % rwbi_num

# Find aoatter in mups dump
aoatter_after_dump_start = zeros(aoatter_num, dtype=bool) # pre-allocate array
aoatter_before_dump_stop = zeros(aoatter_num, dtype=bool) # pre-allocate array
aoatter_in_dump = zeros(aoatter_num, dtype=bool) # pre-allocate array
if filter_mups_dump:
    if (dump_num > 0):
        for n in range(dump_num):
            aoatter_after_dump_start = (aoatter[0,:] >= dump_times[0,n])
            aoatter_before_dump_stop = (aoatter[0,:] <= dump_times[1,n])
            aoatter_in_dump = aoatter_in_dump | (aoatter_after_dump_start & aoatter_before_dump_stop)
        print 'Number of aoatter during dumps = %d' % aoatter_in_dump.sum()

# Find aoatter in RW bias disable
aoatter_after_disa_start = zeros(aoatter_num, dtype=bool) # pre-allocate array
aoatter_before_disa_stop = zeros(aoatter_num, dtype=bool) # pre-allocate array
aoatter_in_disa = zeros(aoatter_num, dtype=bool) # pre-allocate array
if filter_rwbi_disable:
    if (rwbi_num > 0):
        for n in range(rwbi_num):
            aoatter_after_disa_start = (aoatter[0,:] >= rwbi_times[0,n])
            aoatter_before_disa_stop = (aoatter[0,:] <= rwbi_times[1,n])
            aoatter_in_disa = aoatter_in_disa | (aoatter_after_disa_start & aoatter_before_disa_stop)
        print 'Number of aoatter during RW bias disable = %d' % aoatter_in_disa.sum()

# Get filtered aoatter data
aoatter_filter = aoatter_in_npnt & (~ aoatter_in_dump) & (~ aoatter_in_disa)
print 'Number of aoatter in NPNT & not in dump or disable = %d' % aoatter_filter.sum()

figure(3)
clf()
title('Attitude Error Filtered for MUPS Dump & RW Bias DISA')
xlabel('Time')
ylabel('Attitude Err (arcsec)')
plot_cxctime(aoatter[0, aoatter_filter], aoatteryz[aoatter_filter] * 180.0 / pi * 3600, '.r')
grid('on')
draw()
show(block=False)
if save_plots:
    savefig('YZ_Attitude_Error_Filtered_More_'+tstart[0:4]+'-'+tstart[5:8]+'-'+tstart[9:11]+'.png')

# sys.exit(0)

if filter_aoatter_bad_times:
    textin = open('aoatter_bad_times.dat','rU').readlines()
    timelist = [x.split() for x in textin]
    startstr = array([x[0] for x in timelist])
    stopstr = array([x[1] for x in timelist])
    num_bad = startstr.shape[0]
    bad_times = zeros((2,num_bad))
    bad_times[0, :] = Chandra.Time.DateTime(startstr).secs
    bad_times[1, :] = Chandra.Time.DateTime(stopstr).secs
    aoatter_after_bad_start = zeros(aoatter_num, dtype=bool) # pre-allocate array
    aoatter_before_bad_stop = zeros(aoatter_num, dtype=bool) # pre-allocate array
    aoatter_in_bad = zeros(aoatter_num, dtype=bool) # pre-allocate array
    for n in range(num_bad):
        aoatter_after_bad_start = (aoatter[0,:] >= bad_times[0,n])
        aoatter_before_bad_stop = (aoatter[0,:] <= bad_times[1,n])
        aoatter_in_bad = aoatter_in_bad | (aoatter_after_bad_start & aoatter_before_bad_stop)
    print 'Number of aoatter in bad = %d' % aoatter_in_bad.sum()
    aoatter_filter = aoatter_filter & (~ aoatter_in_bad)

figure(4)
clf()
title('Attitude Error Filtered with Bad Times File')
xlabel('Time')
ylabel('Attitude Err (arcsec)')
plot_cxctime(aoatter[0, aoatter_filter], aoatteryz[aoatter_filter] * 180.0 / pi * 3600, '.r')
grid('on')
draw()
show(block=False)
if save_plots:
    savefig('YZ_Attitude_Error_Filtered_Even_More_'+tstart[0:4]+'-'+tstart[5:8]+'-'+tstart[9:11]+'.png')

print 'Done!'

