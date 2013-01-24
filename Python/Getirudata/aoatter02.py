#!/usr/bin/env python
"""
aoatter.py

Describe what this does.
"""
import argparse
import pprint

import matplotlib.pyplot as plt
import numpy as np

import Ska.engarchive.fetch as fetch
from Ska.Matplotlib import plot_cxctime
from Chandra.Time import DateTime

from arraydata import getstrstartstop


def get_opt():
    parser = argparse.ArgumentParser(description='Read attitude errors')
    parser.add_argument('--start', type=str,
                        default='2012:336',
                        help='Processing start date (DateTime format, default=2012:336)')
    parser.add_argument('--stop', type=str,
                        default='2013:021',
                        help='Processing stop date (DateTime format, default=2013:021)')
    parser.add_argument('--bad-aoatter-limit', type=float,
                        default=1.0,
                        help='Bad AOATTER limit (arcsec, default=1.0)')
    parser.add_argument('--settle-time', type=float,
                        default=300.0,
                        help='Kalman settling time (sec, default=300)')
    parser.add_argument('--dump-damp', type=float,
                        default=180.0,
                        help='Dump damping time (sec, default=180)')
    parser.add_argument('--adj-kalm', type=str,
                        default='True',
                        help='Adjust start time of NPNT to KALM time + settle time (default=True)')
    parser.add_argument('--filter-aoatter-bad-times', type=str,
                        default='True',
                        help='Remove AOATTER bad times (default=True)')
    parser.add_argument('--adj-mups-dump', type=str,
                        default='True',
                        help='Adjust  MUPS firing end time to remove damping motion (default=True)')
    parser.add_argument('--filter-mups-dump', type=str,
                        default='True',
                        help='Remove MUPS dumps (default=True)')
    parser.add_argument('--filter-rwbi-disable', type=str,
                        default='True',
                        help='Remove RW bias disable (SCS-107) (default=True)')
    parser.add_argument('--save-plots', type=str,
                        default='False',
                        help='Write plots to png files (default=False)')
    opt = parser.parse_args()
    return opt


def string_to_bool(val):
    val = val.lower()
    if val in ('n', 'f', 'false'):
        out = False
    elif val in ('y', 't', 'true'):
        out = True
    else:
        raise ValueError('Boolean flag value "{}" must be one of Y, N, T, F, True, or False'
                         .format(val))
    return out

# Get run time options and process accordingly
opt = get_opt()

tstart = DateTime(opt.start).date
tstop = DateTime(opt.stop).date

# aoatter bad if > bad_aoatter_limit, after other effects removed
bad_aoatter_limit = np.radians(opt.bad_aoatter_limit / 3600)

adj_kalm = string_to_bool(opt.adj_kalm)
filter_aoatter_bad_times = string_to_bool(opt.filter_aoatter_bad_times)
adj_mups_dump = string_to_bool(opt.adj_mups_dump)
filter_mups_dump = string_to_bool(opt.filter_mups_dump)
filter_rwbi_disable = string_to_bool(opt.filter_rwbi_disable)
save_plots = string_to_bool(opt.save_plots)

# Start processing
print 'Running aoatter.py with options:'
pprint.pprint(vars(opt))

# Read PCAD mode flag and get NPNT & NMAN start and stop times
print 'Requested time interval is %s to %s' % (tstart, tstop)
print 'Fetch AOPCADMD and AOACASEQ'
data = fetch.MSIDset(['AOPCADMD', 'AOACASEQ'], tstart, tstop, filter_bad=True)
aopcadmd_vals = np.array(data['AOPCADMD'].vals)
aopcadmd_times = data['AOPCADMD'].times[0:]
print 'PCAD mode times %s to %s' % (DateTime(aopcadmd_times[0]).date,
                                    DateTime(aopcadmd_times[-1]).date)
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
kalm_time_in_npnt = np.zeros(kalm_num, dtype=bool)
if adj_kalm:
    for n in range(npnt_num):
        kalm_time_in_npnt = ((npnt_times[0, n] < kalm_times[0, :])
                             & (kalm_times[0, :] < npnt_times[1, n]))
        if kalm_time_in_npnt.any():
            npnt_times[0, n] = kalm_times[0, kalm_time_in_npnt].max() + opt.settle_time

# Read PCAD attitude errors and times, compute YZ error
print 'Fetch AOATTER1, AOATTER2, and AOATTER3'
data = fetch.MSIDset(['AOATTER1', 'AOATTER2', 'AOATTER3'], tstart, tstop, filter_bad=True)
aoatter = np.array([data['AOATTER1'].times[0:],  # time of attitude error
                    data['AOATTER1'].vals,  # roll attitude error
                    data['AOATTER2'].vals,  # pitch attitude error
                    data['AOATTER3'].vals])  # yaw attitude error
print 'Attitude error times %s to %s' % (DateTime(aoatter[0, 0]).date,
                                         DateTime(aoatter[0, -1]).date)
aoatter_num = aoatter.shape[1]
print 'AOATTER number = %d' % aoatter_num
aoatteryz = np.sqrt(aoatter[2, :] ** 2 + aoatter[3, :] ** 2)

plt.figure(1)
plt.clf()
plt.title('YZ Attitude Error, NMAN & NPNT')
plt.xlabel('Time')
plt.ylabel('Attitude Err (arcsec)')
plot_cxctime(aoatter[0, :], aoatteryz * 180.0 / np.pi * 3600, '.r')
plt.grid('on')
plt.draw()
plt.show(block=False)
if save_plots:
    plt.savefig('YZ_Attitude_Error_{}-{}-{}.png'
                .format(tstart[0:4], tstart[5:8], tstart[9:11]))

# Find aoatter in NPNT
aoatter_after_npnt_start = np.zeros(aoatter_num, dtype=bool)  # pre-allocate array
aoatter_before_npnt_stop = np.zeros(aoatter_num, dtype=bool)  # pre-allocate array
aoatter_in_npnt = np.zeros(aoatter_num, dtype=bool)  # pre-allocate array
for n in range(npnt_num):
    aoatter_after_npnt_start = (aoatter[0, :] >= npnt_times[0, n])
    aoatter_before_npnt_stop = (aoatter[0, :] <= npnt_times[1, n])
    aoatter_in_npnt = aoatter_in_npnt | (aoatter_after_npnt_start & aoatter_before_npnt_stop)
print 'Number of aoatter in NPNT = %d' % aoatter_in_npnt.sum()

plt.figure(2)
plt.clf()
plt.title('YZ Attitude Error during NPNT')
plt.xlabel('Time')
plt.ylabel('Attitude Err (arcsec)')
plot_cxctime(aoatter[0, aoatter_in_npnt], aoatteryz[aoatter_in_npnt] * 180.0 / np.pi * 3600, '.r')
plt.grid('on')
plt.draw()
plt.show(block=False)
if save_plots:
    plt.savefig('YZ_Attitude_Error_NPNT_{}-{}-{}.png'
                .format(tstart[0:4], tstart[5:8], tstart[9:11]))

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
    dump_times[1, :] = dump_times[1, :] + opt.dump_damp
(rwbi_indices, rwbi_times) = getstrstartstop(aorwbias_vals, aorwbias_times, 'DISA')
rwbi_num = rwbi_indices.shape[1]
print 'Number of RW bias disables = %d' % rwbi_num

# Find aoatter in mups dump
aoatter_after_dump_start = np.zeros(aoatter_num, dtype=bool)  # pre-allocate array
aoatter_before_dump_stop = np.zeros(aoatter_num, dtype=bool)  # pre-allocate array
aoatter_in_dump = np.zeros(aoatter_num, dtype=bool)  # pre-allocate array
if filter_mups_dump:
    if (dump_num > 0):
        for n in range(dump_num):
            aoatter_after_dump_start = (aoatter[0, :] >= dump_times[0, n])
            aoatter_before_dump_stop = (aoatter[0, :] <= dump_times[1, n])
            aoatter_in_dump = aoatter_in_dump | (aoatter_after_dump_start &
                                                 aoatter_before_dump_stop)
        print 'Number of aoatter during dumps = %d' % aoatter_in_dump.sum()

# Find aoatter in RW bias disable
aoatter_after_disa_start = np.zeros(aoatter_num, dtype=bool)  # pre-allocate array
aoatter_before_disa_stop = np.zeros(aoatter_num, dtype=bool)  # pre-allocate array
aoatter_in_disa = np.zeros(aoatter_num, dtype=bool)  # pre-allocate array
if filter_rwbi_disable:
    if (rwbi_num > 0):
        for n in range(rwbi_num):
            aoatter_after_disa_start = (aoatter[0, :] >= rwbi_times[0, n])
            aoatter_before_disa_stop = (aoatter[0, :] <= rwbi_times[1, n])
            aoatter_in_disa = aoatter_in_disa | (aoatter_after_disa_start &
                                                 aoatter_before_disa_stop)
        print 'Number of aoatter during RW bias disable = %d' % aoatter_in_disa.sum()

# Get filtered aoatter data
aoatter_filter = aoatter_in_npnt & (~ aoatter_in_dump) & (~ aoatter_in_disa)
print 'Number of aoatter in NPNT & not in dump or disable = %d' % aoatter_filter.sum()

plt.figure(3)
plt.clf()
plt.title('Attitude Error Filtered for MUPS Dump & RW Bias DISA')
plt.xlabel('Time')
plt.ylabel('Attitude Err (arcsec)')
plot_cxctime(aoatter[0, aoatter_filter], aoatteryz[aoatter_filter] * 180.0 / np.pi * 3600, '.r')
plt.grid('on')
plt.draw()
plt.show(block=False)
if save_plots:
    plt.savefig('YZ_Attitude_Error_Filtered_More_{}-{}-{}.png'
                .format(tstart[0:4], tstart[5:8], tstart[9:11]))

# sys.exit(0)

if filter_aoatter_bad_times:
    textin = open('aoatter_bad_times.dat', 'rU').readlines()
    timelist = [x.split() for x in textin]
    startstr = np.array([x[0] for x in timelist])
    stopstr = np.array([x[1] for x in timelist])
    num_bad = startstr.shape[0]
    bad_times = np.zeros((2, num_bad))
    bad_times[0, :] = DateTime(startstr).secs
    bad_times[1, :] = DateTime(stopstr).secs
    aoatter_after_bad_start = np.zeros(aoatter_num, dtype=bool)  # pre-allocate array
    aoatter_before_bad_stop = np.zeros(aoatter_num, dtype=bool)  # pre-allocate array
    aoatter_in_bad = np.zeros(aoatter_num, dtype=bool)  # pre-allocate array
    for n in range(num_bad):
        aoatter_after_bad_start = (aoatter[0, :] >= bad_times[0, n])
        aoatter_before_bad_stop = (aoatter[0, :] <= bad_times[1, n])
        aoatter_in_bad = aoatter_in_bad | (aoatter_after_bad_start & aoatter_before_bad_stop)
    print 'Number of aoatter in bad = %d' % aoatter_in_bad.sum()
    aoatter_filter = aoatter_filter & (~ aoatter_in_bad)

plt.figure(4)
plt.clf()
plt.title('Attitude Error Filtered with Bad Times File')
plt.xlabel('Time')
plt.ylabel('Attitude Err (arcsec)')
plot_cxctime(aoatter[0, aoatter_filter], aoatteryz[aoatter_filter] * 180.0 / np.pi * 3600, '.r')
plt.grid('on')
plt.draw()
plt.show(block=False)
if save_plots:
    plt.savefig('YZ_Attitude_Error_Filtered_Even_More_{}-{}-{}.png'
                .format(tstart[0:4], tstart[5:8], tstart[9:11]))

print 'Done!'

