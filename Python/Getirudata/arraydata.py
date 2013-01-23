import Ska.Numpy
from pylab import *
import Chandra.Time
def findstrstartstop(data_vals,data_times,strval):
    """Function to find the start and stop times for string mneumonic of specific value.
    Goal is to have equal number of start and stop times with first start time before
    first stop time and last stop time after last start time.  Start time and stop 
    time cannot be equal.  
    Input:
    data_vals  : np.array of strings;
    data_times : np.array of floats; same size as data_vals
    strval     : string;  
    Output:
    sstimes    : float; np.array(num,2) of start and stop times
    Special cases:
    1. (data_vals[0] == strval) & (data_vals[1] == strval)
        data_times[0] is first start time in time interval
    2. (data_vals[0] == strval) & (data_vals[1] != strval)
        data_times[0] is not first start or first stop time
    3. (data_vals[-1] == strval) & (data_vals[-2] == strval)
        data_times[-1] is last stop time in time interval
    4. (data_vals[-1] == strval) & (data_vals[-2] != strval)
        data_times[-1] is not last start time or last stop time
    5. all(data_vals == strval) is True
        one time interval, data_times[0] is start time, data_times[-1] is stop time 
    6. any(data_vals == strval) is False
        There are no start and stop times in time interval."""
    numdata_vals = data_vals.shape[0]
    numdata_times = data_times.shape[0]    
    if (numdata_vals != numdata_times):
        message = 'number of data values (%d) not equal to number of data times (%d)' % (numdata_vals,numdata_times)
        return message
    else:
        numdata = numdata_vals
        
#   put in tests for string and float for input, somehow...
    
#   find start indices for strval
    start_index = find(((data_vals[1:] == strval) & (data_vals[:-1] != strval))) + 1
    num_start = start_index.shape[0]
    
#   find start indices for strval
    stop_index = find((data_vals[:-1] == strval) & (data_vals[1:] != strval))
    num_stop = stop_index.shape[0]
    
#   special case 1. (data_vals[0] == strval) and (data_vals[1] == strval)
    if ((data_vals[0] == strval) & (data_vals[1] == strval)):
        start_index = resize(start_index,(num_start+1,)) # lengthen array by 1
        start_index[1:] = start_index[:-1] # move data by 1 towards end
        start_index[0] = 0 # set initial index
        num_start = start_index.shape[0]
    
#   special case 2. (data_vals[0] == strval) and (data_vals[1] != strval)
    if (data_vals[0] == strval) & (data_vals[1] != strval):
        stop_index[:-1] = stop_index[1:] # move data by 1 towards beginning
        stop_index = resize(stop_index,(num_stop-1,)) # shorten array by 1, removing last stop index
        num_stop = stop_index.shape[0]
    
#   special case 3. (data_vals[-1] == strval) and (data_vals[-2] == strval)
    if ((data_vals[-1] == strval) & (data_vals[-2] == strval)):
        stop_index = resize(stop_index,(num_stop+1,)) # lengthen array by 1
        stop_index[-1] = data_vals.shape[0] - 1 # set final index
        num_stop = stop_index.shape[0]
    
#   special case 4. (data_vals[-1] == strval) and (data_vals[-2] != strval)
    if ((data_vals[-1] == strval) & (data_vals[-2] != strval)):
        start_index = resize(start_index,(num_start-1,)) # shorten array by 1, removing last start index
        num_start = start_index.shape[0]
    
#   obtain start and stop times
    if ((num_start > 0) and (num_start == num_stop)):
        sstimes = zeros((num_start,2))
        sstimes[:,0] = data_times[start_index]
        sstimes[:,1] = data_times[stop_index]
    else:
        sstimes = None

    ssindices = vstack((start_index,stop_index))
    ssindices = ssindices.transpose()

    return (ssindices,sstimes)

# getstrstartstop is an alternate version of findstrstartstop
def getstrstartstop(data_vals,data_times,strval):
    """Function to find the start and stop times for string mneumonic of specific value.
    Goal is to have equal number of start and stop times with first start time before
    first stop time and last stop time after last start time.  Start time and stop 
    time cannot be equal.  
    Input:
    data_vals  : np.array of strings;
    data_times : np.array of floats; same size as data_vals
    strval     : string;  
    Output:
    sstimes    : np.array(2,num) of start and stop times 
				ssindices  : np.array(2,num) of start and stop indices within data_vals
    Special cases:
    1. (data_vals[0] == strval) & (data_vals[1] == strval)
        data_times[0] is first start time in time interval
    2. (data_vals[0] == strval) & (data_vals[1] != strval)
        data_times[0] is not first start or first stop time
    3. (data_vals[-1] == strval) & (data_vals[-2] == strval)
        data_times[-1] is last stop time in time interval
    4. (data_vals[-1] == strval) & (data_vals[-2] != strval)
        data_times[-1] is not last start time or last stop time
    5. all(data_vals == strval) is True
        one time interval, data_times[0] is start time, data_times[-1] is stop time 
    6. any(data_vals == strval) is False
        There are no start and stop times in time interval."""
    numdata_vals = data_vals.shape[0]
    numdata_times = data_times.shape[0]    
    if (numdata_vals != numdata_times):
        message = 'number of data values (%d) not equal to number of data times (%d)' % (numdata_vals,numdata_times)
        return message
    else:
        numdata = numdata_vals
        
#   put in tests for string and float for input, somehow...
    
#   find start indices for strval
    start_index = find(((data_vals[1:] == strval) & (data_vals[:-1] != strval))) + 1
    num_start = start_index.shape[0]
    
#   find start indices for strval
    stop_index = find((data_vals[:-1] == strval) & (data_vals[1:] != strval))
    num_stop = stop_index.shape[0]
    
#   special case 1. (data_vals[0] == strval) and (data_vals[1] == strval)
    if ((data_vals[0] == strval) & (data_vals[1] == strval)):
        start_index = resize(start_index,(num_start+1,)) # lengthen array by 1
        start_index[1:] = start_index[:-1] # move data by 1 towards end
        start_index[0] = 0 # set initial index
        num_start = start_index.shape[0]
    
#   special case 2. (data_vals[0] == strval) and (data_vals[1] != strval)
    if (data_vals[0] == strval) & (data_vals[1] != strval):
        stop_index[:-1] = stop_index[1:] # move data by 1 towards beginning
        stop_index = resize(stop_index,(num_stop-1,)) # shorten array by 1, removing last stop index
        num_stop = stop_index.shape[0]
    
#   special case 3. (data_vals[-1] == strval) and (data_vals[-2] == strval)
    if ((data_vals[-1] == strval) & (data_vals[-2] == strval)):
        stop_index = resize(stop_index,(num_stop+1,)) # lengthen array by 1
        stop_index[-1] = data_vals.shape[0] - 1 # set final index
        num_stop = stop_index.shape[0]
    
#   special case 4. (data_vals[-1] == strval) and (data_vals[-2] != strval)
    if ((data_vals[-1] == strval) & (data_vals[-2] != strval)):
        start_index = resize(start_index,(num_start-1,)) # shorten array by 1, removing last start index
        num_start = start_index.shape[0]
    
#   obtain start and stop times
    if ((num_start > 0) and (num_start == num_stop)):
        sstimes = zeros((2,num_start))
        sstimes[0,:] = data_times[start_index]
        sstimes[1,:] = data_times[stop_index]
    else:
        sstimes = None

    ssindices = vstack((start_index,stop_index))
    ssindices = ssindices

    return (ssindices,sstimes)

