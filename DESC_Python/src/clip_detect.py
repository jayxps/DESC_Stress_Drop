import numpy as np


def clip_detect(data, F=None):
    """
    Simple function to check for clipped waveforms.
    
    Computes the maximum and minimum value in data. Then counts the number 
    of data points within a fraction F of maxval and minval:
        - (data >= (1-F)*maxval) --> NF1max
        - (data <= (1-F)*minval) --> NF1min
    
    If both of these counts are greater than the counts of the next
    fraction down:
        - ((1-2*F)*maxval <= data < (1-F)*maxval) --> NF2max
        - ((1-2*F)*minval >= data > (1-F)*minval) --> NF2min
    
    then flag = 1.0 (likely clip). If only one count is greater, set flag = 0.5
    (possible clipping). If none, set flag = 0.0 (no clipping).
    
    Parameters
    ----------
    data : array-like
        Input waveform data
    F : float, optional
        Fraction parameter (default: 1/3)
    
    Returns
    -------
    clip : float
        Clipping indicator value
    """
    # Default value of F is 1/3
    if F is None or np.abs(F) > 0.5:
        F = 1.0 / 3.0
    
    # Convert to numpy array
    data = np.asarray(data)
    
    # Get max and min values
    maxval = np.max(data)
    minval = np.min(data)
    
    # Indices of data points in different intervals
    I1max = data >= (1 - F) * maxval
    I2max = data >= (1 - 2 * F) * maxval
    I1min = data <= (1 - F) * minval
    I2min = data <= (1 - 2 * F) * minval
    
    # Counts of data points in different intervals
    NF1max = np.sum(I1max)
    NF2max = np.sum(I2max & ~I1max)
    NF2min = np.sum(I2min & ~I1min)
    NF1min = np.sum(I1min)
    
    # Classify data (original commented-out logic)
    # if (NF2max < NF1max) and (NF2min < NF1min):  # likely clipped
    #     clip = 1
    # elif (NF2max < NF1max) or (NF2min < NF1min):  # maybe clipped
    #     clip = 0.5
    # else:  # no worries
    #     clip = 0
    
    # Compute clip_val, modified by Jiewen Zhang on Jun 6, 2020
    clip = (NF1max / NF2max + NF1min / NF2min) / 2.0
    
    return clip