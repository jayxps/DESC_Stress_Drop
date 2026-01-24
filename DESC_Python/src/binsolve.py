import numpy as np


def binsolve(a, npts, index, nbin, itype):
    """
    Computes mean or median within specified bins of data array.
    
    Parameters
    ----------
    a : array-like
        Input data array
    npts : int
        Number of points
    index : array-like
        Bin indices for each point
    nbin : int
        Number of bins
    itype : int
        Type of operation: 1 for median, 2 for mean, 3 for weighted mean
    
    Returns
    -------
    amed : ndarray
        Binned mean/median values
    temp : ndarray
        Temporary array for weighted calculations
    atemp : ndarray
        Temporary weights array
    d : ndarray
        Temporary distance/difference array
    ninbin : ndarray
        Number of points in each bin
    """
    amed = np.zeros(nbin)
    ninbin = np.zeros(npts, dtype=int)
    temp = np.zeros(npts)
    atemp = np.zeros(npts)
    d = np.zeros(npts)
    
    if itype != 3:
        for ibin in range(nbin):
            nin = np.where(index == ibin + 1)[0]  # Convert to 1-based indexing
            ninbin[ibin] = len(nin)
            
            if ninbin[ibin] == 0:
                amed[ibin] = 0
                continue
            
            if itype == 1:  # compute median
                amed[ibin] = np.median(a[nin])
            elif itype == 2:  # compute mean
                amed[ibin] = np.mean(a[nin])
    else:
        # Weighted mean calculation
        xgap = 0.2
        xgapinv = 1.0 / xgap
        nit = 5
        atemp = np.ones(npts)
        
        for it in range(nit):
            for ibin in range(nbin):
                nin = np.where(index == ibin + 1)[0]  # Convert to 1-based indexing
                ninbin[ibin] = len(nin)
                
                if ninbin[ibin] == 0:
                    amed[ibin] = 0
                    continue
                
                temp = atemp[nin] * a[nin]
                amed[ibin] = np.mean(temp)
            
            if it != nit - 1:
                d = np.abs(a - amed[ibin])
                nle = np.where(d <= xgap)[0]
                atemp = 1.0 / d
                atemp[nle] = xgapinv
    
    return amed, temp, atemp, d, ninbin
