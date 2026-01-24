import numpy as np


def specprocess_binsolve(a, index, nbin, itype):
    """
    Computes mean or median within specified bins of data array.
    
    Parameters
    ----------
    a : ndarray
        Input data array of shape (M, N)
    index : array-like
        Bin indices for each point
    nbin : int
        Number of bins
    itype : int
        Type of operation: 1 for median, 2 for mean, 3 for weighted mean
    
    Returns
    -------
    amed : ndarray
        Binned mean/median values of shape (M, nbin)
    ninbin : ndarray
        Number of points in each bin of shape (M, nbin)
    """
    M, N = a.shape
    
    amed = np.zeros((M, nbin))
    ninbin = np.zeros((M, nbin), dtype=int)
    
    if itype != 3:
        for ibin in range(nbin):
            nin = np.where(index == ibin + 1)[0]  # Convert to 1-based indexing
            ninbin[0, ibin] = len(nin)  # MATLAB: ninbin(ibin) = length(nin) - first row only
            
            if len(nin) == 0:
                amed[:, ibin] = 0
                continue
            
            if itype == 1:  # compute median
                amed[:, ibin] = np.median(a[:, nin], axis=1)
            elif itype == 2:  # compute mean
                amed[:, ibin] = np.mean(a[:, nin], axis=1)
    else:
        # Weighted mean calculation
        xgap = 0.5
        xgapinv = 1.0 / xgap
        nit = 5
        atemp = np.ones_like(a)
        
        for it in range(nit):
            atemp_new = np.ones_like(a)
            
            for ibin in range(nbin):
                nin = np.where(index == ibin + 1)[0]  # Convert to 1-based indexing
                ninbin[0, ibin] = len(nin)  # MATLAB: ninbin(ibin) = length(nin) - first row only
                
                if len(nin) == 0:
                    amed[:, ibin] = 0
                    continue
                
                temp = atemp[:, nin] * a[:, nin]
                sumw = np.sum(atemp[:, nin], axis=1)
                amed[:, ibin] = np.sum(temp, axis=1) / sumw
                
                if it != nit - 1:
                    d = np.abs(a[:, nin] - np.tile(amed[:, ibin:ibin+1], (1, len(nin))))
                    # Match MATLAB exactly: atempx = 1./d; atempx(d<=xgap) = xgapinv;
                    # Suppress warnings for divide by zero (MATLAB handles this silently)
                    with np.errstate(divide='ignore', invalid='ignore'):
                        atempx = 1.0 / d
                    atempx[d <= xgap] = xgapinv
                    atemp_new[:, nin] = atempx
            
            atemp = atemp_new
    
    return amed, ninbin
