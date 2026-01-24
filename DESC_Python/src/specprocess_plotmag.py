import numpy as np
import matplotlib.pyplot as plt
from scipy.io import savemat
from scipy import stats


def specprocess_plotmag(evspec, freqlim, nspecmin, savefig=True):
    """
    Calibrate catalog magnitude to spectral amplitude via robust linear fit.

    Parameters
    ----------
    evspec : list of dict
        Each element must have keys: 'freq', 'spec', 'nspec', 'qid', 'qmag'.
    freqlim : tuple or list
        Frequency bounds (low, high) for computing mean spectral amplitude.
    nspecmin : int
        Minimum number of spectra required for inclusion.
    savefig : bool, optional
        Save the fit figure to disk (default: True).

    Returns
    -------
    y0 : float
        Intercept of the magnitude fit.
    slopeline : float
        Slope of the magnitude fit.
    """

    freq = np.asarray(evspec[0]['freq'])

    # Frequency mask
    nfx = (freq >= freqlim[0]) & (freq <= freqlim[1])

    # Indices meeting minimum nspec
    nspec_arr = np.array([ev['nspec'] for ev in evspec])
    ng = np.where(nspec_arr >= nspecmin)[0]

    # Compute mean spectral amplitude in band
    specamp = np.array([np.mean(evspec[i]['spec'][nfx]) for i in ng])
    mag = np.array([evspec[i]['qmag'] for i in ng])

    # All qualifying for plotting background
    ng_all = ng
    specamp_all = np.array([np.mean(evspec[i]['spec'][nfx]) for i in ng_all])

    # Use robust regression (MATLAB's robustfit equivalent)
    # robustfit uses iteratively reweighted least squares with bisquare weights
    # scipy.stats doesn't have robustfit, but we can use HuberRegressor or implement it
    # For now, use a simple robust regression via Theil-Sen estimator as close approximation
    # However, to match MATLAB exactly, we need to implement bisquare robust regression
    
    # MATLAB robustfit with default bisquare weights:
    # Iteratively reweighted least squares
    X = np.column_stack([np.ones_like(specamp), specamp])
    y = mag
    
    # Initialize with OLS
    beta = np.linalg.lstsq(X, y, rcond=None)[0]
    
    # Iterative robust regression (bisquare weights, matching MATLAB's robustfit)
    tune = 4.685  # Tuning constant for bisquare
    max_iter = 50
    tol = 1e-6
    
    for iteration in range(max_iter):
        beta_old = beta.copy()
        
        # Compute residuals
        r = y - X @ beta
        
        # Compute robust scale estimate (MAD)
        mad = np.median(np.abs(r - np.median(r)))
        s = mad / 0.6745  # Convert MAD to standard deviation estimate
        
        if s < 1e-10:
            s = 1e-10
        
        # Compute bisquare weights
        u = r / (tune * s)
        weights = np.where(np.abs(u) < 1, (1 - u**2)**2, 0)
        
        # Weighted least squares
        W = np.diag(weights)
        beta = np.linalg.lstsq(X.T @ W @ X, X.T @ W @ y, rcond=None)[0]
        
        # Check convergence
        if np.max(np.abs(beta - beta_old)) < tol:
            break
    
    y0 = beta[0]
    slopeline = beta[1]

    specampx = np.arange(-6, 20.1, 0.1)
    ypredx = y0 + slopeline * specampx

    # Save specamp and qid for compatibility
    qid = np.array([evspec[i]['qid'] for i in ng])
    savemat('specamp.mat', {'specamp': specamp, 'qid': qid})

    # Plot
    if savefig:
        plt.ioff()
        fig = plt.figure(999)
        plt.plot(specamp_all, [evspec[i]['qmag'] for i in ng_all], '.', markersize=15, color=[0.7, 0.7, 0.7])
        plt.plot(specamp, mag, 'k.', markersize=15)
        plt.plot(specampx, ypredx, linewidth=2)
        x_min = specamp.min() - 0.5
        x_max = specamp.max() + 0.5
        plt.xlim([x_min, x_max])
        titlex = f"freqlim={freqlim[0]:.2f}-{freqlim[1]:.2f} slope={slopeline:.4f} y0={y0:.4f}"
        plt.title(titlex)
        plt.xlabel('relative log magnitude')
        plt.ylabel('catalog magnitude')
        figout = f"magfit-freq-{freqlim[0]}-{freqlim[1]}.jpg"
        plt.savefig(figout)
        plt.close(fig)

    return y0, slopeline
