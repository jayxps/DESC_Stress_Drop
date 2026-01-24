"""
Objective function for differential evolution optimization.

Author: Rainer Storn (original MATLAB), converted to Python3 and modified by Jiewen Zhang
Description: Implements the cost function to be minimized.
"""

import numpy as np


def objfun(FVr_temp, S_struct, S):
    """
    Objective function for differential evolution optimization.
    
    Parameters
    ----------
    FVr_temp : array-like
        Parameter vector (log10 stress drop values)
    S_struct : dict
        Contains a variety of parameters (for compatibility)
    S : dict
        Contains problem-specific data including:
        - nfreqg: frequency indices
        - freqg: frequency values
        - nlow: low frequency indices
        - evspec_DE: event spectra data
        - specall: all spectra array
        - ibin: bin indices
        - fmom: seismic moments
        - nfit: fitting indices
        - min_nspec: minimum number of spectra per bin
    
    Returns
    -------
    S_MSE : dict
        Dictionary containing:
        - I_nc: Number of constraints (0)
        - FVr_ca: Constraint values (0)
        - I_no: Number of objectives (1)
        - FVr_oa: Objective function values
    power : dict
        Dictionary containing:
        - ECSs_med: Median ECS values
        - ECSs: ECS values for all events
        - fc: Corner frequencies
    """
    
    # Extract variables
    nfreqg = S['nfreqg']
    freqg = S['freqg']
    nlow = S['nlow']
    evspec = S['evspec_DE']
    specall = S['specall']
    ibin = S['ibin']
    fmom = S['fmom']
    nfit = S['nfit']
    min_nspec = S['min_nspec']
    
    # Filter events with valid magnitude residuals
    ng = np.where(np.array([ev['qmagresid'] for ev in evspec]) <= 100)[0]
    evspec = [evspec[i] for i in ng]
    specall = specall[ng, :]
    ibin = ibin[ng]
    fmom = fmom[ng]
    
    beta = 3464
    fact = (0.42 * beta) ** 3
    
    # Optimize for matrix computation
    # Convert log10 stress drop to Pa
    dsigtarg = (10.0 ** FVr_temp[ibin]) * 1e6
    
    # Vectorized corner frequency calculation
    fc = ((dsigtarg * fact) / fmom) ** (1.0 / 3.0)
    
    # Compute predicted spectra
    freqg_mat = np.tile(freqg, (len(evspec), 1))
    fc_mat = np.tile(fc[:, np.newaxis], (1, len(freqg)))
    pred = -np.log10(1 + (freqg_mat / fc_mat) ** 2)
    
    # Calculate offset
    yoff = np.mean(specall[:, nfreqg[nlow]] - pred[:, nfreqg[nlow]], axis=1)
    
    # Calculate Empirical Coda Source spectra (ECS)
    ECSs = specall - yoff[:, np.newaxis] - pred
    
    # Calculate median ECS for each bin
    ibin_u = np.unique(ibin[~np.isnan(ibin)])
    ECSs_med = np.full((len(ibin_u), ECSs.shape[1]), np.nan)
    
    for iibin in range(len(ibin_u)):
        if np.sum(ibin == ibin_u[iibin]) >= min_nspec:
            ECSs_med[iibin, :] = np.median(ECSs[ibin == ibin_u[iibin], :], axis=0)
    
    # Calculate standard deviation across bins at fitting frequencies
    ST = np.sum(np.nanstd(ECSs_med[:, nfit], axis=0))
    
    # Store outputs
    power = {
        'ECSs_med': ECSs_med,
        'ECSs': ECSs,
        'fc': fc
    }
    
    FVr_oa = [1 / ST]

    # Output structure
    S_MSE = {
        'I_nc': 0,           # no constraints
        'FVr_ca': 0,         # no constraint array
        'I_no': 1,           # number of objectives (costs)
        'FVr_oa': FVr_oa         # objective function value
    }
    
    return S_MSE, power
