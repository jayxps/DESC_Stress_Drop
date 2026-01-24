"""
Initialization and run script for the differential evolution optimizer.

Author: Rainer Storn (original MATLAB), converted to Python3 and modified by Jiewen Zhang
"""

import numpy as np
from .deopt import deopt


def Rundeopt(S):
    """
    Script for the initialization and run of the differential evolution optimizer.
    
    Parameters
    ----------
    S : dict
        Structure containing:
        - nmagg: magnitude bin indices
        - log10sig: [min, max] log10 stress drop range
        - evspec_DE: event spectra data
        - Other problem-specific parameters
    
    Returns
    -------
    FVr_x : ndarray
        Optimized parameter vector (log10 stress drop values)
    ECSs : dict
        Dictionary containing:
        - ECSs_med: Median ECS values
        - ECSs: ECS values for all events
    """
    
    # Value To Reach (stop when objective function < F_VTR)
    F_VTR = -10000
    
    # Extract variables from struct S
    nmagg = S['nmagg']
    log10sig = S['log10sig']
    
    # I_D: number of parameters of the objective function
    I_D = len(nmagg)
    
    # FVr_minbound, FVr_maxbound: vector of lower and upper bounds of initial population
    # Note: these are bound constraints when I_bnd_constr = 1
    FVr_minbound = log10sig[0] * np.ones(len(nmagg))
    FVr_maxbound = log10sig[1] * np.ones(len(nmagg))
    
    # I_bnd_constr: 1 = use bounds as constraints, 0 = no bound constraints
    I_bnd_constr = 1
    
    # I_NP: number of population members
    I_NP = 25
    
    # I_itermax: maximum number of iterations (generations)
    I_itermax = 1000
    
    # F_weight: DE-stepsize F_weight in [0, 2]
    F_weight = 0.6
    
    # F_CR: crossover probability constant in [0, 1]
    F_CR = 0.9
    
    # I_strategy: DE strategy
    # 1 --> DE/rand/1: classical version
    # 2 --> DE/local-to-best/1: balance between robustness and fast convergence
    # 3 --> DE/best/1 with jitter: for small populations and fast convergence
    # 4 --> DE/rand/1 with per-vector-dither: more robust
    # 5 --> DE/rand/1 with per-generation-dither: more robust (F_weight=0.3 is good)
    # 6 --> DE/rand/1 either-or-algorithm: alternates between mutations
    I_strategy = 5
    
    # I_refresh: intermediate output produced after I_refresh iterations
    # No intermediate output if I_refresh < 1
    I_refresh = 10
    
    # I_plotting: 1 = use plotting, 0 = skip plotting
    I_plotting = 0
    
    # Problem dependent constant values
    FVr_bound = np.array([
        0, 1.2, 1.88, 3.119, 6.06879,
        11.25312, 20.93868, 38.99973,
        72.66066687999998, 135.385869312,
        252.2654194687999, 470.0511374131198,
        875.857310322687, 1632.00640736133,
        3040.958067344504, 5666.29295426548,
        10558.14502289265, 19673.25510067688,
        36657.66721873185, 68305.14622427958,
        127274.6837195391
    ])
    
    # Definition of tolerance scheme
    I_lentol = 50
    FVr_x = np.linspace(-1, 1, I_lentol)
    FVr_lim_up = np.ones(I_lentol)
    FVr_lim_lo = -np.ones(I_lentol)
    
    # Tie all important values to a structure
    S_struct = {
        'FVr_bound': FVr_bound,
        'I_lentol': I_lentol,
        'FVr_x': FVr_x,
        'FVr_lim_up': FVr_lim_up,
        'FVr_lim_lo': FVr_lim_lo,
        'I_NP': I_NP,
        'F_weight': F_weight,
        'F_CR': F_CR,
        'I_D': I_D,
        'FVr_minbound': FVr_minbound,
        'FVr_maxbound': FVr_maxbound,
        'I_bnd_constr': I_bnd_constr,
        'I_itermax': I_itermax,
        'F_VTR': F_VTR,
        'I_strategy': I_strategy,
        'I_refresh': I_refresh,
        'I_plotting': I_plotting,
    }
    
    # Pack S struct variables
    S['ibin'] = np.array([ev['ibin'] for ev in S['evspec_DE']])
    S['fmom'] = np.array([ev['qmomest'] for ev in S['evspec_DE']])
    
    # Values just needed for plotting
    if I_plotting == 1:
        FVr_xplot = np.linspace(-1.3, 1.3, 100)
        # Note: step function implementation needed for plotting
        # FVr_lim_lo_plot = FVr_bound[I_D] - (FVr_bound[I_D] + 1) * step(FVr_xplot + 1.2) + (FVr_bound[I_D] + 1) * step(FVr_xplot - 1.2)
        S_struct['FVr_xplot'] = FVr_xplot
        # S_struct['FVr_lim_lo_plot'] = FVr_lim_lo_plot
    
    # Start of optimization
    best_params, S_y, I_nf, _, ECSs = deopt('objfun', S_struct, S)
    
    return best_params, ECSs