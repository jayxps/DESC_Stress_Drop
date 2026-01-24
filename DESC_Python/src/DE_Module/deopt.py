"""
Differential Evolution Optimization

Python translation of deopt.m - Minimization using differential evolution algorithm.

Authors: Rainer Storn, Ken Price, Arnold Neumaier, Jim Van Zandt (original MATLAB)
         Python translation for DESC project

Description:
    Minimization of a user-supplied function with respect to x(1:I_D),
    using the differential evolution (DE) algorithm.
    
    DE works best if [FVr_minbound,FVr_maxbound] covers the region where the
    global minimum is expected. DE is also somewhat sensitive to
    the choice of the stepsize F_weight. A good initial guess is to
    choose F_weight from interval [0.5, 1], e.g. 0.8. F_CR, the crossover
    probability constant from interval [0, 1] helps to maintain
    the diversity of the population but should be close to 1 for most 
    practical cases. Only separable problems do better with CR close to 0.
    If the parameters are correlated, high values of F_CR work better.
    The reverse is true for no correlation.
    
    The number of population members I_NP is also not very critical. A
    good initial guess is 10*I_D. Depending on the difficulty of the
    problem I_NP can be lower than 10*I_D or must be higher than 10*I_D
    to achieve convergence.

Parameters:
    fname: str or callable
        Name of objective function (if str, imports from objfun module) or callable
    S_struct: dict
        Configuration parameters containing:
        - F_VTR: Value to reach (stopping criterion)
        - FVr_minbound: Lower bounds of initial population
        - FVr_maxbound: Upper bounds of initial population
        - I_D: Number of parameters
        - I_NP: Number of population members
        - I_itermax: Maximum iterations
        - F_weight: DE stepsize [0, 2]
        - F_CR: Crossover probability [0, 1]
        - I_strategy: DE strategy (1-6)
        - I_refresh: Intermediate output frequency
        - I_bnd_constr: Use boundary constraints (1) or not (0)
        - I_plotting: Enable plotting (1) or not (0)
    S: dict
        Problem-specific data passed to objective function

Returns:
    FVr_bestmem: Best parameter vector found
    S_bestval: Best objective value (dict with I_nc, FVr_ca, I_no, FVr_oa)
    I_nfeval: Number of function evaluations
    vp_model: History of best parameters over iterations
    power: Final power/objective info from objective function
"""

import numpy as np
from .objfun import objfun
from .left_win import left_win


def deopt(fname, S_struct, S):
    """
    Differential Evolution Optimization main function.
    
    Args:
        fname: str or callable - objective function name or callable
               (if string 'objfun', uses imported objfun function)
        S_struct: dict of DE configuration parameters
        S: dict of problem-specific data
    
    Returns:
        FVr_bestmem: best parameter vector
        S_bestval: best objective value structure
        I_nfeval: number of function evaluations
        vp_model: parameter history
        power: final power/objective info
    """
    
    # Determine which objective function to use
    # If fname is a string 'objfun', use the imported objfun function
    # Otherwise, if it's callable, use it directly
    if isinstance(fname, str) and fname == 'objfun':
        obj_func = objfun
    elif callable(fname):
        obj_func = fname
    else:
        raise ValueError(f"fname must be 'objfun' string or a callable, got {type(fname)}")
    
    # Extract parameters from S_struct
    I_NP = S_struct['I_NP']
    F_weight = S_struct['F_weight']
    F_CR = S_struct['F_CR']
    I_D = S_struct['I_D']
    FVr_minbound = np.array(S_struct['FVr_minbound'])
    FVr_maxbound = np.array(S_struct['FVr_maxbound'])
    I_bnd_constr = S_struct['I_bnd_constr']
    I_itermax = S_struct['I_itermax']
    F_VTR = S_struct['F_VTR']
    I_strategy = S_struct['I_strategy']
    I_refresh = S_struct['I_refresh']
    I_plotting = S_struct.get('I_plotting', 0)
    
    # Check input variables
    if I_NP < 5:
        I_NP = 5
        print(' I_NP increased to minimal value 5')
    
    if F_CR < 0 or F_CR > 1:
        F_CR = 0.5
        print('F_CR should be from interval [0,1]; set to default value 0.5')
    
    if I_itermax <= 0:
        I_itermax = 200
        print('I_itermax should be > 0; set to default value 200')
    
    I_refresh = int(I_refresh)
    
    # Initialize population
    FM_pop = np.zeros((I_NP, I_D))
    
    # Initialize with random values between min and max bounds
    for k in range(I_NP):
        FM_pop[k, :] = FVr_minbound + np.random.rand(I_D) * (FVr_maxbound - FVr_minbound)
    
    FM_popold = np.zeros_like(FM_pop)
    FVr_bestmem = np.zeros(I_D)
    FVr_bestmemit = np.zeros(I_D)
    I_nfeval = 0
    
    # Evaluate the best member after initialization
    I_best_index = 0
    S_val = [None] * I_NP
    S_val[0], rms_best = obj_func(FM_pop[0, :], S_struct, S)
    
    S_bestval = S_val[0]
    I_nfeval = I_nfeval + 1
    
    for k in range(1, I_NP):
        S_val[k], rms_temp = obj_func(FM_pop[k, :], S_struct, S)
        I_nfeval = I_nfeval + 1
        if left_win(S_val[k], S_bestval) == 1:
            I_best_index = k
            S_bestval = S_val[k]
            rms_best = rms_temp
    
    FVr_bestmemit = FM_pop[I_best_index, :].copy()
    S_bestvalit = S_bestval
    FVr_bestmem = FVr_bestmemit.copy()
    
    # Initialize arrays for DE iterations
    FM_pm1 = np.zeros((I_NP, I_D))
    FM_pm2 = np.zeros((I_NP, I_D))
    FM_pm3 = np.zeros((I_NP, I_D))
    FM_pm4 = np.zeros((I_NP, I_D))
    FM_pm5 = np.zeros((I_NP, I_D))
    FM_bm = np.zeros((I_NP, I_D))
    FM_ui = np.zeros((I_NP, I_D))
    FM_mui = np.zeros((I_NP, I_D))
    FM_mpo = np.zeros((I_NP, I_D))
    FVr_rot = np.arange(I_NP)
    FVr_rotd = np.arange(I_D)
    
    # Initialize model history
    vp_model = np.zeros((I_itermax, I_D))
    
    I_iter = 1
    
    # Main DE loop
    while (I_iter < I_itermax) and (S_bestval['FVr_oa'][0] > F_VTR):
        FM_popold = FM_pop.copy()
        S_struct['FM_pop'] = FM_pop
        S_struct['FVr_bestmem'] = FVr_bestmem
        
        # NOTE: MATLAB's randperm(n) returns values 1 to n
        # NumPy's permutation(n) returns values 0 to n-1
        # Add 1 to match MATLAB behavior
        FVr_ind = np.random.permutation(4) + 1
        
        # Same for population shuffling - add 1 to match MATLAB 1-indexing
        FVr_a1 = np.random.permutation(I_NP) + 1
        FVr_rt = (FVr_rot + FVr_ind[0]) % I_NP
        FVr_a2 = FVr_a1[FVr_rt]
        FVr_rt = (FVr_rot + FVr_ind[1]) % I_NP
        FVr_a3 = FVr_a2[FVr_rt]
        FVr_rt = (FVr_rot + FVr_ind[2]) % I_NP
        FVr_a4 = FVr_a3[FVr_rt]
        FVr_rt = (FVr_rot + FVr_ind[3]) % I_NP
        FVr_a5 = FVr_a4[FVr_rt]
        
        # Subtract 1 to convert back to 0-based indexing for NumPy
        FM_pm1 = FM_popold[FVr_a1 - 1, :]
        FM_pm2 = FM_popold[FVr_a2 - 1, :]
        FM_pm3 = FM_popold[FVr_a3 - 1, :]
        FM_pm4 = FM_popold[FVr_a4 - 1, :]
        FM_pm5 = FM_popold[FVr_a5 - 1, :]
        
        # Population filled with best member
        for k in range(I_NP):
            FM_bm[k, :] = FVr_bestmemit
        
        FM_mui = np.random.rand(I_NP, I_D) < F_CR
        FM_mpo = FM_mui < 0.5
        
        # Apply DE strategy
        if I_strategy == 1:  # DE/rand/1
            FM_ui = FM_pm3 + F_weight * (FM_pm1 - FM_pm2)
            FM_ui = FM_popold * FM_mpo + FM_ui * FM_mui
            FM_origin = FM_pm3
        elif I_strategy == 2:  # DE/local-to-best/1
            FM_ui = FM_popold + F_weight * (FM_bm - FM_popold) + F_weight * (FM_pm1 - FM_pm2)
            FM_ui = FM_popold * FM_mpo + FM_ui * FM_mui
            FM_origin = FM_popold
        elif I_strategy == 3:  # DE/best/1 with jitter
            FM_ui = FM_bm + (FM_pm1 - FM_pm2) * ((1 - 0.9999) * np.random.rand(I_NP, I_D) + F_weight)
            FM_ui = FM_popold * FM_mpo + FM_ui * FM_mui
            FM_origin = FM_bm
        elif I_strategy == 4:  # DE/rand/1 with per-vector-dither
            f1 = ((1 - F_weight) * np.random.rand(I_NP, 1) + F_weight)
            FM_pm5 = np.tile(f1, (1, I_D))
            FM_ui = FM_pm3 + (FM_pm1 - FM_pm2) * FM_pm5
            FM_origin = FM_pm3
            FM_ui = FM_popold * FM_mpo + FM_ui * FM_mui
        elif I_strategy == 5:  # DE/rand/1 with per-generation-dither
            f1 = ((1 - F_weight) * np.random.rand() + F_weight)
            FM_ui = FM_pm3 + (FM_pm1 - FM_pm2) * f1
            FM_origin = FM_pm3
            FM_ui = FM_popold * FM_mpo + FM_ui * FM_mui
        else:  # I_strategy == 6: either-or-algorithm
            if np.random.rand() < 0.5:
                FM_ui = FM_pm3 + F_weight * (FM_pm1 - FM_pm2)
                FM_origin = FM_pm3
            else:
                FM_ui = FM_pm3 + 0.5 * (F_weight + 1.0) * (FM_pm1 + FM_pm2 - 2 * FM_pm3)
            FM_ui = FM_popold * FM_mpo + FM_ui * FM_mui
        
        # Select which vectors are allowed to enter the new population
        for k in range(I_NP):
            # Boundary constraints via bounce back
            if I_bnd_constr == 1:
                for j in range(I_D):
                    if FM_ui[k, j] > FVr_maxbound[j]:
                        FM_ui[k, j] = FVr_maxbound[j] + np.random.rand() * (FM_origin[k, j] - FVr_maxbound[j])
                    if FM_ui[k, j] < FVr_minbound[j]:
                        FM_ui[k, j] = FVr_minbound[j] + np.random.rand() * (FM_origin[k, j] - FVr_minbound[j])
            
            # Check cost of competitor
            S_tempval, rms_temp = obj_func(FM_ui[k, :], S_struct, S)
            I_nfeval = I_nfeval + 1
            
            if left_win(S_tempval, S_val[k]) == 1:
                FM_pop[k, :] = FM_ui[k, :]
                S_val[k] = S_tempval
                
                # Update S_bestval only in case of success
                if left_win(S_tempval, S_bestval) == 1:
                    S_bestval = S_tempval
                    FVr_bestmem = FM_ui[k, :].copy()
                    rms_best = rms_temp
        
        power = rms_best
        vp_model[I_iter - 1, :] = FVr_bestmem
        FVr_bestmemit = FVr_bestmem.copy()
        
        # Output section
        if I_refresh > 0:
            if (I_iter % I_refresh == 0) or (I_iter == 1):
                print(f'Iteration: {I_iter},  Best: {S_bestval["FVr_oa"][0]:.12f},  '
                      f'F_weight: {F_weight},  F_CR: {F_CR},  I_NP: {I_NP}')
                for n in range(I_D):
                    print(f'best({n}) = {FVr_bestmem[n]:.6e}')
                
                if I_plotting == 1:
                    # PlotIt would be called here if implemented
                    pass
        
        I_iter = I_iter + 1
    
    # Trim vp_model to actual iterations
    vp_model = vp_model[:I_iter - 1, :]
    
    return FVr_bestmem, S_bestval, I_nfeval, vp_model, power
