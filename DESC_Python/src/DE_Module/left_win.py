import numpy as np
"""
Function to compare two objective function value structures.

Author: Rainer Storn
Description: left_win(S_x, S_y) takes structures S_x and S_y as arguments.
The function returns 1 if the left structure S_x wins. If the right structure S_y wins,
the result is 0.

Parameters:
    S_x.I_nc  (I)    Number of constraints (must be the same for x and y).
    S_x.I_no  (I)    Number of objectives (must be the same for x and y).
    S_x.FVr_ca (I)   Constraint array containing constraint violation values.
                     If value is 0, constraint is met. If > 0, still violated.
    S_x.FVr_oa (I)   Objective array containing cost values to be minimized.

Return value:
    I_z  (O)  If S_x wins over S_y then I_z=1, else I_z=0.
"""


def left_win(S_x, S_y):
    """
    Compare two objective function value structures.
    
    Returns 1 if S_x wins (is better), 0 if S_y wins.
    Minimization problem: lower values are better.
    
    Parameters
    ----------
    S_x : dict
        Structure with keys:
        - 'I_no': Number of objectives
        - 'FVr_oa': Array of objective values
        - 'I_nc': Number of constraints (optional)
        - 'FVr_ca': Array of constraint values (optional)
    
    S_y : dict
        Same structure as S_x
    
    Returns
    -------
    I_z : int
        1 if S_x wins, 0 if S_y wins
    """
    
    I_z = 1  # start with I_z = 1 (assume S_x wins)
    
    # Deal with objectives
    # For minimization: S_x wins if ALL objectives of S_x are <= corresponding objectives of S_y
    # S_x loses (I_z = 0) if ANY objective of S_x is > corresponding objective of S_y
    
    if 'I_no' in S_x and S_x['I_no'] > 0:
        for k in range(S_x['I_no']):
            # If any objective of S_x is greater (worse) than S_y, S_x loses
            if S_x['FVr_oa'][k] < S_y['FVr_oa'][k]:
                I_z = 0
                break
    
    return I_z
