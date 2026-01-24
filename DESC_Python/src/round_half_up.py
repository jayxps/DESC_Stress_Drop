"""
Rounding utility for consistent cross-platform behavior.

Python and NumPy use "round half to even" (banker's rounding),
while MATLAB uses "round half away from zero". This module provides
a function that matches MATLAB's behavior for numerical compatibility.
"""
import numpy as np


def round_half_up(x):
    """
    Round half away from zero (round .5 up to next integer).
    
    This matches MATLAB's round() behavior, where values ending in .5
    are always rounded away from zero (e.g., 2.5 -> 3, -2.5 -> -3).
    
    Parameters
    ----------
    x : float or array-like
        Value(s) to round
    
    Returns
    -------
    float or ndarray
        Rounded value(s)
    
    Examples
    --------
    >>> round_half_up(2.5)
    3.0
    >>> round_half_up(3.5)
    4.0
    >>> round_half_up(-2.5)
    -3.0
    
    Notes
    -----
    Python's built-in round() and np.round() use "round half to even":
    - round(2.5) -> 2 (rounds to nearest even)
    - round(3.5) -> 4 (rounds to nearest even)
    
    This function always rounds away from zero:
    - round_half_up(2.5) -> 3
    - round_half_up(3.5) -> 4
    """
    return np.sign(x) * np.floor(np.abs(x) + 0.5)
