import numpy as np

from .fget_sac import fget_sac


def fget_sac_cut(filename, tstart=None, tend=None):
    """
    Read SAC and cut the data by time window [tstart, tend].

    Returns (t_cut, data_cut, SAChdr) and updates SAChdr['times']['npts'].
    """
    t, data, SAChdr = fget_sac(filename)
    if tstart is None:
        tstart = -1e6
    if tend is None:
        tend = 1e6

    tmin = float(np.min(t))
    tmax = float(np.max(t))
    if tstart < tmin:
        tstart = tmin
    if tend > tmax:
        tend = tmax

    mask = (t >= tstart) & (t <= tend)
    t_cut = t[mask]
    data_cut = data[mask]

    # Update npts in header struct
    if 'times' in SAChdr:
        SAChdr['times']['npts'] = int(len(data_cut))

    return t_cut, data_cut, SAChdr
