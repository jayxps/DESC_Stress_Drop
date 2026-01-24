import numpy as np


def rd_sac_head(sacFile):
    """
    Read header of SAC format data using ObsPy and return a dict
    with selected fields matching rdSacHead.m.
    """
    try:
        from obspy import read as obspy_read
        from obspy.geodetics import gps2dist_azimuth
    except Exception:
        raise RuntimeError("ObsPy is required: pip install obspy or conda install obspy")

    st = obspy_read(sacFile)
    if len(st) == 0:
        raise ValueError(f"No trace in file: {sacFile}")
    tr = st[0]
    sac = getattr(tr.stats, 'sac', None)
    if sac is None:
        raise ValueError("SAC header missing in trace")

    hd = {
        'delta': float(getattr(sac, 'delta', tr.stats.delta)),
        'b': float(getattr(sac, 'b', 0.0)),
        'stla': getattr(sac, 'stla', np.nan),
        'stlo': getattr(sac, 'stlo', np.nan),
        'evla': getattr(sac, 'evla', np.nan),
        'evlo': getattr(sac, 'evlo', np.nan),
        'dist': getattr(sac, 'dist', np.nan),
        'az': getattr(sac, 'az', np.nan),
        'baz': getattr(sac, 'baz', np.nan),
        'gcarc': getattr(sac, 'gcarc', np.nan),
        'npts': int(tr.stats.npts),
    }

    # If event/station coordinates exist, compute distance/azimuth
    stla, stlo, evla, evlo = hd['stla'], hd['stlo'], hd['evla'], hd['evlo']
    if np.isfinite(stla) and np.isfinite(stlo) and np.isfinite(evla) and np.isfinite(evlo):
        try:
            dist_m, az, baz = gps2dist_azimuth(evla, evlo, stla, stlo)
            hd['dist'] = dist_m / 1000.0
            hd['az'] = az
            hd['baz'] = baz
        except Exception:
            pass

    return hd
