def new_sac_header(n, dt, b):
    """
    Form a new SAC header-like dict.

    Returns a minimal header dict with keys similar to newSacHeader.m:
    - npts, delta, b, e, o, iztype, iftype, leven, version
    """
    e = b + (n - 1) * dt
    hd = {
        'npts': int(n),
        'delta': float(dt),
        'b': float(b),
        'o': 0.0,
        'e': float(e),
        'iztype': 11,   # origin time
        'iftype': 1,    # time series
        'leven': 1,     # evenly spaced
        'version': 6,
    }
    return hd
