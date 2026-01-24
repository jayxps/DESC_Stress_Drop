import numpy as np

from .sachdr import sachdr_from_sac


def fget_sac(filename):
    """
    Read SAC file and return time vector, data, and SAC header struct.

    Returns (t, data, SAChdr) where:
    - t: numpy array from b to e by delta using npts
    - data: numpy array of trace samples
    - SAChdr: dict built by `sachdr_from_sac`
    """
    try:
        from obspy import read as obspy_read
    except Exception:
        raise RuntimeError("ObsPy is required: pip install obspy or conda install obspy")

    st = obspy_read(filename)
    if len(st) == 0:
        raise ValueError(f"No trace in file: {filename}")
    tr = st[0]
    sac = getattr(tr.stats, 'sac', None)
    if sac is None:
        raise ValueError("SAC header missing in trace")

    delta = float(getattr(sac, 'delta', tr.stats.delta))
    b = float(getattr(sac, 'b', 0.0))
    npts = int(tr.stats.npts)
    t = b + np.arange(npts) * delta
    SAChdr = sachdr_from_sac(sac)

    return t, tr.data.astype(np.float64), SAChdr
