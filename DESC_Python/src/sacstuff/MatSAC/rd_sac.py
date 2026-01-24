import numpy as np

from .rd_sac_head import rd_sac_head


def rd_sac(sacFile):
    """
    Read SAC format data using ObsPy and return (data, hd).
    hd is a dict of header fields similar to rdSacHead.m.
    """
    try:
        from obspy import read as obspy_read
    except Exception:
        raise RuntimeError("ObsPy is required: pip install obspy or conda install obspy")

    st = obspy_read(sacFile)
    if len(st) == 0:
        raise ValueError(f"No trace in file: {sacFile}")
    tr = st[0]
    data = tr.data.astype(np.float32)
    hd = rd_sac_head(sacFile)
    return data, hd
