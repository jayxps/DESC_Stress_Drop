import os
import re
import sys
import numpy as np
from scipy.io import savemat
from datetime import datetime, timedelta

sys.path.append(os.path.dirname(os.path.dirname(__file__)))
from date2epoch import date2epoch
from round_half_up import round_half_up


def sac2mat(sacdir_in: str, evdir: str, evid: int, eqinfo: dict, sta: dict):
    """
    Convert SAC files for a given event into a MATLAB-like `.mat` data file.

    Parameters
    ----------
    sacdir_in : str
        Root SAC directory containing per-event subfolders.
    evdir : str
        Output directory for event `.mat` file.
    evid : int
        Event ID.
    eqinfo : dict-like
        Catalog structure with fields like `id`, `qlat`, `qlon`, `qtime`, `mb`, etc. Arrays indexed by event index.
    sta : dict-like
        Station table with arrays: `stid`, `slat`, `slon`, `selev`.

    Returns
    -------
    data : dict
        Data structure saved to `<evdir>/ev<evid>.mat`.
    """
    try:
        from obspy import read as obspy_read
        from obspy.geodetics import gps2dist_azimuth
    except Exception:
        obspy_read = None
        gps2dist_azimuth = None
    if obspy_read is None:
        raise RuntimeError("ObsPy is required: pip install obspy")

    sacdir_in = sacdir_in if sacdir_in.endswith('/') else sacdir_in + '/'
    evdir = evdir if evdir.endswith('/') else evdir + '/'

    idstr = str(evid)
    sacdir = os.path.join(sacdir_in, idstr)
    evfile = os.path.join(evdir, f"ev{evid}.mat")

    if not os.path.isdir(sacdir):
        return {}
    os.makedirs(evdir, exist_ok=True)

    # Event info
    ids = np.asarray(eqinfo['id']).flatten()
    ie_matches = np.where(ids == evid)[0]
    if ie_matches.size == 0:
        print(f"target event not in the catalog {evid}")
        return {}
    ie = ie_matches[0]

    data = {}
    # Copy all fields (scalar/array) at index ie
    for key, arr in eqinfo.items():
        try:
            data[key] = np.asarray(arr).flatten()[ie]
        except Exception:
            data[key] = arr
    data['cuspid'] = data.get('id', evid)
    data['qmag1'] = np.round(data.get('mb', np.nan), 5)
    data['id0'] = idstr

    # Gather SAC files
    saclist = [f for f in os.listdir(sacdir) if not f.startswith('.')]
    saclist.sort()
    numts = 0

    # Preallocate lists
    nets, stas, stids, chans = [], [], [], []
    stime_list = []
    syt, smon, sday, shr, smn, ssc = [], [], [], [], [], []
    sdoy = []
    samprate, dt_list = [], []
    pick1, pick2 = [], []
    comp = []
    dist_list, sazi_list, pazi_list = [], [], []
    seis_list, aseis_list, npts_list = [], [], []
    stlon, stlat, stel = [], [], []

    # Load station lookup arrays
    sta_stid = np.asarray(sta.get('stid', []), dtype=str)
    sta_slat = np.asarray(sta.get('slat', []))
    sta_slon = np.asarray(sta.get('slon', []))
    sta_selev = np.asarray(sta.get('selev', []))

    qlat = float(data.get('qlat', np.nan))
    qlon = float(data.get('qlon', np.nan))

    for fname in saclist:
        fpath = os.path.join(sacdir, fname)
        try:
            st = obspy_read(fpath)
        except Exception:
            continue
        if len(st) == 0:
            continue
        tr = st[0]
        sac = getattr(tr.stats, 'sac', None)
        if sac is None:
            continue

        # Event origin time components from SAC header (reference start)
        nzyear = int(getattr(sac, 'nzyear', 1970))
        nzjday = int(getattr(sac, 'nzjday', 1))
        nzhour = int(getattr(sac, 'nzhour', 0))
        nzmin = int(getattr(sac, 'nzmin', 0))
        nzsec = int(getattr(sac, 'nzsec', 0))
        nzmsec = int(getattr(sac, 'nzmsec', 0))
        b = float(getattr(sac, 'b', 0.0))
        # Convert jday to month/day
        origin = datetime(nzyear, 1, 1) + timedelta(days=nzjday - 1)
        syr, smon0, sday0 = origin.year, origin.month, origin.day
        shr0, smn0, ssc0 = nzhour, nzmin, nzsec + round(nzmsec / 1000.0, 3) + b
        # Use shared utility; expects sequence of vectors, returns array
        stime_epoch = float(date2epoch([(syr, smon0, sday0, shr0, smn0, ssc0)]))

        # Station and channel names
        staname = str(getattr(sac, 'kstnm', ''))
        # strip non-alphanum tail like MATLAB
        m = re.match(r"([A-Za-z0-9]+)", staname)
        staname = m.group(1) if m else staname
        netname = str(getattr(sac, 'knetwk', ''))[:2]
        chname = str(getattr(sac, 'kcmpnm', ''))[:3]
        stid = netname + staname
        npts = int(tr.stats.npts)

        numts += 1
        nets.append(netname)
        stas.append(staname)
        stids.append(stid)
        stime_list.append(stime_epoch)
        syt.append(syr)
        smon.append(smon0)
        sday.append(sday0)
        shr.append(shr0)
        smn.append(smn0)
        ssc.append(ssc0)
        sdoy.append(nzjday)
        sr = round_half_up(1.0 / float(getattr(sac, 'delta', tr.stats.delta)))
        samprate.append(sr)
        dt_list.append(1.0 / sr)
        pick1.append(float(getattr(sac, 'a', -12345.0)))
        pick2.append(float(getattr(sac, 't0', -12345.0)))
        chans.append(chname)

        # Station coordinates: lookup table or SAC header
        # Find station in provided list
        nfsta = np.where(sta_stid == stid)[0]
        if nfsta.size == 0:
            stlat.append(float(getattr(sac, 'stla', np.nan)))
            stlon.append(float(getattr(sac, 'stlo', np.nan)))
            stel.append(float(getattr(sac, 'stel', np.nan)) / 1000.0 if hasattr(sac, 'stel') else np.nan)
        else:
            idx = int(nfsta[0])
            stlat.append(float(sta_slat[idx]))
            stlon.append(float(sta_slon[idx]))
            stel.append(float(sta_selev[idx]))

        # Compute distance and azimuth
        if gps2dist_azimuth is not None and np.isfinite(qlat) and np.isfinite(qlon) and np.isfinite(stlat[-1]) and np.isfinite(stlon[-1]):
            dist_m, azi, back_azi = gps2dist_azimuth(qlat, qlon, stlat[-1], stlon[-1])
            dist_km = dist_m / 1000.0
        else:
            R = 6371.0
            phi1 = np.radians(qlat)
            phi2 = np.radians(stlat)
            dphi = np.radians(stlat - qlat)
            dlambda = np.radians(stlon - qlon)
            a_h = np.sin(dphi / 2) ** 2 + np.cos(phi1) * np.cos(phi2) * np.sin(dlambda / 2) ** 2
            dist_km = 2 * R * np.arcsin(np.sqrt(a_h))
            azi = 0.0
            back_azi = 180.0

        # Component identification
        comp_val = 0
        if len(chname) >= 3:
            if chname[2] in ('E', 'X', '2'):
                comp_val = 3
            elif chname[2] in ('N', 'Y', '1'):
                comp_val = 2
            elif chname[2] == 'Z':
                comp_val = 1
        if netname.upper() in ('BP', 'SF') and len(chname) >= 3 and chname[2].isdigit():
            comp_val = int(chname[2])

        # Displacement vs acceleration channel handling
        acc_flag = 1 if len(chname) >= 2 and chname[1] in ('N', 'Y') else 0
        x = tr.data.astype(np.float64)
        traces = x - np.mean(x)
        if acc_flag:
            tracea = traces
            # integrate acceleration to velocity then to displacement
            traces = np.cumsum(tracea) / sr
        else:
            tracea = np.zeros_like(traces)
            tracea[:-1] = np.diff(traces) * sr

        dist_list.append(dist_km)
        sazi_list.append(azi)
        pazi_list.append(azi + 180.0)

        seis_list.append(traces)
        aseis_list.append(tracea)
        npts_list.append(npts)
        comp.append(comp_val)

    # Build data structure
    data['numts'] = numts
    data['net'] = np.array(nets, dtype=object)
    data['stime'] = np.array(np.round(stime_list, 3), dtype=float)
    data['syr'] = np.array(syt, dtype=int)
    data['smon'] = np.array(smon, dtype=int)
    data['sdy'] = np.array(sday, dtype=int)
    data['shr'] = np.array(shr, dtype=int)
    data['smn'] = np.array(smn, dtype=int)
    data['ssc'] = np.array(np.round(ssc, 3), dtype=float)
    data['sdoy'] = np.array(sdoy, dtype=int)
    data['samprate'] = np.array(samprate, dtype=int)
    data['dt'] = np.array(np.round(dt_list, 5), dtype=float)
    data['pick1'] = np.array(np.round(pick1, 5), dtype=float)
    data['pick2'] = np.array(np.round(pick2, 5), dtype=float)
    data['sta'] = np.array(stas, dtype=object)
    data['stid'] = np.array(stids, dtype=object)
    data['chan'] = np.array(chans, dtype=object)
    data['comp'] = np.array(comp, dtype=int)
    data['dist'] = np.array(np.round(dist_list, 5), dtype=float)
    data['sazi'] = np.array(np.round(sazi_list, 5), dtype=float)
    data['pazi'] = np.array(np.round(pazi_list, 5), dtype=float)
    data['stlat'] = np.array(np.round(stlat, 6), dtype=float)
    data['stlon'] = np.array(np.round(stlon, 6), dtype=float)
    data['stelev'] = np.array(np.round(stel, 5), dtype=float)

    # Build 2D arrays (numts × max_npts) for waveforms, padding with NaN
    max_npts = max(npts_list) if npts_list else 0
    data['npts'] = np.array(npts_list, dtype=int)
    data['seis'] = np.full((numts, max_npts), np.nan, dtype=float)
    data['aseis'] = np.full((numts, max_npts), np.nan, dtype=float)
    for i in range(numts):
        n = len(seis_list[i])
        data['seis'][i, :n] = seis_list[i]
        data['aseis'][i, :n] = aseis_list[i]

    # Save
    savemat(evfile, {'data': data}, do_compression=True)
    return data
