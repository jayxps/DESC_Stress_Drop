import os
import numpy as np
from scipy.io import loadmat

from round_half_up import round_half_up


def specprocess_organize(P):
    """
    Organize spectra by year/month/distance with selection and SNR screening.

    Parameters
    ----------
    P : dict
        Parameter dictionary containing:
        - eqinfo: path to eqinfo.mat
        - evdir: base directory for event .mat files
        - qlat1, qlat2, qlon1, qlon2: spatial bounds
        - minmag, maxmag: magnitude bounds
        - randfrac: random fraction for subsampling
        - targchan: two-letter channel prefix to keep (e.g., 'BH')
        - targnet: network filter or 'all'
        - targsta: station filter or 'all'
        - target: 'P', 'S', or 'C'
        - freq1, freq2, freq3: two-element lists defining SNR bands
        - stnmin: list of minimum SNR thresholds for the three bands
        - fmin_tstar, fmax_tstar: frequency window for t* fit
        - minavgstn: minimum average SNR over t* band
        - dtt: distance bin width
    
    Returns
    -------
    spec : list of dict
        Collected spectra and metadata (bin indices remain 1-based to match MATLAB logic).
    flog : ndarray or None
        Frequency log array from the last processed event (if available).
    """

    # Load eqinfo
    eqinfo_mat = loadmat(P['eqinfo'], squeeze_me=True, struct_as_record=False)
    eqinfo = eqinfo_mat['eqinfo']

    # Random fraction
    nrand = np.random.rand() / len(eqinfo.id)

    mask = (
        (eqinfo.qlat >= P['qlat1']) & (eqinfo.qlat <= P['qlat2']) &
        (eqinfo.qlon >= P['qlon1']) & (eqinfo.qlon <= P['qlon2']) &
        (eqinfo.mb >= P['minmag']) & (eqinfo.mb <= P['maxmag']) &
        (nrand <= P['randfrac'])
    )

    n_select = np.where(mask)[0]
    neq = len(n_select)
    print(f"total number of events within limits: {neq} out of {len(eqinfo.id)} events")

    spec = []
    starray = []
    nsta = 0
    nspec = 0
    flog = None

    for iq in range(neq):
        indi = n_select[iq]
        qyr = str(eqinfo.qyr[indi])
        qmon_val = int(eqinfo.qmon[indi])
        qmon = f"{qmon_val:02d}"
        event_file = os.path.join(P['evdir'], qyr, qmon, f"ev{eqinfo.id[indi]}.mat")
        print(event_file)
        if not os.path.exists(event_file):
            continue

        evdata = loadmat(event_file, squeeze_me=True, struct_as_record=False)
        if 'data' not in evdata:
            continue
        data = evdata['data']
        if getattr(data, 'numts', 0) == 0:
            continue

        ustid = np.unique(np.array(data.stid, dtype=str))
        for st_id in ustid:
            trace_all = np.where(np.array(data.stid, dtype=str) == st_id)[0]
            for itr in trace_all:
                chpref = str(data.chan[itr])[:2]
                if chpref != P['targchan']:
                    continue
                if P['targnet'] != 'all' and str(data.net[itr]) != P['targnet']:
                    continue
                if P['targsta'] != 'all' and str(data.sta[itr]) != P['targsta']:
                    continue
                if P['target'] == 'P' and data.pick1[itr] < 0:
                    continue
                if P['target'] == 'S' and data.pick2[itr] < 0:
                    continue

            itsz = 0
            itse = 0
            itsn = 0
            for itr in trace_all:
                chtemp = str(data.chan[itr])
                if len(chtemp) >= 3 and (chtemp[2] == 'Z' or chtemp[2] == '1'):
                    itsz = itr
                if len(chtemp) >= 3 and (chtemp[2] == 'E' or chtemp[2] == '3'):
                    itse = itr
                if len(chtemp) >= 3 and (chtemp[2] == 'N' or chtemp[2] == '2'):
                    itsn = itr

            if P['target'] == 'P' and itsz == 0:
                continue
            if P['target'] == 'S' and (itse == 0 or itsn == 0):
                continue

            if data.pspecres[itsn, 0] == 0 or data.pspecres[itse, 0] == 0 or data.pspecres[itsz, 0] == 0:
                continue

            if P['target'] == 'P':
                spectemp = np.sqrt(data.pspec[itsn, :] ** 2 + data.pspec[itse, :] ** 2 + data.pspec[itsz, :] ** 2)
                noisetemp = np.sqrt(data.pnspec[itsn, :] ** 2 + data.pnspec[itse, :] ** 2 + data.pnspec[itsz, :] ** 2)
                specres = np.sqrt(data.pspecres[itsn, :] ** 2 + data.pspecres[itse, :] ** 2 + data.pspecres[itsz, :] ** 2)
                freq = data.pfreq[itsz, :]
                ttime = np.mean(data.pick1[trace_all])
            elif P['target'] == 'S':
                spectemp = np.sqrt(data.sspec[itse, :] ** 2 + data.sspec[itsn, :] ** 2)
                noisetemp = np.sqrt(data.pnspec[itsn, :] ** 2 + data.pnspec[itse, :] ** 2)
                specres = np.sqrt(data.sspecres[itse, :] ** 2 + data.sspecres[itsn, :] ** 2)
                freq = data.sfreq[itse, :]
                ttime = np.mean(data.spred[trace_all])
            else:  # Coda (Only a placeholder. Don't use, not ready yet)
                spectemp = data.cspec[itsz, :]
                noisetemp = data.cnspec[itsz, :]
                specres = data.cspecres[itsz, :]
                freq = data.cfreq[itsz, :]
                ttime = np.mean(data.spred[trace_all])

            nbf1 = np.where((freq >= P['freq1'][0]) & (freq <= P['freq1'][1]))[0]
            stn1 = np.min(np.log10(spectemp[nbf1]) - np.log10(noisetemp[nbf1])) if len(nbf1) else -np.inf
            nbf2 = np.where((freq >= P['freq2'][0]) & (freq <= P['freq2'][1]))[0]
            stn2 = np.min(np.log10(spectemp[nbf2]) - np.log10(noisetemp[nbf2])) if len(nbf2) else -np.inf
            nbf3 = np.where((freq >= P['freq3'][0]) & (freq <= P['freq3'][1]))[0]
            stn3 = np.min(np.log10(spectemp[nbf3]) - np.log10(noisetemp[nbf3])) if len(nbf3) else -np.inf

            nf_tstar = np.where((freq >= P['fmin_tstar']) & (freq <= P['fmax_tstar']))[0]
            if len(nf_tstar) == 0:
                continue
            sumnoise = np.sum(np.log10(spectemp[nf_tstar]) - np.log10(noisetemp[nf_tstar]))
            avgstn = sumnoise / len(nf_tstar)

            if (stn1 < P['stnmin'][0] or stn2 < P['stnmin'][1] or stn3 < P['stnmin'][2] or avgstn < P['minavgstn']):
                continue

            specx = np.log10(specres)
            if not np.isfinite(specx).all():
                continue

            nspec += 1
            entry = {}
            entry['spec'] = specx
            entry['inq'] = n_select[iq] + 1
            entry['qlat'] = data.qlat
            entry['qlon'] = data.qlon
            entry['qdep'] = data.qdep
            entry['qmag'] = getattr(data, 'qmag1', getattr(data, 'mb', np.nan))
            entry['qtime'] = data.qtime
            entry['qid'] = getattr(data, 'cuspid', getattr(data, 'id', eqinfo.id[indi]))

            stname = str(data.sta[itsz])
            stid = f"{stname:<6}{str(data.chan[itsz])[:2]:<2}"
            if stid in starray:
                entry['ins'] = starray.index(stid) + 1
            else:
                starray.append(stid)
                nsta += 1
                entry['ins'] = nsta
            entry['stid'] = stid
            entry['stlat'] = data.stlat[itsz]
            entry['stlon'] = data.stlon[itsz]
            entry['stelev'] = data.stelev[itsz]
            entry['dist'] = data.dist[itsz]

            ix = int(round_half_up(ttime / P['dtt']))
            if ix <= 1:
                ix = 1
            entry['inx'] = ix
            entry['ttime'] = ttime

            qdp = data.qdep
            id_bin = int(round_half_up(qdp + 0.6))
            if id_bin <= 1:
                id_bin = 1
            if id_bin >= 30:
                id_bin = 30
            entry['ind'] = id_bin

            qmag_val = entry['qmag']
            im = int(round_half_up(qmag_val * 5.0 + 0.6)) if np.isfinite(qmag_val) else 1
            im = max(1, min(im, 40))
            entry['inm'] = im

            flog = getattr(data, 'flog', None)

            x = freq[nf_tstar]
            y = np.log10(spectemp[nf_tstar])
            xmean = np.mean(x)
            ymean = np.mean(y)
            sumxx = np.sum((x - xmean) ** 2)
            sumxy = np.sum((x - xmean) * (y - ymean))
            b = sumxy / sumxx if sumxx != 0 else np.nan
            entry['amp'] = ymean - b * xmean if np.isfinite(b) else 0
            entry['tstar'] = -b / 1.364 if np.isfinite(b) else 0

            spec.append(entry)

    return spec, flog
