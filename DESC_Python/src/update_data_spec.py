import os
import numpy as np
from scipy.io import loadmat, savemat
import multitaper
from round_half_up import round_half_up


def _matstruct_to_dict(obj):
    """Recursively convert MATLAB structs to Python dicts."""
    if isinstance(obj, np.ndarray) and obj.dtype.names:
        return {name: _matstruct_to_dict(obj[name]) for name in obj.dtype.names}
    if hasattr(obj, "_fieldnames"):
        return {fn: _matstruct_to_dict(getattr(obj, fn)) for fn in obj._fieldnames}
    if isinstance(obj, np.ndarray):
        return np.array(obj)
    return obj


def _multitaper_psd(signal, fs, nfft):
    mt = multitaper.MTSpec(signal, nw=2.5, kspec=3, dt=1/fs, nfft=nfft)
    spec = mt.spec.squeeze()[:nfft // 2 + 1]
    freq = mt.freq.squeeze()[:nfft // 2 + 1]
    if np.isnan(spec[0]):
        spec = np.zeros(spec.shape)
    return np.sqrt(spec), freq


def update_data_spec(evfile, disp_flag=1):
    """
    Compute spectra
    """
    mat = loadmat(evfile, squeeze_me=True, struct_as_record=False)
    if "data" not in mat:
        return
    data = _matstruct_to_dict(mat["data"])
    if data.get("numts", 0) == 0:
        return

    if disp_flag:
        print(f"Working on {evfile}")

    Nfft = 1024
    Nfftc = 4096

    numts = int(data["numts"])
    samprate_arr = np.asarray(data["samprate"]).flatten()

    def zmat(shape):
        return np.zeros(shape)

    data["pseis"] = np.empty(numts, dtype=object)
    data["pnoise"] = np.empty(numts, dtype=object)
    data["sseis"] = np.empty(numts, dtype=object)
    data["snoise"] = np.empty(numts, dtype=object)
    data["cseis"] = np.empty(numts, dtype=object)
    data["cnoise"] = np.empty(numts, dtype=object)
    data["pspec"] = zmat((numts, Nfft // 2 + 1))
    data["sspec"] = zmat((numts, Nfft // 2 + 1))
    data["cspec"] = zmat((numts, Nfftc // 2 + 1))
    data["pnspec"] = zmat((numts, Nfft // 2 + 1))
    data["snspec"] = zmat((numts, Nfft // 2 + 1))
    data["cnspec"] = zmat((numts, Nfftc // 2 + 1))
    data["pfreq"] = zmat((numts, Nfft // 2 + 1))
    data["sfreq"] = zmat((numts, Nfft // 2 + 1))
    data["cfreq"] = zmat((numts, Nfftc // 2 + 1))

    data["np1"] = zmat(numts)
    data["tp1"] = zmat(numts)
    data["ns1"] = zmat(numts)
    data["ts1"] = zmat(numts)
    data["codawin"] = zmat(numts)
    data["tswin"] = zmat(numts)
    data["tpwin"] = zmat(numts)
    data["pwin_mode"] = np.empty(numts, dtype=object)
    data["pnwin"] = zmat(numts)
    data["snwin"] = zmat(numts)
    data["cnwin"] = zmat(numts)
    data["pstn"] = zmat(numts)
    data["sstn"] = zmat(numts)

    samprate = samprate_arr[0]
    N_freq_samples = 1000 if samprate > 100 else 500
    data["flog"] = np.logspace(np.log10(0.5), np.log10(samprate * 0.45), N_freq_samples)

    # local-time adjustment
    data["stime"] += 8 * 3600

    for ns in range(numts):
        ppick = data["pick1"][ns]
        spick = data["pick2"][ns]
        time1 = data["stime"][ns]
        fs = data["samprate"][ns]
        npts = int(data["npts"][ns])

        if ppick > 0:
            data["np1"][ns] = int(np.floor((ppick + data["qtime"] - time1) * fs))
            data["tp1"][ns] = data["np1"][ns] / fs
        else:
            data["np1"][ns] = -12345
            data["tp1"][ns] = -12345

        if spick > 0:
            data["ns1"][ns] = int(np.floor((spick + data["qtime"] - time1) * fs))
            data["ts1"][ns] = data["ns1"][ns] / fs
        else:
            data["ns1"][ns] = -12345
            data["ts1"][ns] = -12345

        data["codawin"][ns] = 10.0
        qmag = data.get("qmag1", data.get("mb", 0.0))
        if qmag >= 4.0:
            data["tswin"][ns] = 10.0
        elif qmag >= 2.5:
            data["tswin"][ns] = 5.0
        else:
            data["tswin"][ns] = 3.0
        data["tpwin"][ns] = 2.0

        if ppick < 0:
            data["tpwin"][ns] = -12345
            data["pwin_mode"][ns] = ""
        elif spick < 0:
            data["pwin_mode"][ns] = "fixed"
        else:
            data["pwin_mode"][ns] = "psdiff"
            data["tpwin"][ns] = min(data["tpwin"][ns], spick - ppick - 0.1)

        # P-wave spectra
        if data["np1"][ns] > 0 and data["tpwin"][ns] >= 0.5:
            np1 = int(data["np1"][ns] - round_half_up(0.1 * fs))
            if np1 < 1:
                np1 = 1
            nwin = int(round_half_up(data["tpwin"][ns] * fs))
            nend = min(np1 + nwin - 1, npts)
            nwin = nend - np1 + 1
            data["pnwin"][ns] = nwin
            seis = np.asarray(data["seis"][ns]).flatten()
            pseg = seis[np1 - 1:nend]
            if np1 - nwin < 1:
                noise = np.zeros(nwin)
                noise[: np1 - 1] = seis[: np1 - 1]
            else:
                noise = seis[np1 - nwin - 1: np1 - 1]
            data["pseis"][ns] = pseg
            data["pnoise"][ns] = noise
            if nwin >= 50:
                spec, f = _multitaper_psd(pseg, fs, Nfft)
                nspec, _ = _multitaper_psd(noise, fs, Nfft)
                omega = 2 * np.pi * f
                nfg = f > 0
                data["pspec"][ns, 0] = spec[0]
                data["pnspec"][ns, 0] = nspec[0]
                data["pspec"][ns, nfg] = spec[nfg] / omega[nfg]
                data["pnspec"][ns, nfg] = nspec[nfg] / omega[nfg]
                data["pfreq"][ns, nfg] = f[nfg]
                specres = np.interp(data["flog"], f[nfg], data["pspec"][ns, nfg])
                data.setdefault("pspecres", np.zeros((numts, len(data["flog"]))))
                data["pspecres"][ns, : len(specres)] = specres

        # S-wave spectra
        if data["ns1"][ns] > 0:
            ns1 = max(int(data["ns1"][ns]), 1)
            nwin = int(round_half_up(data["tswin"][ns] * fs))
            nend = min(ns1 + nwin - 1, npts)
            nwin = nend - ns1 + 1
            data["snwin"][ns] = nwin
            seis = np.asarray(data["seis"][ns]).flatten()
            sseg = seis[ns1 - 1: ns1 - 1 + nwin]
            if ns1 - nwin < 1:
                snoise = np.zeros(nwin)
                snoise[: ns1 - 1] = seis[: ns1 - 1]
            else:
                snoise = seis[ns1 - nwin - 1: ns1 - 1]
            data["sseis"][ns] = sseg
            data["snoise"][ns] = snoise
            if nwin >= 50:
                spec, f = _multitaper_psd(sseg, fs, Nfft)
                nspec, _ = _multitaper_psd(snoise, fs, Nfft)
                omega = 2 * np.pi * f
                nfg = f > 0
                data["sspec"][ns, 0] = spec[0]
                data["snspec"][ns, 0] = nspec[0]
                data["sspec"][ns, nfg] = spec[nfg] / omega[nfg]
                data["snspec"][ns, nfg] = nspec[nfg] / omega[nfg]
                data["sfreq"][ns, nfg] = f[nfg]
                specres = np.interp(data["flog"], f[nfg], data["sspec"][ns, nfg])
                data.setdefault("sspecres", np.zeros((numts, len(data["flog"]))))
                data["sspecres"][ns, : len(specres)] = specres

            # Coda (approximate)
            ns1_c = max(ns1 - int(round_half_up(0.1 * fs)), 1)
            nwin_c = int(round_half_up(data["codawin"][ns] * fs))
            nend_c = min(ns1_c + nwin_c - 1, npts)
            nwin_c = nend_c - ns1_c + 1
            data["cnwin"][ns] = nwin_c
            seis = np.asarray(data["seis"][ns]).flatten()
            cseg = seis[ns1_c - 1: nend_c]
            if ns1_c - nwin_c < 1:
                cnoise = np.zeros(nwin_c)
                cnoise[: ns1_c - 1] = seis[: ns1_c - 1]
            else:
                cnoise = seis[ns1_c - nwin_c - 1: ns1_c - 1]
            spec, f = _multitaper_psd(cseg, fs, Nfftc)
            nspec, _ = _multitaper_psd(cnoise, fs, Nfftc)
            omega = 2 * np.pi * f
            nfg = f > 0
            data["cseis"][ns] = cseg
            data["cnoise"][ns] = cnoise
            data["cspec"][ns, 0] = spec[0]
            data["cnspec"][ns, 0] = nspec[0]
            data["cspec"][ns, nfg] = spec[nfg] / omega[nfg]
            data["cnspec"][ns, nfg] = nspec[nfg] / omega[nfg]
            data["cfreq"][ns, nfg] = f[nfg]
            specres = np.interp(data["flog"], f[nfg], data["cspec"][ns, nfg])
            data.setdefault("cspecres", np.zeros((numts, len(data["flog"]))))
            data["cspecres"][ns, : len(specres)] = specres

    # Ensure unfilled object array elements are set to empty arrays
    for ns in range(numts):
        if data["pseis"][ns] is None:
            data["pseis"][ns] = np.array([])
        if data["pnoise"][ns] is None:
            data["pnoise"][ns] = np.array([])
        if data["sseis"][ns] is None:
            data["sseis"][ns] = np.array([])
        if data["snoise"][ns] is None:
            data["snoise"][ns] = np.array([])
        if data["cseis"][ns] is None:
            data["cseis"][ns] = np.array([])
        if data["cnoise"][ns] is None:
            data["cnoise"][ns] = np.array([])

    savemat(evfile, {"data": data}, do_compression=True)
