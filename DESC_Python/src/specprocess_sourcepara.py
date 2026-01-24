import os
import numpy as np
import matplotlib.pyplot as plt
from sourcepara import sourcepara


def specprocess_sourcepara(rdir, freqfit, evspec, imethod, iwave, EGF, min_nspec, isavefig):
    """
    Estimate source parameters using a constant EGF across events.

    Parameters
    ----------
    rdir : str
        Output directory.
    freqfit : tuple/list
        Frequency range for fitting (low, high).
    evspec : list of dict
        Event spectra with keys: 'freq', 'spec', 'qmagest', 'qmag', 'qid', 'qmomest',
        'qlat', 'qlon', 'qdep', 'qtime', 'nspec'.
    imethod : int
        Source fitting method (passed to sourcepara).
    iwave : int
        Wave type flag for sourcepara.
    EGF : ndarray
        Empirical Green's Function spectrum to subtract.
    min_nspec : int
        Minimum spectra count to include an event.
    isavefig : int/bool
        Whether to save diagnostic plots.

    Returns
    -------
    para_EGF : list of dict
        Fitted parameters per qualifying event.
    """

    if not os.path.isdir(rdir):
        os.makedirs(rdir, exist_ok=True)
    if not rdir.endswith('/'):
        rdir = rdir + '/'

    freq = np.array(evspec[0]['freq'])
    frqlim = [0.5, 50]
    flow = [0.5, 5]
    fhigh = [20, 30]
    nlow = np.where((freq >= flow[0]) & (freq <= flow[1]))[0]
    nhigh = np.where((freq >= fhigh[0]) & (freq <= fhigh[1]))[0]
    nfit = np.where((freq >= freqfit[0]) & (freq <= freqfit[1]))[0]

    out_path = os.path.join(rdir, 'result-constant-egf.txt')
    fmt_header = (
        "    i1   magest magcat              id1       o0             fc1"
        "          r    delsig    asig        G          Es-all   Es-model"
        "   Es-obs   Es-1       Es-2   error fcbound1 fcbound2\n"
    )
    fmt_line = (
        " %5d"  # i1
        "%8.2f"  # magest
        "%8.2f"  # magcat
        "%15d"   # id1
        "%15.2f" # o0
        "%15.2f" # fc1
        "%15.2e" # r
        "%10.3e" # delsig
        "%10.3e" # asig
        "%11.3e" # G
        "%10.3e" # Es-all
        "%10.3e" # Es-model
        "%10.3e" # Es-obs
        "%10.3e" # Es-1
        "%10.3e" # Es-2
        " %10.3e" # error
        " %10.3e" # fcb1
        " %10.3e" # fcb2
        "\n"
    )

    para_EGF = []
    nqsolve = 0

    with open(out_path, 'w') as fid4:
        fid4.write(fmt_header)

        for iq, ev in enumerate(evspec, start=1):
            if ev['nspec'] < min_nspec:
                continue

            dspecc = np.array(ev['spec']) - EGF
            if (dspecc.max() - dspecc.min()) <= 1e-2:
                continue

            nqsolve += 1
            fc, fc2, o0, Es, delsig, asig, r, G, err, dspecfit, iflag, ftest, etest, fcb1, fcb2, _ = \
                sourcepara(dspecc, freq, nfit, nlow, imethod, ev['qmomest'], iwave)

            lhratio = dspecc[nlow].mean() - dspecc[nhigh].mean()
            entry = {
                'specfit': dspecfit,
                'ftest': ftest,
                'etest': etest,
                'exist': 1,
                'fc': fc,
                'fc2': fc2,
                'err': err,
                'Es': Es,
                'delsig': delsig,
                'asig': asig,
                'fcb1': fcb1,
                'fcb2': fcb2,
                'qid': ev['qid'],
                'qlat': ev['qlat'],
                'qlon': ev['qlon'],
                'qdep': ev['qdep'],
                'qmag': ev['qmagest'],
                'qmom': ev['qmomest'],
                'qmagcat': ev['qmag'],
                'qtime': ev['qtime'],
                'lhratio': lhratio,
                'G': G,
                'r': r,
                'nspec': ev['nspec'],
            }
            para_EGF.append(entry)

            if nqsolve % 100 == 0:
                print(f"event: {nqsolve}, id: {ev['qid']}")

            if np.min(Es) > 0 and iflag == 1 and fc < 100:
                fid4.write(fmt_line % (
                    iq, ev['qmagest'], ev['qmag'], ev['qid'],
                    o0, fc, r, delsig, asig,
                    G, Es[0], Es[1], Es[2], Es[3], Es[4],
                    err, fcb1, fcb2
                ))

                if isavefig:
                    try:
                        plt.ioff()
                        fig = plt.figure(1, figsize=(8, 6))
                        plt.subplot(2, 2, 1)
                        plt.semilogx(freq, dspecc)
                        plt.semilogx(freq, dspecfit, '-.')
                        plt.xlabel('frequency')
                        plt.xlim(frqlim)
                        plt.ylabel('amplitude')

                        plt.subplot(2, 2, 2)
                        plt.plot(ftest, etest / np.min(etest))
                        ytest = np.min(etest) * 1.05 / np.min(etest)
                        plt.plot([ftest.min(), ftest.max()], [ytest, ytest], '-.')
                        plt.plot(fc, err / np.min(etest), '*')
                        plt.xlabel('frequency')
                        plt.ylabel('variance')
                        plt.xlim([ftest.min(), ftest.max()])
                        plt.title(f"shearer-stacking-{ev['qid']}")

                        fileout = os.path.join(rdir, f"fig-shearer-{ev['qid']}")
                        plt.savefig(fileout, format='pdf')
                        plt.close(fig)
                    except Exception:
                        # Skip plotting errors
                        pass

    return para_EGF
