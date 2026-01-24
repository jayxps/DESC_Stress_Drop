import numpy as np
import matplotlib.pyplot as plt


def specprocess_stackmag(evspec, nspecmin, maginfo, isavefig, mbin):
    """
    Stack spectra by magnitude bins and compute average spectra and moments.

    Parameters
    ----------
    evspec : list of dict
        Each dict needs: 'freq', 'spec', 'nspec', 'qmag', 'qmagest', 'qmagresid'.
    nspecmin : int
        Minimum number of spectra to include an event in stacking.
    maginfo : dict
        Contains: itype, fmw, y0, slope.
    isavefig : bool
        Whether to save the stacking plot.
    mbin : float
        Magnitude bin width.

    Returns
    -------
    stack_spec : list of dict
        Stacked spectra per magnitude bin with keys: spec, nspec, mag, freq, fmom.
    """
    itype = maginfo['itype']
    fmw = maginfo['fmw']
    y0 = maginfo['y0']
    slopeline = maginfo['slope']

    freq = np.asarray(evspec[0]['freq'])

    ng = np.where(np.array([ev['nspec'] for ev in evspec]) >= nspecmin)[0]
    specmag = np.array([evspec[i]['qmagest'] for i in ng])
    mag = np.array([evspec[i]['qmag'] for i in ng])

    if itype == 1:
        qmag = mag
    else:
        qmag = specmag
        a = (fmw + 10.7) * 1.5 * slopeline - fmw

    resid = np.array([evspec[i]['qmagresid'] for i in ng])

    magbins = np.arange(-1, 8 + mbin, mbin)
    im_total = len(magbins) + 1
    im_ev = np.floor((qmag - magbins.min()) / mbin).astype(int) + 1

    stack_spec = []
    for im in range(im_total):
        stack_spec.append({
            'spec': np.zeros_like(freq),
            'nspec': 0,
            'mag': (im) * mbin - mbin / 2 + magbins.min(),
            'freq': freq,
            'fmom': 0.0,
        })

    for im in range(im_total):
        nin = np.where((im_ev == im) & (np.abs(resid) <= 1.5))[0]
        nspec = len(nin)
        if nspec < 1:
            continue
        spec_sum = np.zeros_like(freq)
        for idx in nin:
            iq = ng[idx]
            spec_sum += evspec[iq]['spec']
        stack_spec[im]['spec'] = spec_sum / nspec
        stack_spec[im]['nspec'] = nspec
        if itype == 1:
            stack_spec[im]['fmom'] = 10 ** ((stack_spec[im]['mag'] + 10.7) * 1.5) / 1e7
        else:
            stack_spec[im]['fmom'] = 10 ** ((stack_spec[im]['mag'] + a) / slopeline) / 1e7

    if isavefig:
        plt.ioff()
        fig = plt.figure(100)
        for im in range(im_total):
            if stack_spec[im]['nspec'] < 1:
                continue
            plt.semilogx(freq, stack_spec[im]['spec'])
            textm = f"M{stack_spec[im]['mag']:.2f} NS {stack_spec[im]['nspec']}"
            plt.text(freq[1], stack_spec[im]['spec'][1], textm)
        plt.xlabel('frequency (Hz)')
        plt.ylabel('amplitude')
        plt.xlim([0.01, freq.max()])
        plt.savefig('stackmag.pdf')
        plt.savefig('stackmag.png', dpi=300)
        plt.close(fig)

    return stack_spec
