import numpy as np
from decimal import Decimal, getcontext


def specprocess_specmag(evspec, freqlim, maginfo):
    """
    Estimate spectral magnitudes and moments for events.

    Parameters
    ----------
    evspec : list of dict
        Each dict must contain keys: 'freq', 'spec', 'qmag', 'qid'.
    freqlim : tuple/list
        Frequency bounds (low, high) used to average spectral amplitude.
    maginfo : dict
        Contains keys:
            - itype: 1 uses catalog magnitude directly; 2 uses spectral scaling
            - y0: intercept of magnitude fit
            - slope: slope of magnitude fit
            - fmw: magnitude pivot used in scaling

    Returns
    -------
    evspec : list of dict
        Updated with qmagest, qmomest, qmagresid
    """
    itype = maginfo['itype']
    y0 = maginfo['y0']
    slopeline = maginfo['slope']
    fmw = maginfo['fmw']

    if itype == 1:
        for ev in evspec:
            ev['qmagest'] = ev['qmag']
            ev['qmomest'] = 10 ** (1.5 * (ev['qmag'] + 10.7)) / 1e7
            ev['qmagresid'] = 0.0
        return evspec

    if itype == 2:
        # Set high precision for Decimal operations (15 significant figures)
        getcontext().prec = 15
        
        a = (fmw + 10.7) * 1.5 * slopeline - fmw

        freq = np.array(evspec[0]['freq'])
        nfx = (freq >= freqlim[0]) & (freq <= freqlim[1])

        with open('tempmom', 'w') as fid:
            for ev in evspec:
                spec = np.array(ev['spec'])
                if spec.size == 0:
                    continue
                qmom = np.mean(spec[nfx]) if np.any(nfx) else 0.0
                if qmom != 0:
                    qmagest = y0 + qmom * slopeline
                else:
                    qmagest = 0.0
                
                # Use Decimal for high-precision exponentiation in moment calculation
                # This mitigates floating-point errors in the 10^x operation
                exponent = Decimal(str(qmagest + a)) / Decimal(str(slopeline))
                qmomest_decimal = Decimal(10) ** exponent
                qmomest = float(qmomest_decimal)
                
                ev['qmomest'] = qmomest / 1e7  # dyne-cm to Nm
                ev['qmagest'] = qmagest
                ev['qmagresid'] = ev['qmag'] - ev['qmagest']
                fid.write(f"{ev['qmag']:10.6f} {qmom:10.6f} {np.log10(qmomest):10.6f} {qmagest:10.6f} {ev['qid']:13d} {a:10.6f}\n")

        return evspec

    # If an unknown itype is provided, just return input unchanged.
    return evspec
