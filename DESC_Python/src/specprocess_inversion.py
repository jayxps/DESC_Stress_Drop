import numpy as np
from specprocess_binsolve import specprocess_binsolve


def specprocess_inversion(spec, P, freq):
    """
    Invert spectra to separate event, station, distance, depth, and magnitude terms.

    Parameters
    ----------
    spec : list of dict
        Each dict contains keys: inq, ind, ins, inx, inm, tstar, spec (array),
        qid, qtime, qlat, qlon, qdep, qmag, stlat, stlon, stelev, stid.
        Bin indices (inq, ind, ins, inx, inm) are expected to be 1-based integers.
    P : dict
        Inversion parameters with keys: nitmax (iterations), itype (1 median, 2 mean, 3 weighted), dtt (distance bin width).
    freq : ndarray
        Frequency vector.

    Returns
    -------
    evspec : list of dict
        Event spectra and metadata.
    stspec : list of dict
        Station spectra and metadata.
    distspec : list of dict
        Distance spectra and metadata.
    """

    nf0 = len(freq)

    nq = max(s['inq'] for s in spec)
    ev_spec = np.zeros((nf0, nq))

    ndep = max(s['ind'] for s in spec)
    dep_spec = np.zeros((nf0, ndep))

    nsta = max(s['ins'] for s in spec)
    st_spec = np.zeros((nf0, nsta))

    ndmax = max(s['inx'] for s in spec)
    dist_spec = np.zeros((nf0, ndmax))

    nmag = max(s['inm'] for s in spec)
    mag_spec = np.zeros((nf0, nmag))

    depterm = np.zeros(ndep)
    distterm = np.zeros(ndmax)
    magterm = np.zeros(nmag)
    stterm = np.zeros(nsta)
    evterm = np.zeros(nq)

    npts = len(spec)

    # Pre-extract index arrays (1-based as in MATLAB)
    idx_inx = np.array([s['inx'] for s in spec], dtype=int)
    idx_ins = np.array([s['ins'] for s in spec], dtype=int)
    idx_inq = np.array([s['inq'] for s in spec], dtype=int)
    idx_ind = np.array([s['ind'] for s in spec], dtype=int)
    idx_inm = np.array([s['inm'] for s in spec], dtype=int)

    for it in range(P['nitmax']):
        y = np.zeros(npts)
        yspec = np.zeros((nf0, npts))

        # Match MATLAB evaluation order exactly (no Kahan summation)
        for ispec, sp in enumerate(spec):
            is_ = sp['ins'] - 1
            iq = sp['inq'] - 1
            im = sp['inm'] - 1
            ix = sp['inx'] - 1
            id_ = sp['ind'] - 1

            # MATLAB: spec(ispec).tstar - depterm(id) - magterm(im) - stterm(is) - evterm(iq)
            y[ispec] = sp['tstar'] - depterm[id_] - magterm[im] - stterm[is_] - evterm[iq]
            
            # MATLAB: spec(ispec).spec' - dep_spec(:,id) - mag_spec(:,im) - st_spec(:,is) - ev_spec(:,iq)
            yspec[:, ispec] = sp['spec'] - dep_spec[:, id_] - mag_spec[:, im] - st_spec[:, is_] - ev_spec[:, iq]

        distterm, _ = specprocess_binsolve(y.reshape(1, -1), idx_inx, ndmax, P['itype'])
        distterm = distterm.flatten()
        dist_spec, _ = specprocess_binsolve(yspec, idx_inx, ndmax, P['itype'])

        for ispec, sp in enumerate(spec):
            is_ = sp['ins'] - 1
            iq = sp['inq'] - 1
            im = sp['inm'] - 1
            ix = sp['inx'] - 1
            id_ = sp['ind'] - 1

            # MATLAB: spec(ispec).tstar - depterm(id) - magterm(im) - distterm(ix) - evterm(iq)
            y[ispec] = sp['tstar'] - depterm[id_] - magterm[im] - distterm[ix] - evterm[iq]
            
            # MATLAB: spec(ispec).spec' - dep_spec(:,id) - mag_spec(:,im) - dist_spec(:,ix) - ev_spec(:,iq)
            yspec[:, ispec] = sp['spec'] - dep_spec[:, id_] - mag_spec[:, im] - dist_spec[:, ix] - ev_spec[:, iq]

        stterm, _ = specprocess_binsolve(y.reshape(1, -1), idx_ins, nsta, P['itype'])
        stterm = stterm.flatten()
        st_spec, _ = specprocess_binsolve(yspec, idx_ins, nsta, P['itype'])

        for ispec, sp in enumerate(spec):
            is_ = sp['ins'] - 1
            iq = sp['inq'] - 1
            im = sp['inm'] - 1
            ix = sp['inx'] - 1
            id_ = sp['ind'] - 1

            # MATLAB: spec(ispec).tstar - depterm(id) - magterm(im) - distterm(ix) - stterm(is)
            y[ispec] = sp['tstar'] - depterm[id_] - magterm[im] - distterm[ix] - stterm[is_]
            
            # MATLAB: spec(ispec).spec' - dep_spec(:,id) - mag_spec(:,im) - dist_spec(:,ix) - st_spec(:,is)
            yspec[:, ispec] = sp['spec'] - dep_spec[:, id_] - mag_spec[:, im] - dist_spec[:, ix] - st_spec[:, is_]

        evterm, _ = specprocess_binsolve(y.reshape(1, -1), idx_inq, nq, P['itype'])
        evterm = evterm.flatten()
        ev_spec, _ = specprocess_binsolve(yspec, idx_inq, nq, P['itype'])

        # RMS computation - match MATLAB exactly (simple accumulation, no Kahan)
        sum2 = 0.0
        sumy = 0.0
        
        for ispec, sp in enumerate(spec):
            is_ = sp['ins'] - 1
            iq = sp['inq'] - 1
            im = sp['inm'] - 1
            ix = sp['inx'] - 1
            id_ = sp['ind'] - 1

            # MATLAB: y(ispec) = spec(ispec).tstar - depterm(id) - magterm(im) - distterm(ix) - stterm(is) - evterm(iq)
            y[ispec] = sp['tstar'] - depterm[id_] - magterm[im] - distterm[ix] - stterm[is_] - evterm[iq]
            sumy = sumy + y[ispec] ** 2

            # MATLAB: yspec(:,ispec) = spec(ispec).spec' - dep_spec(:,id) - mag_spec(:,im) - dist_spec(:,ix) - st_spec(:,is) - ev_spec(:,iq)
            yspec[:, ispec] = sp['spec'] - dep_spec[:, id_] - mag_spec[:, im] - dist_spec[:, ix] - st_spec[:, is_] - ev_spec[:, iq]
            sum2 = sum2 + np.sum(yspec[:, ispec] ** 2)

        rms2 = np.sqrt(sum2) / (npts * nf0)
        rms = np.sqrt(sumy) / npts
        print(f"it,rms,rms2= {it+1} {rms:12.8f} {rms2:12.8f}")

    # Build outputs
    evspec = []
    iq_out = 0
    for i in range(1, nq + 1):
        nin = np.where(idx_inq == i)[0]
        if len(nin) == 0:
            continue
        iq_out += 1
        first = spec[nin[0]]
        evspec.append({
            'spec': ev_spec[:, i - 1],
            'nspec': len(nin),
            'freq': freq,
            'qid': first['qid'],
            'qtime': first['qtime'],
            'qlat': first['qlat'],
            'qlon': first['qlon'],
            'qdep': first['qdep'],
            'qmag': first['qmag'],
            'tstar': evterm[i - 1]
        })

    stspec = []
    is_out = 0
    for i in range(1, nsta + 1):
        nin = np.where(idx_ins == i)[0]
        if len(nin) == 0:
            continue
        is_out += 1
        first = spec[nin[0]]
        stspec.append({
            'spec': st_spec[:, i - 1],
            'nspec': len(nin),
            'freq': freq,
            'stlat': first['stlat'],
            'stlon': first['stlon'],
            'stelev': first['stelev'],
            'stid': first['stid'],
            'tstar': stterm[i - 1]
        })

    distspec = []
    ix_out = 0
    for i in range(1, ndmax + 1):
        nin = np.where(idx_inx == i)[0]
        if len(nin) == 0:
            continue
        ix_out += 1
        distspec.append({
            'spec': dist_spec[:, i - 1],
            'nspec': len(nin),
            'freq': freq,
            'tt': i * P['dtt'] - P['dtt'] / 2,
            'tstar': distterm[i - 1]
        })

    return evspec, stspec, distspec
