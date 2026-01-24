import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from scipy.stats import gmean
from scipy.ndimage import uniform_filter1d

from DE_Module.Rundeopt import Rundeopt
from sourcepara import sourcepara


def specprocess_floategf_allEQs(evspec, stack_spec, freqfit, ddsig, maglim, 
                                  imethod, min_nspec, nlow):
    """
    Process stacked spectra with floating EGF correction for all earthquakes.
    
    Input:
        evspec: event spectra data
        stack_spec: stacked-spec
        freqfit: frequency range for fitting, e.g., [1, 20]
        ddsig: stress drop range
        maglim: magnitude limit for fitting
        min_nspec: for each stacked bin, required at least min_nspec events
        imethod: fitting method
        nlow: low frequency index
    
    Output:
        EGF: Empirical Green's Function
        stack_spec: updated with correction
        ECSs_info: ECS information
        log10dsigbins: log10 of stress drop bins
        nmagg: magnitude indices
        evspec_DE: filtered event spectra
    """
    global falloffin
    falloffin = 2
    
    isavefig = 1
    
    beta = 3464
    # fact = (0.42 * beta) ** 3
    
    freq = stack_spec[0]['freq']
    nfreqg = np.where(freq >= 0)[0]
    freqg = freq[nfreqg]
    
    n1 = np.where(freqg <= freqfit[0])[0]
    n2 = np.where(freqg >= freqfit[1])[0]
    nfit = np.arange(n1[-1], n2[0] + 1)
    
    xmag = np.array([s['mag'] for s in stack_spec])
    dmag = xmag[1] - xmag[0]
    xnspec = np.array([s['nspec'] for s in stack_spec])
    
    nmaggx = np.where((xmag >= maglim[0] - 0.05) & 
                      (xmag <= maglim[1] + 0.05) & 
                      (xnspec >= min_nspec))[0]
    
    if len(nmaggx) == 0:
        EGF = np.zeros_like(freq)
        stack_spec = np.zeros_like(freq)
        return EGF, stack_spec, None, None, None, None
    
    nmagg = np.arange(nmaggx.min(), nmaggx.max() + 1)
    
    # Initialize new fields in stack_spec
    for i in range(len(stack_spec)):
        stack_spec[i]['fc'] = np.nan
        stack_spec[i]['dsig'] = np.nan
        stack_spec[i]['fc2'] = np.nan
        stack_spec[i]['allspecs'] = []
        stack_spec[i]['evspecind'] = []
        stack_spec[i]['fmom'] = []
    
    ifit = np.zeros(len(stack_spec), dtype=int)
    ifit[nmagg] = 1
    
    # Filter evspec_DE
    evspec_DE = [ev for ev in evspec 
                 if ev['qmagest'] >= xmag[nmagg[0]] - dmag/2 and 
                    ev['qmagest'] < xmag[nmagg[-1]] + dmag/2]
    
    # Assign events to bins
    for iq, ev in enumerate(evspec_DE):
        qmagest = ev['qmagest']
        ind = np.argmin(np.abs(qmagest - xmag[nmagg]))
        evspec_DE[iq]['ibin'] = ind
        
        stack_spec[nmagg[ind]]['allspecs'].append(ev['spec'])
        stack_spec[nmagg[ind]]['evspecind'].append(iq)
        stack_spec[nmagg[ind]]['fmom'].append(ev['qmomest'])
        evspec_DE[iq]['magbin'] = xmag[nmagg[ind]]
    
    # Convert allspecs to arrays
    for i in nmagg:
        if stack_spec[i]['allspecs']:
            stack_spec[i]['allspecs'] = np.array(stack_spec[i]['allspecs'])
            stack_spec[i]['fmom'] = np.array(stack_spec[i]['fmom'])
    
    # Calculate fmom_all for each event
    for iq, ev in enumerate(evspec_DE):
        ibin = ev['ibin']
        evspec_DE[iq]['fmom_all'] = gmean(stack_spec[nmagg[ibin]]['fmom'])
    
    # Organize inputs for objfun
    log10sig = [np.log10(ddsig[0]), np.log10(ddsig[1])]
    specall = np.array([ev['spec'] for ev in evspec_DE])
    
    S = {
        'nfreqg': nfreqg,
        'freqg': freqg,
        'nlow': nlow,
        'stack_spec': stack_spec,
        'nmagg': nmagg,
        'log10sig': log10sig,
        'evspec_DE': evspec_DE,
        'specall': specall,
        'nfit': nfit,
        'min_nspec': min_nspec
    }
    
    # Run DE optimization
    log10dsigbins, outputs = Rundeopt(S)
    
    ECSs_allEQs = outputs['ECSs_med']
    ECSs = outputs['ECSs']
    ECSs_info = {
        'ECSs_allEQs': ECSs_allEQs,
        'ECSs': ECSs
    }
    
    # Update stack_spec with stress drops
    for ibin in range(len(nmagg)):
        im = nmagg[ibin]
        stack_spec[im]['dsig'] = 10 ** log10dsigbins[ibin]
    
    # Calculate EGF
    EGF = np.nanmedian(ECSs_allEQs, axis=0)
    EGF = uniform_filter1d(EGF, size=10, mode='nearest')
    
    # Process each magnitude bin
    specdiff = []
    fcpred = np.zeros(len(stack_spec))
    o0pred = np.zeros(len(stack_spec))
    mompred = np.zeros(len(stack_spec))
    
    for img in range(len(nmagg)):
        im = nmagg[img]
        ifit[im] = 1
        fmom = gmean(stack_spec[im]['fmom'])
        dspec = stack_spec[im]['spec'] - ECSs_allEQs[img, :]
        
        # Call sourcepara function (placeholder - needs actual implementation)
        fc, fc2, o0, _, _, _, _, _, err, dspecfit, iflag, ftest, etest, fcb1, fcb2, _ = \
            sourcepara(dspec, freqg, nfit, nlow, imethod, fmom, 1)
        
        if imethod == 11:
            if fc2 < fc:
                fc, fc2 = fc2, fc
        
        stack_spec[im]['fc'] = fc
        stack_spec[im]['o0'] = o0
        stack_spec[im]['ftest'] = ftest
        stack_spec[im]['etest'] = etest
        stack_spec[im]['specfit'] = dspecfit
        stack_spec[im]['speccorr'] = dspec
        stack_spec[im]['fc2'] = fc2
        specdiff.append(dspecfit - dspec)
        
        if img == 0:
            fcref = fc
            fmomref = gmean(stack_spec[im]['fmom'])
        
        fcpred[im] = fmomref ** (1/3) / gmean(stack_spec[im]['fmom']) ** (1/3) * fcref
        o0pred[im] = np.interp(fcpred[im], freqg, dspec)
        mompred[im] = gmean(stack_spec[im]['fmom'])
    
    nfcfit = np.where(fcpred > 0)[0]
    
    # Save figures
    if isavefig == 1:
        plt.close('all')
        
        fig = plt.figure(100, figsize=(10, 8))
        
        ng = [i for i in range(len(stack_spec)) if stack_spec[i]['nspec'] > 1]
        
        for img in range(len(ng)):
            im = ng[img]
            if ifit[im] == 1:
                specx = stack_spec[im]['speccorr']
                plt.semilogx(freqg, stack_spec[im]['speccorr'], 'k')
                plt.semilogx(freqg, stack_spec[im]['specfit'], 'r:')
                corr_spec = stack_spec[im]['spec'] - stack_spec[im]['speccorr'] + np.mean(specx[nlow])
                plt.semilogx(freqg, corr_spec, 'b--')
                
                textm = f"M{stack_spec[im]['mag']:5.2f} NEQ={stack_spec[im]['nspec']} Δσ={stack_spec[im]['dsig']:5.2f}MPa"
                plt.text(freqg[1], specx[1] + 0.1, textm)
                
                specamp = np.interp(stack_spec[im]['fc'], freqg, stack_spec[im]['specfit'])
                plt.semilogx(stack_spec[im]['fc'], specamp, 'ro')
                
                if stack_spec[im]['fc2'] < 50:
                    specamp2 = np.interp(stack_spec[im]['fc2'], freqg, stack_spec[im]['specfit'])
                    plt.semilogx(stack_spec[im]['fc2'], specamp2, 'bd')
            
            elif stack_spec[im]['nspec'] > 1 and len(stack_spec[im]['fmom']) > 0:
                specx = stack_spec[im]['speccorr'][nfreqg]
                plt.semilogx(freqg, specx, 'k-.')
                textm = f"M{stack_spec[im]['mag']:5.2f} NEQ={stack_spec[im]['nspec']} Δσ={stack_spec[im]['dsig']:5.2f}MPa"
                plt.text(freqg[1], specx[1] + 0.1, textm)
        
        plt.semilogx(fcpred[nfcfit], o0pred[nfcfit], 'k--')
        plt.xlabel('frequency (Hz)')
        plt.ylabel('amplitude')
        
        median_dsig = np.median([stack_spec[i]['dsig'] for i in nmagg])
        titlex = f'ECS corrected spectra, median Δσ={median_dsig:5.2f}MPa'
        plt.title(titlex)
        
        plt.savefig('stacked-spectra-after-EGF.jpg')
        plt.close(100)
    
    return EGF, stack_spec, ECSs_info, log10dsigbins, nmagg, evspec_DE