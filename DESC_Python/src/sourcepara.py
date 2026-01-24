import numpy as np
from scipy.optimize import fminbound, minimize
from scipy.interpolate import UnivariateSpline, interp1d
from scipy.integrate import trapezoid

# Global variable for falloff rate
falloffin = 2


def sourcepara(dspec, f, nfrq, nlow, imethod, moment, iwave):
    """
    Estimate source parameters from displacement spectra.
    
    Parameters
    ----------
    dspec : array-like
        Displacement spectrum
    f : array-like
        Frequency array
    nfrq : array-like
        Indices for fitting frequency range
    nlow : array-like
        Indices for low frequency range
    imethod : int
        Fitting method (1-11, see method descriptions)
    moment : float
        Seismic moment
    iwave : int
        Wave type (1 for P-wave, 2 for S-wave)
    
    Returns
    -------
    fc : float
        Corner frequency for target event
    fc2 : float
        Second corner frequency (for EGF or double-corner models)
    o0 : float
        Low-frequency spectral level
    Es : array
        Radiated energy [total, modeled, measured, low_extrap, high_extrap]
    delsig : float
        Stress drop in MPa
    asig : float
        Apparent stress in MPa
    r : float
        Source radius
    G : float
        Fracture energy
    err : float
        Fitting error
    dspecfit : array
        Fitted spectrum
    iflag : int
        Optimization flag
    ftest : array
        Test frequencies for corner frequency uncertainty
    etest : array
        Error at test frequencies
    fcb1 : float
        Lower bound of corner frequency
    fcb2 : float
        Upper bound of corner frequency
    falloff : float
        Falloff rate
    """
    global falloffin
    
    # Physical constants
    beta = 3464
    pho = 2700
    mu = pho * beta ** 2
    kfact = 0.26
    
    # Wave-dependent parameters
    if iwave == 1:  # P-wave
        fact = (0.42 * beta) ** 3
        alpha = 1.732 * beta
        scale = 2.0 / (15 * np.pi * pho * alpha ** 5) * 14.4
    else:  # S-wave
        kfact = 0.26
        fact = (0.2766 * beta) ** 3
        scale = 1.07 / (5 * np.pi * pho * beta ** 5)
    
    # Resample to log-spaced frequencies if needed
    fmin = np.min(f[nfrq])
    fmax = np.max(f[nfrq])
    dfx = np.diff(f)
    dfxmean = np.abs(dfx - np.mean(dfx))
    
    if np.max(dfxmean) < 1e-3:  # Linear spaced
        t = np.logspace(np.log10(fmin), np.log10(fmax), 100)
        spl = UnivariateSpline(f[nfrq], dspec[nfrq], s=0)
        y = spl(t)
    else:  # Already log spaced
        t = f[nfrq]
        y = dspec[nfrq]
    
    # Set up method parameters
    im = imethod
    npar = 2
    falloff = falloffin
    
    if imethod == 3:
        im = 5
        npar = 3
        falloff = falloffin
    elif imethod == 2:
        im = 1
        npar = 3
        falloff = falloffin
    elif imethod == 4:
        im = 6
        npar = 4
    elif imethod == 8:
        im = 7
        npar = 4
    elif imethod == 9:
        im = 10
        npar = 4
    elif imethod in [10, 11]:
        npar = 3
    
    # Initial parameter guess
    nfx = nlow
    x0 = np.mean(dspec[nfx])
    
    if npar == 3:
        pstart = [x0, 10, 50]
    elif npar == 2:
        pstart = [x0, 10]
    elif npar == 4:
        pstart = [x0, 10, 50, 2]
    
    if imethod in [6, 7]:
        npar = 3
        pstart = [x0, 10, 2]
    
    if imethod == 10:
        pstart = [x0, 10, 2]
    
    if imethod == 11:
        pstart = [x0, 2, 20]
    
    # Optimize parameters
    if imethod != 11:
        result = minimize(
            lambda x: method(x, t, y, imethod),
            pstart,
            method='Nelder-Mead',
            options={'xatol': 1e-10, 'maxfev': 1500, 'maxiter': 500}
        )
        p_method3 = result.x
        iflag = 1 if result.success else 0
    else:
        # Grid search for double-corner frequency model
        errmin = 99999
        fc_grid = np.arange(0.1, 100.1, 0.1)
        p_method3 = np.zeros(3)
        
        for ifc1 in range(len(fc_grid)):
            for ifc2 in range(ifc1 + 1, len(fc_grid)):
                fc1 = fc_grid[ifc1]
                fc2 = fc_grid[ifc2]
                upred = -0.5 * np.log10(1 + (t / fc1) ** 2) - 0.5 * np.log10(1 + (t / fc2) ** 2)
                udiff = y - upred
                omega0 = np.mean(udiff)
                ndiff = udiff - omega0
                err = np.linalg.norm(ndiff)
                
                if err < errmin:
                    errmin = err
                    p_method3[1] = fc1
                    p_method3[2] = fc2
                    p_method3[0] = omega0
        
        iflag = 1
    
    err = method(p_method3, t, y, imethod)
    u3 = fitmethod(p_method3, f, imethod)
    
    dspecfit = u3
    o0 = p_method3[0]
    fc = np.abs(p_method3[1])
    
    if npar == 3:
        fc2 = np.abs(p_method3[2])
    elif npar == 2:
        fc2 = 1000
    elif npar == 4:
        fc2 = np.abs(p_method3[2])
        falloff = p_method3[3]
    
    if imethod in [6, 7]:
        fc2 = 1000
    
    if imethod == 10:
        fc2 = 1000
        falloff = p_method3[2]
    
    # Test corner frequency uncertainty
    fc1 = fc - 5
    if fc1 <= 0:
        fc1 = 0.05
    ftest = np.linspace(fc1, fc + 5, 50)
    etest = np.zeros(len(ftest))
    
    for i in range(len(ftest)):
        p_test = p_method3.copy()
        result = minimize(
            lambda x: testmethod(x, ftest[i], t, y, imethod),
            p_test,
            method='Nelder-Mead',
            options={'xatol': 1e-10, 'maxfev': 1500, 'maxiter': 500}
        )
        etest[i] = result.fun
    
    nlow_idx = np.where(etest <= 1.05 * np.min(etest))[0]
    
    if len(nlow_idx) > 0:
        fcb1 = ftest[nlow_idx[0]]
        fcb2 = ftest[nlow_idx[-1]]
    else:
        fcb1 = -1
        fcb2 = -1
    
    # Prepare parameters for energy calculation
    p = np.zeros(3)
    p[0] = np.log10(moment)
    p[1] = p_method3[1]
    
    if imethod == 4:
        p[2] = p_method3[3]
        falloff = p[2]
    elif imethod == 8:
        p[2] = p_method3[3]
        falloff = p[2]
    elif imethod == 6:
        falloff = p_method3[2]
        p[2] = falloff
    elif imethod == 7:
        falloff = p_method3[2]
        p[2] = falloff
    elif imethod == 9:
        falloff = p_method3[3]
        p[2] = falloff
    elif imethod == 10:
        falloff = p_method3[2]
        p[2] = falloff
    elif imethod == 11:
        p[2] = p_method3[2]
    
    dspecx = dspec - o0 + p[0]
    nfrqx = nfrq
    
    Es, fex, dspecex, dspecfitx = energy(f, dspecx, nfrqx, p, im)
    
    # Calculate source parameters
    delsig = (fc ** 3) * moment / fact / 1e6  # in MPa
    
    Es = Es * scale
    asig = mu * Es[0] / moment / 1e6  # in MPa
    r = kfact * beta / fc
    A = np.pi * r ** 2
    D = moment / (mu * A)
    G = 0.5 * (delsig - 2 * asig) * D * 1e6
    
    return fc, fc2, o0, Es, delsig, asig, r, G, err, dspecfit, iflag, ftest, etest, fcb1, fcb2, falloff


def testmethod(p, p2, t, y, im):
    """
    Test method for corner frequency uncertainty estimation.
    
    Parameters
    ----------
    p : array-like
        Parameter vector
    p2 : float
        Fixed corner frequency
    t : array-like
        Frequency array
    y : array-like
        Observed spectrum
    im : int
        Method index
    
    Returns
    -------
    err : float
        Fitting error
    """
    global falloffin
    
    if im == 1:
        u = p[0] - np.log10(1 + (t / p2) ** falloffin)
    elif im == 2:
        u = p[0] - np.log10(1 + (t / p2) ** falloffin) + np.log10(1 + (t / p[2]) ** falloffin)
    elif im == 3:
        u = p[0] - np.log10(1 + (t / p2) ** (2 * falloffin)) / 2 + np.log10(1 + (t / p[2]) ** (2 * falloffin)) / 2
    elif im == 4:
        u = p[0] - (1.0 / p[3]) * np.log10(1 + (t / p2) ** (2 * p[3])) + (1.0 / p[3]) * np.log10(1 + (t / p[2]) ** (2 * p[3]))
    elif im == 5:
        u = p[0] - np.log10((1 + (t / p2) ** (2 * falloffin)) ** 0.5)
    elif im == 6:
        u = p[0] - np.log10((1 + (t / p2) ** (2 * p[2])) ** (1.0 / p[2]))
    elif im == 7:
        u = p[0] - np.log10(1 + (t / p2) ** p[2])
    elif im == 8:
        u = p[0] - np.log10(1 + (t / p2) ** p[3]) + np.log10(1 + (t / p[2]) ** p[3])
    elif im == 9:
        u = p[0] - 0.5 * np.log10(1 + (t / p2) ** (2 * p[3])) + 0.5 * np.log10(1 + (t / p[2]) ** (2 * p[3]))
    elif im == 10:
        u = p[0] - 0.5 * np.log10(1 + (t / p2) ** (2 * p[2]))
    elif im == 11:
        u = p[0] - 0.5 * np.log10(1 + (t / p[1]) ** 2) - 0.5 * np.log10(1 + (t / p[2]) ** 2)
    
    err = np.linalg.norm(u - y)
    return err


def method(p, t, y, im):
    """
    Fit various spectral models to observed spectra.
    
    Model descriptions:
    - Brune family: 1, 2, 7, 8, 11
    - Boatwright family: 3, 4, 5, 6, 9, 10
    
    Parameters
    ----------
    p : array-like
        Parameter vector
    t : array-like
        Frequency array
    y : array-like
        Observed spectrum
    im : int
        Method index
    
    Returns
    -------
    err : float
        Fitting error
    """
    global falloffin
    
    if im == 1:  # Pure Brune model
        u = p[0] - np.log10(1 + (t / p[1]) ** falloffin)
    elif im == 2:  # Brune model with ratio
        u = p[0] - np.log10(1 + (t / p[1]) ** falloffin) + np.log10(1 + (t / p[2]) ** falloffin)
    elif im == 3:  # Boatwright with ratio
        u = p[0] - np.log10(1 + (t / p[1]) ** (2 * falloffin)) / 2 + np.log10(1 + (t / p[2]) ** (2 * falloffin)) / 2
    elif im == 4:  # Boatwright ratio with gamma float
        u = p[0] - (1.0 / p[3]) * np.log10(1 + (t / p[1]) ** (2 * p[3])) + (1.0 / p[3]) * np.log10(1 + (t / p[2]) ** (2 * p[3]))
    elif im == 5:  # Boatwright (no ratio)
        u = p[0] - np.log10((1 + (t / p[1]) ** (2 * falloffin)) ** 0.5)
    elif im == 6:  # Boatwright with gamma float (no ratio)
        u = p[0] - np.log10((1 + (t / p[1]) ** (2 * p[2])) ** (1.0 / p[2]))
    elif im == 7:  # Brune with n float (no ratio)
        u = p[0] - np.log10(1 + (t / p[1]) ** p[2])
    elif im == 8:  # Brune ratio with n float
        u = p[0] - np.log10(1 + (t / p[1]) ** p[3]) + np.log10(1 + (t / p[2]) ** p[3])
    elif im == 9:  # Boatwright ratio with n float
        u = p[0] - 0.5 * np.log10(1 + (t / p[1]) ** (2 * p[3])) + 0.5 * np.log10(1 + (t / p[2]) ** (2 * p[3]))
    elif im == 10:  # Boatwright with n float
        u = p[0] - 0.5 * np.log10(1 + (t / p[1]) ** (2 * p[2]))
    elif im == 11:  # Double-corner frequency model
        u = p[0] - 0.5 * np.log10(1 + (t / p[1]) ** 2) - 0.5 * np.log10(1 + (t / p[2]) ** 2)
    
    err = np.linalg.norm(u - y)
    return err


def fitmethod(p, t, im):
    """
    Generate fitted spectrum using specified model.
    
    Parameters
    ----------
    p : array-like
        Parameter vector
    t : array-like
        Frequency array
    im : int
        Method index
    
    Returns
    -------
    u : ndarray
        Fitted spectrum
    """
    global falloffin
    
    if im == 1:
        u = p[0] - np.log10(1 + (t / p[1]) ** falloffin)
    elif im == 2:
        u = p[0] - np.log10(1 + (t / p[1]) ** falloffin) + np.log10(1 + (t / p[2]) ** falloffin)
    elif im == 3:
        u = p[0] - np.log10(1 + (t / p[1]) ** (2 * falloffin)) / 2 + np.log10(1 + (t / p[2]) ** (2 * falloffin)) / 2
    elif im == 4:
        u = p[0] - (1.0 / p[3]) * np.log10(1 + (t / p[1]) ** (2 * p[3])) + (1.0 / p[3]) * np.log10(1 + (t / p[2]) ** (2 * p[3]))
    elif im == 5:
        u = p[0] - np.log10((1 + (t / p[1]) ** (2 * falloffin)) ** 0.5)
    elif im == 6:
        u = p[0] - np.log10((1 + (t / p[1]) ** (2 * p[2])) ** (1.0 / p[2]))
    elif im == 7:
        u = p[0] - np.log10(1 + (t / p[1]) ** p[2])
    elif im == 8:
        u = p[0] - np.log10(1 + (t / p[1]) ** p[3]) + np.log10(1 + (t / p[2]) ** p[3])
    elif im == 9:
        u = p[0] - 0.5 * np.log10(1 + (t / p[1]) ** (2 * p[3])) + 0.5 * np.log10(1 + (t / p[2]) ** (2 * p[3]))
    elif im == 10:
        u = p[0] - 0.5 * np.log10(1 + (t / p[1]) ** (2 * p[2]))
    elif im == 11:
        u = p[0] - 0.5 * np.log10(1 + (t / p[1]) ** 2) - 0.5 * np.log10(1 + (t / p[2]) ** 2)
    
    return u


def energy(f, dspec, nfrq, p, im):
    """
    Estimate radiated energy from spectra.
    
    Method: Direct integration for frequencies in nfrq, then extrapolate
    spectra using fitted model.
    
    Parameters
    ----------
    f : array-like
        Frequency array
    dspec : array-like
        Displacement spectrum
    nfrq : array-like
        Frequency indices for integration
    p : array-like
        Model parameters
    im : int
        Method index
    
    Returns
    -------
    Esx : ndarray
        Energy array [total, modeled, measured, low_extrap, high_extrap]
    fex : ndarray
        Extended frequency array
    dspecex : ndarray
        Extended spectrum
    dspecfit : ndarray
        Fitted spectrum on extended frequency array
    """
    # Direct integration in measured frequency range
    w0 = f[nfrq] * 2 * np.pi
    m0 = (w0 * 10 ** dspec[nfrq]) ** 2
    Es0 = trapezoid(m0, f[nfrq])
    
    # Extrapolate below minimum frequency
    df = 1e-3
    fmin = 1e-10
    fmax = 2000
    
    f1 = np.arange(fmin, f[nfrq[0]], df)
    u1 = fitmethod(p, f1, im)
    w1 = f1 * 2 * np.pi
    m1 = (w1 * 10 ** u1) ** 2
    
    # Extrapolate above maximum frequency
    f2 = np.arange(f[nfrq[-1]] + df, fmax + df, df)
    u2 = fitmethod(p, f2, im)
    w2 = f2 * 2 * np.pi
    m2 = (w2 * 10 ** u2) ** 2
    
    # Full frequency range
    f3 = np.arange(fmin, fmax + df, df)
    u3 = fitmethod(p, f3, im)
    w3 = f3 * 2 * np.pi
    m3 = (w3 * 10 ** u3) ** 2
    
    Es1 = trapezoid(m1, f1)
    Es2 = trapezoid(m2, f2)
    Es3 = trapezoid(m3, f3)
    
    fex = np.concatenate([f1, f[nfrq], f2])
    dspecex = np.concatenate([u1, dspec[nfrq], u2])
    dspecfit = fitmethod(p, fex, im)
    
    Es = Es0 + Es1 + Es2
    
    Esx = np.array([Es, Es3, Es0, Es1, Es2])
    # total energy, modeled energy, measured energy, extrapolated low, extrapolated high
    
    return Esx, fex, dspecex, dspecfit
