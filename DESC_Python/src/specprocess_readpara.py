def specprocess_readpara(filename):
    """
    Read parameter file and return parameter dictionary P.

    The file format mirrors the MATLAB version; lines starting with '*' are comments.
    Parameters are expected in order, one per non-comment line.
    """
    P = {}
    with open(filename, 'r') as fid:
        npara = 0
        for line in fid:
            line = line.strip('\n')
            if not line or line.startswith('*'):
                continue
            npara += 1
            if npara == 1:  # event directory
                P['evdir'] = line
            elif npara == 2:  # target phase
                P['target'] = line
            elif npara == 3:  # maximum distance
                x = list(map(float, line.split()))
                P['delmax'] = x[0]
            elif npara == 4:  # magnitude range
                x = list(map(float, line.split()))
                P['minmag'], P['maxmag'] = x[0], x[1]
            elif npara == 5:  # three frequency bands
                x = list(map(float, line.split()))
                P['freq1'] = [x[0], x[1]]
                P['freq2'] = [x[2], x[3]]
                P['freq3'] = [x[4], x[5]]
            elif npara == 6:  # STN for the three bins
                x = list(map(float, line.split()))
                P['stnmin'] = [float(__import__('math').log10(v)) for v in x[:3]]
            elif npara == 7:  # fraction of events to check
                x = list(map(float, line.split()))
                P['randfrac'] = x[0]
            elif npara == 8:  # area limits
                x = list(map(float, line.split()))
                P['qlat1'], P['qlat2'], P['qlon1'], P['qlon2'] = x[0], x[1], x[2], x[3]
            elif npara == 9:  # error type
                x = list(map(float, line.split()))
                P['itype'] = int(x[0])
            elif npara == 10:  # number of iterations
                x = list(map(float, line.split()))
                P['nitmax'] = int(x[0])
            elif npara == 11:  # time interval for travel-time bin stack
                x = list(map(float, line.split()))
                P['dtt'] = x[0]
            elif npara == 12:  # tstar frequency range
                x = list(map(float, line.split()))
                P['fmin_tstar'], P['fmax_tstar'] = x[0], x[1]
            elif npara == 13:  # minimum average stn for tstar
                x = list(map(float, line.split()))
                import math
                P['minavgstn'] = math.log10(x[0])
            elif npara == 14:  # target network(s)
                P['targnet'] = line.split()
            elif npara == 15:  # target station(s)
                P['targsta'] = line.split()
            elif npara == 16:  # target channel(s)
                P['targchan'] = line.split()
            elif npara == 17:  # eqinfo path
                P['eqinfo'] = line
    # Default eqinfo if not provided explicitly
    if 'eqinfo' not in P and 'evdir' in P:
        P['eqinfo'] = f"{P['evdir']}/eqinfo.mat"
    return P
