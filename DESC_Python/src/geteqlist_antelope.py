import numpy as np


def _loadtxt_2d(filename):
    """Load numeric text and ensure a 2-D shape even for single-line files."""
    mdata = np.loadtxt(filename)
    if mdata.ndim == 1:
        mdata = mdata.reshape(1, -1)
    return mdata


def geteqlist_antelope(filename, itype):
    """
    Read earthquake catalog in various formats.
    
    Parameters
    ----------
    filename : str
        Path to catalog file.
    itype : int
        Format type:
        1 = SCEC, 2 = HK, 3 = EFS, 4 = NCSN, 5 = JMA, 
        6 = hypoDD, 7 = ANSS, 8 = XYZ
    
    Returns
    -------
    eqlist : ndarray
        Array with columns: [year, month, day, hour, minute, second, 
                             lat, lon, depth, mag, evid]
    """
    if itype == 1:
        eqlist = _get_scec(filename)
    elif itype == 2:
        eqlist = _get_hk(filename)
    elif itype == 3:
        eqlist = _get_efs(filename)
    elif itype == 4:
        eqlist = _get_ncsn(filename)
    elif itype == 5:
        eqlist = _get_jma(filename)
    elif itype == 6:
        eqlist = _get_dd(filename)
    elif itype == 7:
        eqlist = _get_anss(filename)
    elif itype == 8:
        eqlist = _get_xyz(filename)
    else:
        raise ValueError(f"Unknown catalog type: {itype}")
    
    return eqlist


def _get_jma(filename):
    """JMA format: load numeric data from file."""
    mData = _loadtxt_2d(filename)
    M, N = mData.shape
    uOutput = np.zeros((M, 11))
    
    if N == 10:
        uOutput[:, 10] = np.arange(1, M + 1)
    elif N >= 11:
        uOutput[:, 10] = mData[:, 10]
    
    uOutput[:, :10] = mData[:, :10]
    return uOutput


def _get_xyz(filename):
    """XYZ format: columns rearranged [yr mon day hr min sec, mag, lat lon dep, evid]."""
    mData = _loadtxt_2d(filename)
    M, N = mData.shape
    uOutput = np.zeros((M, 11))
    
    if N == 10:
        uOutput[:, 10] = np.arange(1, M + 1)
    elif N >= 11:
        uOutput[:, 10] = mData[:, 10]
    
    uOutput[:, :6] = mData[:, :6]
    uOutput[:, 9] = mData[:, 6]      # mag
    uOutput[:, 6:9] = mData[:, 7:10] # lat, lon, dep
    return uOutput


def _get_scec(filename):
    """SCEC format: fixed-width text parser."""
    uOutput = []
    
    with open(filename, 'r') as f:
        for line_no, line in enumerate(f, 1):
            line = line.replace(' ', '0')
            try:
                row = _read_values_scec(line)
                if not np.isnan(row[0]):
                    uOutput.append(row)
            except Exception as e:
                print(f"Import: Problem in line {line_no} of {filename}. Line ignored.")
                continue
    
    uOutput = np.array(uOutput)
    # Filter invalid latitudes
    valid = uOutput[:, 6] <= 100
    return uOutput[valid]


def _read_values_scec(line):
    """Parse SCEC fixed-width format."""
    row = np.zeros(11)
    row[0] = float(line[0:4])       # year
    row[1] = float(line[5:7])       # month
    row[2] = float(line[8:10])      # day
    row[3] = float(line[11:13])     # hour
    row[4] = float(line[14:16])     # minute
    row[5] = float(line[17:23])     # second
    row[6] = float(line[35:43])     # latitude
    row[7] = float(line[43:51])     # longitude
    row[8] = float(line[53:58])     # depth
    row[9] = float(line[26:30])     # magnitude
    row[10] = float(line[59:67])    # event ID
    return row


def _get_efs(filename):
    """EFS format: fixed-width text parser."""
    uOutput = []
    
    with open(filename, 'r') as f:
        for line_no, line in enumerate(f, 1):
            line = line.replace(' ', '0')
            try:
                row = _read_values_efs(line)
                if not np.isnan(row[0]):
                    uOutput.append(row)
            except Exception as e:
                print(f"Import: Problem in line {line_no} of {filename}. Line ignored.")
                continue
    
    uOutput = np.array(uOutput)
    valid = uOutput[:, 6] <= 100
    return uOutput[valid]


def _read_values_efs(line):
    """Parse EFS fixed-width format."""
    row = np.zeros(11)
    row[0] = float(line[9:13])      # year
    row[1] = float(line[14:16])     # month
    row[2] = float(line[17:19])     # day
    row[3] = float(line[20:22])     # hour
    row[4] = float(line[23:25])     # minute
    row[5] = float(line[26:32])     # second
    row[6] = float(line[45:53])     # latitude
    row[7] = float(line[54:64])     # longitude
    row[8] = float(line[67:72])     # depth
    row[9] = float(line[36:41])     # magnitude
    row[10] = float(line[0:8])      # scecid
    return row


def _get_hk(filename):
    """HK format: fixed-width text parser."""
    uOutput = []
    
    with open(filename, 'r') as f:
        for line_no, line in enumerate(f, 1):
            line = line.replace(' ', '0')
            try:
                row = _read_values_hk(line)
                if not np.isnan(row[0]):
                    uOutput.append(row)
            except Exception as e:
                print(f"Import: Problem in line {line_no} of {filename}. Line ignored.")
                continue
    
    return np.array(uOutput)


def _read_values_hk(line):
    """Parse HK fixed-width format."""
    row = np.zeros(11)
    row[0] = float(line[0:4])       # year
    row[1] = float(line[5:7])       # month
    row[2] = float(line[8:10])      # day
    row[3] = float(line[11:13])     # hour
    row[4] = float(line[14:16])     # minute
    row[5] = float(line[17:23])     # second
    row[6] = float(line[34:42])     # latitude
    row[7] = float(line[43:53])     # longitude
    row[8] = float(line[55:60])     # depth
    row[9] = float(line[62:68])     # magnitude
    row[10] = float(line[24:33])    # cuspid
    return row


def _get_ncsn(filename):
    """NCSN format: fixed-width text parser."""
    uOutput = []
    
    with open(filename, 'r') as f:
        for line_no, line in enumerate(f, 1):
            line = line.replace(' ', '0')
            try:
                row = _read_values_ncsn(line)
                if not np.isnan(row[0]):
                    uOutput.append(row)
            except Exception as e:
                print(f"Import: Problem in line {line_no} of {filename}. Line ignored.")
                continue
    
    uOutput = np.array(uOutput)
    valid = uOutput[:, 6] <= 100
    return uOutput[valid]


def _read_values_ncsn(line):
    """Parse NCSN fixed-width format."""
    row = np.zeros(11)
    row[0] = float(line[0:4])       # year
    row[1] = float(line[5:7])       # month
    row[2] = float(line[8:10])      # day
    row[3] = float(line[11:13])     # hour
    row[4] = float(line[14:16])     # minute
    row[5] = float(line[17:22])     # second
    row[6] = float(line[24:33])     # latitude
    row[7] = float(line[33:44])     # longitude
    row[8] = float(line[45:52])     # depth
    row[9] = float(line[53:58])     # magnitude
    row[10] = float(line[89:97])    # ncsn id
    return row


def _get_dd(filename):
    """hypoDD format: load and rearrange columns."""
    mData = _loadtxt_2d(filename)
    uOutput = np.zeros((mData.shape[0], 11))
    uOutput[:, 0] = mData[:, 10]    # year
    uOutput[:, 1] = mData[:, 11]    # month
    uOutput[:, 2] = mData[:, 12]    # day
    uOutput[:, 3] = mData[:, 13]    # hour
    uOutput[:, 4] = mData[:, 14]    # minute
    uOutput[:, 5] = mData[:, 15]    # second
    uOutput[:, 6] = mData[:, 1]     # lat
    uOutput[:, 7] = mData[:, 2]     # lon
    uOutput[:, 8] = mData[:, 3]     # depth
    uOutput[:, 9] = mData[:, 16]    # mag
    uOutput[:, 10] = mData[:, 0]    # id
    return uOutput


def _get_anss(filename):
    """ANSS format: fixed-width text parser."""
    uOutput = []
    
    with open(filename, 'r') as f:
        for line_no, line in enumerate(f, 1):
            line = line.replace(' ', '0')
            try:
                row = _read_values_anss(line)
                if not np.isnan(row[0]):
                    uOutput.append(row)
            except Exception as e:
                print(f"Import: Problem in line {line_no} of {filename}. Line ignored.")
                continue
    
    uOutput = np.array(uOutput)
    valid = uOutput[:, 6] <= 100
    return uOutput[valid]


def _read_values_anss(line):
    """Parse ANSS fixed-width format."""
    row = np.zeros(11)
    row[0] = float(line[0:4])       # year
    row[1] = float(line[5:7])       # month
    row[2] = float(line[8:10])      # day
    row[3] = float(line[11:13])     # hour
    row[4] = float(line[14:16])     # minute
    row[5] = float(line[17:22])     # second
    row[6] = float(line[22:31])     # latitude
    row[7] = float(line[31:41])     # longitude
    row[8] = float(line[41:48])     # depth
    row[9] = float(line[48:54])     # magnitude
    row[10] = float(line[83:96])    # ncsn id
    return row
