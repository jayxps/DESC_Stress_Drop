import sys
import os
import numpy as np
from datetime import datetime

# Add parent to path for imports
sys.path.insert(0, os.path.dirname(__file__))
from geteqlist_antelope import geteqlist_antelope
from date2epoch import date2epoch


def update_eqinfo(catfile, itype):
    """
    Generate eqinfo dict from raw catalog file.
    
    Parameters
    ----------
    catfile : str
        Path to catalog file.
    itype : int
        Catalog format type (see geteqlist_antelope for details).
    
    Returns
    -------
    eqinfo : dict
        Earthquake catalog with fields:
        - id, qlat, qlon, qdep, mb, qtime (epoch seconds)
        - qyr, qmon, qdy, qdoy, qhr, qmn, qsc
        - index (initialized to 0)
    """
    eqlist = geteqlist_antelope(catfile, itype)
    neq = len(eqlist)
    
    # Initialize arrays
    eqinfo = {
        'id': np.zeros(neq, dtype=int),
        'qlat': np.zeros(neq),
        'qlon': np.zeros(neq),
        'qdep': np.zeros(neq),
        'mb': np.zeros(neq),
        'qtime': np.zeros(neq),
        'qyr': np.zeros(neq, dtype=int),
        'qmon': np.zeros(neq, dtype=int),
        'qdy': np.zeros(neq, dtype=int),
        'qdoy': np.zeros(neq, dtype=int),
        'qhr': np.zeros(neq, dtype=int),
        'qmn': np.zeros(neq, dtype=int),
        'qsc': np.zeros(neq),
        'index': np.zeros(neq, dtype=int),
    }
    
    # eqlist columns: [year, month, day, hour, minute, second, lat, lon, depth, mag, evid]
    datevectors = eqlist[:, :6]  # year, month, day, hour, minute, second
    qtime = date2epoch(datevectors)
    
    for ie in range(neq):
        # Day of year based on calendar date
        date_obj = datetime(int(eqlist[ie, 0]), int(eqlist[ie, 1]), int(eqlist[ie, 2]))
        jan1 = datetime(int(eqlist[ie, 0]), 1, 1)
        doy = (date_obj - jan1).days + 1
        
        eqinfo['id'][ie] = int(eqlist[ie, 10])
        eqinfo['qlat'][ie] = eqlist[ie, 6]
        eqinfo['qlon'][ie] = eqlist[ie, 7]
        eqinfo['qdep'][ie] = eqlist[ie, 8]
        eqinfo['mb'][ie] = eqlist[ie, 9]
        eqinfo['qtime'][ie] = qtime[ie]
        eqinfo['qyr'][ie] = int(eqlist[ie, 0])
        eqinfo['qmon'][ie] = int(eqlist[ie, 1])
        eqinfo['qdy'][ie] = int(eqlist[ie, 2])
        eqinfo['qdoy'][ie] = doy
        eqinfo['qhr'][ie] = int(eqlist[ie, 3])
        eqinfo['qmn'][ie] = int(eqlist[ie, 4])
        eqinfo['qsc'][ie] = eqlist[ie, 5]
    
    return eqinfo
