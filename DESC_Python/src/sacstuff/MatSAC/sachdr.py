import numpy as np

def _get(sac, name, default=None):
    return getattr(sac, name, default)

def sachdr_from_sac(sac):
    """
    Build a Python dict resembling MATLAB `sachdr.m` structure from an ObsPy SAC header.

    Returns keys: 'times', 'station', 'event', 'user', 'descrip', 'evsta', 'llnl', 'response', 'data'.
    Missing fields are filled with None or -12345 where appropriate.
    """
    times = {
        'delta': _get(sac, 'delta', np.nan),
        'b': _get(sac, 'b', 0.0),
        'e': _get(sac, 'e', None),
        'o': _get(sac, 'o', None),
        'a': _get(sac, 'a', None),
        't0': _get(sac, 't0', None), 't1': _get(sac, 't1', None), 't2': _get(sac, 't2', None), 't3': _get(sac, 't3', None), 't4': _get(sac, 't4', None),
        't5': _get(sac, 't5', None), 't6': _get(sac, 't6', None), 't7': _get(sac, 't7', None), 't8': _get(sac, 't8', None), 't9': _get(sac, 't9', None),
        'k0': _get(sac, 'k0', ''), 'ka': _get(sac, 'ka', ''), 'kt0': _get(sac, 'kt0', ''), 'kt1': _get(sac, 'kt1', ''), 'kt2': _get(sac, 'kt2', ''),
        'kt3': _get(sac, 'kt3', ''), 'kt4': _get(sac, 'kt4', ''), 'kt5': _get(sac, 'kt5', ''), 'kt6': _get(sac, 'kt6', ''), 'kt7': _get(sac, 'kt7', ''),
        'kt8': _get(sac, 'kt8', ''), 'kt9': _get(sac, 'kt9', ''), 'kf': _get(sac, 'kf', ''),
        'npts': _get(sac, 'npts', None),
    }

    station = {
        'stla': _get(sac, 'stla', None),
        'stlo': _get(sac, 'stlo', None),
        'stel': _get(sac, 'stel', None),
        'stdp': _get(sac, 'stdp', None),
        'cmpaz': _get(sac, 'cmpaz', None),
        'cmpinc': _get(sac, 'cmpinc', None),
        'kstnm': _get(sac, 'kstnm', ''),
        'kcmpnm': _get(sac, 'kcmpnm', ''),
        'knetwk': _get(sac, 'knetwk', ''),
    }

    event = {
        'evla': _get(sac, 'evla', None), 'evlo': _get(sac, 'evlo', None), 'evel': _get(sac, 'evel', None), 'evdp': _get(sac, 'evdp', None),
        'nzyear': _get(sac, 'nzyear', None), 'nzjday': _get(sac, 'nzjday', None), 'nzhour': _get(sac, 'nzhour', None), 'nzmin': _get(sac, 'nzmin', None),
        'nzsec': _get(sac, 'nzsec', None), 'nzmsec': _get(sac, 'nzmsec', None), 'kevnm': _get(sac, 'kevnm', ''), 'mag': _get(sac, 'mag', None),
        'imagtyp': None, 'imagsrc': None,
    }

    user = {
        'data': [
            _get(sac, 'user0', None), _get(sac, 'user1', None), _get(sac, 'user2', None), _get(sac, 'user3', None), _get(sac, 'user4', None),
            _get(sac, 'user5', None), _get(sac, 'user6', None), _get(sac, 'user7', None), _get(sac, 'user8', None), _get(sac, 'user9', None),
        ],
        'label': [_get(sac, 'kuser0', ''), _get(sac, 'kuser1', ''), _get(sac, 'kuser2', '')],
    }

    descrip = {
        'iftype': _get(sac, 'iftype', None), 'idep': _get(sac, 'idep', None), 'iztype': _get(sac, 'iztype', None), 'iinst': _get(sac, 'iinst', None),
        'istreg': _get(sac, 'istreg', None), 'ievreg': _get(sac, 'ievreg', None), 'ievtyp': _get(sac, 'ievtyp', None), 'iqual': _get(sac, 'iqual', None),
        'isynth': _get(sac, 'isynth', None),
    }

    evsta = {
        'dist': _get(sac, 'dist', None), 'az': _get(sac, 'az', None), 'baz': _get(sac, 'baz', None), 'gcarc': _get(sac, 'gcarc', None),
    }

    llnl = {
        'xminimum': _get(sac, 'xminimum', None), 'xmaximum': _get(sac, 'xmaximum', None), 'yminimum': _get(sac, 'yminimum', None), 'ymaximum': _get(sac, 'ymaximum', None),
        'norid': None, 'nevid': None, 'nxsize': None, 'nysize': None,
    }

    response = [
        _get(sac, 'resp0', None), _get(sac, 'resp1', None), _get(sac, 'resp2', None), _get(sac, 'resp3', None), _get(sac, 'resp4', None),
        _get(sac, 'resp5', None), _get(sac, 'resp6', None), _get(sac, 'resp7', None), _get(sac, 'resp8', None), _get(sac, 'resp9', None),
    ]

    data = {
        'trcLen': _get(sac, 'npts', None), 'scale': _get(sac, 'scale', None),
    }

    return {
        'times': times,
        'station': station,
        'event': event,
        'user': user,
        'descrip': descrip,
        'evsta': evsta,
        'llnl': llnl,
        'response': response,
        'data': data,
    }
