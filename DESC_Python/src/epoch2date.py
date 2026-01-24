import numpy as np
from datetime import datetime, timezone

def epoch2date(epoch_times):
    # Determine if epoch_times is a single value or an array
    if np.isscalar(epoch_times):
        epoch_times = np.array([epoch_times])
        
    datevectors = np.zeros((len(epoch_times), 6))
    for i, epoch in enumerate(epoch_times):
        msecond = epoch - int(epoch)
        dt = datetime.fromtimestamp(int(epoch), tz=timezone.utc)
        datevectors[i, 0] = dt.year
        datevectors[i, 1] = dt.month
        datevectors[i, 2] = dt.day
        datevectors[i, 3] = dt.hour
        datevectors[i, 4] = dt.minute
        datevectors[i, 5] = dt.second + msecond
    return datevectors