import numpy as np
from datetime import datetime, timezone

def date2epoch(datevector):
    epochs = np.zeros(len(datevector))
    for i, date in enumerate(datevector):
        year = int(date[0])
        month = int(date[1])
        day = int(date[2])
        hour = int(date[3])
        minute = int(date[4])
        second = date[5]
        time = datetime(year, month, day, hour, minute, 0, tzinfo=timezone.utc)
        epochs[i] = time.timestamp() + second
        
    return epochs