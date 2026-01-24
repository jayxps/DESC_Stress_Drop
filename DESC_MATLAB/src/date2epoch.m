function tepoch = date2epoch(datevector)

t = datenum(datevector);
t0 = datenum([1970 1 1 0 0 0]);

tepoch = ((t-t0)*86400);
