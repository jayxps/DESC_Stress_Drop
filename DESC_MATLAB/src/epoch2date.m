function qdate = epoch2date(qtime)

t0 = datenum([1970 1 1 0 0 0]);
%for i = 1:nlen
%timestring=epoch2str(qtime(i),'%Y/%m/%d,%H:%M:%S.%s');
%timex=sscanf(timestring,'%d/%d/%d,%d:%d:%f');
timex = datevec(t0+qtime/86400);
qdate = datenum(timex);
%end
