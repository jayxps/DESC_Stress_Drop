function yeartime = epoch2year(epochtime)

for i = 1:length(epochtime)
    daten = epoch2date(epochtime(i));
    %yearn = year(daten);
    datev = datevec(daten);
    yearn = datev(1);
    datev(1) = 0;
    yeartime(i) = yearn + datenum(datev)/yeardays(yearn);
end