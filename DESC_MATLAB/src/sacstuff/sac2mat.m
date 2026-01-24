function data = sac2mat(sacdir_in, evdir, evid, eqinfo, sta)
% input parameters:
% 	sacdir_in: sac directory
%	evdir: output event file name.
%	evid: event ID
%	eqinfo: catalog
%	sta: station list
% output: data file

if (sacdir_in(end) ~= '/')
    sacdir_in(end+1) = '/';
end

if (evdir(end) ~= '/')
    evdir(end+1) = '/';
end

idstr = num2str(evid);
sacdir = [sacdir_in, idstr];
evfile = [evdir 'ev' num2str(evid) '.mat'];

if (~exist(sacdir,'dir'))
    data = [];
    return;
end

if (~exist(evdir,'dir'))
    mkdir(evdir);
end

%% first set event information in the datafile.

ie = find(eqinfo.id == evid);
if (isempty(ie))
    disp(['target event not in the catalog',num2str(evid)]);
end

fields = fieldnames(eqinfo);
for i = 1:length(fields)
    data.(fields{i}) = eqinfo.(fields{i})(ie);
end
data.cuspid = eqinfo.id(ie);
data.qmag1 = eqinfo.mb(ie);
data.id0 = idstr;

% event information updated.

%% get sacfile list
saclist = dir(sacdir);
saclist = saclist(3:end);

nsac = length(saclist);         % sac file number

numts = 0;

for isac = 1:nsac
    [~,Zdata,ZSAChdr] = fget_sac([sacdir,'/',saclist(isac).name]);
    [syr,smon,sday] = datevec(datenum(ZSAChdr.event.nzyear,1,1) + ZSAChdr.event.nzjday-1);
    shr = ZSAChdr.event.nzhour;
    smn = ZSAChdr.event.nzmin;
    ssc = ZSAChdr.event.nzsec + ZSAChdr.event.nzmsec/1000 + ZSAChdr.times.b;
    stime = [syr, smon, sday, shr, smn, ssc];

    xx = isstrprop(ZSAChdr.station.kstnm, 'alphanum');
    nbad = find(xx == 0);
    nstop = nbad(1) - 1;
    staname = ZSAChdr.station.kstnm(1:nstop);
    netname = ZSAChdr.stations.knetwk(1:2);
    chname  = ZSAChdr.stations.kcmpnm(1:3);
    stid = [netname,staname];
    npts = ZSAChdr.data.trcLen;
    numts = numts + 1;
    data.net{numts} = netname;
    data.stime(numts) = date2epoch(stime);
    data.syr(numts) = syr;
    data.smon(numts) = smon;
    data.sdy(numts) = sday;
    data.shr(numts) = shr;
    data.smn(numts) = smn;
    data.ssc(numts) = ssc;
    data.sdoy(numts) = ZSAChdr.event.nzjday;
    data.samprate(numts) = 1 ./ ZSAChdr.times.delta;
    data.dt(numts) = ZSAChdr.times.delta;
    data.pick1(numts) = ZSAChdr.times.a;    % Change according to your SAC header
    data.pick2(numts) = ZSAChdr.times.t0;   % Change according to your SAC header
    data.sta{numts} = staname;
    data.stid{numts} = [netname, staname];
    nfsta = find(ismember(sta.stid, stid));
    if (isempty(nfsta))
        disp([netname ':' staname]);
        data.stlat(numts) = ZSAChdr.station.stla;
        data.stlon(numts) = ZSAChdr.station.stlo;
        data.stelev(numts) = ZSAChdr.station.stel / 1000;
    else
        %if (length(nfsta)>1)
        %	disp([netname,':',staname]);
        %	sta.sta(nfsta)
        %end
        data.stlat(numts) = sta.slat(nfsta(1));
        data.stlon(numts) = sta.slon(nfsta(1));
        data.stelev(numts) = sta.selev(nfsta(1));
    end
    data.chan{numts} = chname;
    data.comp(numts) = 0;
    if (chname(3:3)=='E' || chname(3:3)=='X' || chname(3:3)=='2')
        data.comp(numts) = 3;
    elseif (chname(3:3)=='N' || chname(3:3)=='Y' || chname(3:3)=='1')
        data.comp(numts) = 2;
    elseif (chname(3:3)=='Z')
        data.comp(numts) = 1;
    end
    if (strcmpi(netname,'BP') || strcmpi(netname,'SF'))
        data.comp(numts) = str2num(chname(3:3));
    end
    if (chname(2:2)=='N' || chname(2:2)=='Y')
        data.acc(numts) = 1;
        tracea = Zdata - mean(Zdata);
        traces = cumsum(tracea) / data.samprate(numts);
    else
        traces = Zdata - mean(Zdata);
        tracea = zeros(size(Zdata));
        tracea(1:length(Zdata) - 1) = diff(traces) * data.samprate(numts);
    end

    [dist,azi] = distance(data.qlat, data.qlon, data.stlat(numts), data.stlon(numts));
    dist = dist * 111.19;

    data.del(numts) = dist;
    data.sazi(numts) = azi;
    data.pazi(numts) = azi + 180;

    data.seis(numts, 1:npts) = traces;
    data.aseis(numts, 1:npts) = tracea;
    data.npts(numts) = npts;
end
data.numts = numts;

save(evfile,'data','-v7.3');

