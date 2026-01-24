clear;
clc;
addpath('../src/')
addpath(genpath('../src/sacstuff/'));

%% get event and station catalogs
catfile = 'catalog';
stafile = 'station.txt';
sacdir_in = '../data/sac';
evdir = '../data/matfile';
sactype = 'scec';
mode = 'new';

if(~exist(evdir))
    mkdir(evdir);
end

eqinfo = update_eqinfo(catfile,8);
save([evdir,'/eqinfo.mat'],'eqinfo');

fid = fopen(stafile,'r');
for i = 1:(2^31-1)
    if (feof(fid))
        break;
    end
    a = fgetl(fid);
    xx = textscan(a,'%s %f %f %f');
    sta.stid{i} = xx{1}{1};
    sta.slat(i) = xx{2};
    sta.slon(i) = xx{3};
    sta.selev(i) = xx{4};
end
fclose(fid);

neq = length(eqinfo.id);
if (sacdir_in(end) ~= '/')
    sacdir_in(end+1) = '/';
end
if (evdir(end) ~= '/')
    evdir(end+1) = '/';
end

%% get input parameters
n1 = find(eqinfo.qyr<=9999);

for ii = 1:length(eqinfo.id)
    evid  = eqinfo.id(ii);
    idstr = num2str(evid);
    year  = eqinfo.qyr(ii);
    month = eqinfo.qmon(ii);
    day   = eqinfo.qdy(ii);

    % Organize your own sacdir and outevdir to fit in your database structure
    sacdir = sacdir_in;
    outevdir = sprintf('%s%4.4d%s%2.2d', evdir, year, '/', month);

    disp(['Working on event:',num2str(evid)]);
    data = sac2mat(sacdir, outevdir, evid, eqinfo, sta);
end

