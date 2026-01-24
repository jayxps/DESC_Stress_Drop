function P = specprocess_readpara(filename)

fid = fopen(filename,'r');

npara = 0;
while (1 == 1)
    a = fgetl(fid);
    if (a == -1) break; end
    if (a(1:1) ~= '*')
        npara = npara+1;
        if (npara == 1)    % this is for event directory
            P.evdir = a;
        elseif (npara==2)  % this is for target Phase.
            P.target = a;
        elseif (npara==3)
            x = sscanf(a,'%f'); % this is for maximum distance
            P.delmax = x(1);
        elseif (npara == 4)
            x = sscanf(a,'%f'); % this is for magnitude range
            P.minmag = x(1);
            P.maxmag = x(2);
        elseif (npara == 5)
            x = sscanf(a,'%f'); % this is for the three frequency band
            P.freq1(1) = x(1);
            P.freq1(2) = x(2);
            P.freq2(1) = x(3);
            P.freq2(2) = x(4);
            P.freq3(1) = x(5);
            P.freq3(2) = x(6);
        elseif (npara == 6)
            x = sscanf(a,'%f'); % this is for STN for the three bins
            P.stnmin = log10(x(1:3));
        elseif (npara == 7)
            x = sscanf(a,'%f'); % this is for fraction of events to check, 1 for all
            P.randfrac = x(1);
        elseif (npara == 8)
            x = sscanf(a,'%f'); % this is for limits of area
            P.qlat1 = x(1);
            P.qlat2 = x(2);
            P.qlon1 = x(3);
            P.qlon2 = x(4);
        elseif (npara == 9)
            x = sscanf(a,'%f'); % this is for error type.
            P.itype = x(1);
        elseif (npara == 10)
            x = sscanf(a,'%f'); % this is for number of iterations
            P.nitmax = x(1);
        elseif (npara == 11)
            x = sscanf(a,'%f');  % this is for time interval for travel-time bin stack
            P.dtt = x(1);
        elseif (npara == 12)
            x = sscanf(a,'%f'); % this is for tstar frequency range
            P.fmin_tstar = x(1);
            P.fmax_tstar = x(2);
        elseif (npara == 13)
            x = sscanf(a,'%f'); % this is for minium average stn for tstar
            P.minavgstn = log10(x(1));
        elseif (npara == 14)  % * means all networks
            x = strsplit(a,' ');
            P.targnet = x;
        elseif (npara == 15)  % * means all stations
            x = strsplit(a,' ');
            P.targsta = x;
        elseif (npara == 16)
            x = strsplit(a,' ');
            P.targchan = x;
        elseif (npara == 17)
            P.eqinfo = a;
        end
    end
end

if (npara == 15)
    P.eqinfo = [P.evdir,'/eqinfo.mat'];
end

