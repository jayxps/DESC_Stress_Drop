% this is to recompute spectrum for data files.
function [] = update_data_spec(evfile,disp_flag)

if(nargin == 1)
    disp_flag = 1; % when disp_flag=1, the function will show matfile path
end

load(evfile);
if (~exist('data')); return; end

if (data.numts == 0)
    return;
end

if(disp_flag == 1)
    disp(strcat('Working on',evfile));
end
Nfft = 1024;
Nfftc = 4096;

data.pspec = zeros(data.numts, Nfft/2+1);
data.sspec = zeros(data.numts, Nfft/2+1);
data.cspec = zeros(data.numts, Nfftc/2+1);
data.pnspec = zeros(data.numts, Nfft/2+1);
data.snspec = zeros(data.numts, Nfft/2+1);
data.cnspec = zeros(data.numts, Nfft/2+1);
data.pfreq = zeros(data.numts, Nfft/2+1);
data.sfreq = zeros(data.numts, Nfft/2+1);
data.cfreq = zeros(data.numts, Nfft/2+1);

data.np1 = zeros(1, data.numts);
data.tp1 = zeros(1, data.numts);
data.ns1 = zeros(1, data.numts);
data.ts1 = zeros(1, data.numts);
data.codawin = zeros(1, data.numts);
data.tswin = zeros(1, data.numts);
data.tpwin = zeros(1, data.numts);
data.pwin_mode = cell(1, data.numts);
data.pnwin = zeros(1, data.numts);
data.snwin = zeros(1, data.numts);
data.cnwin = zeros(1, data.numts);
data.pstn = zeros(1, data.numts);
data.sstn = zeros(1, data.numts);

samprate = data.samprate(1);
N_freq_samples = 500;
if (round(samprate) > 100)
    N_freq_samples = 1000;
end
data.flog = logspace(log10(0.5), log10(samprate * 0.45), N_freq_samples);

data.pspec = zeros(data.numts, Nfft / 2 + 1);
data.sspec = zeros(data.numts, Nfft / 2 + 1);
data.cspec = zeros(data.numts, Nfftc / 2 + 1);
data.pnspec = zeros(data.numts, Nfft / 2 + 1);
data.snspec = zeros(data.numts, Nfft / 2 + 1);
data.cnspec = zeros(data.numts, Nfftc / 2 + 1);
data.pfreq = zeros(data.numts, Nfft / 2 + 1);
data.sfreq = zeros(data.numts, Nfft / 2 + 1);
data.cfreq = zeros(data.numts, Nfftc / 2 + 1);

%%% Only for local-time cases, adjust UTC to Beijing Time
data.stime = data.stime + 8*3600;

for ns = 1:data.numts
    ppick = data.pick1(ns);
    spick = data.pick2(ns);

    time1 = data.stime(ns);
    samprate = data.samprate(ns);

    if (ppick > 0)
        data.np1(ns) = floor((ppick + data.qtime - time1)*samprate);
        data.tp1(ns) = data.np1(ns) / samprate;
    else
        data.np1(ns) = -12345;
        data.tp1(ns) = -12345;
    end


    if (spick > 0)
        data.ns1(ns) = floor((spick + data.qtime - time1)*samprate);
        data.ts1(ns) = data.ns1(ns)/samprate;
    else
        data.ns1(ns) = -12345;
        data.ts1(ns) = -12345;
    end

    % Define windows for spectra calculation
    data.codawin(ns) = 10;
    try
        qmag = data.qmag1;
    catch
        qmag = data.mb;
    end
    %%% Define your own S window lengths for different magnitudes
    if (qmag >= 4.0)
        data.tswin(ns) = 10.0;
    elseif (qmag >= 2.5)
        data.tswin(ns) = 5.0;
    else
        data.tswin(ns) = 3.0;
    end
    %%% Define your own P window lengths according to your waveforms
    data.tpwin(ns) = 2.0;
    % if (qmag >= 4.0)
    %     data.tpwin(ns) = 2.0;
    % elseif (qmag >= 2.5)
    %     data.tpwin(ns) = 2.0;
    % else
    %     data.tpwin(ns) = 2.0;
    % end
    if (ppick < 0)
        data.tpwin(ns) = -12345;
        data.pwin_mode{ns} = '';
    elseif (spick < 0)
        data.pwin_mode{ns} = 'fixed';
    else
        data.pwin_mode{ns} = 'psdiff';
    end
    if(strcmp(data.pwin_mode{ns},'psdiff'))     % Use S-P as P window
        data.tpwin(ns) = min(data.tpwin(ns), spick - ppick - 0.1);
    end


    %########### now will compute P and S signal to noise ratio with 10-40 bandpass
    data.pstn(ns) = 0;
    data.sstn(ns) = 0;

    if (data.np1(ns) > 0 && data.tpwin(ns) >= 0.5)
        %### This is for P wave.
        np1 = data.np1(ns) - round(0.1 * data.samprate(ns)); %allow 0.1 second before
        if (np1 < 1)
            np1 = 1;
        end
        nwin = round(data.tpwin(ns) * data.samprate(ns));
        nend = np1 + nwin - 1;
        if (nend > data.npts(ns))
            nend = data.npts(ns);
            nwin = data.npts(ns) - np1 + 1;
        end
        data.pnwin(ns) = nwin;
        data.pseis(ns,1:nwin) = data.seis(ns,np1:nend);
        if (np1 - nwin < 1)
            data.pnoise(ns,1:nwin) = zeros(1,nwin);
            data.pnoise(ns,1:np1-1) = data.seis(ns,1:np1-1);
        else
            data.pnoise(ns,1:nwin) = data.seis(ns,np1-nwin:np1-1);
        end
        if (nwin >= 50)
            nx = Nfft;
            [spec,~] = pmtm(data.pseis(ns,1:nwin),3,nx,data.samprate(ns));
            [nspec,f] = pmtm(data.pnoise(ns,1:nwin),3,nx,data.samprate(ns));
            omega = 2 * pi * f;
            nfg = find(f > 0);
            spec = sqrt(spec);
            nspec = sqrt(nspec);
            data.pspec(ns,1) = spec(1);
            data.pnspec(ns,1) = nspec(1);
            data.pspec(ns,nfg) = spec(nfg) ./ omega(nfg);
            data.pnspec(ns,nfg) = nspec(nfg) ./ omega(nfg);
            data.pfreq(ns,nfg) = f(nfg);
            % nband = find(f(nfg)>=4 & f(nfg)<=10);
            % data.pstn(ns) = mean(log10(data.pspec(ns,nband))-log10(data.pnspec(ns,nband)));
            specres = interp1(f(nfg),data.pspec(ns,nfg),data.flog,'linear','extrap');
            data.pspecres(ns,1:length(data.flog)) = specres;
        end
    end

    if (data.ns1(ns) > 0)
        %#### this is for S-wave.
        ns1 = data.ns1(ns);
        if (ns1 < 1)
            ns1 = 1;
        end
        nwin = round(data.tswin(ns)*data.samprate(ns));
        nend = ns1+nwin-1;
        if (nend > data.npts(ns))
            nend = data.npts(ns);
            nwin = data.npts(ns)-ns1+1;
        end
        data.snwin(ns) = nwin;
        if (nwin >= 50)

            data.sseis(ns,1:nwin) = data.seis(ns,ns1:ns1+nwin-1);
            if (ns1-nwin<1)
                data.snoise(ns,1:nwin) = zeros(1,nwin);
                data.snoise(ns,1:ns1-1) = data.seis(ns,1:ns1-1);
            else
                data.snoise(ns,1:nwin) = data.seis(ns,ns1-nwin:ns1-1);
            end
            nx = Nfft;

            [spec,~] = pmtm(data.sseis(ns,1:nwin),3,nx,data.samprate(ns));
            [nspec,f] = pmtm(data.snoise(ns,1:nwin),3,nx,data.samprate(ns));
            omega = 2*pi*f;
            nfg = find(f>0);
            spec = sqrt(spec);
            nspec = sqrt(nspec);
            data.sspec(ns,1) = spec(1);
            data.snspec(ns,1) = nspec(1);
            data.sspec(ns,nfg) = spec(nfg)./omega(nfg);
            data.snspec(ns,nfg) = nspec(nfg)./omega(nfg);
            data.sfreq(ns,nfg) = f(nfg);
            % nband = find(f(nfg)>=4 & f(nfg)<=10);
            % data.sstn(ns) = mean(log10(data.sspec(ns,nband))-log10(data.snspec(ns,nband)));
            specres = interp1(f(nfg),data.sspec(ns,nfg),data.flog,'linear','extrap');
            data.sspecres(ns,1:length(data.flog)) = specres;
        end

        % this is for coda wave (Don't use yet. Not ready).
        ns1 = data.ns1(ns)-round(0.1*data.samprate(ns)); % allow 0.1 second before
        if (ns1<1) ns1 = 1; end
        nwin = round(data.codawin(ns)*data.samprate(ns));
        nend = ns1+nwin-1;
        if (nend > data.npts(ns)) nend = data.npts(ns); nwin = data.npts(ns)-ns1+1; end
        data.cnwin(ns) = nwin;
        data.cseis(ns,1:nwin) = data.seis(ns,ns1:nend);
        nxx1 = ns1-nwin;
        nwin1 = nwin;

        %if (ns1-nwin<1)
        %    nxx1 = 5;
        %    nwin1 = ns1-nxx1;
        %end
        if (ns1-nwin<1)
            data.cnoise(ns,1:nwin) = zeros(1,nwin);
            data.cnoise(ns,1:ns1-1) = data.seis(ns,1:ns1-1);
        else
            data.cnoise(ns,1:nwin) = data.seis(ns,ns1-nwin:ns1-1);
        end
        %data.cnoise(ns,1:nwin1) = data.seis(ns,nxx1:ns1-1);
        %nx = 2^(nextpow2(nwin));
        nx = Nfftc;
        [spec,~] = pmtm(data.cseis(ns,1:nwin),3,nx,data.samprate(ns));
        [nspec,f] = pmtm(data.cnoise(ns,1:nwin1),3,nx,data.samprate(ns));
        %factor = nx*samprate(iw)*2*pi;
        omega = 2*pi*f;
        nfg = find(f>0);
        spec = sqrt(spec);
        nspec = sqrt(nspec);
        data.pspec(ns,1) = spec(1);
        data.pnspec(ns,1) = nspec(1);
        data.cspec(ns,nfg) = spec(nfg)./omega(nfg);
        data.cnspec(ns,nfg) = nspec(nfg)./omega(nfg);
        data.cfreq(ns,nfg) = f(nfg);
        specres = interp1(f(nfg),data.cspec(ns,nfg),data.flog,'linear','extrap');
        data.cspecres(ns,1:length(data.flog)) = specres;
    end

end

try
    data.pseis;
catch
    data.pseis = [];
    data.pnoise = [];
    data.pspecres = [];
end
try
    data.sseis;
catch
    data.sseis = [];
    data.snoise = [];
    data.sspecres = [];
end
try
    data.cseis;
catch
    data.cseis = [];
    data.cnoise = [];
    data.cspecres = [];
end
numTraceAvail = size(data.pseis, 1);
trace_diff = data.numts - numTraceAvail;

data.pseis = [data.pseis; zeros(trace_diff, size(data.pseis, 2))];
data.pnoise = [data.pnoise; zeros(trace_diff, size(data.pnoise, 2))];
data.pspecres = [data.pspecres; zeros(trace_diff, size(data.pspecres, 2))];

data.sseis = [data.sseis; zeros(trace_diff, size(data.sseis, 2))];
data.snoise = [data.snoise; zeros(trace_diff, size(data.snoise, 2))];
data.sspecres = [data.sspecres; zeros(trace_diff, size(data.sspecres, 2))];
data.cseis = [data.cseis; zeros(trace_diff, size(data.cseis, 2))];
data.cnoise = [data.cnoise; zeros(trace_diff, size(data.cnoise, 2))];
data.cspecres = [data.cspecres; zeros(trace_diff, size(data.cspecres, 2))];

save(evfile,'data','-v7.3');

