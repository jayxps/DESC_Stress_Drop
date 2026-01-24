clear
clc

%% Add library paths
addpath('../src');
addpath('../src/sacstuff')
addpath('../src/matlabtool/export_fig');

%% Read in parmeters, see file 'specinp'
P = specprocess_readpara('specinp');
save P.mat P
load P.mat

%% Rename dataset by study cases
name_suffix = 'python';
[spec,freq] = specprocess_organize(P);
save spec.mat spec freq
movefile('spec.mat',['spec_' name_suffix '.mat']);
load(['spec_' name_suffix '.mat']);

%% Separate observed spectra into evspec (event term), stspec (site term)
% and distspec (path term, not including attenuation)
[evspec,stspec,distspec] = specprocess_inversion(spec,P,freq);
save(['specsolve_' name_suffix '.mat'],'evspec','stspec','distspec');
load(['specsolve_' name_suffix '.mat']);


%% At least nspecmin stations record all the earthquakes, see Shearer et al. (2006), Figure 10
nspecmin = 4;
evspec_save = evspec;
ng = find([evspec_save.nspec] >= nspecmin);
evspec = evspec_save(ng);

%% Magnitude calibration, convert other magnitude scales to moment magnitude
% using low-frequency amplitudes, see Shearer et al. (2006), Equation 3
freqlim = [2, 3]; % Frequency band to calculate low-freq amplitude
[y0,slopeline] = specprocess_plotmag(evspec,freqlim,nspecmin);
itype = 2;
fmw = 3.0;
maginfo.y0 = y0;
maginfo.slope = slopeline;
maginfo.fmw = fmw;
maginfo.itype = itype;
evspec = specprocess_specmag(evspec,freqlim,maginfo);
save maginfo.mat maginfo

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%----- Now evspec (event term) can be used for further analysis -----%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

save(['evspec_' name_suffix '.mat'],'freq','evspec','y0','slopeline','distspec','stspec','freqlim')