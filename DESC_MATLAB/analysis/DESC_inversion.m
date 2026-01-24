clear

%% Add paths
addpath('../src');
addpath('../src/matlabtool/export_fig');
addpath('../src/DE_Module');

%% Predefine parameters and load necessary dataset
prefix_all = {'brune','boatwright','doublecorner'}; % Don't change
method_all = [1 5 11]; % Don't change
name_suffix = 'test';
load(['evspec_' name_suffix '.mat']);
load P.mat
nspecmin = 4;

%% Initialize stress drop inversion parameters
falloff = 2:0.1:2; % Falloff range, don't need to change
maglim = [0.9 4.0]; % Mw bin range, usually wider is better
d_mag = 0.2; % Mw bin size, usually 0.2-0.3
ddsig = [0.01, 100]; % Stress drop searching range. in log10
min_nspec = 10; % Min number of events in each Mw bin
iwave = 1; % 1=P, 2=S
freqfit = [2, 40]; % Frequency band for spectral fitting

%% Build magnitude bins information
load maginfo.mat
itype = 2; % Don't change
[stack_spec] = specprocess_stackmag(evspec,nspecmin,maginfo,1,d_mag);  %% use the same stack for all source models

%% Stress drop inversion
isource = 1; % Use Brune's spectrum model
disp(['Solve for ' prefix_all{method_all==isource} ' model:'])
flow = freqfit(1);
fhigh = freqfit(2);
disp(['working on freqfit=' num2str(flow) '-' num2str(fhigh)])
tic
freqfit = [flow fhigh];
imethod = method_all(isource);
prefix = prefix_all{isource};
nlow = find(freq >= freqlim(1) & freq <= freqlim(2));

[ECS,stack_spec,EGF_all,~,~,~] = specprocess_floategf_allEQs(evspec,stack_spec,freqfit,...
    ddsig,maglim,imethod,min_nspec,nlow);

outdir = ['result-indieq-',prefix,'-f',num2str(freqfit(1)),'-',num2str(fhigh)];
outpara = ['para_ECS_',prefix,'.mat'];
outfig = ['result-all-',prefix,'-f',num2str(freqfit(1)),'-',num2str(fhigh)];
if (~exist(outfig,'dir'))
    mkdir(outfig);
end

para_EGF = specprocess_sourcepara(outdir,freqfit,evspec,...
    imethod,iwave, ECS, 1, 0);
save(outpara,'para_EGF','stack_spec','ECS','distspec','stspec','evspec');

movefile('*.pdf',outfig);
movefile(outpara, outfig);
toc
