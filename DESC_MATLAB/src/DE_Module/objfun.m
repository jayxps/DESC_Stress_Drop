%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function:         S_MSE= objfun(FVr_temp, S_struct)
% Author:           Rainer Storn
% Description:      Implements the cost function to be minimized.
% Parameters:       FVr_temp     (I)    Paramter vector
%                   S_Struct     (I)    Contains a variety of parameters.
%                                       For details see Rundeopt.m
% Return value:     S_MSE.I_nc   (O)    Number of constraints
%                   S_MSE.FVr_ca (O)    Constraint values. 0 means the constraints
%                                       are met. Values > 0 measure the distance
%                                       to a particular constraint.
%                   S_MSE.I_no   (O)    Number of objectives.
%                   S_MSE.FVr_oa (O)    Objective function values.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [S_MSE, power]= objfun(FVr_temp, S_struct, S)

%---Objective function----------------------------------------------

%% Extract variables
nfreqg     = S.nfreqg;
freqg      = S.freqg;
nlow       = S.nlow;
evspec     = S.evspec_DE;
specall    = S.specall;
ibin       = S.ibin;
fmom       = S.fmom;
nfit       = S.nfit;
min_nspec  = S.min_nspec;

%% Organize variables

%% Compute objective function values
% Remember: here FVr_temp is log10 stress drop!!!
ng = find([evspec.qmagresid]<=100);
evspec = evspec(ng);
specall = specall(ng,:);
ibin = ibin(ng);
fmom = fmom(ng);
beta = 3464;
fact = (0.42*beta)^3;

%%% Optimize for matrix computation
tic
dsigtarg = (10.^FVr_temp(ibin)) * 1e6; % convert to Pa.
fc = (dsigtarg.*fact./fmom).^(1/3);

pred = -log10(1+(repmat(freqg,length(evspec),1)./repmat(fc',1,length(freqg))).^2);
yoff = mean(specall(:,nfreqg(nlow))-pred(:,nfreqg(nlow)),2);
ECSs = specall - repmat(yoff,1,length(freqg)) - pred;

ibin_u = unique(ibin(~isnan(ibin)));
ECSs_med = nan*ones(length(ibin_u),size(ECSs,2));
for iibin = 1:length(ibin_u)
    if(sum(ibin==ibin_u(iibin))>=min_nspec)
        ECSs_med(iibin,:) = median(ECSs(ibin==ibin_u(iibin),:),1);
    end
end

ST = sum(nanstd(ECSs_med(:,nfit),[],1));

power.ECSs_med = ECSs_med;
power.ECSs = ECSs;
power.fc = fc;

%% Output
% ----strategy to put everything into a cost function------------
S_MSE.I_nc      = 0;%no constraints
S_MSE.FVr_ca    = 0;%no constraint array
S_MSE.I_no      = 1;%number of objectives (costs)
S_MSE.FVr_oa(1) = 1/ST;
