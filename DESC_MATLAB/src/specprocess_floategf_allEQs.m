% input:
% 	stacked-spec.
%   frequency range for fitting, e.g., [1 20]
%   falltarg: falloff rate
%   dsigtarg: stress drop range
%   maglim: magnitude limit for fitting
%   min_nspec: for each stacked bin, required at least min_nspec events.
%   isavefig: save the stacked corrected figures
%
% output:
% 	EGF, dsigbest, fcc, stack_spec (updated with correction).
%
function [EGF, stack_spec, ECSs_info, log10dsigbins, nmagg, evspec_DE] = specprocess_floategf_allEQs(evspec,stack_spec,freqfit,...
    ddsig,maglim,imethod,min_nspec,nlow)
global falloffin

falloffin = 2;

isavefig = 1;

% if (ismember(imethod,brunemethod))
%     isource = 1;
% else
%     isource = 2;
% end

beta = 3464;
% fact = (0.42*beta)^3;

freq = stack_spec(1).freq;
nfreqg = find(freq>=0);
freqg = freq(nfreqg);

n1 = find(freqg<=freqfit(1));
n2 = find(freqg>=freqfit(2));
nfit = n1(end):n2(1);

xmag = [stack_spec.mag];
dmag = xmag(2) - xmag(1);
xnspec = [stack_spec.nspec];

nmaggx = find(xmag >=maglim(1)-0.05 & xmag <=maglim(2)+0.05 ...
    & xnspec >= min_nspec);
nmagg = min(nmaggx):max(nmaggx);
if (isempty(nmagg))
    EGF = zeros(size(freq));
    stack_spec = zeros(size(freq));
    return;
end

for i = 1:length([stack_spec.mag])
    stack_spec(i).fc = nan;
    stack_spec(i).dsig = nan;
    stack_spec(i).fc2 = nan;
    stack_spec(i).allspecs = [];
    stack_spec(i).evspecind = [];
    stack_spec(i).fmom = [];
end
% nx=find(xmag == 3.5);
% nsum = length(nmagg);
ifit = zeros(size(stack_spec));
ifit(nmagg) = 1;

evspec_DE = evspec([evspec.qmagest]>=xmag(nmagg(1))-dmag/2 & [evspec.qmagest]<xmag(nmagg(end))+dmag/2);

for iq = 1:length(evspec_DE)
    qmagest = evspec_DE(iq).qmagest;
    [~,ind] = min(abs(qmagest-xmag(nmagg)));
    evspec_DE(iq).ibin = ind;
    stack_spec(nmagg(ind)).allspecs = [stack_spec(nmagg(ind)).allspecs;evspec_DE(iq).spec];
    stack_spec(nmagg(ind)).evspecind = [stack_spec(nmagg(ind)).evspecind iq];
    stack_spec(nmagg(ind)).fmom = [stack_spec(nmagg(ind)).fmom evspec_DE(iq).qmomest];
    evspec_DE(iq).magbin = xmag(nmagg(ind));
end

for iq = 1:length(evspec_DE)
    evspec_DE(iq).fmom_all = geomean([stack_spec(nmagg(evspec_DE(iq).ibin)).fmom]);
end

% Organize Inputs for objfun
log10sig = [log10(ddsig(1)),log10(ddsig(2))];
specall = reshape([evspec_DE.spec]',length(freq),length(evspec_DE))';

S.nfreqg = nfreqg;
S.freqg = freqg;
S.nlow = nlow;
S.stack_spec = stack_spec;
S.nmagg = nmagg;
S.log10sig = log10sig;
S.evspec_DE = evspec_DE;
S.specall = specall;
S.nfit = nfit;
S.min_nspec = min_nspec;

[log10dsigbins,outputs] = Rundeopt(S);

ECSs_allEQs = outputs.ECSs_med;
ECSs = outputs.ECSs;
ECSs_info.ECSs_allEQs = ECSs_allEQs;
ECSs_info.ECSs = ECSs;

for ibin = 1:length(nmagg)
    im = nmagg(ibin);
    stack_spec(im).dsig = 10^log10dsigbins(ibin);
end

EGF = nanmedian(ECSs_allEQs,1);
EGF = smooth(EGF,10)';

clear specdiff dsig
for img = 1:length(nmagg)
    im = nmagg(img);
    ifit(im) = 1;
    fmom = geomean(stack_spec(im).fmom);
    dspec = stack_spec(im).spec - ECSs_allEQs(img,:);

    [fc,fc2,o0,~,~,~,~,~,err,dspecfit,iflag,ftest,etest,fcb1,fcb2,~] = ...
        sourcepara(dspec,freqg,nfit,nlow,imethod,fmom,1);
    if (imethod == 11)
        if (fc2<fc)
            x = fc;
            fc = fc2;
            fc2 = x;
        end
    end

    stack_spec(im).fc = fc;
    stack_spec(im).o0 = o0;
    stack_spec(im).ftest = ftest;
    stack_spec(im).etest = etest;
    stack_spec(im).specfit = dspecfit;
    stack_spec(im).speccorr = dspec;
    %stack_spec(im).dsig = fc^3*fmom/fact/1e6;
    stack_spec(im).fc2 = fc2;
    specdiff(img,:) = dspecfit - dspec;

    if (img == 1)
        fcref = fc;
        fmomref = geomean(stack_spec(im).fmom);
    end
    fcpred(im) =  fmomref^(1/3) /geomean(stack_spec(im).fmom)^(1./3)  * fcref;
    o0pred(im) = interp1(freqg,dspec,fcpred(im));
    mompred(im) = geomean(stack_spec(im).fmom);

end

nfcfit = find(fcpred>0);

if (isavefig == 1)
    close all

    f2=figure(100);
    set(f2,'visible','off');

    ng = find([stack_spec.nspec]>1);
    % nfreqx = find(freqg>=1 & freqg<=5);
    for img = 1:length(ng)
        im = ng(img);
        if (ifit(im) == 1)
            specx = stack_spec(im).speccorr;
            semilogx(freqg,stack_spec(im).speccorr,'k');
            hold on
            semilogx(freqg,stack_spec(im).specfit,'r:')
            semilogx(freqg,stack_spec(im).spec-stack_spec(im).speccorr+mean(specx(nlow)),'b--');
            textm = [sprintf('M%5.2f NEQ=%d',stack_spec(im).mag,stack_spec(im).nspec) ' ' char(916) char(963) '=' num2str(stack_spec(im).dsig,'%5.2f') 'MPa'];
            text(freqg(2),specx(2)+0.1,textm);
            hold on
            specamp = interp1(freqg,stack_spec(im).specfit,stack_spec(im).fc);
            semilogx(stack_spec(im).fc,specamp,'ro');
            if (stack_spec(im).fc2<50)
                specamp2 = interp1(freqg,stack_spec(im).specfit,stack_spec(im).fc2);
                semilogx(stack_spec(im).fc2,specamp2,'bd');
            end
        elseif (stack_spec(im).nspec>1 & ~isempty(stack_spec(im).fmom))
            specx = stack_spec(im).speccorr(nfreqg);
            semilogx(freqg,specx,'k-.');
            hold on
            textm = [sprintf('M%5.2f NEQ=%d',stack_spec(im).mag,stack_spec(im).nspec) ' ' char(916) char(963) '=' num2str(stack_spec(im).dsig,'%5.2f') 'MPa'];
            text(freqg(2),specx(2)+0.1,textm);
        end
    end
    semilogx(fcpred(nfcfit),o0pred(nfcfit),'k--');
    xlabel('frequency (Hz)');
    ylabel('amplitude');
    titlex = sprintf('%s%8.3f%s','ECS corrected spectra');
    title([titlex ', median \Delta\sigma=' num2str(median([stack_spec(nmagg).dsig]),'%5.2f') 'MPa']);
    saveas(gcf, 'stacked-spectra-after-EGF.jpg');
    close(100);

end
