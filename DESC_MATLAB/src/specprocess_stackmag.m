% this is the function to process stacking.
% input:
% 	evspec structure file
%   frequency limit for moment calibration
%	moment calibration choices
% output:
% 	stack
function [stack_spec] = specprocess_stackmag(evspec,nspecmin,maginfo,isavefig,mbin)
itype = maginfo.itype;
fmw = maginfo.fmw;
y0 = maginfo.y0;
slopeline = maginfo.slope;

freq = evspec(1).freq;

%***** now will compute linear scaling between different frequency bands
%and magnitude..

%nfx = find(freq >= freqlim(1) & freq<=freqlim(2));
ng = find([evspec.nspec]>=nspecmin);
specmag = [evspec(ng).qmagest];
mag = [evspec(ng).qmag];

%specmag = zeros(size(ng));
%mag = zeros(size(ng));

%for k = 1:length(ng)
%    	specmag(k) = y0 + slopeline * mean(evspec(ng(k)).spec(nfx));
%	mag(k) = evspec(ng(k)).qmag;
%end
clear a
if (itype == 1);
    qmag = mag;
else
    qmag = specmag;
    a = (fmw+10.7)*1.5*slopeline-fmw;
    %    resid = mag - specmag;
end

%qmag = [evspec(ng).qmagest];
resid = [evspec(ng).qmagresid];

% im_min = round(min(qmag(abs(resid)<=1.5))*5. + 0.501);
% im_max = round(max(qmag(abs(resid)<=1.5))*5. + 0.501);
% [min(qmag(abs(resid)<=1.5)) max(qmag(abs(resid)<=1.5))]
%
% im_corr = 1-im_min;
% im_total = im_max - im_min + 1;
%
% im_ev = round(qmag*5.+0.501) + im_corr;

%mbin = 0.2;
magbins = -1:mbin:8;
im_min = 1;
im_max = length(magbins);
im_corr = 0;
im_total = length(magbins)+1;
for i = 1:length(qmag)
    im_ev(i) = floor((qmag(i)-min(magbins))/mbin)+1;
end

%for k = 1:length(ng)
%	fprintf(1,'%d %f %d\n',[ng(k) qmag(k) im_ev(k)]);
%end

for im = 1: im_total
    stack_spec(im).spec = zeros(size(freq));
    stack_spec(im).nspec = 0;
    stack_spec(im).mag = (im-im_corr)*mbin-mbin/2 + min(magbins);
    stack_spec(im).freq = freq;
    stack_spec(im).fmom = 0;
    %stack_spec(im).fmw = 0;
    %if (itype == 1)
    %	nin = find(im_ev == im);
    %elseif (itype == 2)
    nin = find(im_ev == im & abs(resid)<=1.5);
    %end

    nspec = length(nin);
    if (nspec < 1)
        continue;
    end

    for ii = 1:nspec
        iq = ng(nin(ii));
        stack_spec(im).spec = stack_spec(im).spec + evspec(iq).spec;
    end
    stack_spec(im).spec = stack_spec(im).spec / nspec;
    stack_spec(im).nspec = nspec;
    %stack_spec(im).fmom = mean([evspec(ng(nin)).qmomest]);
    if (itype == 1)
        stack_spec(im).fmom = 10.^((stack_spec(im).mag+10.7)*1.5)/1e7;
    else
        stack_spec(im).fmom = 10.^((stack_spec(im).mag+a)/slopeline)/1e7;
    end

    %stack_spec(im).fmw = log10(stack_spec(im).fmom)/1.5-10.7;
end
if (isavefig == 1)
    f1=figure(100);
    set(f1,'visible','off');
    for im = 1: im_total
        if (stack_spec(im).nspec < 1) continue; end
        semilogx(freq,stack_spec(im).spec);
        hold on
        textm = sprintf('M%5.2f NS %d',stack_spec(im).mag,stack_spec(im).nspec);
        text(freq(2),stack_spec(im).spec(2),textm);
    end
    hold off
    xlabel('frequency (Hz)');
    ylabel('amplitude');
    xlim([0.01 max(freq)]);
    saveas(gcf,'stackmag.pdf');
    close(100);
end


