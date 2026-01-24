function [y0,slopeline] = specprocess_plotmag(varargin)
% this program will calibrate local magnitude to spectrum magnitude
% input: evspec, frequence limits to compute spectral amplitude
% output: y0, slopeline and a figure saved in local directory.
% input variable list:
% 1. evspec
% 2. freqlim
% 3. nspecmin
% 4. savefig (optional)

evspec = varargin{1};
freqlim = varargin{2};
nspecmin = varargin{3};

if(nargin == 3)
    savefig = 1;
elseif(nargin == 4)
    savefig = varargin{4};
end

freq = evspec(1).freq;

%***** now will compute linear scaling between different frequency bands
%and magnitude..

magcut = 0;

nfx = freq >= freqlim(1) & freq<=freqlim(2);
ng = find([evspec.nspec]>=nspecmin);
qid = [evspec(ng).qid];

for k = 1:length(ng)
    specamp(k) = mean(evspec(ng(k)).spec(nfx));
    mag(k) = evspec(ng(k)).qmag;% + (rand(1)-0.5)/10.;
end

ng_all = find([evspec.nspec]>=nspecmin);
for kk = 1:length(ng_all)
    specamp_all(kk) = mean(evspec(ng_all(kk)).spec(nfx));
end

bfit = robustfit(specamp,mag);

save specamp.mat specamp qid
ypred = bfit(1)+specamp*bfit(2);
specampx = -6:0.1:20;
ypredx = bfit(1)+specampx*bfit(2);
slopeline = bfit(2);
y0 = bfit(1);

if(savefig == 1)
    f1=figure(999);
    set(f1,'visible','off');
    plot(specamp_all,[evspec(ng_all).qmag],'.','MarkerSize',15,'color',[0.7,0.7,0.7]);
    hold on
    plot(specamp,mag,'k.','MarkerSize',15);
    plot(specampx,ypredx,'LineWidth',2);
    xlim([min(specamp)-0.5, max(specamp)+0.5]);
    hold off
    titlex = sprintf('%s%4.2f%s%4.2f%s%4.2f%s%4.2f','freqlim=',freqlim(1),'-',freqlim(2),...
        ' slope=',bfit(2),' y0=',bfit(1));
    title(titlex);
    xlabel('relative log magnitude');
    ylabel('catalog magnitude');
    figout = strcat('magfit-freq-',num2str(freqlim(1)),'-',num2str(freqlim(2)),'.jpg');
    saveas(gcf, figout)
    close(f1);
end
