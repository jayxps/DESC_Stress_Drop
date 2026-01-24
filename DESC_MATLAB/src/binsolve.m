function [amed, temp,atemp,d,ninbin] = binsolve(a,npts,index,nbin,itype)
% computes mean or median within specified bins of data array
amed = zeros(1,nbin);
ninbin = zeros(1,npts);
temp = zeros(1,npts);
atemp= zeros(1,npts);
d= zeros(1,npts);
if (itype ~=3)
    for ibin = 1:nbin
        nin = find(index == ibin);
        ninbin(ibin) = length(nin);
        if (ninbin(ibin)==0)
            amed(ibin) = 0;
            continue;
        end
        if (itype == 1) % compute median
            amed(ibin) = median(a(nin));
        elseif (itype == 2) % compute mean
            amed(ibin) = mean(a(nin));
        end
    end
else
    xgap = 0.2;
    xgapinv = 1./xgap;
    nit = 5;
    atemp = ones(1,npts);
    for it = 1:nit
        for ibin = 1:nbin
            nin = find(index == ibin);
            ninbin(ibin) = length(nin);
            if (ninbin(ibin)==0)
                amed(ibin) = 0;
                continue;
            end
            temp = atemp(nin).*a(nin);
            amed(ibin) = mean(temp);
        end
        if (it~=nit)
            d = abs(a-amed(ibin));
            nle = d<=xgap;
            atemp = 1./d;
            atemp(nle) = xgapinv;
        end
    end
end

