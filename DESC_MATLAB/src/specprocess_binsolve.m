function [amed, ninbin] = specprocess_binsolve(a,index,nbin,itype)
% computes mean or median within specified bins of data array
[M, N] = size(a);

amed = zeros(M,nbin);
ninbin = zeros(M,nbin);

if (itype ~=3)
    for ibin = 1:nbin
        nin = find(index == ibin);
        ninbin(ibin) = length(nin);
        if (ninbin(ibin)==0)
            amed(1:M,ibin) = 0;
            continue;
        end
        if (itype == 1) % compute median
            amed(:,ibin) = median(a(:,nin),2);
        elseif (itype == 2) % compute mean
            amed(:,ibin) = mean(a(:,nin),2);
        end
    end
else
    xgap = 0.5;
    xgapinv = 1./xgap;
    nit = 5;
    atemp = ones(size(a));
    for it = 1:nit
        atemp_new = ones(size(a));
        for ibin = 1:nbin
            nin = find(index == ibin);
            ninbin(ibin) = length(nin);
            if (ninbin(ibin)==0)
                amed(1:M,ibin) = 0;
                continue;
            end
            temp = atemp(:,nin).*a(:,nin);
            sumw = sum(atemp(:,nin),2);
            amed(:,ibin) = sum(temp,2)./sumw;

            if (it ~= nit)
                d = abs(a(1:M,nin) - repmat(amed(1:M,ibin),1,length(nin)));
                atempx = 1./d;
                atempx(d<=xgap) = xgapinv;
                atemp_new(:,nin) = atempx;
            end
        end
        atemp = atemp_new;
    end
end
