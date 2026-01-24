function evspec = specprocess_specmag(evspec,freqlim,maginfo)
itype = maginfo.itype;
y0 = maginfo.y0;
slopeline = maginfo.slope;   % slopeline(k0) and y0 indicates the relationship between catalog magnitude and spectra amplitude (Mc=k0log10(A)+y0)
fmw = maginfo.fmw;

if (itype == 1)
    for i = 1:length(evspec)
        evspec(i).qmagest = evspec(i).qmag;
        evspec(i).qmomest = 10.^(1.5*(evspec(i).qmag+10.7))/1e7;
        evspec(i).qmagresid = 0;
    end
    return;
end

if (itype == 2)
    a=(fmw+10.7)*1.5*slopeline-fmw;

    freq = evspec(1).freq;
    nfx = freq >= freqlim(1) & freq<=freqlim(2);
    fid = fopen('tempmom','w');

    for i = 1:length(evspec)
        if (isempty(evspec(i).spec)) continue; end
        qmom = mean(evspec(i).spec(nfx));
        if(qmom~=0)
            qmagest = y0+qmom*slopeline;  % Use spectral amplitude and k0,y0 to estimate spectral magnitude
        else
            qmagest = 0;
        end
        qmomest = 10^((qmagest+a)/slopeline);  % Use spectral amplitude and new scaling to calculate estimated moments
        evspec(i).qmomest = qmomest/1.e7;  % Turn dyne.cm into Nm
        evspec(i).qmagest = qmagest;
        evspec(i).qmagresid = evspec(i).qmag - evspec(i).qmagest;
        fprintf(fid,'%10.6f %10.6f %10.6f %10.6f %13d %10.6f\n',evspec(i).qmag,...
        	qmom,log10(qmomest),qmagest,evspec(i).qid,a);
    end
    fclose(fid);
end
