function [evspec,stspec,distspec] ...
    = specprocess_inversion(spec,P,freq)

nf0 = length(freq);

nq = max([spec.inq]);
ev_spec = zeros(nf0,nq);

ndep = max([spec.ind]);
dep_spec = zeros(nf0,ndep);

nsta = max([spec.ins]);
st_spec = zeros(nf0,nsta);

ndmax = max([spec.inx]);
dist_spec = zeros(nf0,ndmax);

nmag = max([spec.inm]);
mag_spec = zeros(nf0,nmag);

depterm = zeros(1,ndep);
distterm = zeros(1,ndmax);
magterm = zeros(1,nmag);
stterm = zeros(1,nsta);
evterm = zeros(1,nq);

npts = length(spec);

for it = 1:P.nitmax
    tic
    for ispec = 1:npts
        is = spec(ispec).ins;
        iq = spec(ispec).inq;
        im = spec(ispec).inm;
        ix = spec(ispec).inx;
        id = spec(ispec).ind;

        y(ispec) = spec(ispec).tstar ...
            - depterm(id) - magterm(im) ...
            - stterm(is) - evterm(iq);

        yspec(:,ispec) = spec(ispec).spec' - dep_spec(:,id) ...
            - mag_spec(:,im) - st_spec(:,is) ...
            - ev_spec(:,iq);
    end
    [distterm,~] = specprocess_binsolve(y,[spec(:).inx],ndmax,P.itype);
    [dist_spec,~] = specprocess_binsolve(yspec,[spec(:).inx],ndmax,P.itype);
    % ****** 'done with dist_spec';
    for ispec = 1:npts
        is = spec(ispec).ins;
        iq = spec(ispec).inq;
        im = spec(ispec).inm;
        ix = spec(ispec).inx;
        id = spec(ispec).ind;

        y(ispec) = spec(ispec).tstar ...
            - depterm(id) - magterm(im) ...
            - distterm(ix) - evterm(iq);

        yspec(:,ispec) = spec(ispec).spec' - dep_spec(:,id) ...
            - mag_spec(:,im) - dist_spec(:,ix) ...
            - ev_spec(:,iq);
    end

    [stterm,~] = specprocess_binsolve(y,[spec(:).ins],nsta,P.itype);
    [st_spec,~] = specprocess_binsolve(yspec,[spec(:).ins],nsta,P.itype);
    % **** 'done with station terms'

    for ispec = 1:npts
        is = spec(ispec).ins;
        iq = spec(ispec).inq;
        im = spec(ispec).inm;
        ix = spec(ispec).inx;
        id = spec(ispec).ind;

        y(ispec) = spec(ispec).tstar ...
            - depterm(id) - magterm(im) ...
            - distterm(ix) - stterm(is);

        yspec(:,ispec) = spec(ispec).spec' - dep_spec(:,id) ...
            - mag_spec(:,im) - dist_spec(:,ix) ...
            - st_spec(:,is);
    end

    [evterm,~] = specprocess_binsolve(y,[spec(:).inq],nq,P.itype);
    [ev_spec,~] = specprocess_binsolve(yspec,[spec(:).inq],nq,P.itype);
    % **** 'done with event terms'

    disp(strcat('iteration:',num2str(it)));
    toc

    % **** 'compute rms'
    disp('computing RMS');
    sum2 = 0;
    sumy = 0;
    for ispec = 1:npts
        is = spec(ispec).ins;
        iq = spec(ispec).inq;
        im = spec(ispec).inm;
        ix = spec(ispec).inx;
        id = spec(ispec).ind;

        y(ispec) = spec(ispec).tstar ...
            - depterm(id) - magterm(im) ...
            - distterm(ix) - stterm(is) - evterm(iq);
        sumy = sumy + y(ispec)^2;

        yspec(:,ispec) = spec(ispec).spec' - dep_spec(:,id) ...
            - mag_spec(:,im) - dist_spec(:,ix) ...
            - st_spec(:,is) - ev_spec(:,iq);
        sum2 = sum2 + sum(yspec(:,ispec).^2);
    end

    rms2 = sqrt(sum2)/(npts*nf0);
    rms = sqrt(sumy)/npts;
    fprintf(1,'it,rms,rms2= %d %12.8f %12.8f\n',it,rms,rms2);
end

% now done with inversion, will organize data, and output.
% first is event_term
% will output: event spectrum, event time, id, location, number of spectrum.
iq = 0;
for i = 1:nq
    nin = find([spec.inq]==i);
    if ( isempty(nin));
        continue;
    end
    iq = iq + 1;
    evspec(iq).spec = ev_spec(:,i)';
    evspec(iq).nspec = length(nin);
    evspec(iq).freq = freq;
    evspec(iq).qid = spec(nin(1)).qid;
    evspec(iq).qtime = spec(nin(1)).qtime;
    evspec(iq).qlat = spec(nin(1)).qlat;
    evspec(iq).qlon = spec(nin(1)).qlon;
    evspec(iq).qdep = spec(nin(1)).qdep;
    evspec(iq).qmag = spec(nin(1)).qmag;
    evspec(iq).tstar = evterm(iq);
end

% now is station terms.
is = 0;
for i = 1:nsta
    nin = find([spec.ins] == i);
    if ( isempty(nin));
        continue;
    end
    is = is + 1;
    stspec(is).spec = st_spec(:,i)';
    nin = find([spec.ins] == is);
    if (isempty(nin)) continue; end
    stspec(is).nspec = length(nin);
    stspec(is).freq = freq;
    stspec(is).stlat = spec(nin(1)).stlat;
    stspec(is).stlon = spec(nin(1)).stlon;
    stspec(is).stelev = spec(nin(1)).stelev;
    stspec(is).stid = spec(nin(1)).stid;
    stspec(is).tstar = stterm(is);
end

% now is distance terms.
ix = 0;
for i = 1:ndmax
    nin = find([spec.inx] == i);
    if ( isempty(nin))
        continue;
    end
    ix = ix + 1;
    distspec(ix).spec = dist_spec(:,i)';
    distspec(ix).nspec = length(nin);
    distspec(ix).freq = freq;
    distspec(ix).tt = i*P.dtt-P.dtt/2;
    distspec(ix).tstar = distterm(ix);
end

