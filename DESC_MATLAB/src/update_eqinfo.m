% this is used to generate eqinfo.mat from raw catalog
function eqinfo = update_eqinfo(catfile,itype)

%run('/opt/antelope/5.5/setup.m');
eqlist = geteqlist_antelope(catfile,itype);

neq = length(eqlist(:,1));

for ie = 1:neq
    qtime(ie) = eqlist(ie,5);
    timenum = epoch2date(qtime(ie));
    timex = datevec(timenum);
    %timestring=epoch2str(qtime(ie),'%Y/%m/%d,%H:%M:%S.%s');
    %timex=sscanf(timestring,'%d/%d/%d,%d:%d:%f');
    doy = datenum(timex(1:3)) - datenum([timex(1),1,1])+1;

    eqinfo.id(ie) = eqlist(ie,6);
    eqinfo.qlat(ie) = eqlist(ie,1);
    eqinfo.qlon(ie) = eqlist(ie,2);
    eqinfo.qdep(ie) = eqlist(ie,3);
    eqinfo.mb(ie) = eqlist(ie,4);
    eqinfo.qtime(ie) = qtime(ie);   %this is in epoch time. need to convert to regular time.
    eqinfo.qyr(ie) = timex(1);
    eqinfo.qmon(ie) = timex(2);
    eqinfo.qdy(ie) = timex(3);
    eqinfo.qdoy(ie) = doy;
    eqinfo.qhr(ie) = timex(4);
    eqinfo.qmn(ie) = timex(5);
    eqinfo.qsc(ie) = timex(6);
    eqinfo.index(ie) = 0;
end
