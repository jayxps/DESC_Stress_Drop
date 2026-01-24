function [spec,flog] = specprocess_organize(P)

%eqfile = strcat(P.evdir,'eqinfo.mat');
eqfile = P.eqinfo;
load(eqfile);

nrand = randi(length(eqinfo.id),1)/length(eqinfo.id);

n_select = find(eqinfo.qlat >= P.qlat1 & ...
    eqinfo.qlat <= P.qlat2 & ...
    eqinfo.qlon >= P.qlon1 & ...
    eqinfo.qlon <= P.qlon2 & ...
    eqinfo.mb >= P.minmag & ...
    eqinfo.mb <= P.maxmag & ...
    nrand <= P.randfrac);

neq = length(n_select);
fprintf(1,'total number of events within limits: %d out of %d events\n',[neq,length(eqinfo.id)]);

nsta = 0;
nspec = 0;
for iq = 1:neq
    indi = n_select(iq);
    qyr = num2str(eqinfo.qyr(indi));
    if(eqinfo.qmon(indi)<10)
        qmon = ['0' num2str(eqinfo.qmon(indi))];
    else
        qmon = num2str(eqinfo.qmon(indi));
    end
    event_file = strcat(P.evdir, qyr, '/', qmon, '/ev', num2str(eqinfo.id(indi)), '.mat');
    % if(eqinfo.id(indi) ~= 653);continue;end  % For debug use

    disp(event_file);
    if (~exist(event_file,'file'))
        continue;
    end

    load(event_file);
    if (data.numts == 0)
        continue;
    end
    ustid = unique(data.stid);
    NumOfStn = length(ustid);
    for istn = 1:NumOfStn
        trace_all = strmatch(ustid(istn),data.stid);
        % % % [~,trace_nomultiple_ind,~] = unique(data.chan(trace_all));
        % % % trace_nomultiple = trace_all(trace_nomultiple_ind);
        %        	its = trace_final(1);        Modified by Jiewen Zhang,3/30/2017
        %if (strmatch(data.chan{its}(1:2),'EH') == 0 | ...
        %	strmatch(data.chan{its}(1:2),'HH') == 0 );
        % % % trace_acc = 1;
        % % % trace_final = [];
        % % % for itr = 1:length(trace_nomultiple)
        for itr = 1:numel(trace_all)
            its = trace_all(itr);
            % its = trace_nomultiple(itr);
            if (~strcmp(data.chan{its}(1:2),P.targchan))
                continue;
            end
            if (~strcmp(P.targnet,'all') && ~strcmp(data.net{its},P.targnet))  % network not matching, skip
                continue;
            end
            if (~strcmp(P.targsta,'all') && ~strcmp(data.sta{its},P.targsta))  % station not matching, skip
                continue;
            end
            % if (strcmp(data.net{its},'SF')==0 && sum(isstrprop(data.sta{its},'digit')) > 0)   % Station name includes digits like MH009, will exclude
            %     continue;
            % end
            if (strcmp(P.target,'P') && data.pick1(its) < 0)
                continue;
            end

            if (strcmp(P.target,'S') && data.pick2(its) < 0)
                continue;
            end
            % trace_final(trace_acc) = trace_nomultiple(itr);
            % trace_acc = trace_acc+1;
        end
        % if(length(trace_final) == 0); continue; end
        % for Ntrace = 1:length(trace_final)
        itsz = 0;
        itse = 0;
        itsn = 0;
        for Ntrace = 1:numel(trace_all)
            iTrace = trace_all(Ntrace);
            chtemp = data.chan{iTrace};
            if (strcmp(chtemp(3:3),'Z') == 1 || strcmp(chtemp(3:3),'1') == 1)
                itsz = iTrace;
            end
            if (strcmp(chtemp(3:3),'E') == 1 || strcmp(chtemp(3:3),'3') == 1)
                itse = iTrace;
            end
            if (strcmp(chtemp(3:3),'N') == 1 || strcmp(chtemp(3:3),'2') == 1)%& data.clip(iTrace) == 0)
                itsn = iTrace;
            end
        end
        % if (strcmp(P.target,'P') & itsz == 0)
        %     continue;
        % end
        % if (strmatch(P.target,'S')==1 & itse*itsn == 0)
        %     continue;
        % end
        % % 		[itsz data.numts]
        % %         	data
        % %if (data.samprate(itsz)<=90) continue; end

        % if (data.samprate(iTrace) < 150)

        if (data.pspecres(itsn, 1) == 0 || data.pspecres(itse, 1) == 0 || data.pspecres(itsz, 1) == 0)
            continue
        end
        if (strcmp(P.target,'P'))
            spectemp = sqrt(data.pspec(itsn,:).^2 + data.pspec(itse,:).^2 + data.pspec(itsz,:).^2);
            noisetemp = sqrt(data.pnspec(itsn,:).^2 + data.pnspec(itse,:).^2 + data.pnspec(itsz,:).^2);
            specres = sqrt(data.pspecres(itsn,:).^2 + data.pspecres(itse,:).^2 + data.pspecres(itsz,:).^2);
            freq = data.pfreq(itsz,:);
            ttime = mean(data.pick1(trace_all));
        elseif (strcmp(P.target,'S'))
            spectemp = sqrt(data.sspec(itse,:).^2 + data.sspec(itsn,:).^2);
            noisetemp = sqrt(data.pnspec(itsn,:).^2 + data.pnspec(itse,:).^2);
            specres = sqrt(data.sspecres(itse,:).^2 + data.sspecres(itsn,:).^2);
            %spectemp = data.sspec(itse,:);
            %noisetemp = data.snspec(itse,:);
            %specres = data.sspecres(itse,:);
            freq = data.sfreq(itse,:);
            ttime = mean(data.spred(trace_all));
        elseif (strcmp(P.target,'C'))  % for coda wave, DON'T USE YET
            spectemp = data.cspec(its,:);
            noisetemp = data.cnspec(its,:);
            specres = data.cspecres(its,:);
            freq = data.cfreq(itse,:);
            ttime = mean(data.spred(trace_all));
        end
        % else
        %     if (strmatch(P.target,'P')==1)  % done with P-wave, simple, just use vertical channel
        %         spectemp = (data.pspec(itsz,:));
        %         noisetemp = (data.pnspec(itsz,:));
        %         specres = (data.pspecres100(itsz,:))*D;
        %         freq = data.pfreq(itsz,:);
        %         ttime = mean(data.pick1(trace_final));
        %     elseif (strmatch(P.target,'S')==1 ) % will do rotation or use average??
        %         spectemp = (data.sspec(itse,:)+data.sspec(itsn,:))/2;
        %         noisetemp = (data.snspec(itse,:)+data.snspec(itsn,:))/2*D;
        %         specres = (data.sspecres100(itse,:)+data.sspecres100(itsn,:))/2;
        %         %spectemp = data.sspec(itse,:);
        %         %noisetemp = data.snspec(itse,:);
        %         %specres = data.sspecres(itse,:);
        %         freq = data.sfreq(itse,:);
        %         ttime = mean(data.spred(trace_final));
        %     elseif (strmatch(P.target,'C')==1)  % for coda wave
        %         spectemp = data.cspec(its,:);
        %         noisetemp = data.cnspec(its,:);
        %         specres = data.cspecres100(its,:)*D;
        %         freq = data.cfreq(its,:);
        %         ttime = mean(data.spred(trace_final));
        %     end
        % end

        %nnan = find(isnan(specres));
        %if (~isempty(nnan)) continue; end
        % now will test signal_to_noise ratio
        nbf1 = find(freq>=P.freq1(1) & freq<=P.freq1(2));
        stn1 = min(log10(spectemp(nbf1))-log10(noisetemp(nbf1)));

        nbf2 = find(freq>=P.freq2(1) & freq<=P.freq2(2));
        stn2 = min(log10(spectemp(nbf2))-log10(noisetemp(nbf2)));

        nbf3 = find(freq>=P.freq3(1) & freq<=P.freq3(2));
        stn3 = min(log10(spectemp(nbf3))-log10(noisetemp(nbf3)));

        nf_tstar = find(freq>=P.fmin_tstar & freq<=P.fmax_tstar);
        sumnoise = sum(log10(spectemp(nf_tstar))-log10(noisetemp(nf_tstar)));
        avgstn = sumnoise/length(nf_tstar);

        if (stn1<P.stnmin(1) || stn2<P.stnmin(2) || stn3<P.stnmin(3) ...
                || avgstn < P.minavgstn)
            continue;
        end
        % if (ttime > 50 || data.del(iTrace)>P.delmax)
        %     continue;
        % end

        % now will flaten the low frequenc part below 2 Hz, added on
        % 4/10/2017
        %{
            n1 = data.flog100<=1.2;
            n12 = data.flog100>1.2 & data.flog100<=2.;
            ratio = mean(specres(n12))/mean(specres(n1));
            if (ratio>=1.5 | ratio <= 0.66)
                n2 = find(data.flog100<=2);
                for iff = length(n2):-1:1
                    ii =  n2(iff);
                    i1 = ii - 2;
                    i2 = ii + 2;
                    if (i1<1) i1 = 1; i2 = i2+5; end
                    specres(ii) = mean(specres(i1:i2));
                end
                specres(n1) = mean(specres(n12));
            end
        %}
        specx = log10(specres);
        ninf = find(isinf(specx), 1);
        if (~isempty(ninf) || ~isreal(specx))
            continue;
        end
        % work on individual earthquake match
        nspec = nspec+1;

        spec(nspec).spec = specx;  % here I only use log spectrum
        spec(nspec).inq = n_select(iq);
        spec(nspec).qlat = data.qlat;
        spec(nspec).qlon = data.qlon;
        spec(nspec).qdep = data.qdep;
        try
            spec(nspec).qmag = data.qmag1;
        catch
            spec(nspec).qmag = data.mb;
        end
        spec(nspec).qtime = data.qtime;
        try
            spec(nspec).qid = data.cuspid;
        catch
            spec(nspec).qid = data.id;
        end
        % work on individual station match
        stname = data.sta{itsz};
        stid = blanks(10);
        stid(1:length(stname)) = stname;
        chanx = data.chan{itsz}(1:2);
        disp(['channel name:' chanx])
        stid(7:8) = chanx;
        if (nsta == 0)
            nsta = nsta+1;
            spec(nspec).ins = nsta;
            spec(nspec).stid = stid;
            starray{nsta} = stid;
        else

            nfound = strmatch(stid,starray);
            if (isempty(nfound)) % a new station
                nsta = nsta+1;
                spec(nspec).ins = nsta;
                spec(nspec).stid = stid;
                starray{nsta} = stid;
            else
                spec(nspec).ins = nfound(1);
                spec(nspec).stid = stid;
            end
        end
        spec(nspec).stlat = data.stlat(itsz);
        spec(nspec).stlon = data.stlon(itsz);
        spec(nspec).stelev = data.stelev(itsz);
        try
            spec(nspec).del = data.del(itsz);
        catch
            spec(nspec).del = data.dist(itsz);
        end

        % now work on individual distance bin match
        ix = round(ttime/P.dtt);
        if (ix<=1) ix = 1; end
        % if (ix>=30) ix = 30; end
        spec(nspec).inx = ix; % this is for distance bin
        spec(nspec).ttime = ttime;

        % depth index
        qdp = data.qdep;
        id = round(qdp+0.6);
        if (id<=1) id = 1; end
        if (id>=30) id = 30; end
        spec(nspec).ind = id;


        % mag index
        try
            qmag = data.qmag1;
        catch
            qmag = data.mb;
        end
        im = round(qmag*5.+0.6);
        if (im<1) im = 1; end
        if (im>40) im = 40; end
        spec(nspec).inm = im;

        flog = data.flog;

        % now do tstar for each spectrum
        x = freq(nf_tstar);
        y = log10(spectemp(nf_tstar));
        sumx = sum(x);
        sumy = sum(y);

        xmean = sumx/length(nf_tstar);
        ymean = sumy/length(nf_tstar);

        sumxx = sum((x-xmean).^2);
        sumxy = sum((x-xmean).*(y-ymean));
        %[isreal(sumxx) isreal(sumxy) xmean ymean sumx sumy length(nf_tstar)]
        %eqinfo.id(indi)
        b = sumxy/sumxx;
        spec(nspec).amp = ymean-b*xmean;
        spec(nspec).tstar = -b/1.364;
        if (isnan(b))
            spec(nspec).amp = 0;
            spec(nspec).tstar = 0;
        end
    end

end


