function  para_EGF = specprocess_sourcepara(rdir,freqfit,evspec,...
    imethod,iwave, EGF,min_nspec,isavefig)
% this is only for constant EGF!
% individual pair will  use a different program to use individual stations

if (~exist(rdir,'dir'))
    mkdir(rdir);
end
if (rdir(end:end)~='/') rdir(end+1) = '/'; end

freq = evspec(1).freq;
frqlim=[0.5 50];
flow = [0.5 5];
fhigh= [20 30];
nlow = find(freq>=flow(1) & freq<=flow(2));
nhigh= find(freq>=fhigh(1) & freq<=fhigh(2));
nfit = find(freq>=freqfit(1) & freq<=freqfit(2));

fid4 = fopen(strcat(rdir,'result-constant-egf.txt'),'w');
format4 = [' %5d','%8.2f','%8.2f','%15d',...
    '%15.2f','%15.2f','%15.2e','%10.3e','%10.3e',...
    '%11.3e','%10.3e','%10.3e','%10.3e','%10.3e','%10.3e %10.3e' ...
    ' %15.2f %15.2f \n'];

header4 = ['    i1   magest magcat              id1  ',...
    '     o0             fc1          r    delsig    asig    ',...
    '     G          Es-all   Es-model   ',...
    'Es-obs   Es-1       Es-2   error fcbound1 fcbound2 \n'];
fprintf(fid4,header4);

im = imethod;
%if (imethod == 2) im = 1; end
%if (imethod == 3) im = 5; end
nq = length(evspec);

nqsolve = 0;
for iq = 1:nq
    if (evspec(iq).nspec >= min_nspec) % then will solve using constant-EGF
        dspecc = evspec(iq).spec - EGF;

        % compute a few quality parameters
        if ((max(dspecc)-min(dspecc))>1e-2) % exist
            iexist = 1;
            nqsolve = nqsolve+1;
            [fc,fc2,o0,Es,delsig,asig,r,G,err,dspecfit,iflag,ftest,etest,fcb1,fcb2] = ...
                sourcepara (dspecc,freq,nfit,nlow,im,evspec(iq).qmomest,iwave);
            lhratio = mean(dspecc(nlow))-mean(dspecc(nhigh));
            para_EGF(nqsolve).specfit = dspecfit;
            para_EGF(nqsolve).ftest = ftest;
            para_EGF(nqsolve).etest = etest;
            para_EGF(nqsolve).exist = 1;
            para_EGF(nqsolve).fc = fc;
            para_EGF(nqsolve).fc2 = fc2;
            para_EGF(nqsolve).err = err;
            para_EGF(nqsolve).Es = Es;
            para_EGF(nqsolve).delsig = delsig;
            para_EGF(nqsolve).asig = asig;
            para_EGF(nqsolve).fcb1 = fcb1;
            para_EGF(nqsolve).fcb2 = fcb2;
            para_EGF(nqsolve).qid = evspec(iq).qid;
            para_EGF(nqsolve).qlat = evspec(iq).qlat;
            para_EGF(nqsolve).qlon = evspec(iq).qlon;
            para_EGF(nqsolve).qdep = evspec(iq).qdep;
            para_EGF(nqsolve).qmag = evspec(iq).qmagest;
            para_EGF(nqsolve).qmom = evspec(iq).qmomest;
            para_EGF(nqsolve).qmagcat = evspec(iq).qmag;
            para_EGF(nqsolve).qtime = evspec(iq).qtime;
            para_EGF(nqsolve).lhratio = lhratio;
            para_EGF(nqsolve).G = G;
            para_EGF(nqsolve).r = r;
            para_EGF(nqsolve).nspec = evspec(iq).nspec;

  	      if (mod(nqsolve,100) == 0)
  	      	fprintf(1,'event: %d, id: %d \n',nqsolve,evspec(iq).qid);
          end
          if (min(Es)>0 && iflag == 1 && fc<100)
              fprintf(fid4,format4,iq,evspec(iq).qmagest, evspec(iq).qmag,evspec(iq).qid,...
                  o0,fc,r,delsig,asig,...
                  G,Es(1:5),err, fcb1,fcb2);
              if (iexist == 1 & isavefig == 1 ) % just plot egf-stack fit
                  h1=figure(1);
                  set(h1,'visible','off');
                  subplot(2,2,1)
                  semilogx(freq,dspecc);
                  hold on
                  semilogx(freq,dspecfit,'-.');
                  hold off
                  xlabel('frequency');
                  xlim(frqlim);
                  ylabel('amplitude');

                  subplot(2,2,2);
                  plot(ftest,etest/min(etest));
                  ytest = min(etest)*1.05*[1 1]/min(etest);
                  hold on
                  plot([min(ftest) max(ftest)],ytest,'-.');
                  plot(fc,err/min(etest),'*');
                  hold off
                  xlabel('frequency');
                  ylabel('variance');
                  xlim([min(ftest) max(ftest)]);
                  title(strcat('shearer-stacking-',num2str(evspec(iq).qid)));
                  fileout = strcat(rdir,'fig-shearer-',num2str(evspec(iq).qid));
      			try
                    export_fig(fileout,'-pdf','-painters','-transparent');
    			catch
    				continue;
    			end
                close(1);
              end
          end
        end
    end
end
