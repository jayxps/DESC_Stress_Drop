clear variables
clc

% Specify number of threads to use for parallel processing
% You need to confirm how many cores are available on your computer
threadCnt = 8;
try
    parpool(threadCnt);
end

addpath('../src');

eqfile = '../data/matfile/eqinfo.mat';
evdir = '../data/matfile';
evdir_target = [evdir '_spec'];

load(eqfile);

mkdir(evdir_target);
% Sort events into different parallel sessions
EQ_par_sections = round(linspace(0,length(eqinfo.id),threadCnt+1));
bound_labs = EQ_par_sections;

spmd

    lab_Nmat = bound_labs(spmdIndex+1) - bound_labs(spmdIndex);
    matcount = 1;
    matcorrupt = [];
    for iev = bound_labs(spmdIndex)+1:bound_labs(spmdIndex+1)
        tic
        % Parse your own matfile0, matpath and matfile based off your database
        % structure.
        evid = num2str(eqinfo.id(iev));
        evyr = num2str(eqinfo.qyr(iev));
        if(eqinfo.qmon(iev) >= 10)
            evmon = num2str(eqinfo.qmon(iev));
        else
            evmon = ['0' num2str(eqinfo.qmon(iev))];
        end
        matfile0 = [evdir '/' evyr '/' evmon '/ev' evid '.mat'];
        matpath = [evdir_target '/' evyr '/' evmon];
        matfile = [evdir_target '/' evyr '/' evmon '/ev' evid '.mat'];
        if(~exist(matpath,'dir'))
            mkdir(matpath);
        end
        if(exist(matfile0,'file'))
            copyfile(matfile0,matfile);
            if(exist(matfile,'file'))
                disp(['Computing spectra for event ' evid '... (' num2str(matcount) '/' num2str(lab_Nmat) ')'])
                update_data_spec(matfile,0);
            else
                disp([evid ' does not exist! skipped (' num2str(matcount) '/' num2str(lab_Nmat) ')'])
            end
            matcount = matcount + 1;
        end
        toc
    end

end

