clear variables
clc

addpath('../src');

% Change eqfile, evdir and evdir_target accordingly in your case
eqfile = '../data/matfile/eqinfo.mat';
evdir = '../data/matfile';
evdir_target = [evdir '_spec'];
load(eqfile);
mkdir(evdir_target);

matcount = 1;

for iev = 1:length(eqinfo.id)
    % Parse your own matfile0, matpath, matfile accordingly to fit in your
    % database structure
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
    if(~exist(matpath, 'dir'))
        mkdir(matpath);
    end
    if(exist(matfile0, 'file'))
        % if (~exist(matfile, 'file'))
        copyfile(matfile0,matfile);
        % end
        disp(['Processing event ' evid ' ...'])
        update_data_spec(matfile, 0);
    else
        disp([evid ' does not exist! skipped (' num2str(matcount) '/' num2str(length(eqinfo.id)) ')'])

    end
    matcount = matcount + 1;
end


