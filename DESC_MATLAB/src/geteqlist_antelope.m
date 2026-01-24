function [eqinfo] = geteqlist_antelope(filename,itype)

if (itype==1) %scec format
	eqlist = getscec(filename);
elseif (itype==2)
	eqlist = gethk(filename);
elseif (itype==3)
    eqlist = getefs(filename);
elseif (itype==4)
    eqlist = getncsn(filename);
elseif (itype==5)
    eqlist = getjma(filename);
elseif (itype==8)
    eqlist = getxyz(filename);
elseif (itype==6) % hypoDD loc file
    eqlist = getdd(filename);
elseif (itype==7) % anss file
    eqlist = getanss(filename);
end

eqinfo(:,1) = eqlist(:,7); %lat
eqinfo(:,2) = eqlist(:,8); %lon
eqinfo(:,3) = eqlist(:,9); %dep
eqinfo(:,4) = eqlist(:,10);%mag
datex = datenum(eqlist(:,1:6));
eqinfo(:,5) = date2epoch(datex);
%% to verify
%	fprintf(1,'%d %d %d %d %d %f\n',eqlist(1,1:6))
%	timenum = epoch2date(eqinfo(1,5))
%	timex = datevec(timenum)
%	fprintf(1,'%d %d %d %d %d %f\n',datevec(datenum(eqlist(1,1:6))))
%for i = 1:neq
%   eqinfo(i,5) = str2epoch(strcat([num2str(eqlist(i,1)),'/',...
%       num2str(eqlist(i,2)),'/',num2str(eqlist(i,3)),'/',...
%       num2str(eqlist(i,4)),':',num2str(eqlist(i,5)),':',...
%       num2str(eqlist(i,6))]));  
%end
eqinfo(:,6) = eqlist(:,11);

function uOutput = getjma(sFilename)
   mData = load(sFilename);
   [M, N ] = size(mData);
   if (N == 10) % then evid exist
       uOutput(:,11) = 1:M;
   elseif (N >= 11)
       uOutput(:,11) = mData(:,11);
   end

   uOutput(:,1:10) = mData(:,1:10);

function uOutput = getxyz(sFilename)
   mData = load(sFilename);
   [M, N ] = size(mData)
   if (N == 10) % then evid exist
       uOutput(:,11) = 1:M;
   elseif (N >= 11)
       uOutput(:,11) = mData(:,11);
   end
   uOutput(:,1:6) = mData(:,1:6);
   uOutput(:,10) = mData(:,7);
   uOutput(:,7:9) = mData(:,8:10);
   %uOutput(:,1:10) = mData(:,1:10);

function uOutput = getscec(sFilename)
    mData = textread(sFilename, '%s', 'delimiter', '\n', 'whitespace', '');
    % Create empty catalog
    uOutput = zeros(length(mData), 10); 
    % Loop through all lines of catalog and convert them
     mData = char(mData);
     l = find( mData == ' ' );
     mData(l) = '0';
     errc = 1;
    try % this is the fast method .... 
        i = 1:length(mData(:,1));   
        [uOutput] = readvaluesscec(mData,i,i);
         
    catch  % this is the line by line method .. 
        for i = 1:length(mData(:,1))
%             if rem(i,100) == 0 ; disp([ num2str(i) ' of ' num2str(length(mData)) ' events processed ']); end          
            try      
                uOutput(i,:) = readvaluesscec(mData,i,1); 
            catch
                errc = errc + 1; 
                if errc == 100
                    if stoploop 
                        return
                    end 
                end      
                disp(['Import: Problem in line ' num2str(i) ' of ' sFilename '. Line ignored.']);
                uOutput(i,:) = uOutput(i,:)*nan;
            end
        end
        
    end
    l = isnan(uOutput(:,1));
    uOutput(l,:) = [];
    nx = find(uOutput(:,7)>100);


function [uOutput] = readvaluesscec(mData,i,k)

% the following is the new SCEC format, but when there is negative depth,
% this will be wrong
uOutput(k,1) = str2num(mData(i,1:4));    % Year
uOutput(k,2) = str2num(mData(i,6:7));    % Month
uOutput(k,3) = str2num(mData(i,9:10));   % Day
uOutput(k,4) = str2num(mData(i,12:13));  % Hour
uOutput(k,5) = str2num(mData(i,15:16));  % Minute
uOutput(k,6) = str2num(mData(i,18:23));  % seconds
uOutput(k,7) = str2num(mData(i,36:43));  % Latitude
uOutput(k,8) = str2num(mData(i,44:51)) ;  % Longitude
uOutput(k,9) = str2num(mData(i,54:58));  % Depth
uOutput(k,10) = str2num(mData(i,27:30));  % Magnitude
uOutput(k,11) = str2num(mData(i,60:67));  % Minute

function uOutput = getefs(sFilename)
    mData = textread(sFilename, '%s', 'delimiter', '\n', 'whitespace', '');
    % Create empty catalog
    uOutput = zeros(length(mData), 10); 
    % Loop through all lines of catalog and convert them
     mData = char(mData);
     l = find( mData == ' ' );
     mData(l) = '0';
     errc = 1;
    try % this is the fast method .... 
        i = 1:length(mData(:,1));   
        [uOutput] = readvaluesefs(mData,i,i);
         
    catch  % this is the line by line method .. 
        for i = 1:length(mData(:,1))
%             if rem(i,100) == 0 ; disp([ num2str(i) ' of ' num2str(length(mData)) ' events processed ']); end          
            try      
                uOutput(i,:) = readvaluesefs(mData,i,1); 
            catch
                errc = errc + 1; 
                if errc == 100
                    if stoploop 
                        return
                    end 
                end      
                disp(['Import: Problem in line ' num2str(i) ' of ' sFilename '. Line ignored.']);
                uOutput(i,:) = uOutput(i,:)*nan;
            end
        end
        
    end
    l = isnan(uOutput(:,1));
    uOutput(l,:) = [];
    nx = find(uOutput(:,7)>100);


function [uOutput] = readvaluesefs(mData,i,k)

% the following is the new SCEC format, but when there is negative depth,
% this will be wrong
uOutput(k,1) = str2num(mData(i,10:13));    % Year
uOutput(k,2) = str2num(mData(i,15:16));    % Month
uOutput(k,3) = str2num(mData(i,18:19));   % Day
uOutput(k,4) = str2num(mData(i,21:22));  % Hour
uOutput(k,5) = str2num(mData(i,24:25));  % Minute
uOutput(k,6) = str2num(mData(i,27:32));  % seconds
uOutput(k,7) = str2num(mData(i,46:53));  % Latitude
uOutput(k,8) = str2num(mData(i,55:64)) ;  % Longitude
uOutput(k,9) = str2num(mData(i,68:72));  % Depth
uOutput(k,10) = str2num(mData(i,37:41));  % Magnitude
uOutput(k,11) = str2num(mData(i,1:8));  % scecid


function uOutput = gethk(sFilename)
    mData = textread(sFilename, '%s', 'delimiter', '\n', 'whitespace', '');
    % Create empty catalog
    uOutput = zeros(length(mData),11);
    Error_lines = []; 
    % Loop through all lines of catalog and convert them
    mData = char(mData);
    l = find( mData == ' ' );
    mData(l) = '0';
    errc = 1;
    try % this is the fast method .... 
        i = 1:length(mData(:,1));   
        [uOutput] = readvalueshk(mData,i,i);
         
    catch  % this is the line by line method .. 
        for i = 1:length(mData(:,1))
            if rem(i,100) == 0 ; disp([ num2str(i) ' of ' num2str(length(mData)) ' events processed ']); end          
            try      
                uOutput(i,:) = readvalueshk(mData,i,1); 
            catch
                errc = errc + 1; 
                if errc == 100
                    if stoploop 
                        return
                    end 
                end      
                disp(['Import: Problem in line ' num2str(i) ' of ' sFilename '. Line ignored.']);
                uOutput(i,:) = uOutput(i,:)*nan;
                Error_lines = [Error_lines; i];
            end
        end
        
    end
    l = isnan(uOutput(:,1));
    uOutput(l,:) = [];
    
function [uOutput] = readvalueshk(mData,i,k)

uOutput(k,1) = str2num(mData(i,(1:4)));    % Year
uOutput(k,2) = str2num(mData(i,(6:7)));    % Month
uOutput(k,3) = str2num(mData(i,(9:10)));   % Day
uOutput(k,4) = str2num(mData(i,(12:13)));  % Hour
uOutput(k,5) = str2num(mData(i,(15:16)));  % Minute
uOutput(k,6) = str2num(mData(i,(18:23)));  % Second
uOutput(k,7) = str2num(mData(i,(35:42)));  % Latitude [deg]
uOutput(k,8) = str2num(mData(i,(44:53)));  % Longitude [deg]
uOutput(k,9) = str2num(mData(i,(56:60)));  % Depth [km]
uOutput(k,10) = str2num(mData(i,(63:68)));  % Magnitude
uOutput(k,11) = str2num(mData(i,(25:33))); % SCSN cuspid (up to 9 digits)


function uOutput = getncsn(sFilename)
    mData = textread(sFilename, '%s', 'delimiter', '\n', 'whitespace', '');
    % Create empty catalog
    uOutput = zeros(length(mData), 10); 
    % Loop through all lines of catalog and convert them
     mData = char(mData);
     l = find( mData == ' ' );
     mData(l) = '0';
     errc = 1;
    try % this is the fast method .... 
        i = 1:length(mData(:,1));   
        [uOutput] = readvaluesncsn(mData,i,i);
         
    catch  % this is the line by line method .. 
        for i = 1:length(mData(:,1))
%             if rem(i,100) == 0 ; disp([ num2str(i) ' of ' num2str(length(mData)) ' events processed ']); end          
            try      
                uOutput(i,:) = readvaluesncsn(mData,i,1); 
            catch
                errc = errc + 1; 
                if errc == 100
                    if stoploop 
                        return
                    end 
                end      
                disp(['Import: Problem in line ' num2str(i) ' of ' sFilename '. Line ignored.']);
                uOutput(i,:) = uOutput(i,:)*nan;
            end
        end
        
    end
    l = isnan(uOutput(:,1));
    uOutput(l,:) = [];
    nx = find(uOutput(:,7)>100);


function [uOutput] = readvaluesncsn(mData,i,k)

% the following is the new SCEC format, but when there is negative depth,
% this will be wrong
uOutput(k,1) = str2num(mData(i,1:4));    % Year
uOutput(k,2) = str2num(mData(i,6:7));    % Month
uOutput(k,3) = str2num(mData(i,9:10));   % Day
uOutput(k,4) = str2num(mData(i,12:13));  % Hour
uOutput(k,5) = str2num(mData(i,15:16));  % Minute
uOutput(k,6) = str2num(mData(i,18:22));  % seconds
uOutput(k,7) = str2num(mData(i,25:33));  % Latitude
uOutput(k,8) = str2num(mData(i,34:44)) ;  % Longitude
uOutput(k,9) = str2num(mData(i,46:52));  % Depth
uOutput(k,10) = str2num(mData(i,54:58));  % Magnitude
uOutput(k,11) = str2num(mData(i,90:97));  % ncsn id

function uOutput = getdd(sFilename)
    mData = load(sFilename);
    uOutput(:,1) = mData(:,11); %year
    uOutput(:,2) = mData(:,12);
    uOutput(:,3) = mData(:,13);
    uOutput(:,4) = mData(:,14);
    uOutput(:,5) = mData(:,15);
    uOutput(:,6) = mData(:,16); %second
    uOutput(:,7) = mData(:,2); %lat
    uOutput(:,8) = mData(:,3);
    uOutput(:,9) = mData(:,4); %dep
    uOutput(:,10) = mData(:,17); %mag
    uOutput(:,11) = mData(:,1); %id

    
function uOutput = getanss(sFilename)
mData = textread(sFilename, '%s', 'delimiter', '\n', 'whitespace', '');
    % Create empty catalog
    uOutput = zeros(length(mData), 11); 
    % Loop through all lines of catalog and convert them
     mData = char(mData);
     l = find( mData == ' ' );
     mData(l) = '0';
     errc = 1;
    try % this is the fast method .... 
        i = 1:length(mData(:,1));   
        [uOutput] = readvaluesanss(mData,i,i);
         
    catch  % this is the line by line method .. 
        for i = 1:length(mData(:,1))
%             if rem(i,100) == 0 ; disp([ num2str(i) ' of ' num2str(length(mData)) ' events processed ']); end          
            try      
                uOutput(i,:) = readvaluesanss(mData,i,1); 
            catch
                errc = errc + 1; 
                if errc == 100 
                    if stoploop 
                        return
                    end 
                end      
                disp(['Import: Problem in line ' num2str(i) ' of ' sFilename '. Line ignored.']);
                uOutput(i,:) = uOutput(i,:)*nan;
            end
        end
        
    end
    l = isnan(uOutput(:,1));
    uOutput(l,:) = [];
    nx = find(uOutput(:,7)>100);


function [uOutput] = readvaluesanss(mData,i,k)

% the following is the new SCEC format, but when there is negative depth,
% this will be wrong
uOutput(k,1) = str2num(mData(i,1:4));    % Year
uOutput(k,2) = str2num(mData(i,6:7));    % Month
uOutput(k,3) = str2num(mData(i,9:10));   % Day
uOutput(k,4) = str2num(mData(i,12:13));  % Hour
uOutput(k,5) = str2num(mData(i,15:16));  % Minute
uOutput(k,6) = str2num(mData(i,18:22));  % seconds
uOutput(k,7) = str2num(mData(i,23:31));  % Latitude
uOutput(k,8) = str2num(mData(i,32:41)) ;  % Longitude
uOutput(k,9) = str2num(mData(i,42:48));  % Depth
uOutput(k,10) = str2num(mData(i,49:54));  % Magnitude
uOutput(k,11) = str2num(mData(i,84:96));  % ncsn id

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%%    DO NOT CHANGE %%%%%%%%%%%

function [mystop] = stoploop()

ButtonName=questdlg('More than 100 lines could not be read. Continue?', ...
    'Interrupt?', ...
    'Yes','No','Nope');

switch ButtonName
case 'Yes'
    disp('going on');
    mystop = 0;
case 'No'
    mystop = 1;
end % switch
