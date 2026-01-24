function filtdata = eqfiltfilt(data,ch1,ch2,dt,ih)
%filtdata = eqfiltfilt(data,ch1,ch2,dt,ih)
%	Function to perform butterworth	filter on "data"
%	vector,	and return filtered data vector	"filtdata"
%
%	filter corners (above "ch1" and	below "ch2" can	pass, in Hz)
%	  if ch1 < 0,	    % LOW-PASS FILTER  (-inf,ch2)
%	  elseif ch2 < 0,   % HIGH-PASS	FILTER (ch1,inf)
%	  elseif ch1 > ch2, % BAND-STOP	FILTER (-inf,ch2) & (ch1,inf)
%	  else		    % BAND-PASS	FILTER (ch1,ch2)
%	sampling rate (dt, in sec/sample)
%	butterworth filter order (ih, default =	4)

%
% Last update 5/2/97, Cheng-Ping Li
% Modified by Zhigang Peng, Mon Oct 14 17:04:30 PDT 2002

if nargin <= 4,	ih=4; end
if nargin <= 3,	dt=1; end
if nargin <= 2
    error(['Not enough input argument, (at least 3,	data ch1, and ch2)']);
end
samp = 1/dt;		% sampling frequency
ch1n = ch1*2/samp;
ch2n = ch2*2/samp;

if ch1 < 0			% LOW-PASS FILTER
    [b,a]	= butter(ih,ch2n);
elseif ch2 < 0			% HIGH-PASS FILTER
    [b,a]	= butter(ih,ch1n,'high');
elseif ch1 > ch2		% BAND-STOP FILTER
    [b,a]	= butter(ih,[ch2n ch1n],'stop');
else				% BAND-PASS FILTER
    [b,a]	= butter(ih,[ch1n ch2n]);
end

filtdata = filtfilt(b,a,data); % zero phase filter the data
