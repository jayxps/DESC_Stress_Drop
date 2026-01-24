function [t_cut,data_cut,SAChdr] = fget_sac_cut(filename,tstart,tend)
%[t,data,SAChdr] = fget_sac(filename)

% read sac into matlab
% written by Zhigang Peng
% program called
% [head1, head2, head3, data]=sac(filename);
% [SAChdr]=sachdr(head1, head2, head3);

% Updated Mon Jul 30 11:21:24 PDT 2001
% This version differs from fget_sac in that in ask
% you to input start and end time to cut the data.

if nargin <3, tend = 1000000; end
if nargin <2, tstart = -1000000; end
if nargin <1, error('ERROR!! No input file name'); end

[head1, head2, head3, data]=sac(filename);
[SAChdr]=sachdr(head1, head2, head3);
t = [SAChdr.times.b:SAChdr.times.delta:(SAChdr.data.trcLen-1)*SAChdr.times.delta+SAChdr.times.b]';
tmin = min(t);
tmax = max(t);
if tstart<tmin, tstart=tmin; end
if tend>tmax, tstart=tmin; end

t_cut = t(t>=tstart & t<=tend);
data_cut = data(t>=tstart & t<=tend);
SAChdr.times.npts = length(data_cut);
