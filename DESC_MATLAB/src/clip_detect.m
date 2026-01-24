function [clip] = clip_detect(data, F)
% % Simple function to check for clipped waveforms.
% % Computes the maximum and minimum value in data. Then counts the number
% % of data points within a fraction F of maxval and minval:
% % (data >= (1-F)*maxval) --> NF1max
% % (data <= (1-F)*minval) --> NF1min
% % If both of these are counts are greater than the the counts of the next
% % fraction down:
% % ( (1-2*F)*maxval <= data(i) < (1-F)*maxval ) --> NF2max
% % ( (1-2*F)*minval >= data(i) > (1-F)*minval ) --> NF2min
% % then flag = 1.0 (likely clip). If only one count is greater, set flag = 0.5
% % (possible clipping). If none, set flag = 0.0 (no clipping).

% default value of F is 1/3;
if nargin < 2 || abs(F) > 0.5
    F = 1.0/3.0;
end

% get max and min values
maxval = max(data); minval = min(data);

% indices of data points in different intervals
I1max = data >= (1-F)*maxval;
I2max = data >= (1-2*F)*maxval;
I1min = data <= (1-F)*minval;
I2min = data <= (1-2*F)*minval;

% counts of data points in different intervals
NF1max = sum(I1max);
NF2max = sum(I2max & ~I1max);
NF2min = sum(I2min & ~I1min);
NF1min = sum(I1min);

% classify data
%     if (NF2max < NF1max) && (NF2min < NF1min) % likely clipped
%         clip = 1;
%     elseif (NF2max < NF1max) || (NF2min < NF1min) % maybe clipped
%         clip = 0.5;
%     else % no worries
%         clip = 0;
%     end

% compute clip_val, modified by Jiewen Zhang on Jun 6, 2020
clip = (NF1max/NF2max+NF1min/NF2min)/2;

end