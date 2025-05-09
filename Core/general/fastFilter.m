function [Y] = fastFilter(B,X)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Y = fastFilter(B,X);
% 
% B = the impulse response of the filter. It is assumed that the
%          filter order is even. If not, one extra zero is added to the
%          end of the impulse response so that group delay can be
%          corrrected for.
% X = the waveform(s) to be filtered.  If data are a matrix, the data are
%          assumed to be in columns.
% Y = filtered data.
%
% This function implements the Matlab function fftfilt to perform
% efficient filtering of long sequences with relatively short filters.
% Used alone, fftfilt truncates the filtered sequence to ensure the same
% length as the original sequence.  However, this does not properly correct
% for group delay.  
%
% This function (fastFilter) assumes an FIR filter with a symmetric impulse response.
% The group delay is then constant and equals half the filter order divided by 
% the sampling rate.  The group delay can be accounted for by
% the proper application of zero padding before applying fftfilt, and then 
% removal of the padding afterwards.  Note that the filter order must be
% even in order for this to work.  
% 
% Author: Shawn Goodman
% Date: July 24, 2008
% Updated: September 11, 2011; ssg
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Y = [];
M = size(B,1) -1;
if mod(M,2) ~= 0
    B = [B;zeros(1,size(B,2))];
end
pad = zeros(M/2,size(X,2));
X = [X;pad];
Y = fftfilt(B,X);
Y = Y(size(pad,1)+1:end,:);
