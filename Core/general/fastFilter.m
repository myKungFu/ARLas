function [Y] = fastFilter(B,X)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Y = fastFilter(B,X);
%
% This function implements the Matlab function fftfilt to perform
% efficient filtering of long sequences with relatively short FIR filters.
% Used alone, the Matlab function fftfilt truncates the filtered sequence to 
% ensure the same output length as the original sequence.  However, this does
% not properly correct for group delay.  
%
% This function (fastFilter) assumes an FIR filter with a symmetric impulse 
% response (in other words, linear phase).
% The group delay is then constant and equals half the filter order divided by 
% the sampling rate.  The group delay can be accounted for by
% the proper application of zero padding before applying fftfilt, and then 
% removal of the padding afterwards.  
%
% IMPORTANT: The filter order must be even in order for this to work!
% 
% B = the impulse response of the filter. Should be a column vector.
%          It is assumed that the filter order is even. If not, one extra zero is added to the
%          end of the impulse response so that group delay can be
%          corrrected for. However, this "cheat" results in some inaccuracy
%          and should be avoided.
% X = the waveform(s) to be filtered. Data should be in columns.
%          If data are a matrix, the data are assumed to be in columns and
%          the filter will be applied to each column separately
% Y = the filtered data, corrected for group delay.
%
% 
% Author: Shawn Goodman, PhD.
% Date: July 24, 2008
% Updated: September 11, 2011; ssg
% Updated: November 9, 2020; ssg
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    Y = [];
    M = size(B,1) -1; % FIR filter order is filter length -1
    if mod(M,2) ~= 0  % This is a hack in the case that milter order is odd.
        B = [B;zeros(1,size(B,2))]; % should avoid this if at all possible!
    end
    pad = zeros(M/2,size(X,2)); % create a zero padded signal
    X = [X;pad]; 
    Y = fftfilt(B,X); % apply Matlab's fast fftfilt function
    Y = Y(size(pad,1)+1:end,:); % strip the zero padding off the end.

end