function [data] = ARLas_hpFilter(data,fs,cutoff)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% data = ARLas_hpFilter(data,fs,cutoff);
%
% Highpass filtering routine for use with ARLas (Auditory Research Laboratory auditory software)
% data = a matrix of recorded data
% fs = sampling rate in Hz
% cutoff = highpass cutoff in Hz
%
% Auditory Research Lab, The University of Iowa
% Deptartment of Communication Sciences & Disorders
% The University of Iowa
% Author: Shawn S. Goodman, PhD
% Date: October 27, 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% The data are recorded with no software filters. As a result, there is
% usually a lot of low-frequency noise, and HP filtering is suggested.
% A few HP FIR filters with commonly-used cutoffs are saved in 
% ARLas/Core/general/, or create and apply your own custom filters. 
% Note that if you use the provided function fastFilter.m, group delay
% of the filter will be corrected for automatically; however, it is
% important that the filter must be FIR, must be saved as a column vector, 
% and group delay in samples must be an even number (odd length).
% 
% To avoid edge discontinuities, the following filter method is suggested

[R,C] = size(data); % get the original size
b = getHPfilter(fs,cutoff);
data = fastFilter(b,data(:)); % filter as one long vector
data = reshape(data,R,C); % put back into matrix form

% subfunctions ------------------------------------------------------------
function b = getHPfilter(fs,cutoff)
Fstop = 0;             % Stopband Frequency
Fpass = cutoff;             % Passband Frequency
Dstop = 0.0001;          % Stopband Attenuation
Dpass = 0.05;  % Passband Ripple
flag  = 'scale';         % Sampling Flag
% Calculate the order from the parameters using KAISERORD.
[N,Wn,BETA,TYPE] = kaiserord([Fstop Fpass]/(fs/2), [0 1], [Dpass Dstop]);
if mod(N,2)~=0 % force order to be even (length odd)
    N = N+1;
end
b  = fir1(N, Wn, TYPE, kaiser(N+1, BETA), flag)'; % create and transpose into column

