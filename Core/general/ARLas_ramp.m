function [data] = ARLas_ramp(data,fs,len)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% data = ARLas_ramp(data,fs,len);
%
% Onset/offset ramp routine for use with ARLas (Auditory Research Laboratory auditory software)
% data = a matrix of recorded data
% fs = sampling rate in Hz
% cutoff = desired length of the ramps, in seconds
%
% Auditory Research Lab, The University of Iowa
% Deptartment of Communication Sciences & Disorders
% The University of Iowa
% Author: Shawn S. Goodman, PhD
% Date: October 27, 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N = round(fs * len);
h = hann(N * 2);
ramp = h(1:N);
[~,m] = size(data);
Ramp = repmat(ramp,1,m);
data(1:N,:) = data(1:N,:) .* Ramp;

h = hann(N * 2);
ramp = h(N+1:end);
Ramp = repmat(ramp,1,m);
data(end-N+1:end,:) = data(end-N+1:end,:) .* Ramp;
