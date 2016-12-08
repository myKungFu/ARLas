function [frequency,signal,noiseFloor] = ARLas_fda(data,fs,ref)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [frequency,signal,noiseFloor] = ARLas_fda(data,fs,ref);
%
% Frequency-domain analysis routine for use with ARLas (Auditory Research Laboratory auditory software)
% data = a matrix of recorded data
% fs = sampling rate in Hz
% ref = reference: usually 0.00002 Pa or 1 V
%
% Auditory Research Lab, The University of Iowa
% Deptartment of Communication Sciences & Disorders
% The University of Iowa
% Author: Shawn S. Goodman, PhD
% Date: October 27, 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x = data;
[N,K] = size(x); % Input data are arranged into a matrix of N samples by K buffers.
Xk = fft(x,N,1); % take the DFT of each buffer in the matrix (down columns)
Nyquist = [floor(N/2 + 1)]-1; % calculate the location of the Nyquist frequency
deltaF = fs / N; % DFT bandwidth
frequency = (0:deltaF:fs-deltaF)'; % frequency vector
frequency = frequency(1:Nyquist);
Xk = Xk(1:Nyquist,:); % discard aliased portion above Nyquist rate
Xbar = (1/K) * sum(Xk,2); % signal is the coherent mean
Xbar2 = abs(Xbar) .^2; % signal energy
phase = angle(Xbar); % signal phase
if K > 1,
    XBAR = repmat(Xbar,1,K); % replicate to matrix (to avoid running a loop)
    S2 = (1/(K-1)) * sum((Xk - XBAR) .* conj(Xk - XBAR),2); % variance (this is the HA noise floor)
        % Se = sqrt((1/K) * S2); % standard error of mean (calculation not necessary)
    Se2 = (1/K) * S2; % energy of the standard error
else
    S2 = [];
    Se2 = [];
end
% scale dc and Nyquist properly
Xbar2(1) = Xbar2(1) * 0.5;
Xbar2(end) = Xbar2(end) * 0.5;
if ~isempty(S2),
    S2(1) = S2(1) * 0.5;
    S2(end) = S2(end) * 0.5;
    Se2(1) = Se2(1) * 0.5;
    Se2(end) = Se2(end) * 0.5;
end
Xbar2 = Xbar2 * (2/N^2);
if ~isempty(S2),
    S2 = S2 * (2/N^2);
    Se2 = Se2 * (2/N^2);
end
signal = 10*log10(Xbar2/(ref^2));
noise = 10*log10(S2/(ref^2));    
noiseFloor = 10*log10(Se2/(ref^2));    


