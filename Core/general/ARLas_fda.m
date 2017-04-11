function [frequency,signal,noiseFloor,phase] = ARLas_fda(data,fs,ref,nfft)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [frequency,signal,noiseFloor] = ARLas_fda(data,fs,ref,nfft);
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
% Updated: March 27, 2017 ssg -- added phase as an output argument
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x = data;
[N,K] = size(x); % Input data are arranged into a matrix of N samples by K buffers.
if nargin == 3 % if fft size is not specified, use N
    nfft = N;
else % if nfft > N, fft will be zero padded to specified nfft
    if nfft < N
        error('nfft must be >= number of rows in data matrix')
    end
end
Xk = fft(x,nfft,1); % take the DFT of each buffer in the matrix (down columns)
Nyquist = [floor(nfft/2 + 1)]-1; % calculate the location of the Nyquist frequency
deltaF = fs / nfft; % DFT bandwidth
frequency = (0:deltaF:fs-deltaF)'; % frequency vector
frequency = frequency(1:Nyquist);
Xk = Xk(1:Nyquist,:); % discard aliased portion above Nyquist rate
Xbar = (1/K) * sum(Xk,2); % signal is the coherent mean
Xbar2 = abs(Xbar) .^2; % signal energy
phase = angle(Xbar); % signal phase
if K > 1,
    XBAR = repmat(Xbar,1,K); % replicate to matrix (to avoid running a loop)
    S2 = (1/(K-1)) * sum((Xk - XBAR) .* conj(Xk - XBAR),2); % variance (this is the noise floor)
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
Xbar2 = Xbar2 * (2/N^2); % scale magnitude by N (NOT by nfft; do not include zero padding when scaling)
if ~isempty(S2),
    S2 = S2 * (2/N^2); % scale noise
    Se2 = Se2 * (2/N^2);
end
signal = 10*log10(Xbar2/(ref^2));
noise = 10*log10(S2/(ref^2));    
noiseFloor = 10*log10(Se2/(ref^2));    


