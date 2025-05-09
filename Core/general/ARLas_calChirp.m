function [stimulus] = ARLas_calChirp(samplingRate,fmin,fmax)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [stimulus] = ARLas_calChirp(samplingRate,fmin,fmax);
%
% Calibration chirp for use with ARLas (Auditory Research Laboratory auditory software)
% samplingRate = sampling rate (Hz)
% fmin = the low frequency (starting frequency) of the chirp
% fmax = the high frequency (ending frequency) of the chirp
%
% Auditory Research Lab, The University of Iowa
% Deptartment of Communication Sciences & Disorders
% The University of Iowa
% Author: Shawn S. Goodman, PhD
% Date: October 31, 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

desiredGain_dB = -30; % output = 30 db down from full out
% sweep rate is constant at 100 Hz/ms

bw = fmax - fmin; % total bandwidth in Hz
len = bw / 50 / 1000; % stimulus length in seconds
N = round(len * samplingRate); % total number of samples
if mod(N,2)~=0 % force to be an even number of samples
    N = N + 1;
end

rampCycles = 4; % ramps this number of cycles of the stimulus frequency at onset and offset
rampN_on = round((rampCycles / fmin) * samplingRate);
if mod(rampN_on,2)~=0
    rampN_on = rampN_on + 1;
end
rampN_off = round((rampCycles / fmax) * samplingRate);
if mod(rampN_off,2)~=0
    rampN_off = rampN_off + 1;
end

% edgeExtra = 0.05; % extra for onset and offset portions of chirp
% fmin = fmin-(edgeExtra*fmin);
% fmax = fmax+(edgeExtra*fmax);
f = [fmin,fmax]; % frequency minimum and maximum (Hz)

 stepSize = (fmax - fmin)/(N-1); % frequency step size (Hz)
F = (fmin:stepSize:fmax)';


F = 2.^((log2(fmin):log2(stepSize):log2(fmax)))';
F = linspace(log2(fmin),log2(fmax),N);
F = 2.^(F);

% fOn = ones(rampN_on,1)*f(1);
% fOff = ones(rampN_off,1)*f(end);
% F = [fOn;F;fOff];
% N = length(F);

phi = 0;
for jj = 1:N
    p(jj,1) = cos(phi);
    phi = phi + ((2*pi) / (1/F(jj) * samplingRate));
end

rampOn = hann(rampN_on*2);
rampOn = rampOn(1:rampN_on);
rampOff = hann(rampN_off*2);
rampOff = rampOn(1:rampN_off);
p(1:rampN_on) = p(1:rampN_on) .* rampOn;
p = flipud(p);
p(1:rampN_off) = p(1:rampN_off) .* rampOff;
p = flipud(p);

multiplier = 10.^(desiredGain_dB / 20);
p = p * multiplier;

padLen = 0.005; % add 5 ms of zero padding to beginning and end
padN = round(padLen * samplingRate);
pad = zeros(padN,1);
stimulus = [pad;p;pad];

ref = 0.00002;
[frequency,signal,noiseFloor] = ARLas_fda(stimulus,samplingRate,ref);
figure
plot(frequency,signal)



keyboard
