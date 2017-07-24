function [chirp] = ARLas_chirp(fmin,fmax,linear,sweepRate,fs)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [chirp] = ARLas_chirp(fmin,fmax,linear,sweepRate,fs);
%
% Create a (upward moving) frequency glide.
% fmin = minimum frequency (Hz)
% fmax = maximum frequency (Hz)
% linear = 1 (linear sweep) or 0 (logarithmic, octave-based sweep)
% sweepRate = if linear sweep, frequency change in Hz/s
%             if log sweep, frequency change in octaves/s
% fs = sampling rate (Hz)
%
% Examples: 
%           linear chirp from 1-8 kHz, half a second long
%               linChirp = ARLas_chirp(1000,8000,1,3500,44100);
%           log chirp from 1-8 kHz, 250 ms per octave
%               logChirp = ARLas_chirp(1000,8000,0,0.25,44100);
%
% Author: Shawn Goodman
% Date: July 7, 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if fmax > (fs/2)
    error('Maximum frequency (fmax) cannot exceed the Nyquist rate.')
end
if fmax <= fmin
    error('Maximum frequency (fmax) must be greater than minimum frequency (fmin).')
end
    
if linear == 1 % if a linear sweep is desired
    sweepLen = (fmax-fmin) / sweepRate; % length of sweep in seconds
    sweepN = round(sweepLen * fs); % number of samples in sweep
        if mod(sweepN,2)~=0 % force number of samples to be even
            sweepN = sweepN + 1;
        end
    fStepSize = (fmax-fmin) / sweepN; % frequency step size (Hz);
    F = (fmin:fStepSize:fmax)'; % frequency vector
    F = F(1:sweepN);
elseif linear == 0 % if a log sweep is desired
    nOctaves = log2(fmax)-log2(fmin); % number of requested octaves
    sweepLen = sweepRate * nOctaves; % length of sweep in seconds
    sweepN = round(sweepLen * fs); % number of samples in sweep
        if mod(sweepN,2)~=0 % force number of samples to be even
            sweepN = sweepN + 1;
        end
    oStepSize = nOctaves / sweepN; % frequency step size (fractional octaves);
    F = 2.^((log2(fmin):oStepSize:log2(fmax))'); % frequency vector
    F = F(1:sweepN);
else
   error('Unrecognized value: linear must be 1 or 0.') 
end
% extend ends to account for onset and offset ramps
rampLen = 0.002; % ramp length (s)
rampN = round(rampLen * fs); % number of samples in each ramp
    if mod(rampN,2)~=0 % force to even number
        rampN = rampN + 1;
    end
F = [F(1)*ones(rampN,1);F;F(end)*ones(rampN,1)]; % extend the edges with constant freq

% Rotate the phase:
% At each time sample, time has advanced by the sampling period (1/fs).
% At each time sample, the phase has rotated by the frequency (f) times
% the sampling period: phase advance (in cycles) = f * (1/fs), or f/fs.
% To make this a radian change, multiply by 2pi: phi = (2*pi*f)/fs;
% Finally, add the newly accumulated phase to the previous phase value.
phi = 0; % starting phase is zero
chirp = zeros(sweepN,1); % initialize amplitude vector
for jj = 1:sweepN
    chirp(jj,1) = cos(phi);
    phi = phi + ((2*pi*F(jj))/fs);
end

% ramp the chirp on and off
h = hann(rampN*2);
h = h(1:rampN);
chirp(1:rampN) = chirp(1:rampN) .* h;
chirp = flipud(chirp);
chirp(1:rampN) = chirp(1:rampN) .* h;
chirp = flipud(chirp);

padLen = 0.002; % zero-pad the chirp with 2 ms
padN = round(padLen * fs);
pad = zeros(padN,1);
chirp = [pad;chirp;pad];
chirp = chirp / max(abs(chirp));  % rescale so 1 is max out
chirp = chirp * .99;

end