function [stimClick,stimNoise,durSamp,Time] = memrStim_v8(fs)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [stimClick,stimNoise,durSamp,Time] = memrStim_v8(fs);
%
% Create stimuli for MEMR test.
% The vectors are then folded into matrices for efficiency presenting
% through ARLas.
% This function replicates the experiment paradigm of Mepani et al.
% This function is called by MEMR_v5.m
%
% Author: Shawn S. Goodman
% Date: July 7, 2021
% Last Updated: July 15, 2021 -- ssg
% Last Updated: August 23, 2021 -- ssg -- no major changes; bring version
%                       number in line with other current versions
% Last Updated: February 27, 2022 -- ssg -- bring in line with newest version.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % zeroPadding click1 silence1 noise silence2 click2 silence3
    durSec = [0.005, 0.005,    0.03,   0.5,   0.005,  0.005,  1.0-0.0005]; % duration of each piece in seconds
    durSamp = round(durSec*fs); % duration of each piece in samples
    zPad = zeros(durSamp(1),1);
    click1 = zeros(durSamp(2),1);
      click1(1) = 0.99;
    silence1 = zeros(durSamp(3),1);
    noise = randn(durSamp(4),1); % random, white noise with mean 0 and std(1)
    
    % assume the noise should be flat from 100 Hz to 16 kHz
    b = getbpf; % bandpass filter from 1-16 kHz
    noise = fastFilter(b,noise);
        noise = noise * 0.5; % scale down so doesn't overdrive the output (max = -1 to 1)
        % here we maked "low-noise noise"
        noiseMax = 1;
        indx = find(abs(noise) >= noiseMax); % find the samples that exceed the max
        replacement = (rand(length(indx),1)*2)-1;
        noise(indx) = replacement;
        noise = noise * .99; 
    
      winLen = 0.0025;
      winN = round(winLen * fs);
      h = hann(winN*2);
      h = h(1:winN);
      noise(1:winN) = noise(1:winN).*h;
      noise = flipud(noise);
      noise(1:winN) = noise(1:winN).*h;
      noise = flipud(noise);

    silence2 = zeros(durSamp(5),1);
    click2 = zeros(durSamp(6),1);
      click2(1) = 0.99;
    silence3 = zeros(durSamp(7),1);
    stimClick = [zPad;click1;silence1;noise*0;silence2;click2;silence3];
    stimNoise = [zPad;click1*0;silence1;noise;silence2;click2*0;silence3];
    
    foldLength = 0.1; % fold into a matrix with columns this length (sec)
    foldN = round(foldLength * fs); % number of samples in each column
    NN = length(stimClick);
    nFolds = ceil(NN / foldN);
    extra = mod(NN,foldN);
    if extra ~= 0
        pad = zeros(foldN - extra,1);
        stimClick = [stimClick;pad];
        stimNoise = [stimNoise;pad];
    end
    stimClick = reshape(stimClick,foldN,nFolds);
    stimNoise = reshape(stimNoise,foldN,nFolds);
    NN = length(stimClick(:));
    time = (0:1:NN-1)'/fs;
    Time = reshape(time,foldN,nFolds);
end

% Internal Functions ------------------------------------------------------
function b = getbpf()
    Fs = 96000;
    Fstop1 = 0;             % First Stopband Frequency
    Fpass1 = 100;            % First Passband Frequency
    Fpass2 = 16000;           % Second Passband Frequency
    Fstop2 = 22000;           % Second Stopband Frequency
    Dstop1 = 0.001;           % First Stopband Attenuation
    Dpass  = 0.057501127785;  % Passband Ripple
    Dstop2 = 0.0001;          % Second Stopband Attenuation
    flag   = 'scale';         % Sampling Flag
    % Calculate the order from the parameters using KAISERORD.
    [N,Wn,BETA,TYPE] = kaiserord([Fstop1 Fpass1 Fpass2 Fstop2]/(Fs/2), [0 ...
                                 1 0], [Dstop1 Dpass Dstop2]);
    if mod(N,2)~= 0
        N = N + 1;
    end
    b  = fir1(N, Wn, TYPE, kaiser(N+1, BETA), flag)';
end
