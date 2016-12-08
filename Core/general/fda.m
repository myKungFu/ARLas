function [signal,noiseFloor,sem,rms,phase,frequency] = fda(x,fs,parameters,nfft)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [signal,noiseFloor,sem,rms,phase,frequency] = fda(x,fs,parameters,nfft);
% Frequency Domain Analysis
% x = input time-domain waveform.  Data in columns, matrix or vector. 
% fs = sampling rate in Hz
%
% parameters is an optional input structure that allows the user to control certain
% options when calling the function.  Each value should be set to 1 (on) or 0 (off).
% Usage:  parameters = struct('norm1HzBandwidth',1,'doOctaves',1,'doEqSpectLvl',1); 
%       norm1HzBandwidth controls whether output is normalized to a 1-Hz
%           wide bandwidth.
%       doOctaves controls whether output is on a normal linear scale or
%           integrated over 1/3 octaves.
%       doEqSpectLvl controls conversion of 1/3 octave levels to eqivalent
%           spectum levels.  Only used if doOctaves == 1.
%
% For accuracy, observe the following guidelines:
%   Sinusoids should normally be measured in linear bands; use 
%       parameters = struct('norm1HzBandwidth',0,'doOctaves',0,'doEqSpectLvl',0);
%   If for some resaon you need to measure sinusoids in 1/3 octave bands, use 
%       parameters = struct('norm1HzBandwidth',0,'doOctaves',1,'doEqSpectLvl',0);
%   When measuring noise in 1/3 octave bands, use 
%       parameters = struct('norm1HzBandwidth',1,'doOctaves',1,'doEqSpectLvl',1);
%           or
%       parameters = struct('norm1HzBandwidth',1,'doOctaves',1,'doEqSpectLvl',0);
%   When measuring noise in normal linear bands, use 
%       parameters = struct('norm1HzBandwidth',1,'doOctaves',0,'doEqSpectLvl',0);
%
% Note: Each of the parameters values may also be controlled manually by setting their
% values in the "switchbox" and omitting parameters as an input argument.
%
% Author: Shawn Goodman
% Date: March 9, 2007
% Modified: May 18, 2007, October 17, 2008
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%--------------------- switchbox -------------------------------------
norm1HzBandwidth = 0; % normalize to a 1-Hz bandwidth
doOctaves = 0; % integrate over 1/3 octaves 
doEqSpectLvl = 0;  % convert from 1/3 octaves to equivalent spectrum level

Pref = 0.00002; %.00002; % pressure reference (Pa)
doPlot = 0; % turn plotting on and off
doWindow = 0; % window entire signal & normalize
%---------------------------------------------------------------------
signal=[]; noiseFloor=[]; sem=[]; rms=[]; phase=[]; frequency=[]; nZeros = 0;

if nargin >= 3, % if a parameters input argument exists, its values overide the switchbox settings
    norm1HzBandwidth = parameters.norm1HzBandwidth;
    doOctaves = parameters.doOctaves;
    doEqSpectLvl = parameters.doEqSpectLvl;
end

[originalN,K] = size(x); % Input data are arranged into a matrix of N samples by K buffers.
if nargin == 4, % if number of samples in the fft is specified,
    nZeros = nfft - size(x,1); % pad with zeros to get to the correct size
    if nZeros < 0,
        disp('ERROR: Requested FFT size is smaller than the data to be transformed.')
        nZeros = 0;
    end
end

if nZeros > 0,
    padding = zeros(nZeros,K);
    x = [x;padding];
    N = size(x,1);
else
    N = originalN;
end

% if mod(N,2) ~= 0, % force N to be even
%     N = N + 1;
%     disp('WARNING: Requested FFT size is odd. Adding one additional sample.')
% end


%originalN = N; % code to be used when zero padding is used...
%originalN = 111;

%rmsCheck = 20*log10(sqrt(mean(mean(x,2).^2))/Pref); % time-domain check for correct freq-domain rms calculations
rmsCheck = mean(x,2); % take average across buffers
rmsCheck = sqrt((1/originalN)*sum(rmsCheck.^2));
rmsCheck = 20*log10(rmsCheck / Pref);

if doWindow == 1,
    w = hann(N);
    w = w ./ (sqrt((1/N)*sum(w.^2)));
    w = repmat(w,1,K);
    x = x .* w;
    disp('Window ON')
end

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
%    antimask = ~(S2 == 0); % used later to set variance to back to floating point precision
%    mask = (S2 == 0) * eps; % set variance of zero to floating point precision to avoid division by zero
%    S2 = S2 + mask;
    Se2(1) = Se2(1) * 0.5;
    Se2(end) = Se2(end) * 0.5;
%    Se2 = Se2 + mask;
end

%Xbar2 = Xbar2 * (2/N^2); % added ssg may 16
Xbar2 = Xbar2 * (2/originalN^2); % added ssg Oct 17

if ~isempty(S2),
%     S2 = S2 * (2/N^2);
%     Se2 = Se2 * (2/N^2);
     S2 = S2 * (2/originalN^2);
     Se2 = Se2 * (2/originalN^2);
end

%rms = 10*log10(sum(Xbar2/Pref^2));
rms = 10*log10(sum(Xbar2*(originalN/N)) /Pref^2); % put into dB SPL, properly scaled--use this if added zero padding
if round(rms) ~= round(rmsCheck),
    disp('Warning: time and frequency domain rms calculations do not match!')
    disp(['Difference = ',num2str(rmsCheck - rms),' dB.'])
    disp('Using time-domain calculation.')
    %keyboard
    rms = rmsCheck;
end

if doOctaves == 1, % calculate the rms energy in fractional octave bands
    % ----> [freq,Xbar2] = thirdOctaves(frequency,Xbar2,doEqSpectLvl);
    [freq,Xbar2] = thirdOctaves(frequency,Xbar2*(originalN/N),doEqSpectLvl);
    %[freq,Xbar2] = thirdOctaves(frequency,Xbar2/(N/fs));
    %[freq,Xbar2] = thirdOctaves(frequency,Xbar2/(originalN/N));
    if K > 1,
        [freq,S2] = thirdOctaves(frequency,S2,doEqSpectLvl);
        [freq,Se2] = thirdOctaves(frequency,Se2,doEqSpectLvl);
    end
    frequency = freq;
end

warning off all
if norm1HzBandwidth == 1,
    deltaF = fs / originalN; % DFT bandwidth % ad
    signal = 10*log10(Xbar2/deltaF/(Pref^2));
	noiseFloor = 10*log10(S2/deltaF/(Pref^2));
	sem = 10*log10(Se2/deltaF/(Pref^2));
else
	signal = 10*log10(Xbar2/(Pref^2));
	noiseFloor = 10*log10(S2/(Pref^2));    
	sem = 10*log10(Se2/(Pref^2));    
end
warning on all

if doPlot == 1,
	plot(frequency,signal)
	hold on
    if ~isempty(noiseFloor),
    	plot(frequency,noiseFloor,'r')
    	plot(frequency,sem,'r:')
    end
	title(['rms = ',num2str(rms),' dB SPL'])
end

function [fc,y] = thirdOctaves(frequency,x,doEqSpectLvl)
% partition the signal energy into 1/3 octaves.  
% frequency = frequency vector (Hz)
% x = energy vector
% Center frequencies are derived from ANSI S1.11-1986, 
% "Specification for Octave-Band and Fractional Octave Band Analog and Digital Filters", 
% and from ANSI S1.6-1984, "Preferred Frequencies, Frequency Levels, and Band Numbers for
% Acoustical Measurements".  The defining formulas yeild some values that
% are not integral numbers.
fc = [10, 12.5, 16, 20, 25, 31.5, 40, 50, 63, 80, 100, 125, 160, 200, 250,...
    315, 400, 500, 630, 800, 1000, 1250, 1600, 2000, 2500, 3150, 4000, 5000,...
    6300, 8000, 10000, 12500, 16000, 20000]; % filter center frequencies
fL = [9.2, 10.9, 14.3, 17.9, 22.4, 28, 33.5, 45, 56, 71, 90, 112, 140, 180, 224,...
    280, 355, 450, 560, 710, 900, 1120, 1400, 1800, 2240, 2800, 3550, 4500,...
    5600, 7100, 9000, 11200, 14000, 18000]; % lower filter cutoff
fu = [10.9, 14.3, 17.9, 22.4, 28, 33.5, 45, 56, 71, 90, 112, 140, 180, 224,...
    280, 355, 450, 560, 710, 900, 1120, 1400, 1800, 2240, 2800, 3550, 4500,...
    5600, 7100, 9000, 11200, 14000, 18000, 22400]; % upper filter cutoff
if fu(end) > frequency(end),
    ii=1;
    while fu(ii) <= frequency(end),
        ii = ii + 1;
    end
    fc(ii+1:end) = [];
    fL(ii+1:end) = [];
    fu(ii+1:end) = [];
end
%indxfc = findClosest(frequency,fc); % find the indices of the center frequencies
%indxfL = findClosest(frequency,fL); % find the indices of the lower cutoff frequencies
%indxfu = findClosest(frequency,fu); % find the indices of the upper cutoff frequencies
for mm=1:length(fc)
    [dummy,indxfc(mm,1)] = min(abs(frequency - fc(mm)));
    [dummy,indxfL(mm,1)] = min(abs(frequency - fL(mm)));
    [dummy,indxfu(mm,1)] = min(abs(frequency - fu(mm)));
end
% you can't do 1/3 octave ave where the fft-bin spacing is too far apart...
% added May 17, 2007  by ssg
maskL = indxfL == indxfc; 
masku = indxfu == indxfc;
cut = max(find(maskL(1:end-1)+masku(1:end-1))) + 1;
if ~isempty(cut),
    indxfc = indxfc(cut:end);
    indxfL = indxfL(cut:end);
    indxfu = indxfu(cut:end);
end
y = zeros(size(indxfc)); % initialize output vector
for ii=1:length(indxfc),
    y(ii,1) = (x(indxfL(ii))*.5)   +   sum(x(indxfL(ii)+1:indxfu(ii)-1))  +   (x(indxfu(ii))*.5);
end
if doEqSpectLvl == 1,
    % Convert levels expressed in 1/3 octave band levels to equivalent pressure spectrum levels.
    % Theoretical bandwidth:  BW = fc * (2^(1/N)-1), with N=3 for 3rd octave filters.
    % The first and last octaves are only half the width of the others, so N=6 (1/6 octaves).
    % Can be applied by dividing energy by BW or by subtracting 10*log10(BW) from the SPL.
    % At this point in the code we are still in energy units, so divide.
    % HOWEVER, the theoretical bandwidths are not the actual bandwidths, so here we should
    % use either 1) calculated rather than theoretical bandwidths, or 2)
    % simply divide by the number of bins in the sum.  This implementation uses the latter.
    correction = indxfu - indxfL;
    y = y ./ correction;
end
if ~isempty(cut),
    fc = fc(cut:end);
    %fu = fu(cut:end);
    %fL = fL(cut:end);
end

% eof


% Extra code-----------------------------------------------------------------------
% Xbar2 = (1/K^2) * abs(sum(Xk,2)).^2; % signal energy--alternate calculation
%rms = sum((2/N^2) * Xbar2); % rms energy in signal
%rms = 10*log10(rms/(Pref^2)); 

% TEST SIGNALS
%
% for testing sinusoids
% fs = 10000;
% nReps = 8;
% t = (0:1/fs:1-(1/fs))';
% f = 4000;
% sig = sin(2*pi*f*t);
% sig = repmat(sig,1,nReps);
% variance = randn(size(sig)) .* .001;
% sig = sig + variance;
% parametersOct = struct('norm1HzBandwidth',0,'doOctaves',1,'doEqSpectLvl',0);
% parametersSin = struct('norm1HzBandwidth',0,'doOctaves',0,'doEqSpectLvl',0);
% [signal,noiseFloor,sem,rms,phase,frequency] = fda(sig,fs,parametersSin);
% [signalB,noiseFloorB,semB,rmsB,phaseB,frequencyB] = fda(sig,fs,parametersOct);
% plot(frequency,signal,'b'); hold on; stem(frequencyB,signalB,'ro-')
% plot(frequency,noiseFloor,'b:'); hold on; plot(frequencyB,noiseFloorB,'r:')
% max(signal)-max(signalB)
% 
% for testing impulses
% t = .25; % time in seconds
% fs = 10000;
% nSamples = round(fs*t);
% nReps = 8;
% sig = zeros(nSamples,nReps);
% variance = randn(1,nReps) .* .1;
% sig(1,:) = 1 + variance;
% parametersOct = struct('norm1HzBandwidth',1,'doOctaves',1,'doEqSpectLvl',1);
% parametersSin = struct('norm1HzBandwidth',1,'doOctaves',0,'doEqSpectLvl',0);
% [signal,noiseFloor,sem,rms,phase,frequency] = fda(sig,fs,parametersSin);
% [signalB,noiseFloorB,semB,rmsB,phaseB,frequencyB] = fda(sig,fs,parametersOct);
% plot(frequency,signal,'b'); hold on; stem(frequencyB,signalB,'ro-')
% plot(frequency,noiseFloor,'b:'); hold on; plot(frequencyB,noiseFloorB,'r:')
% 
% 
% % for testing white noise
% t = .25; % time in seconds
% fs = 10000;
% nSamples = round(fs*t);
% sig = randn(nSamples,1);
% parametersOct = struct('norm1HzBandwidth',1,'doOctaves',1,'doEqSpectLvl',1);
% parametersSin = struct('norm1HzBandwidth',1,'doOctaves',0,'doEqSpectLvl',0);
% N = 1000;
% [signalB,noiseFloor,sem,rms,phase,frequencyB] = fda(sig,fs,parametersOct);
% signal = zeros(nSamples/2 + 1,N);
% signalB = zeros(size(signalB,1),N);
% for ii=1:N,
%     ii
%     sig = randn(nSamples,1);
%     [signal(:,ii),noiseFloor,sem,rms,phase,frequency] = fda(sig,fs,parametersSin);
%     [signalB(:,ii),noiseFloor,sem,rmsB,phase,frequencyB] = fda(sig,fs,parametersOct);
% end
% q = 10.^(signal/10); % you can't average dB, so put into linear units and then back
% qB = 10.^(signalB/10);
% qq = mean(q,2);
% qqB = mean(qB,2);
% qq = 10*log10(qq);
% qqB = 10*log10(qqB);
% plot(frequency,qq,'r-')
% hold on
% plot(frequencyB,qqB,'go-')
% hold off
% mean(qq)
% mean(qqB)
% mean(qq)-mean(qqB)


% old code had this before my May 17 change to thirdOctaves.
% if doOctaves == 1,
%     if N < fs*2,
%         N = fs * 2; % ensure 1/2-Hz bin size--need this for my implementation of 1/3 octave filters...
%     end
% end


% Old (incorrect) way of calculating RMS.  Didn't work with zero padding...
% Replaced May 16, 2007 -- ssg
%rms = sum(Xbar2); % rms energy in signal
%rms = 10*log10((2/N^2) * (rms /(Pref^2))); % put into dB SPL, properly
%scaled


% fResponse = microphone frequency response (dB). Column 1 is frequencies;
%               column 2 is amplitudes. If given, amplitudes will be
%               corrected using these values.
% if nargin < 3,
%     fResponse = [];
% end
% if ~isempty(fResponse),
%     amps = interp1(fResponse(:,1),fResponse(:,2),frequency,'pchip'); % interpolate correction values
%     Xbar2 = 10.^((10*log10(Xbar2) - amps)/10); % subtract the mic frquency response to obtain correct energy values
%     S2 = 10.^((10*log10(S2) - amps)/10);
%     Se2 = 10.^((10*log10(Se2) - amps)/10);
% end

