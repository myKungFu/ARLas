function [] = DW_calibrateAA_full(varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DW_calibrateAA_full(varargin)
%
% Calibrate acoustic assembly. Includes both pure tones and chirp
% calibration, for verficication and comparison.
%
% Auditory Research Lab, The University of Iowa
% Deptartment of Communication Sciences & Disorders
% The University of Iowa
% Author: Shawn S. Goodman, PhD
% Date: January 11, 2017
% Last Updated: January 11, 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

obj = varargin{1}; % get the arlas object
fs = obj.fs; % get the system sampling rate

% USER CAN ADUST PARAMETERS HERE: -----------------------------------------
startFreq = 125; % lowest frequency
finishFreq = 16000; % highest frequency

% Load Input:
label = 'GRAS8'; % GRAS 1/8-inch pressure mic, SN:266921
micSens = 0.00092; % sensitivity in V/Pa
gain = -0.25; % quarter inch GRAS pre-amp gain
ch = 5; % CHANES THIS, IF NEEDED
obj.setRecList(ch,label,micSens,gain);
% make a structure to save
input1.label = label;
input1.micSens = micSens;
input1.gain = gain;
input1.ch = ch;

% Load Input:
label = 'ER10C_AOS';
micSens = 0.05;
gain = 20;
ch = 1; % CHANGE THIS, IF NEEDED!
obj.setRecList(ch,label,micSens,gain);
% make a structure to save
input2.label = label;
input2.micSens = micSens;
input2.gain = gain;
input2.ch = ch;

% Specify output:
% make a structure to save
output1.label = 'ER10C_AOS';
output1.gain = 20;
output1.ch = 1; % CHANGE THIS< IF NEEDED (channel on the sound card)
output1.receiver = 1; % which receiver on the ER10C is being used
% -------------------------------------------------------------------------

% pure tones
fracOct = 1/2; % fractional octave spacing
freqs = 2.^((log2(startFreq):fracOct:log2(finishFreq)))';
freqs = round(freqs);
nFreqs = length(freqs); % number of frequencies to test
amp_pt = 0.05; % set initial output amplitude re: full output (.1 is -20 dB)
obj.objPlayrec.nReps = 4; % number of times to play stimulus

receiver = ['Rec',num2str(output1.receiver)];
rmsG8 = zeros(nFreqs,1); % initializing output variables
rms10C = zeros(nFreqs,1);
disp(' ')
for ii=1:nFreqs % Pure tones: loop acros frequencies...
    disp(['  Testing ',num2str(freqs(ii)),' Hz;  ',num2str(ii),' of ',num2str(nFreqs),'.'])
    stimulus = makeSinusoid(freqs(ii),fs); % Get stimulus (note: is full output (+1 to -1))
    % Load Output:
    obj.setPlayList(stimulus*amp_pt,output1.ch); % load the output
    obj.objPlayrec.run % playback and record
    ok = obj.checkForErrors;
    if ~ok
       return
    end
    % get the recorded data
    [headerG8_pt,dataG8_pt] = obj.retrieveData('GRAS8');
    [header10C_pt,data10C_pt] = obj.retrieveData('ER10C');
    rmsG8(ii,1) = getRMS(dataG8_pt,fs); % analyze
    rms10C(ii,1) = getRMS(data10C_pt,fs);
end
fullOutSPL_pt = rmsG8 - 20*log10(amp_pt/1); % full output of this transducer in this configuration
% denominator is 1, because ARLas/playrec only allows +/-1 as the full
frequency_pt = freqs; % rename for clarity


% chirp: no loop required
disp(['Testing chirp'])
stimulus = getChirp(startFreq,finishFreq,fs);
amp = 0.05; % set initial output amplitude re: full output
obj.objPlayrec.nReps = 50; % number of times to play stimulus
% Load Output:
obj.setPlayList(stimulus*amp,output1.ch); % load the output
obj.objPlayrec.run % playback and record
ok = obj.checkForErrors;
if ~ok
   return
end
% get the recorded data
[headerG8_,dataG8] = obj.retrieveData('GRAS8');
[header10C,data10C] = obj.retrieveData('ER10C');
pRef = 0.00002; % pressure reference (uPa)
[frequency,signal,noiseFloor] = ARLas_fda(dataG8,fs,pRef);
[frequency10,signal10,noiseFloor10] = ARLas_fda(data10C,fs,pRef);        
BW = finishFreq - startFreq; % signal bandwidth
correction = 10*log10(BW); % bandwidth correction
fudgeFactor = 4;
fullOutSPL = signal + correction - 20*log10(amp/1) - fudgeFactor; % full output of this transducer in this configuration
noiseFloor = noiseFloor + correction - 20*log10(amp/1) - fudgeFactor;  
    
figure
plot(frequency/1000,fullOutSPL,'r')
hold on
plot(frequency/1000,noiseFloor,'k')
plot(freqs/1000,fullOutSPL_pt,'b*-')
xlabel('Frequency (kHz)','FontSize',12)
ylabel('Amplitude (dB SPL)','FontSize',12)
legend('PureTone','Chirp')
title('Output Calibration','FontSize',12)
xlim([(startFreq-startFreq*.2)/1000,(finishFreq+finishFreq*.05)/1000])


timeStamp = datestr(now);
parentDir = 'C:\myWork\ARLas\Peripheral\calibrations\';
calDir = 'DW_AAcal_FULL\';
pathName = [parentDir,calDir];
if ~exist(pathName,'dir')
    mkdir(parentDir,calDir)
    addpath(pathName)
end
fileName = ['DW_AAcal_',receiver,'_',datestr(now,'yyyy_mm_dd_hh'),datestr(now,'MM'),'.mat'];
disp(' Saving calibration file...')
try
    save([pathName,fileName],'frequency','fullOutSPL','noiseFloor','frequency_pt','fullOutSPL_pt','fs','input1','input2','output1','timeStamp')
    disp(' Successful save!')
catch ME
    disp(' ERROR: unable to save file!')
    keyboard
end
end % end of experiment file

% internal functions ------------------------------------------------------
function rms = getRMS(data,fs)
    pRef = 0.00002; % pressure reference (uPa)
    [frequency,signal,noiseFloor] = ARLas_fda(data,fs,pRef);
    rms = max(signal);
end
function [stimulus] = makeSinusoid(f,fs)
% 1) CREATE THE STIMULUS --------------------------------------------------
len = 1; % 0.05; % desired stimulus length (s)
nSamples = round(len * fs); % number of samples in stimulus
if mod(nSamples,2) ~= 0 % force to be even number
    nSamples = nSamples + 1;
end
% make pure tones
t = (0:1:nSamples-1)'/fs;
stimulus = sin(2*pi*f*t);
end

function [stimulus] = getChirp(fmin,fmax,samplingRate)
len = .2; %0.075; % stimulus length in seconds
N = round(len * samplingRate); % total number of samples
if mod(N,2)~=0
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
edgeExtra = 0.1;
fmin = fmin-(edgeExtra*fmin);
fmax = fmax+(edgeExtra*fmax);
f = [fmin,fmax]; % frequency minimum and maximum (Hz)
stepSize = (f(2)-f(1))/(N-1); % frequency step size (Hz)
F = (f(1):stepSize:f(end))';
fOn = ones(rampN_on,1)*f(1);
fOff = ones(rampN_off,1)*f(end);
F = [fOn;F;fOff];
N = length(F);
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
padLen = 0.005;
padN = round(padLen * samplingRate);
pad = zeros(padN,1);
stimulus = [pad;p;pad];
end
