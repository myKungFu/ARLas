function [] = ARLas_teoae_linear(varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ARLas_teoae_linear(varargin)
%
% Measure transient otoacoustic emissions, using a linear extration
% paradigm.
%
% Auditory Research Lab, The University of Iowa
% Deptartment of Communication Sciences & Disorders
% The University of Iowa
% Author: Shawn S. Goodman, PhD
% Date: October 27, 2016
% Last Updated: November 29, 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

obj = varargin{1}; % get the arlas object

% 1) CREATE THE STIMULUS --------------------------------------------------
% Each column is one output channel. The stimulus must be a matrix the same
% size as the number of initialized output channels
fs = obj.fs; % get the system sampling rate
len = 0.05; % desired stimulus length (s)
    nSamples = round(len * fs); % number of samples in stimulus
    if mod(nSamples,2) ~= 0 % force to be even number
        nSamples = nSamples + 1;
    end
    
stimulus = zeros(nSamples,1); % create a vector of zeros
offset = 0.003; % time delay of the click re: time zero (s)
indx = round(offset *fs); % number of samples offset of click

stimulus1 = stimulus;
%stimulus2 = stimulus;
stimulus1(indx,1) = .75; % make a single impulse
%stimulus2(indx,1) = .25;
%stimulus2(indx+50,1) = .25; % make a single impulse
stimulus2 = randn(nSamples,100) * .01;


% 2) LOAD THE STIMULUS ----------------------------------------------------

% Load Output:
%   Load one vector at a time. Each vector is a channel of output.
%        Use vector of zeros if desire an output channel with no output.
%        Specify channel number for each (1 through maxN, where maxN is the maximum for the sound card)
obj.setPlayList(stimulus1,2);
obj.setPlayList(stimulus2,3);

% Load Input:
%   Specify channel number for each (1 through maxN, where maxN is the maximum for the sound card)
%        For each, specify a label, mic sensitivity, and gain.
%        Label is a string for idenfification and file saving purposes
%        Mic sens is in V/Pa. If no microphone is used, set = 1.
%        Gain refers to hardware gain applied prior to ADC by the sound card. Specify in dB.
label = 'test6';
micSens = 0.06;
gain = 6;
ch = 6;
obj.setRecList(ch,label,micSens,gain);

label = 'test5';
micSens = 0.05;
gain = 5;
ch = 5;
obj.setRecList(ch,label,micSens,gain);


obj.objPlayrec.nReps = 100; % number of times to play stimulus

% 3) PLAYBACK & RECORD ----------------------------------------------------
obj.objPlayrec.run % run the stimulus
ok = obj.checkForErrors;
if ~ok
   return
end

% channel = 1; % which channel to retrieve data from
% [header,data] = obj.retrieveData(channel);
% peaks(1) = max(mean(data,2))
% 
% 
% obj.setPlayList(stimulus*1.2,ch);
% obj.objPlayrec.run % run the stimulus
% ok = obj.checkForErrors;
% if ~ok
%    return
% end    





% 4) RETRIEVE DATA ----------------------------------------------------
[header6,data6] = obj.retrieveData('Ch6');

[header5,data5] = obj.retrieveData('Ch5');

keyboard

%obj.clearRecList
obj.objPlayrec.playChanList
obj.objPlayrec.loadingDock

obj.clearPlayList(0)

obj.objPlayrec.playChanList
obj.objPlayrec.loadingDock


keyboard


obj.clearRecList([1,3])
obj.objPlayrec.recChanList
obj.objPlayrec.label
obj.objPlayrec.ampGain
obj.objPlayrec.micSens

        


keyboard
return
keyboard
% Alternate way, if you want to get a list of the file names and then load
% them later in the the experiment, instead of immediately.
% fileName = obj.objPlayrec.savedFiles{1}; % the name of the file saved to input channel 1
% pathName = obj.objPlayrec.savedFilesPath; % the location of the saved files

% 5) ANALYZE & PLOT DATA ----------------------------------------------------
% Common Data Dperations: (code is available in ARLas/Core/general/)
%   HP Filter
    cutoff = 500;
    data = ARLas_hpFilter(data,fs,cutoff);
%   Artifact Rejection
    data = ARLas_artifactReject(data);
%   Time Window the Emission
    time = (0:1:nSamples-1)'/fs*1000; % time in s
    [~,indx] = max(abs(mean(data,2))); % location of the peak is time zero
    data = data(indx:end,:); % cut off time before zero
    time = time(indx:end);
    offset = 0.004; % portion to discard (s)
    indx = round(offset *fs); % number of samples offset of click
    data = data(indx:end,:);
    time = time(indx:end);
%   Ramp Edges
    len = 0.002; % use 2 ms ramps
    data = ARLas_ramp(data,fs,len);
%   Mean waveform
    d = mean(data,2); 
    
figure
plot(time,d)
xlim([time(1) time(end)])
xlabel('Time (ms)','FontSize',12)
ylabel('Amplitude (Pa)','FontSize',12)
title('TEOAE WAVEFORM','FontSize',12)
    
%   Frequency Domain Analysis
    ref = 0.00002; % reference is 20 uPa
    [frequency,signal,noiseFloor] = ARLas_fda(data,fs,ref);
	
figure
plot(frequency/1000,signal)
hold on
plot(frequency/1000,noiseFloor,'r')
xlim([0 10])
xlabel('Frequency (kHz)','FontSize',12)
ylabel('Magnitude (dB SPL)','FontSize',12)
title('linear TEOAE SPECTRUM','FontSize',12)

end % end of experiment file