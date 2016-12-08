function [] = ARLas_teoae_2e(varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ARLas_teoae_2e(varargin)
%
% Measure transient otoacoustic emissions, using a double-evoked (2e) extraction
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
Stimulus = zeros(nSamples,obj.chansOut); % create a matrix of zeros
offset = 0.003; % time delay of the click re: time zero (s)
indx = round(offset *fs); % number of samples offset of click
S1 = Stimulus;
S2 = Stimulus;
S1(indx,1) = 0.2512; % make a single impulse in channel 1 (probe click)
S2(indx,2) = 0.99; % make a single impulse in channel 2 (suppressor click)
silence = zeros(size(S1));  
stimulusCh1 = [S1;silence;S1];
stimulusCh2 = [silence;S2;S2];

% 2) LOAD THE STIMULUS ----------------------------------------------------
obj.objPlayrec.stimTrain.Ch1 = stimulusCh1; % load the stimulus
obj.objPlayrec.stimTrain.Ch2 = stimulusCh2; % load the stimulus
obj.objPlayrec.nReps = 500; % number of times to play stimulus

% 3) PLAYBACK & RECORD ----------------------------------------------------
obj.objPlayrec.run % run the stimulus
if obj.killRun
   return
end    


% 4) RETRIEVE DATA ----------------------------------------------------
channel = 1; % which channel to retrieve data from
[header,data] = obj.retrieveData(channel);
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
    
%   Separate into the three parts
    start = 1;
    finish = nSamples;
    P1 = data(start:finish,:);
    Probe = P1;
    probe = mean(Probe,2);
    start = start + nSamples;
    finish = finish + nSamples;
    P2 = data(start:finish,:);
    start = start + nSamples;
    finish = finish + nSamples;
    P12 = data(start:finish,:);
    data = P1 + P2 - P12;
    
%   Time Window the Emission
    [~,indx] = max(abs(probe)); % location of the peak is time zero
    data = data(indx:end,:); % cut off time before zero
    offset = 0.001; % portion to discard (s)
    indx = round(offset *fs); % number of samples offset of click
    data = data(indx:end,:);
%   Ramp Edges
    len = 0.001; % use 1 ms ramps
    data = ARLas_ramp(data,fs,len);
%   Mean waveform
    d = mean(data,2); 
    time = (0:1:length(d)-1)'/fs*1000; % time in s
    
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
title('2e TEOAE SPECTRUM','FontSize',12)

end % end of experiment file