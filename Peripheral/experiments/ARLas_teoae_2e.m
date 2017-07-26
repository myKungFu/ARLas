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
% Date: July 26, 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

obj = varargin{1}; % get the arlas object

% ------------------- USER ADJUSTABLE PARAMETERS --------------------------
nReps = 1000;   % number of times to play each paired click condition
suppAmp = .75; % suppressor amplitude re: full out
probeAmp = suppAmp * 10^(-12/20); % probe amplitude re: suppressor. (here, 12 dB down)

% Specify Inputs and Outputs:
% INPUTS -----
    % ER10x probe A
    inputA.label = '0141'; % change this to be whatever label you want
    inputA.micSens = 0.05; % sensitivity in V/Pa
    inputA.gain = 20; % amplifier gain in dB
    inputA.ch = 7; % 1 % input channel on the sound card

    % ER10x probe B
    inputB.label = '0142'; % change this to be whatever label you want
    inputB.micSens = 0.05; % sensitivity in V/Pa
    inputB.gain = 20; % amplifier gain in dB
    inputB.ch = 8; % 2 % input channel on the sound card
    
% OUTPUTS -----    
    % ER10x probe A
    outputA.label = 'ER10xA';
    outputA.ch = [7,8]; %[1,2]; % output channels on the sound card

    % ER10x probe B
    outputB.label = 'ER10xB';
    outputB.ch = [5,6]; %[3,4]; % output channels on the sound card
%--------------------------------------------------------------------------

% 1) CREATE THE STIMULI --------------------------------------------------
% Each column is one output channel. The stimulus must be a matrix the same
% size as the number of initialized output channels
fs = obj.fs; % get the system sampling rate
len = .05; % desired stimulus length (s)
    nSamples = round(len * fs); % number of samples in stimulus
    if mod(nSamples,2) ~= 0 % force to be even number
        nSamples = nSamples + 1;
    end
probe = zeros(nSamples,1); % create a vector of zeros
supp = zeros(nSamples,1);
zpad = zeros(nSamples,1); % zero padd (buffer of zeros)
offset = 0.003; % time delay of the click re: time zero (s)
indx = round(offset *fs); % number of samples offset of click
probe(indx,1) = probeAmp; % make a single impulse for probe click
supp(indx,1) = suppAmp; % make a single impulse for suppressor click
% set the double evoked sequence
stimCh1 = [probe,zpad,probe]; 
stimCh1 = stimCh1(:);
stimCh2 = [zpad,supp,supp];
stimCh2 = stimCh2(:);

% LOAD THE STIMULUS ----------------------------------------------------
% Load Output:
%   Load one vector at a time. Each vector is a channel of output.
%        Use vector of zeros if desire an output channel with no output.
%        Specify channel number for each (1 through maxN, where maxN is the maximum for the sound card)
obj.setPlayList(stimCh1,outputA.ch(1));
obj.setPlayList(stimCh2,outputA.ch(2));

% Load Input:
%   Specify channel number for each (1 through maxN, where maxN is the maximum for the sound card)
%        For each, specify a label, mic sensitivity, and gain.
obj.setRecList(inputA.ch,inputA.label,inputA.micSens,inputA.gain);

% set desired user-specified info
obj.objPlayrec.userInfo = []; % clear out previous, non-structure array info
obj.objPlayrec.userInfo.fs = fs;

% PLAYBACK & RECORD ----------------------------------------------------
obj.setFilter(1); % filter is on (set to 0 if want no filtering)
obj.setNReps(nReps); % number of times to play stimulus
obj.run % run the stimulus
if obj.killRun
   return
end

% RETRIEVE RECORDED DATA
[header,Data] = obj.retrieveData(['Ch',num2str(inputA.ch)]); % retrive recorded data

% ANALYZE DATA
[time,teoae,frequency,Teoae,nf] = analyzeData(Data,fs,nSamples);
% PLOT DATA
plotData(time,teoae,frequency,Teoae,nf)

end % end of experiment file

% INTERAL FUNCTIONS -------------------------------------------------------
function [time,teoae,frequency,Teoae,nf] = analyzeData(Data,fs,nSamples)
    cutoff = 500;
    Data = ARLas_hpFilter(Data,fs,cutoff);
    start = 1;
    finish = nSamples;
    p1 = Data(start:finish,:);
    start = start + nSamples;
    finish = finish + nSamples;
    p2 = Data(start:finish,:);
    start = start + nSamples;
    finish = finish + nSamples;
    p12 = Data(start:finish,:);
    teoae = p1 + p2 - p12;
    teoae = ARLas_artifactReject(teoae);
    [~,indx0] = max(abs(mean(p2,2))); % time zero
    time = (0:1:nSamples-1)'/fs;
    time = time - (indx0/fs);
    time = time(indx0:end);
    teoae = teoae(indx0:end,:);
    finish = 0.02;
    [~,indx] = min(abs(finish - time));
    time = time(1:indx);
    teoae = teoae(1:indx,:);
    time = time * 1000;
    len = 0.003;
    teoae = ARLas_ramp(teoae,fs,len);
    ref = 0.00002;
    [frequency,Teoae,nf] = ARLas_fda(teoae,fs,ref);
    frequency = frequency / 1000;
    teoae = mean(teoae,2);
    teoae = teoae * 1000;
end

function plotData(time,teoae,frequency,Teoae,nf)
    figure
    plot(time,teoae,'r')
    xlim([time(1),time(end)])
    xlabel('Time (ms)')
    ylabel('Amplitude (mPa)')
    title('TEOAE WAVEFORM')

    figure
    plot(frequency,Teoae,'r')
    hold on
    plot(frequency,nf,'k')
    plot(frequency,nf+6,'k:')
    xlim([.5,10])
    xlabel('Frequency (kHz)')
    ylabel('Magnitude (dB SPL)')
    title('TEOAE SPECTRUM')
    legend('teoae','nf','nf+6')
end