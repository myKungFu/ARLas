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
% Last Updated: July 24, 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

obj = varargin{1}; % get the arlas object

% ------------------- USER ADJUSTABLE PARAMETERS --------------------------
nReps = 100; % number of times to play each paired click condition
clickAmp = .25; % click amplitude re: full out

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
stimulus = zeros(nSamples,1); % create a vector of zeros
offset = 0.003; % time delay of the click re: time zero (s)
indx = round(offset *fs); % number of samples offset of click
stimA = stimulus;
stimA(indx,1) = clickAmp; % make a single impulse


% LOAD THE STIMULUS ----------------------------------------------------
% Load Output:
%   Load one vector at a time. Each vector is a channel of output.
%        Use vector of zeros if desire an output channel with no output.
%        Specify channel number for each (1 through maxN, where maxN is the maximum for the sound card)
obj.setPlayList(stimA,outputA.ch(1));

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
obj.objPlayrec.run % run the stimulus
if obj.killRun
   return
end

% RETRIEVE RECORDED DATA
[headerA,DataA] = obj.retrieveData(['Ch',num2str(inputA.ch)]); % retrive recorded data

figure
plot(DataA)

end % end of experiment file