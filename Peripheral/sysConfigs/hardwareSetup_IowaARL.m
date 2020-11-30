function [inputs,outputs] = hardwareSetup()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [inputs,outputs] = hardwareSetup;
%
% Get hardware inputs and outputs for use with ARLas.
%
% Auditory Research Lab, The University of Iowa
% Deptartment of Communication Sciences & Disorders
% The University of Iowa
% Author: Shawn Goodman
% Date: November 4,2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Specify Inputs ----------
    inputs{1}.label = 'ER10xA'; % data will be saved using this string
    inputs{1}.micSens = 0.05;   % sensitivity in V/Pa
    inputs{1}.gain = 20;        % amplifier gain in dB
    inputs{1}.ch = 3;           % input channel being used on the sound card

    inputs{2}.label = 'ER10xB'; % data will be saved using this string
    inputs{2}.micSens = 0.05;   % sensitivity in V/Pa
    inputs{2}.gain = 20;        % amplifier gain in dB
    inputs{2}.ch = 4;           % input channel on the sound card
    
    inputs{3}.label = '250118'; % G.R.A.S 1/8" Pressure Mic, Type 40DP
    inputs{3}.micSens = 0.00088; % sensitivity in V/Pa
    inputs{3}.gain = 0;         % amplifier gain in dB
    inputs{3}.ch = 5;           % input channel on the sound card    

    inputs{4}.label = 'ER10C1'; % data will be saved using this string
    inputs{4}.micSens = 0.05;   % sensitivity in V/Pa
    inputs{4}.gain = 20;        % amplifier gain in dB
    inputs{4}.ch = 1;           % input channel being used on the sound card
    
    inputs{5}.label = 'ER10C2'; % data will be saved using this string
    inputs{5}.micSens = 0.05;   % sensitivity in V/Pa
    inputs{5}.gain = 20;        % amplifier gain in dB
    inputs{5}.ch = 2;           % input channel being used on the sound card

    inputs{6}.label = 'IEC711'; % data will be saved using this string
    inputs{6}.micSens = 0.01102;% sensitivity in V/Pa
    inputs{6}.gain = 0;         % amplifier gain in dB
    inputs{6}.ch = 5;           % input channel being used on the sound card
    
    
    % Specify Outputs ----------
    outputs{1}.label = 'ER10xA';% not used, but saved for reference
    outputs{1}.ch = [3,4];      % output channel(s) on the sound card
    
    outputs{2}.label = 'ER10xB';% not used, but saved for reference
    outputs{2}.ch = [5,6];      % output channel(s) on the sound card

    outputs{3}.label = 'lights';% blinking yellow light for speech experiment
    outputs{3}.ch = [7,8];      % output channel(s) on the sound card
    
    outputs{4}.label = 'ER10C1';% blinking yellow light for speech experiment
    outputs{4}.ch = [1];      % output channel(s) on the sound card

    outputs{5}.label = 'ER10C2';% blinking yellow light for speech experiment
    outputs{5}.ch = [2];      % output channel(s) on the sound card

%    inputs{3}.label = '109338'; % G.R.A.S 1/8" Pressure Mic, Type 40DP
%     inputs{3}.micSens = 0.0012; % sensitivity in V/Pa
%     inputs{3}.gain = 0;         % amplifier gain in dB
%     inputs{3}.ch = 5;           % input channel on the sound card        
    