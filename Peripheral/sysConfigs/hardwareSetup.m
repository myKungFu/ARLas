function [inputs,outputs,extra] = hardwareSetup()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [inputs,outputs] = hardwareSetup;
%
% Get hardware inputs and outputs for use with ARLas.
% This version is for setupt at USF for FX study.
%
% Auditory Research Lab, The University of Iowa
% Deptartment of Communication Sciences & Disorders
% The University of Iowa
% Author: Shawn Goodman
% Date: November 4,2017
% Last Updated: June 14, 2021 -- ssg -- updated for USF FX study
% Last Updated: August 4, 2021 -- ssg -- added info for coding ER10X
%       settings; used with updated ARLas version to display reminder message
%       about output limiter.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% specify additional info, if desired
    extra.outputLimiter = 0;    % =1 if output limiter being used; else 0

% Specify Inputs -------------------------------------------------------------------------
    inputs{1}.label = 'ER10xA'; % data will be saved using this string
    inputs{1}.micSens = 0.05;   % sensitivity in V/Pa
    inputs{1}.gain = 20;        % amplifier gain in dB
    inputs{1}.ch = 3;           % input channel being used on the sound card

    inputs{2}.label = 'ER10xB'; % data will be saved using this string
    inputs{2}.micSens = 0.05;   % sensitivity in V/Pa
    inputs{2}.gain = 20;        % amplifier gain in dB
    inputs{2}.ch = 4;           % input channel on the sound card

    % Note: Optiamp allows the following gain values 20*log10([5,10,50,100,200,300]*1000)'
        % ans in dB =
        % 5     73.9794
        % 10    80.0000
        % 50    93.9794
        % 100  100.0000
        % 200  106.0206
        % 300  109.5424
    inputs{3}.label = 'optiAmp1'; % IHS bioamplifier for ABR/CAP
    inputs{3}.micSens = 1;      % sensitivity in V/Pa
    inputs{3}.gain = 93.9794;   % amplifier gain in dB
    inputs{3}.ch = 8;           % input channel on the sound card  
    
    inputs{4}.label = 'optiAmp2'; % IHS bioamplifier for ABR/CAP
    inputs{4}.micSens = 1;      % sensitivity in V/Pa
    inputs{4}.gain = 93.9794;   % amplifier gain in dB
    inputs{4}.ch = 7;           % input channel on the sound card  

    inputs{5}.label = '250118'; % G.R.A.S 1/8" Pressure Mic, Type 40DP
    inputs{5}.micSens =0.00088; % sensitivity in V/Pa (shawns)
    inputs{5}.gain = 0;         % amplifier gain in dB
    inputs{5}.ch = 5;           % input channel on the sound card     

    inputs{6}.label = 'clicker'; % data will be saved using this string
    inputs{6}.micSens = 1;   % sensitivity in V/Pa
    inputs{6}.gain = 0;        % amplifier gain in dB
    inputs{6}.ch = 1;           % input channel being used on the sound card
    
    inputs{7}.label = 'IHS'; % data will be saved using this string
    inputs{7}.micSens = 0.5;   % sensitivity in V/Pa (500 IHS vs 50 for the ER stuff!)
    inputs{7}.gain = 0;        % amplifier gain in dB (0 for IHS vs 20 for the ER stuff!)
    inputs{7}.ch = 2;           % input channel being used on the sound card







% Specify Outputs  -------------------------------------------------------------------------
    outputs{1}.label = 'ER10xA';% not used, but saved for reference
    outputs{1}.ch = [3,4];      % output channel(s) on the sound card
    
    outputs{2}.label = 'ER10xB';% not used, but saved for reference
    outputs{2}.ch = [5,6];      % output channel(s) on the sound card

    outputs{3}.label = 'IHS';% not used, but saved for reference
    outputs{3}.ch = [1,2];      % output channel(s) on the sound card

% OLD CODE ----------------------------------------------------------------


%    inputs{3}.label = '109338'; % G.R.A.S 1/8" Pressure Mic, Type 40DP
%     inputs{3}.micSens = 0.0012; % sensitivity in V/Pa
%     inputs{3}.gain = 0;         % amplifier gain in dB
%     inputs{3}.ch = 5;           % input channel on the sound card     
    % 
    % inputs{3}.label = 'GrasBio'; % G.R.A.S bioamplifier for ABR/CAP
    % inputs{3}.micSens = 1; % sensitivity in V/Pa
    % inputs{3}.gain = 80;         % amplifier gain in dB
    % inputs{3}.ch = 8;           % input channel on the sound card    

%
    % inputs{4}.label = '250118'; % G.R.A.S 1/8" Pressure Mic, Type 40DP
    % %inputs{4}.micSens = 0.00088; % sensitivity in V/Pa (shawns)
    % inputs{4}.micSens = 0.00092; % sensitivity in V/Pa (jeffs)
    % inputs{4}.gain = 0;         % amplifier gain in dB
    % inputs{4}.ch = 5;           % input channel on the sound card    
    % 
    % inputs{5}.label = 'IEC711'; % data will be saved using this string
    % inputs{5}.micSens = 0.01102;% sensitivity in V/Pa
    % inputs{5}.gain = 0;         % amplifier gain in dB
    % inputs{5}.ch = 5;           % input channel being used on the sound card
    % 
    % inputs{6}.label = 'clicker'; % clicker for behavioral audiometry
    % inputs{6}.micSens = 1;      % sensitivity in V/Pa
    % inputs{6}.gain = 0;         % amplifier gain in dB
    % inputs{6}.ch = 7;           % input channel on the sound card    
