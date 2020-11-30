function [inputs,outputs] = hardwareSetup_DW2()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [inputs,outputs] = hardwareSetup_DW2;
%
% Get hardware inputs and outputs for use with ARLas.
%
% Authors: Shawn S Goodman, PhD & Jeffery T Lichtenhan PhD
% Auditory Research Lab, Dept. of Comm. Sciences & Disorders, The University of Iowa
% Auditory Physiology Lab, Dept. of Otolaryngology, Washington University School of Medicine, St. Louis
% Date: January 8, 2019
% Updated: October 10, 2019 -- ssg -- This is a clone of Jeff's setup at WU.
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
    
    inputs{3}.label = 'GRASSbio'; % GRASS CP511
    inputs{3}.micSens = 1;      % sensitivity in V/Pa
    inputs{3}.gain = 80;        % amplifier gain in dB
    inputs{3}.ch = 5;           % input channel on the sound card    

    % NOTE: This one is different from Jeff's!
    % This one is shawn's GRAS 1/8-inch 40DP, sn 250118.
    % Has a different sensitivity value.
    inputs{4}.label = 'GRASS8'; % data will be saved using this string
    inputs{4}.micSens = 0.00088;   % sensitivity in V/Pa
    inputs{4}.gain = 0;        % amplifier gain in dB
    inputs{4}.ch = 6;           % input channel being used on the sound card
    
% Specify Outputs ----------
    outputs{1}.label = 'ER10xA';% not used, but saved for reference
    outputs{1}.ch = [3,4];      % output channel(s) on the sound card
    
    outputs{2}.label = 'ER10xB';% not used, but saved for reference
    outputs{2}.ch = [5,6];      % output channel(s) on the sound card

end
