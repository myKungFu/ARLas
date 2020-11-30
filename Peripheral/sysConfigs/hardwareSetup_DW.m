function [inputs,outputs] = hardwareSetup_DW()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [inputs,outputs] = hardwareSetup_DW;
%
% Get hardware inputs and outputs for use with ARLas in Jeff's lab at Wash U.
%
% Authors: Shawn S Goodman, PhD & Jeffery T Lichtenhan PhD
% Auditory Research Lab, Dept. of Comm. Sciences & Disorders, The University of Iowa
% Auditory Physiology Lab, Dept. of Otolaryngology, Washington University School of Medicine, St. Louis
% Date: January 8, 2019
% Updated: October 10, 2019 -- ssg -- Changed the name to include _DW
%                                     Added a version _DW2 as a clone of the WU setup for use at Iowa
%                                     Changed old name to call the new name, for sake of backwards compatibility.
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

    inputs{4}.label = 'GRASS8'; % data will be saved using this string
    inputs{4}.micSens = 0.00092045;   % sensitivity in V/Pa
    inputs{4}.gain = 0;        % amplifier gain in dB
    inputs{4}.ch = 6;           % input channel being used on the sound card
    
% Specify Outputs ----------
    outputs{1}.label = 'ER10xA';% not used, but saved for reference
    outputs{1}.ch = [3,4];      % output channel(s) on the sound card
    
    outputs{2}.label = 'ER10xB';% not used, but saved for reference
    outputs{2}.ch = [5,6];      % output channel(s) on the sound card

end
