function [fileName] = fileNamesSetup(expName)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [fileName] = fileNamesSetup(expName);
%
% Get preferred experiment file name to populate ID window.
% For use with ARLas.
% This version is for set up for FX study.
% Save this in the folder
%       ARLas\Peripheral\sysConfigs\
%
% Auditory Research Lab, The University of Iowa
% Deptartment of Communication Sciences & Disorders
% The University of Iowa
% Author: Shawn Goodman
% Date: July 31, 2021
% Last Updated: July 31, 2021 -- ssg -- updated for USF FX study protocol
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % List the known experiment names and their preferred associated file names
    try
        if contains(expName,'fixedRatioDiscrete')
            fileName = 'dpoaeDiscrete';
        elseif contains(expName,'fixedL1Sweep')
            fileName = 'dpoaeLGFsweep';
        elseif contains(expName,'fixedRatioSweep')
            fileName = 'dpoaeGramSweep';
        elseif contains(expName,'memr')
            fileName = 'memr';
        elseif contains(expName,'audio')
            fileName = 'audio';
        elseif contains(expName,'caps')
            fileName = 'caps';
        elseif contains(expName,'teoae')
            fileName = 'teoae';
        elseif contains(expName,'sfoae')
            fileName = 'sfoae';
        elseif contains(expName,'cal')
            fileName = 'cal';
        elseif contains(expName,'Cal')
            fileName = 'cal';
        else
            fileName = [];
        end
    catch
        fileName = [];
    end
end