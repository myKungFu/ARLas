function [micCorrection] = getMicCorrection(probeInput,doMicCorrection,obj)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%[micCorrection] = getMicCorrection(probeInput,doMicCorrection);
%
% Get the most recent microphone correction.
%
% probeInput = structure describing probe. Obtained from hardwareSetup.m
% doMicCorrection = 0 (do not use a mic correction) or 
%                   1 (do use the most recent correction available.
%                   Note: if 0, mic correction is the scalar value 1.
%
% Author: Shawn Goodman
% Date: August 23, 2021
% Last Updaed: January 6, 2025 -- ssg -- added obj as a third input
%               argument. Old code may break because of this. I added an
%               error message to help.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if isempty(probeInput)
        micCorrection = [];
        return
    end
    if nargin < 3
        warning('ERROR: Input to getMicCorretion must now contain 3 arguments, including obj as a third argument.')
        micCorrection = [];
        return
    end

    [pathName_mic,folderName_mic,fileName_mic] = mostRecentCalibration('mic',probeInput.label,[],obj);
    if doMicCorrection == 1
        if ~isempty(fileName_mic)
            dummy = load([pathName_mic,folderName_mic,fileName_mic]);
            micCorrection = dummy.micCorrection; % microphone correction
        else
            micCorrection = 1; % convolving by this changes nothing
        end
    else
        micCorrection = 1; % convolving by this changes nothing
    end
end
