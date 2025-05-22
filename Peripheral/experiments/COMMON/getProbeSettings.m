function [probeInputL,probeOutputL,probeInputR,probeOutputR,refInput] = getProbeSettings(probeL,probeR,validateStimLevels)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [probeInputL,probeOutputL,probeInputR,probeOutputR,refInput] = getProbeSettings(probeL,probeR,validateStimLevels);
%
% get the probe settings (general) and return them.
% probeL = probe that is in LEFT ear; 'A', 'B', or []
% probeR = probe that is in LEFT ear; 'A', 'B', or []
% validateStimLevels = validate stimulus levels for left ('L') or right ('R') ear, else use empty set ([]).
%
% Author: Shawn Goodman, PhD
% Date: September 16,, 2021
% Last Updated: September 16 2021 -- ssg
% Last Updated: January 29, 2025 -- ssg -- added probe C for IHS system
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    [inputs,outputs] = hardwareSetup; % read in the hardware setup
    
    % For left probe
    if strcmp(probeL,'A')
        probeInputL = inputs{1};         % ER10xA microphone
        probeOutputL = outputs{1};       % ER10xA loudspeakers
    elseif strcmp(probeL,'B')
        probeInputL = inputs{1};         % ER10xB microphone
        probeOutputL = outputs{1};       % ER10xB loudspeakers
    elseif isempty(probeL)
        probeInputL = [];
        probeOutputL = [];        
    else
        error('Unrecognized probe value. Must be A or B.')
    end
    
    % for right probe
    if strcmp(probeR,'A')
        probeInputR = inputs{1};         % ER10xA microphone
        probeOutputR = outputs{1};       % ER10xA loudspeakers
    elseif strcmp(probeR,'B')
        probeInputR = inputs{1};         % ER10xB microphone
        probeOutputR = outputs{1};       % ER10xB loudspeakers
    elseif isempty(probeR)
        probeInputR = [];
        probeOutputR = [];
    else
        error('Unrecognized probe value. Must be A or B.')
    end

    if ~isempty(validateStimLevels)
        refInput = inputs{4};         % Reference condensor microphone
    else
        refInput = [];
    end

end