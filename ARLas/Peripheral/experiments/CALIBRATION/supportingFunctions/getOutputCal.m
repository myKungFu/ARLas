function [C1,C2,calPath1,calPath2] = getOutputCal(calType,probeOutput,map)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%[C1,C2,calPath1,calPath2] = getOutputCal(calType,probeOutput,map);
%
% Get the most recent output calibration files.
%
% calType = 'thev'; % use 'thev' or 'long' for Thevenin source or long tube, respectively
% probeOutput = structure describing probe. Obtained from hardwareSetup.m;
%
% Author: Shawn Goodman
% Date: August 23, 2021
% Last Updated: September 13, 2021 -- ssg -- fixed error in LTC; typo
%                                       applied calPath1 to calPath2
% Last Updated: March 23, 2022 -- ssg -- added map as an input argument.
%                                   This allows users to map MAC and other
%                                   non-standard directories for ARLas.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if nargin < 3 % for backwards compatibility (3/23/2022) -- ssg
        map = [];
    end
    if isempty(probeOutput)
        C1 = [];
        C2 = [];
        calPath1.pathName = [];
        calPath1.folderName = [];
        calPath1.fileName = [];
        calPath2.pathName = [];
        calPath2.folderName = [];
        calPath2.fileName = [];
        return
    end
    if strcmp(calType,'thev')
        [pathName,folderName,fileName] = mostRecentCalibration('thev',probeOutput.label,1,map);
        calPath1.pathName = pathName;
        calPath1.folderName = folderName;
        calPath1.fileName = fileName;
        [pathName,folderName,fileName] = mostRecentCalibration('thev',probeOutput.label,2,map);
        calPath2.pathName = pathName;
        calPath2.folderName = folderName;
        calPath2.fileName = fileName;
    elseif strcmp(calType,'long')
        [pathName,folderName,fileName] = mostRecentCalibrationLT('longTube',probeOutput.label,1,map);
        calPath1.pathName = pathName;
        calPath1.folderName = folderName;
        calPath1.fileName = fileName;
        [pathName,folderName,fileName] = mostRecentCalibrationLT('longTube',probeOutput.label,2,map);
        calPath2.pathName = pathName;
        calPath2.folderName = folderName;
        calPath2.fileName = fileName;
    else
        error('Unrecognized calibration type. Must be thev or long.')
    end
    % load the calibration files
    dummy = load([calPath1.pathName,calPath1.folderName,calPath1.fileName]);
    if strcmp(calType,'thev')
        C1 = dummy.t;
    elseif strcmp(calType,'long')
        C1 = dummy.LTC;
    end
    dummy = load([calPath2.pathName,calPath2.folderName,calPath2.fileName]);
    if strcmp(calType,'thev')
        C2 = dummy.t;
    elseif strcmp(calType,'long')
        C2 = dummy.LTC;
    end
end
