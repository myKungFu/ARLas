function [LC] = getMostRecentLoadCalibration(varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LC = getMostRecentLoadCalibration(obj);
%
% Perform Thevenin-based load calibration.
% Save the results as well as return them as output arguments
%
% Author: Shawn Goodman, PhD
% Auditory Research Lab, the University of Iowa
% Original Date: May 20, 2025
% Last Updated: May 22, 2025 -- ssg
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    obj = varargin{1}; % get the arlas object
    V = '22MAY2025'; % this is the current version number of this program

%------ USER MODIFIABLE PARAMETERS ----------------------------------------
%--------------------------------------------------------------------------

%------ END USER MODIFIABLE PARAMETERS ------------------------------------
%--------------------------------------------------------------------------
    
    LC = [];
    % to get the most recent file, do this:
    pathName = [obj.map.data,obj.subjectID,'\','loadCals\'];
    d = dir([pathName,'*.mat']);
    nFiles = length(d);
    %Q = [];
    for ii=1:nFiles
        Q{ii,:} = d(ii).date; 
    end
    if isempty(Q)
        disp(['No valid load calibration found for this subject.'])
        disp(['Run the program "loadCalibration.m" first.'])
        return
    end
    [pick,iPick] = max(datetime(Q));
    recall = load([pathName,d(iPick).name]);
    LC = recall.LC;
    

end

% INTERNAL FUNCTIONS ------------------------------------------------------

% OLD CODE ----------------------------------------------------------------

    % % specify calibration to use to set stimulus levels
    % calType = 'thev'; % use 'thev' for Thevenin source
    % targetCalType = 'fpl'; % 'fpl', 'spl', 'ipl'.

    % % Define which PROBE is being used (regardless of which ear it is placed in)
    % if ~(isempty(probeA) | isempty(probeB)) % if both ears are being tested
    %     testProbe = 'Both';
    % elseif ~(isempty(probeA))
    %     testProbe = 'Left';
    % elseif ~(isempty(probeB))
    %     testProbe = 'Right';
    % end

    % Define which EAR the probe is placed in (regardless of probe used).
    % Note that if both probes are used, it is assumed that A=Left and
    % B=Right
    % ISCA.(['Rec',num2str(1)]).(['S1']) = iscSA1;
    % ISCA.(['Rec',num2str(1)]).(['S2']) = iscSA2;
    % ISCB.(['Rec',num2str(1)]).(['S1']) = iscSB1;
    % ISCB.(['Rec',num2str(1)]).(['S2']) = iscSB2;            
