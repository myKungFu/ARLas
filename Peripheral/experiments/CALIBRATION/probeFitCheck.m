function [iscCheckHandleA,iscCheckHandleB,OK] = probeFitCheck(varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [iscCheckHandleA,iscCheckHandleB,OK] = probeFitCheck(obj,LC,iscCheckHandleA,iscCheckHandleB)
%
% Check the probe fit. Usually called by loadCalibration.m
% 
%
% Author: Shawn Goodman, PhD
% Auditory Research Lab, the University of Iowa
% Original Date: May 20, 2025
% Last Updated: May 20, 2025 -- ssg
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    obj = varargin{1}; % get the arlas object
    LC = varargin{2};
    iscCheckHandleA = varargin{3};
    iscCheckHandleB = varargin{4};    
    V = '20MAY2025'; % this is the current version number of this program

%------ USER MODIFIABLE PARAMETERS ----------------------------------------
%--------------------------------------------------------------------------

%------ END USER MODIFIABLE PARAMETERS ------------------------------------
%--------------------------------------------------------------------------

    probeOutputA = LC.probeOutputA;
    probeOutputB = LC.probeOutputB;
    OK = 0; % initialize 
    if isempty(iscCheckHandleA) & isempty(iscCheckHandleB)
        questionMe = 1;
    else
        questionMe = 0;
    end

    try
        if ~isempty(probeOutputA)
            labelA = probeOutputA.label;
            dA = dir([obj.savingGrace,'*',labelA,'_inSituCal_*']);
        else
            dA = [];
        end
        if ~isempty(probeOutputB)
            labelB = probeOutputB.label;
            dB = dir([obj.savingGrace,'*',labelB,'_inSituCal_*']);
        else
            dB = [];
        end
        %iscCheckHandleA = [];
        if ~isempty(probeOutputA)
            dummy = load([obj.savingGrace,dA(1).name]); % first calibration chirp (loudspeaker 1)
            iscCheckHandleA = inSituCheck(probeOutputA.label,probeOutputA.ch(1),dummy.header,dummy.data,obj.fs,1,iscCheckHandleA);
            dummy = load([obj.savingGrace,dA(2).name]); % first calibration chirp (loudspeaker 2)
            iscCheckHandleA = inSituCheck(probeOutputA.label,probeOutputA.ch(2),dummy.header,dummy.data,obj.fs,2,iscCheckHandleA);
        end
        %iscCheckHandleB = [];
        if ~isempty(probeOutputB)
            dummy = load([obj.savingGrace,dB(1).name]); % first calibration chirp (loudspeaker 1)
            iscCheckHandleB = inSituCheck(probeOutputR.label,probeOutputR.ch(1),dummy.header,dummy.data,obj.fs,1,iscCheckHandleB);
            dummy = load([obj.savingGrace,dB(2).name]); % first calibration chirp (loudspeaker 2)
            iscCheckHandleB = inSituCheck(probeOutputR.label,probeOutputR.ch(2),dummy.header,dummy.data,obj.fs,2,iscCheckHandleB);
        end
    catch ME
        disp('Start of measurement in-situ calibration check failed!')
        disp('probeFitCheck.m')
    end

    if questionMe == 1
        answer = questdlg('Peaks at low and high frequencies should be between 1 and 2 Pa. Do NOT close figures!',...
            'In-situ Calibration Check',...
            'Continue','Abort','Continue');
        switch answer
            case 'Continue'
                OK = 1;
            case 'Abort'
                OK = 0;
        end
    end
    
end

