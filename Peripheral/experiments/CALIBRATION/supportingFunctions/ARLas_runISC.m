function [iscS1,iscS2,iscS12] = ARLas_runISC(obj,probeInput,probeOutput,calPath1,calPath2,calType,fmin,fmax,micCorrection,nReps,doIndividual,doSimultaneous)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [iscS1,iscS2,iscS12] = ARLas_runISC(obj,probeInput,probeOutput,calPath1,calPath2,calType,fmin,fmax,micCorrection,nReps,doIndividual,doSimultaneous)
%
% Make recordings from ear canal using ER10x. Thse recordings used
% along with Thevenin source parameters to calculate FPL, RPL, and WR.
%
% Input Arguments:
% obj = arlas object
% probeInput = structure describing probe input. Get from hardwareSetup
% probeOutput = structure describing probe output. Get from hardwareSetup
% calPath1 = location of Thevenin calibration file for reciever 1
% calPath2 = location of Thevenin calibration file for reciever 2
% calType = string describing calibration type. Should be 'thev'
% fmin = minimum frequency over which to calibrate (Hz)
% fmax = maximum frequency over which to calibrate (Hz)
% micCorrection = 0 or 1, turns on and off microphone correction
% nReps = number of repetitions to include in average
% doIndividual = 0 or 1, turns on and off calibration of separate recievers
% doSimultaneous = 0 or 1, turns on and off calibration of both receivers simultaneously
%
% OUTPUT ARGUMENTS:
% iscS1 = structure containing in-situ calibration information for receiver 1
% iscS2 = structure containing in-situ calibration information for receiver 2
% iscS12 = structure containing in-situ calibration information for both receivers simultaneously
%
% Auditory Research Lab, The University of Iowa
% Deptartment of Communication Sciences & Disorders
% The University of Iowa
% Author: Shawn Goodman
% Date: July 7, 2021
% Last Updated: July 7, 2021 -- ssg 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    iscS1 = [];
    iscS2 = [];
    iscS12 = [];
    if isempty(probeInput)
        return
    end
    if strcmp(calType,'thev')
        % ARLas_earRecordings10x(obj,input,      output,      calPathName,calFileName);
        params.fmin = fmin;
        params.fmax = fmax;
        params.micCorrection = micCorrection;
        params.nReps = nReps;
        
        if ~isempty(probeOutput)
            if doIndividual == 1 % calibrate the output of each channel individually
                po = probeOutput; % hack to make this work each output channel separately
                po.ch = probeOutput.ch(1);
                iscS1 = ARLas_earRecordings_ER10x(obj,probeInput,po,[calPath1.pathName,calPath1.folderName],calPath1.fileName,params); 
                po.ch = probeOutput.ch(2);
                iscS2 = ARLas_earRecordings_ER10x(obj,probeInput,po,[calPath2.pathName,calPath2.folderName],calPath2.fileName,params); 
            end
            if doSimultaneous == 1 % calibrate the output for both channels presented simultaneously
                po = probeOutput; 
                iscS12 = ARLas_earRecordings_ER10x(obj,probeInput,po,[calPath1.pathName,calPath1.folderName],calPath1.fileName,params); 
            end                
        else
            iscS1 = [];
            iscS2 = [];
            iscS12 = [];
        end
    end
    
end
