function [isc] = ARLas_earRecordings_IHS(varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ARLas_earRecordings_IHS(varargin);
% ARLas_earRecordings_IHS(obj,input,output,calPathName,calFileName);
%
% Make recordings from ear canal using the IHS systgem. Thse recordings used
% along with Thevenin source parameters to calculate FPL, RPL, and WR.
% Run this code separately for each output channel that is being tested.
% It is assumed that this function is being called by another experiment file.
% 
% obj = ARLas object
% input = Structure specifying inputs.
%          Must specify a single input channel.
% output = Structure specifying outputs.
%          Must specify a single output channel.
% calPathName = String specifying the directory location of the Thevenin
%               calibration files.
% calFileName = Structure containing strings specifying the Thevenin 
%               calibration file names.
% params = structure specifying calibration parameters:
%          params.fmin = minimum frequency to calibrate (Hz)
%          params.fmax = maximum frequency to calibrate (Hz)
%          params.nReps = number of stimulus reps (48 recommended)
%
% isc = in-situ calibration structure.
%
% Auditory Research Lab, The University of Iowa
% Deptartment of Communication Sciences & Disorders
% The University of Iowa
% Author: Shawn Goodman
% Date: December 1, 2017
% Updated: November 4, 2018 -- ssg
% Updated: July 2, 2021 -- ssg -- added mic correction
% Updated: July 7, 2021 -- ssg -- name change from ARLas_earRecordings10x.m
%                                 to ARLas_earRecordings_ER10x.m.
%                                 Made handle two recievers simultaneously.
% Updated: May 21, 2025 -- ssg -- name change to ARLas_earRecordings_IHS.m.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    try
        obj = varargin{1};
        input = varargin{2};
        output = varargin{3};
        calPathName = varargin{4};
        calFileName = varargin{5};
        params = varargin{6};
    catch
        error('Unexpected input arguments in ARLas_earRecordings_IHS.')
    end

    fmin = params.fmin;
    fmax = params.fmax;
    nReps = params.nReps;
    micCorrection = params.micCorrection;
    isc = [];
    
    % LOAD AND CHECK THE THEVENIN CALIBRATION FILE
    dummy = load([calPathName,calFileName]); % load the Thevenin object
    t = dummy.t;
    if t.fmin > fmin
        error('Requested value of fmin is less than the calibrated value.')
    end
    if t.fmax < fmax
        error('Requested value of fmax is greater than the calibrated value.')
    end
    if obj.fs ~= t.fs
        error('Current sampling rate is different from the calibrated value.')
    end
    stimulus = t.stimOrig; % get the logarithmic chirp originally used to calibrate
    cut1 = t.cut1;         % get the original cut samples
    cut2 = t.cut2;
    
    % LOAD THE STIMULUS ---------------------------------------------------
    obj.clearRecList % clear the previously used recording list
    input.label = [input.label,'_inSituCal'];
    obj.setRecList(input.ch,input.label,input.micSens,input.gain); % load the recording info for ARLas to use
 
    obj.setNReps(nReps); % number of times to play stimulus
    obj.setFilter(1); % turn on default highpass filter

    obj.clearPlayList % clear the previously used playback list
    obj.setPlayList(stimulus,output.ch(1)); % load the currently tested output channel
    if length(output.ch) == 2
        obj.setPlayList(stimulus,output.ch(2)); % load the currently tested output channel
    end
    
    % PLAYBACK & RECORD ---------------------------------------------------
    obj.objPlayrec.run % run the stimulus
    if obj.killRun
       return
    end    

    % RETRIEVE DATA ----------------------------------------------------
    [header,Data] = obj.retrieveData(input.label); % get chirp recorded in the ear canal

    if ~isempty(micCorrection)
        Data = applyMicCorrection(Data,micCorrection);
    end    
    
    [earRecording,stimf] = cleanData(Data,obj.fs,stimulus,header,cut1,cut2);
    
    isc = ARLas_inSituCal(t,earRecording,stimf); % calculate the load impedance and load characteristic impedance

end

% internal functions ------------------------------------------------------
