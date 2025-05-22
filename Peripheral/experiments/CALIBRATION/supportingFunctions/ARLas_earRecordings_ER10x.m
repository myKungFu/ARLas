function [isc] = ARLas_earRecordings_ER10x(varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ARLas_earRecordings_ER10x(varargin);
% ARLas_earRecordings_ER10x(obj,input,output,calPathName,calFileName);
%
% Make recordings from ear canal using ER10x. Thse recordings used
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    try
        obj = varargin{1};
        input = varargin{2};
        output = varargin{3};
        calPathName = varargin{4};
        calFileName = varargin{5};
    catch
        error('Unexpected input arguments in ARLas_earRecordings_ER10x.')
    end

    try
        params = varargin{6};
        fmin = params.fmin;
        fmax = params.fmax;
        nReps = params.nReps;
        micCorrection = params.micCorrection;
    catch
        fmin = 200;
        fmax = 18000;
        nReps = 48;
        micCorrection = [];
    end
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
 
% use reference mic for comparison
%inputs = hardwareSetup;
%obj.setRecList(inputs{3}.ch,inputs{3}.label,inputs{3}.micSens,inputs{3}.gain);    
    
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

% % compare reference to probe mic
% [headerRef,DataRef] = obj.retrieveData(inputs{3}.label);
% d = mean(Data,2);
% dref = mean(DataRef,2);
% D = fft(d);
% Dref = fft(dref);
% fff = (0:1:length(d)-1)'*(obj.fs/length(d))/1000;
% plt(fff,20*log10(abs(D./Dref)))
% xlim([0 20])


    if ~isempty(micCorrection)
        Data = applyMicCorrection(Data,micCorrection);
    end    
    
    [earRecording,stimf] = cleanData(Data,obj.fs,stimulus,header,cut1,cut2);
    
    isc = ARLas_inSituCal(t,earRecording,stimf); % calculate the load impedance and load characteristic impedance

end

% internal functions ------------------------------------------------------
