function [caParams] = ARLas_earCanalRecordings_DW10x(varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [caParams] = ARLas_earCanalRecordings_DW10x(varargin);
%
% This code is to be called by other experiment files, and is used to make 
% recordings from the ear canal for the purposes of an in-situ calibration. 
% This version (DW10x) is for the ER10X in Jeff Lichtenhan's lab. 
% It assumes that the 10x probe is attached to an ear bar is being used 
% to deliver the sound to guinea pigs.
%
%
% Authors: Shawn S Goodman, PhD & Jeffery T Lichtenhan PhD
% Auditory Research Lab, Dept. of Comm. Sciences & Disorders, The University of Iowa
% Auditory Physiology Lab, Dept. of Otolaryngology, Washington University School of Medicine, St. Louis
% Date: October 25, 2019
% Updated: October 25, 2019 -- ssg
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    obj = varargin{1};
    caParams = varargin{2};
    
% ---  Adjustable Parameteters --------------------------------------------
    nRepsChirp = 48;
% -------------------------------------------------------------------------
    
    % unpack the paramter structure
    fmin = caParams.fmin;
    fmax = caParams.fmax;
    testProbe = caParams.testProbe;
    usingChannels = caParams.usingChannels;
    usingChannels = sort(usingChannels);
    nChannels = length(usingChannels);

    fs = obj.fs; % get the system sampling rate

    [inputs,outputs] = hardwareSetup; %hardwareSetup; % read in the saved hardware setup
    if strcmp(testProbe,'A')
        input = inputs{1};        % ER10xA microphone
        output = outputs{1};      % ER10xA loudspeakers
    elseif strcmp(testProbe,'B')
        input = inputs{2};        % ER10xB microphone
        output = outputs{2};      % ER10xB loudspeakers
    end
    %inputRef = inputs{4};         % GRAS 1/8" reference microphone
    caParams.input = input;
    caParams.output = output;
    
    [pathName_mic,folderName_mic,fileName_mic] = mostRecentCalibration('mic',input.label,[]);
    [pathName_thev1,folderName_thev1,fileName_thev1] = mostRecentCalibration('thev',output.label,1);
    [pathName_thev2,folderName_thev2,fileName_thev2] = mostRecentCalibration('thev',output.label,2);
    caParams.pathName_mic = pathName_mic;
    caParams.folderName_mic = folderName_mic;
    caParams.fileName_mic = fileName_mic;
    
    dummy = load([pathName_mic,folderName_mic,fileName_mic]);
    micCorrection = dummy.micCorrection; % microphone correction
    caParams.micCorrection = micCorrection;
    
    for jj=1:nChannels % loop across both output channels
    
        % LOAD AND CHECK THE THEVENIN CALIBRATION FILE
        stimGain_dB = -5; % reduce so you don't overdrive the system; for click
        if usingChannels(jj) == 1
            [t,cut1,cut2,stimulus] = loadCal(pathName_thev1,folderName_thev1,fileName_thev1,stimGain_dB,fmin,fmax,obj);
        elseif usingChannels(jj) == 2
            [t,cut1,cut2,stimulus] = loadCal(pathName_thev2,folderName_thev2,fileName_thev2,stimGain_dB,fmin,fmax,obj);
        else
            warning('Unrecognized channel number.')
            keyboard
        end

        % LOAD THE STIMULUS ---------------------------------------------------
        obj.clearRecList % clear the previously used recording list
        obj.setRecList(input.ch,input.label,input.micSens,input.gain); % load the recording info for ER10X
        obj.setNReps(nRepsChirp); % number of times to play stimulus
        obj.setFilter(1); % turn on default highpass filter

        % PLAYBACK & RECORD ---------------------------------------------------
        obj.clearPlayList % clear the previously used playback list
        obj.setPlayList(stimulus,output.ch(jj)); % load the currently tested output channel
        obj.objPlayrec.run % run the stimulus
        if obj.killRun
           return
        end    

        % RETRIEVE DATA -------------------------------------------------------
        [header,data] = obj.retrieveData(input.label); % get chirp recorded in the ear canal
        data = fastFilter(micCorrection,data);
        [recording,stimf] = cleanDataDW(data,fs,stimulus,header,cut1,cut2);

        % IN-SITU CALIBRATION -------------------------------------------------
        %isc = ARLas_inSituCal(t,recording,stimf); % perform an in-situ calibration
        isc = ARLas_inSituCal_DW(t,recording,stimf); % perform an in-situ calibration (this version forces use of 0.3 cm diameter
        caParams.isc = isc;
    end
end

% Internal Functions ------------------------------------------------------
function [t,cut1,cut2,stimulus] = loadCal(pathName,folderName,fileName,stimGain_dB,fmin,fmax,obj)
    dummy = load([pathName,folderName,fileName]);
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
    stimulus = stimulus * 10^(stimGain_dB/20); % scale down so don't overdrive output
    cut1 = t.cut1;         % get the original cut samples
    cut2 = t.cut2;
end
