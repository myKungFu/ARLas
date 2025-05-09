function [] = ARLas_calCheck(varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ARLas_calCheck(varargin)
%
% Perform calibration check.
% This code is for periodic checks to make sure that calibration has not changed.
% These checks should be made using the same probe tips placed in the same 
% calibration cavities. 
% If the recordings are nearly identical to previous recordings, this indicates
% no change in either the microphone and the loudspeakers.
%
% If any change is observed, then a more thorough "diagnostic"
% recalibration is necessary. Changes in calCheck do NOT necessarily tell 
% you whether the microphone has changed or whether the loudspeakers, or both.
% 
% Saves results into a structure called Checks,
% located in 'C:\myWork\ARLas\Peripheral\calibrations\calChecks\' in a
% sub-folder under the date of the calibration.
%
% Author: Shawn Goodman, PhD
% Date: March 15, 2022
% Last Updated: March 15, 2022 -- ssg
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    obj = varargin{1}; % get the arlas object
    V = '15MAR2022'; % this is the current version number of this program
    
    %------ USER MODIFIABLE PARAMETERS ----------------------------------------
    probe = 'B'; % which ER10X probe to use (A = left; B = right)
    cavityName= 'longTube'; % Name of calibration cavity. 
                     % This is for reference only. Be sure that you always
                     % use the same cavity and probe tip!!!
    tipColor = 'red';
    
    %--------------------------------------------------------------------------

    fmin = 100; % minimum frequency to include
    fmax = 20000; % maximum frequency to include
    applyMicCorrection = 1; % turn on (1) and off (0) mic correction
    nReps = 64; % Number of stimulus repetitions to record
    doPlot = 1; % turn on (1) and off (0) plotting

    %pathName_calChecks = 'C:\myWork\ARLas\Peripheral\calibrations\calChecks\';
    pathName_calChecks = [obj.map.calibrations,'calChecks',obj.sep];
    if exist(pathName_calChecks,'dir') ~= 7 % check to make sure that the save path exists
        success = mkdir(pathName_calChecks); % if not, try to create it
        if success ~= 1
            warning('Specified folder for saving data does not exist. Entering debug mode.')
            keyboard
        else
            addpath(genpath(pathName_calChecks)) % add all directories and subdirectories
        end
    end    
    
    [inputs,outputs] = hardwareSetup;
    if strcmp(probe,'A')
        probeInput = inputs{1};         % ER10xA microphone
        probeOutput = outputs{1};       % ER10xA loudspeakers
    elseif strcmp(probe,'B')
        probeInput = inputs{2};         % ER10xB microphone
        probeOutput = outputs{2};       % ER10xB loudspeakers
    end
    
    disp(' '); disp(' '); disp(' '); disp(' ')
    alertTxt = {'Starting ARLas Calibration Check.'
         ['  Calibrating probe ',probe,'.']
         ['  Place ER10X probe with ',tipColor,' tip in cavity: ',cavityName,'.']
        };
    nn = size(alertTxt,1);
    for ii=1:nn
        cprintf([0,0,.4],[alertTxt{ii},'\n']);
    end

    % -------------------------------------------------------------------------
    % get most recent calibration files (mic only)
    [pathName_mic,folderName_mic,fileName_mic] = mostRecentCalibration('mic',probeInput.label,[]);
    if applyMicCorrection == 1
        if ~isempty(fileName_mic)
            dummy = load([pathName_mic,folderName_mic,fileName_mic]);
            micCorrection = dummy.micCorrection; % microphone correction
        else
            micCorrection = 1; % convolving by this changes nothing
        end
    else
        micCorrection = 1; % convolving by this changes nothing
    end
    %--------------------------------------------------------------------------
    
    % create stimuli ----------------------------------------------------------
    fs = obj.fs; % get the system sampling rate
    len = 0.05; % stimulus length (sec) -->use 100 ms
    stimN = round(len * fs); % number of samples in stimulus
    stim = zeros(stimN,1); % initialze to all zeros
    stimOffset_sec = 0.003; % stimulus offset from time zero (sec)
    stimOffset_samples = round(stimOffset_sec * fs);
    stim(stimOffset_samples,1) = 1; % impulse click, full output
    multiplier = 0.5;
    stim = stim * multiplier; % scale output to desired level
    
    nChannelsOut = length(probeOutput.ch); % number of output channnels being calibrated (typically 2)
    for kk=1:nChannelsOut % loop over output channels; test output channels one at a time
        alertTxt = {['     Calibrating loudspeaker ',num2str(kk)]};
        nn = size(alertTxt,1);
        for ii=1:nn
            cprintf([0,0,.4],[alertTxt{ii},'\n']);
        end
    
        % COLLECT THE DATA ------------------------------------------------
        obj.clearRecList % clear out whatever was used previously 
        obj.setRecList(probeInput.ch,probeInput.label,probeInput.micSens,probeInput.gain);
        obj.setNReps(nReps); % number of times to play stimulus
        obj.setFilter(0); % note: this is a highpass filter with a 75 Hz cutoff frequency.

        obj.clearPlayList % clear out whatever was used previously
        obj.setPlayList(stim,probeOutput.ch(kk)); % load the currently tested output channel

        obj.objPlayrec.userInfo = []; % clear out previous, non-structure array info
        obj.objPlayrec.userInfo.multiplier = multiplier;
        obj.objPlayrec.userInfo.fs = fs;
        obj.objPlayrec.userInfo.stim = stim;
        obj.objPlayrec.userInfo.caps_version = V;

        obj.objPlayrec.run % playback and record
        if obj.killRun
           return
        end

        % RETRIEVE DATA ---------------------------------------------------
        [header,Data] = obj.retrieveData(probeInput.label); % get raw data
        if applyMicCorrection == 1
            Data = fastFilter(micCorrection,Data);
        end
        
        [Recording,stimulusFiltered,Data,time] = cleanData(Data,fs,stim,header);
        
        originalN = size(Data,1);
        nfft = 1024;
        ref = 0.00002;
        [frequency,signal,noiseFloor] = ARLas_fda(Data,fs,ref,nfft,originalN);
        
        fileName = [probeOutput.label,'_Ch',num2str(kk),'_1.mat'];
        out = exist([pathName_calChecks,fileName]);
        if out == 0 % if file does not exist, create a new one
            counter = 1;
        elseif out == 2 % othewise, load the existing file
            dummy = load([pathName_calChecks,fileName]);
            Check = dummy.Check;
            counter = Check.counter + 1;
        else
            warning('File Not Found!')
            keyboard
        end
        
        field = ['recording',num2str(counter)];
        matrix = [frequency,signal,noiseFloor];
        Check.(field).('data') = matrix;
        Check.(field).('timeStamp') = header.timeStamp;
        Check.('counter') = counter;
        
        save([pathName_calChecks,fileName],'Check')
        
        if doPlot == 1
            figure; hold on
            for jj=1:Check.counter-1
                field = ['recording',num2str(jj)];
                matrix = Check.(field).('data');
                frequency = matrix(:,1);
                signal = matrix(:,2);
                noiseFloor = matrix(:,3);
                plt(frequency/1000,signal,'Color',[1 .7 .7])
                hold on
                plt(frequency/1000,noiseFloor,'Color',[.7 .7 .7])       
            end
            field = ['recording',num2str(counter)];
            matrix = Check.(field).('data');
            frequency = matrix(:,1);
            signal = matrix(:,2);
            noiseFloor = matrix(:,3);            
            plt(frequency/1000,signal,'r')
            plt(frequency/1000,noiseFloor,'k')
            xlim([fmin/1000,fmax/1000])
            xlabel('Frequency (kHz)')
            ylabel('Magnitude (dB SPL)')
            title(['Calibration Check  ',probeOutput.label,' Ch',num2str(kk),'   ',header.timeStamp])
        
        end
    end
end

% INTERNAL FUNCTIONS ------------------------------------------------------
function [data,stimulus,Data,time] = cleanData(X,fs,stim,header)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [data,stimulus,Data,time] = cleanData(X,fs,stim,header);
%
% For use with ARLas (Auditory Research Laboratory auditory software)
% Used in applying long tube (LT) calibration to cavity and/or ear canal recordings.
%
% X (input argument) = matrix of input waveforms
% fs = sampling rate (Hz)
% stim (input argument) = waveform vector; the original electrical stimulus that was used to obtain X
% cut1 = cutoff value to remove excess zero padding at the start of the signal
% cut2 = cutoff value to remove excess zero padding at the end of the signal
%
% X (output argument) = the lowpass Butterworth filtered (22 kHz cutoff), 
%       truncated waveform w/ artifact rejection and 3 ms onset/offset windowing
% stim (output argument) the stimulus with the same filtering and windowing applied.
%
% Auditory Research Lab, The University of Iowa
% Deptartment of Communication Sciences & Disorders
% The University of Iowa
% Author: Shawn Goodman
% Original Date: October 16, 2018; cleanData.m
% Adapted June 7, 2021 -- ssg
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    [rows,cols] = size(X); % orignal size of the recording
    % it is important to do the same processing to the stimulus as to the
    % recorded waveforms.
    stim = repmat(stim,1,cols+2); % replicate stimulus to matrix the same size as recording
    stim = stim(:); % reshape into a long vector
    X = X(:);
    stim = filtfilt(header.SOS,header.G,stim); % apply the playrec highpass filter to stim (same as was for recording)
    [SOS,G] = getLPiir(fs); % get coefficients for a lowpass filter
    stim = filtfilt(SOS,G,stim); % apply the lowpass filter
    X = filtfilt(SOS,G,X);
    X = reshape(X,rows,cols); % reshape into matrices
    stim = reshape(stim,rows,cols+2);
    stimulus = stim(:,2:end-1);
    stimulus = mean(stimulus,2);

    Data = ARLas_artifactReject(X); % perform basic artifact reject
    data = mean(Data,2); % averaged data
    [~,indx] = max(abs(data)); % location of time zero
    preLen = 0.001; % amount of time to include before time zero
    postLen = 0.004; % amount of time to include after time zero
    preSamples = round(preLen * fs); % number of samples before indx zero
    postSamples = round(postLen * fs); % number of samples after indx zero
    
    data = data(indx-preSamples:indx+postSamples); % truncate data; averaged vector
    Data = Data(indx-preSamples:indx+postSamples,:); % truncate data; full data set
    stimulus = stimulus(indx-preSamples:indx+postSamples);
    N = length(data); % number of samples in truncated vector
    if mod(N,2) ~= 0 % force length to be even
        data = data(1:end-1);
        Data = Data(1:end-1,:);
        stimulus = stimulus(1:end-1);
        N = length(data);
    end
    
    data = ARLas_ramp(data,fs,0.001); % ramp onset and offset
    Data = ARLas_ramp(Data,fs,0.001);
    n = (0:1:N-1)';
    time = n/fs*1000;
    time = time - (preLen*1000);
end

% OLD CODE ----------------------------------------------------------------
%         
%         
%     t = datetime('now'); % create the folder name for saving calibrations
%     folderName = datestr(t,'mm_dd_yyyy');
%     folderName = [folderName,'\']; 
%     if exist([pathName_LTC,folderName],'dir') ~= 7 % check to make sure that the save path exists
%         success = mkdir([pathName_LTC,folderName]); % if not, try to create it
%         if success ~= 1
%             warning('Specified folder for saving data does not exist. Entering debug mode.')
%             keyboard
%         else
%             addpath(genpath([pathName_LTC,folderName])) % add all directories and subdirectories
%         end
%     end
%     pathName = [pathName_LTC,folderName];
%     
%     channel = 1;
%     fileName = [LTC1.probeOutput.label,'_Ch',num2str(channel),'.mat'];
%     fileName = ARLas_saveName(pathName,fileName);
%     Check = LTC1;
%     save([pathName,fileName],'Check');
%     
%     if nChannelsOut == 2
%         channel = 2;
%         fileName = [LTC2.probeOutput.label,'_Ch',num2str(channel),'.mat'];
%         fileName = ARLas_saveName(pathName,fileName);
%         Check = LTC2;
%         save([pathName,fileName],'Check');
%     end

% keyboard    
%     % Save Check here -----------------------------
%     alertTxt = {['  Saving calibrations']};
%     nn = size(alertTxt,1);
%     for ii=1:nn
%         cprintf([0,0,.4],[alertTxt{ii},'\n']);
%     end
% 
%     pathName_LTC = 'C:\myWork\ARLas\Peripheral\calibrations\LTCals\';
%     if exist(pathName_LTC,'dir') ~= 7 % check to make sure that the save path exists
%         success = mkdir(pathName_LTC); % if not, try to create it
%         if success ~= 1
%             warning('Specified folder for saving data does not exist. Entering debug mode.')
%             keyboard
%         else
%             addpath(genpath(pathName_LTC)) % add all directories and subdirectories
%         end
%     end
%     
%     t = datetime('now'); % create the folder name for saving calibrations
%     folderName = datestr(t,'mm_dd_yyyy');
%     folderName = [folderName,'\']; 
%     if exist([pathName_LTC,folderName],'dir') ~= 7 % check to make sure that the save path exists
%         success = mkdir([pathName_LTC,folderName]); % if not, try to create it
%         if success ~= 1
%             warning('Specified folder for saving data does not exist. Entering debug mode.')
%             keyboard
%         else
%             addpath(genpath([pathName_LTC,folderName])) % add all directories and subdirectories
%         end
%     end
%     pathName = [pathName_LTC,folderName];
%     
%     channel = 1;
%     fileName = [LTC1.probeOutput.label,'_Ch',num2str(channel),'.mat'];
%     fileName = ARLas_saveName(pathName,fileName);
%     Check = LTC1;
%     save([pathName,fileName],'Check');
%     
%     if nChannelsOut == 2
%         channel = 2;
%         fileName = [LTC2.probeOutput.label,'_Ch',num2str(channel),'.mat'];
%         fileName = ARLas_saveName(pathName,fileName);
%         Check = LTC2;
%         save([pathName,fileName],'Check');
%     end
%     
%     % do plotting here ------------------------------------------
%     if doPlot == 1
%         ch = 1;
%         h1 = makePlots(LTC1,ch);
%         if nChannelsOut == 2
%             ch = 2;
%             h2 = makePlots(LTC2,ch);
%         end
%     end
%     
%     alertTxt = {'Long Tube calibration finished.'};
%     nn = size(alertTxt,1);
%     for ii=1:nn
%         cprintf([0,0,.4],[alertTxt{ii},'\n']);
%     end
%     disp(' ')    
%       
%         %------------------------------------------------------------------
%         if kk == 1
%             [fo_spl,phi,freq,Amax,pRef,nf] = XferFunction(Data,stimulusFiltered,fs,fmin,fmax);
%             Check1 = buildStruct(stim,time,stimulusFiltered,Recording,Data,probeInput,probeOutput,...
%                 Amax,pRef,nReps,len,fs,level,fmin,fmax,micCorrection,cavityDiameter,...
%                 cavityTemperature,z0,fo_spl,phi,freq,nf,V);
%         elseif kk == 2
%             [fo_spl,phi,freq,Amax,pRef,nf] = XferFunction(Data,stimulusFiltered,fs,fmin,fmax);
%             Check2 = buildStruct(stim,time,stimulusFiltered,LTRecording,Data,probeInput,probeOutput,...
%                 Amax,pRef,nReps,len,fs,level,fmin,fmax,micCorrection,cavityDiameter,...
%                 cavityTemperature,z0,fo_spl,phi,freq,nf,V);
%         end
%     end
