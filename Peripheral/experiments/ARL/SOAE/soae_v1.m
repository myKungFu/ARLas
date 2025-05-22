function [] = soae_v1(varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% soae_v1(varargin);
%
% Measure Spontaneous Otoacoustic Emissions
%
% Auditory Research Lab, The University of Iowa
% Deptartment of Communication Sciences & Disorders
% The University of Iowa
% Author: Shawn Goodman 
% Date: June 5, 2023
% Last Updated: June 5, 2023
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    obj = varargin{1}; % get the arlas object
    V = '05JUN2023'; % this is the current version number of this program

    %------ USER MODIFIABLE PARAMETERS ----------------------------------------
    %--------------------------------------------------------------------------

    totalLen = 60; % total recording length (sec)
    chunkLen = 0.5; % chunk length (sec)
    
    % specify which ER10X probes to use    
    probeL = 'A'; % probe that is in LEFT ear; 'A' or []
    probeR = []; % probe that is in RIGHT ear; 'B' or []
                  % Note: Either one or both ears can be tested at the same time.
                  % If either probe is empty set, will not be tested. 

    
    %------ END USER MODIFIABLE PARAMETERS ------------------------------------
    %--------------------------------------------------------------------------

    fmin = 100; % min frequency for calibration (Hz)
    fmax = 20000; % max frequency for calibration (Hz)
    % specify calibration to use to set stimulus levels
    calType = 'thev'; % use 'thev' for Thevenin source
    targetCalType = 'fpl'; % 'fpl', 'spl', 'ipl'
    % additional parameters that user usually doesn't need to mess with
    doMicCorrection = 1; % turn on and off microphone correction
    
    % get initial stimuli
    fs = obj.fs; % get the system sampling rate
    chunkN = round(fs * chunkLen);
    nReps = ceil(totalLen / chunkLen);
    stim = zeros(chunkN,1);
    
    %----------------------------------------------------------------------
    if isempty(obj)
        tempobj = struct;
        tempobj.map = [];
    else
        tempobj = obj;
    end    
    % get probe settings
    [probeInputL,probeOutputL,probeInputR,probeOutputR,refInput] = getProbeSettings(probeL,probeR,[]);
    
    % get most recent calibration files -----
    micCorrectionL = getMicCorrection(probeInputL,doMicCorrection,tempobj.map);
    micCorrectionR = getMicCorrection(probeInputR,doMicCorrection,tempobj.map);
    [C1L,C2L,calPath1L,calPath2L] = getOutputCal(calType,probeOutputL,tempobj.map);
    [C1R,C2R,calPath1R,calPath2R] = getOutputCal(calType,probeOutputR,tempobj.map);
    
    % PERFORM IN-SITU CALIBRATION -----
    disp('----- Running in-situ calibration -----')
    inSituReps = 6;
    doIndividual = 1;
    doSimultaneous = 0;
    [iscS1L,iscS2L,iscS12L] = ARLas_runISC(obj,probeInputL,probeOutputL,calPath1L,calPath2L,calType,fmin,fmax,micCorrectionL,inSituReps,doIndividual,doSimultaneous);
    [iscS1R,iscS2R,iscS12R] = ARLas_runISC(obj,probeInputR,probeOutputR,calPath1R,calPath2R,calType,fmin,fmax,micCorrectionR,inSituReps,doIndividual,doSimultaneous);

   % COLLECT THE DATA --------------------------------------------------------
    disp('----- Collecting data -----')
    nBlocks = 2;
    for ii=1:nBlocks % loop over blocks
        
        % Get the corrected click stimuli ---------------------------------------------------------
        
        % Apply the in-situ calibration to the stimuli ----------------
        % tell ARLas what to record and what to play
        obj.clearPlayList % clear out whatever was played previously
        if ~isempty(probeInputL) % load new signals to play
            obj.setPlayList(stim,probeOutputL.ch(1));
            obj.setPlayList(stim,probeOutputL.ch(2));
        end
        if ~isempty(probeInputR)
           obj.setPlayList(stim,probeOutputR.ch(1));
           obj.setPlayList(stim,probeOutputR.ch(2));
        end

        obj.clearRecList % clear out whatever was used previously 
        if ~isempty(probeInputL)
            if ii==1
                probeInputL.label = [probeInputL.label,'_SOAE'];
            end
            obj.setRecList(probeInputL.ch,probeInputL.label,probeInputL.micSens,probeInputL.gain);    
        end
        if ~isempty(probeInputR)
           if ii==1
               probeInputR.label = [probeInputR.label,'_SOAE'];
           end
           obj.setRecList(probeInputR.ch,probeInputR.label,probeInputR.micSens,probeInputR.gain);    
        end

        obj.setNReps(nReps); % number of times to play stimulus
        obj.setFilter(1); % note: this is a highpass filter with a 75 Hz cutoff frequency.

        obj.objPlayrec.userInfo = []; % clear out previous, non-structure array info
        obj.objPlayrec.userInfo.fs = fs;
        obj.objPlayrec.userInfo.nReps = nReps;
        obj.objPlayrec.userInfo.fmin = fmin;
        obj.objPlayrec.userInfo.fmax = fmax;
        obj.objPlayrec.userInfo.moc_version = V;
        obj.objPlayrec.userInfo.iscS1L = iscS1L;
        obj.objPlayrec.userInfo.iscS2L = iscS2L;
        obj.objPlayrec.userInfo.iscS1R = iscS1R;
        obj.objPlayrec.userInfo.iscS2R = iscS2R;
        obj.objPlayrec.userInfo.C1L = C1L;
        obj.objPlayrec.userInfo.C2L = C2L;
        obj.objPlayrec.userInfo.C1R = C1R;
        obj.objPlayrec.userInfo.C2R = C2R;
        obj.objPlayrec.userInfo.calType = calType;
        obj.objPlayrec.userInfo.targetCalType = targetCalType;
        obj.objPlayrec.userInfo.probeR = probeR;
        obj.objPlayrec.userInfo.probeL = probeL;
        obj.objPlayrec.userInfo.micCorrectionL = micCorrectionL; 
        obj.objPlayrec.userInfo.micCorrectionR = micCorrectionR; 


        disp('----- Recording SOAE data: ')
        obj.objPlayrec.run % playback and record -----
        if obj.killRun
           return
        end

        % extract the recordings -----
        if ~isempty(probeInputL)
            [headerL,DataL] = obj.retrieveData(['Ch',num2str(probeInputL.ch)]); % retrive recorded data
            if doMicCorrection == 1
                DataL = applyMicCorrection(DataL,micCorrectionL);
            end
            %[plL,PlL,phiL,otherL,wfL] = ARLas_convertPL(DataL,iscS1L);
        end


    % Analyze the data here:
    dataL = DataL(:);
    DataL = fft(dataL);
    nfft = length(dataL);
    DataL = abs(DataL)/(nfft/2);
    DataL = 20*log10(DataL/.00002);
    Frequency = (0:1:nfft-1)'*(fs/nfft);
    fmin = 500;
    fmax = 8000;
    [~,indxMin] = min(abs(fmin - Frequency));
    [~,indxMax] = min(abs(fmax - Frequency));
    DataL = DataL(indxMin:indxMax);
    Frequency = Frequency(indxMin:indxMax);
    
    figure(110)
    plt(Frequency,DataL)
    hold on
    xlabel('Frequency (Hz)')
    ylabel('Magnitude (db SPL)')
    pause(0.1)
    
keyboard    
    
    
    

    end

    
    disp('----- Finished with data collection -----')
    disp(' ')
    disp(' ')

end


% INTERNAL FUNCTIONS ------------------------------------------------------
function [obj,nLevels,fs,dL,dR,iscS1L,iscS2L,iscS12L,iscS1R,iscS2R,iscS12R,postPath,postPathSave] = getSavedData(basePath,SUBJ,DATE,EXPRUN,probeInputL,probeInputR)      
    os = '\';
    one = [basePath,SUBJ,os];
    two = [SUBJ,'_',DATE,os];
    three = [SUBJ,'_',DATE,'_',EXPRUN,os];
    four = [SUBJ,'_',DATE,'_',EXPRUN,'_analysis',os];
    postPath = [one,two,three];
    postPathSave = [one,two,three,four];

    % Make sure the analysis directory exists -----
    if exist(postPathSave,'dir') ~= 7 % check to make sure that the backup directory exists
        success = mkdir(postPathSave); % if not, try to create it
        if success ~= 1
            error('Backup directory does not exist and could not be created. Aborting program.')
        end
        addpath(genpath(postPathSave))
    end
    dL = dir([postPath,'Ch',num2str(probeInputL.ch),'_',probeInputL.label,'_moc_*.mat']);
    dR = dir([postPath,'Ch',num2str(probeInputR.ch),'_',probeInputR.label,'_moc_*.mat']);
    nLevels = size(dL,1); % will be used to set iterations of the loop
    fs = 96000;

    % get the saved in-situ calibrations
    dummy = load([postPath,dL(1).name]);
    q = dummy.header;
    try % this is the way it should be. 
        iscS1L = q.userInfo.iscS1L;
        iscS2L = q.userInfo.iscS2L;
    catch
        iscS1L = q.userInfo.iscS1; % work around: found error in saving on 7/10 (only saved right ear)
        iscS2L = q.userInfo.iscS2;
    end
    iscS12L = [];
    C1L = q.userInfo.C2L;
    C2L = q.userInfo.C2L;

    dummy = load([postPath,dR(1).name]);
    q = dummy.header;
    try
        iscS1R = q.userInfo.iscS1L;
        iscS2R = q.userInfo.iscS2L;
    catch
        iscS1R = q.userInfo.iscS1;
        iscS2R = q.userInfo.iscS2;
    end
    iscS12R = [];
    C1R = q.userInfo.C2R;
    C2R = q.userInfo.C2R;

    obj.subjectID = SUBJ;
    obj.timeStamp = q.timeStamp;
end
function b = bpf(Fs)
    Fstop1 = 0;               % First Stopband Frequency
    Fpass1 = 250;             % First Passband Frequency
    Fpass2 = 1000;            % Second Passband Frequency
    Fstop2 = 2000;            % Second Stopband Frequency
    Dstop1 = 0.001;           % First Stopband Attenuation
    Dpass  = 0.057501127785;  % Passband Ripple
    Dstop2 = 0.0001;          % Second Stopband Attenuation
    flag   = 'scale';         % Sampling Flag
    [N,Wn,BETA,TYPE] = kaiserord([Fstop1 Fpass1 Fpass2 Fstop2]/(Fs/2), [0 ...
                                 1 0], [Dstop1 Dpass Dstop2]);
    if mod(N,2)~=0
        N = N + 1;
    end
    b  = fir1(N, Wn, TYPE, kaiser(N+1, BETA), flag)';
end
function [micCorrection] = getMicCorrection(probeInput,doMicCorrection,map)
    if isempty(probeInput)
        micCorrection = [];
        return
    end
    [pathName_mic,folderName_mic,fileName_mic] = mostRecentCalibration('mic',probeInput.label,[],map);
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

% OLD CODE ----------------------------------------------------------------
