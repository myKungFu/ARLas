function [] = DPgram(varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DPgram(varargin)
%
% Measure DPOAE grams using a continuous frequency sweep.
%
% Tones are swept as described in Adbala, Luo, and Shera, 2015,
% "optimizing swept-tone protocols for recording DPOAEs in adults and newborns."
%
% Author: Shawn Goodman, PhD
% Auditory Research Lab, the University of Iowa
% Date: Original Date: March 7, 2022
% Last Updated: March 7, 2022 -- ssg
% Last Updated: July 19, 2024 -- ssg -- New version _vRT. Old version was vSSG November 17, 2023 
%                   With 24 averages, 6 levels, and just over 1 octave
%                   sweep rate of 0.5 sec/oct
%                   Test time was 6.16 minutes, which works out to 
%                   6.16 / 6 = 1 minute per level * nOctaves
%                   nMinutes = 1 minute * nOctaves * nLevels
% Last Updated: July 26, 2024 -- ssg -- updated pre and post cal checks
% Last Updated: January 3, 2025 -- ssg -- name changed from dpoae_fixedRatioSweep_vRT
%               to dp_sweepFrequency_vRT3.m.
% Last Updated: March 27, 2025 -- new version _vJR
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    obj = varargin{1}; % get the arlas object
    V = '27MAR2025'; % this is the current version number of this program

%--------------------------------------------------------------------------
%------ USER MODIFIABLE PARAMETERS ----------------------------------------
    nSweeps = 24; % 24 % Number of sweeps recorded at each level.
                    % Note: At sweep rate of 0.5 and 1-16 kHz, each sweep takes 8 sec.
                    %       Recommended nSweeps is 24, which should get noise floor to -10 dB SPL. This will take 24*8 = 192 sec = 3.2 min per
    
    % specify nominal stimulus parameters:
    %   Note: these will be modified slightly at low and high frequencies
    targetL1 = 65; %[65 65 65 65 65]; % target levels for primaries (dB FPL); Can be scalar or vector for more than one level
    targetL2 = 55; %[25 35 45 55 65]; % length of targetL2 must be same size as targetL1  

    fmin = 750; % min f2 value is 750 to ensure 1000 Hz (8kHz is a 1/2 octave below 10 kHz)
    fmax = 5000; % max f2 value is 18000 to ensure 16000 Hz (14000 is a 1/2 octave above 10 kHz)
    fRatio = 1.22; % ratio of f2/f1. Ratio will remain fixed at this value
    sweepRate = 0.5; % sweep rate in octaves/s; should be 0.5.
                    
% specify which ER10X probes to use    
    probeL = 'A'; %'A'; % probe that is in LEFT ear; 'A', 'B', or []
    probeR = []; %'B'; % probe that is in RIGHT ear; 'A', 'B', or []
                  % Note: Either one or both ears can be tested at the same time.
                  % If either probe is empty set, will not be tested. 
  
% additional parameters that user usually doesn't need to mess with
    doMicCorrection = 1; % turn on and off microphone correction
    doPreCalCheck = 1; % turn on and off the pre-test cal check
    doPostCalCheck = 1; % turn on and off the post-test cal check
    
% specify data location for post-hoc analysis, if using
        basePath = 'C:\myWork\ARLas\Data\'; % <-- be sure to end this one with a \ character, but not the others below it!
        SUBJ = 'dumdum';        % subject ID
        DATE = '19JUL2024';    % date data collected
        EXPRUN = 'dpoae_run7'; % experiment name and run number
        
%------ END USER MODIFIABLE PARAMETERS ------------------------------------
%--------------------------------------------------------------------------

% specify calibration to use to set stimulus levels
    calType = 'thev'; % use 'thev' for Thevenin source
    targetCalType = 'fpl'; % 'fpl', 'spl', 'ipl'.

    if doPostCalCheck==1 && doPreCalCheck==0
        error('If doPostCalCheck=1, doPreCalcheck must also =1.')
    end

    disp(' ')
    disp('----- Starting DPOAE Fixed Ratio Sweep experiment -----')
tic
    ISCL = struct;
    ISCR = struct;
    
    % check to make sure all parameters are valid
%    [targetL1,nSweeps,targetL2,sweepRate,probeR,probeL,calType,targetCalType] = checkInputs(targetL1,nSweeps,targetL2,sweepRate,probeR,probeL,calType,targetCalType);
    
    if doPostCalCheck==1 && doPreCalCheck==0
        error('If doPostCalCheck=1, doPreCalcheck must also =1.')
    end

    if ~isempty(obj) % if running post-hoc analysis, skip this part
        if ~(isempty(probeL) | isempty(probeR)) % if both ears are being tested
            testEar = 'Both';
        else
            % Get stimulus ear:
            prompt = {'Enter the test ear (Left or Right)'};
            title = 'Test Ear'; 
            defaultAns = {'Right'};
            numlines = 1;
            answer = inputdlg(prompt,title,numlines,defaultAns);
            if isempty(answer) % user hit cancel
                return
            elseif strcmp(answer{1},'') % if user left field blank
                return
            else % user put something in the field
                % check for legal input
                testEar = answer{1}; % voltage rms must be positive
                if strcmp(testEar,'left')
                    testEar = 'Left';
                end
                if strcmp(testEar,'right')
                    testEar = 'Right';
                end
                if ~(strcmp(testEar,'Left') | strcmp(testEar,'Right'))
                    disp('Error: Invalid input. Must be Left or Right.')
                    return
                end
            end
        end
    end
    
    % If presenting X dB FPL, this is 3-6 dB higher in dB SPL.
    % In order to make more comparable with other studies that calibrate in
    % dB SPL, can adjust the level accordingly:
    if strcmp(calType,'thev') 
        if strcmp(targetCalType,'fpl')
            adjust = -3; % if 55 dB fpl, this is more like 58 dB SPL. Therefore, subtract 3 to make the level 52 dB FPL ~= 55 dB SPL.
        else
            adjust = 0;
        end
    end
        
    % get probe settings
    validateStimLevels = [];
    [probeInputL,probeOutputL,probeInputR,probeOutputR,~] = getProbeSettings(probeL,probeR,validateStimLevels);
    
    % get most recent calibration files -----
    if isempty(obj)
        tempobj = struct;
        tempobj.map = [];
    else
        tempobj = obj;
    end    

    micCorrectionL = getMicCorrection(probeInputL,doMicCorrection,tempobj.map);
    micCorrectionR = getMicCorrection(probeInputR,doMicCorrection,tempobj.map);
    [C1L,C2L,calPath1L,calPath2L] = getOutputCal(calType,probeOutputL,tempobj.map);
    [C1R,C2R,calPath1R,calPath2R] = getOutputCal(calType,probeOutputR,tempobj.map);
    
    if isempty(obj) % decide whether collecting new data or running post-hoc analysis
        runPostAnalysis = 1;
        disp('----- Starting dpoae post-hoc analysis -----')
    else
        runPostAnalysis = 0;
    end
    if runPostAnalysis == 1 % running post-hoc analyses on previously recorded data
        [obj,nLevels,fs,dL,dR,iscS1L,iscS2L,iscS12L,iscS1R,iscS2R,iscS12R,postPath,postPathSave,testEar] = getSavedData(basePath,SUBJ,DATE,EXPRUN,probeInputL,probeInputR);
    else % actually collecting data

        % PERFORM IN-SITU CALIBRATION -----
        % Only perform in-situ calibration if using Thevenin source;
        % otherwise using long tube and constant voltage paradigm.
        if strcmp(calType,'thev')
            disp('----- Running in-situ calibration -----')
            inSituReps = 6;
            doIndividual = 1;
            doSimultaneous = 0;
            [iscS1L,iscS2L,iscS12L] = ARLas_runISC(obj,probeInputL,probeOutputL,calPath1L,calPath2L,calType,fmin,fmax,micCorrectionL,inSituReps,doIndividual,doSimultaneous);
            [iscS1R,iscS2R,iscS12R] = ARLas_runISC(obj,probeInputR,probeOutputR,calPath1R,calPath2R,calType,fmin,fmax,micCorrectionR,inSituReps,doIndividual,doSimultaneous);
            ISCL.(['Rec',num2str(1)]).(['S1']) = iscS1L;
            ISCL.(['Rec',num2str(1)]).(['S2']) = iscS2L;
            ISCR.(['Rec',num2str(1)]).(['S1']) = iscS1R;
            ISCR.(['Rec',num2str(1)]).(['S2']) = iscS2R;            
        else
            error('Unrecognized calibration type. Aborting DPOAE program.')
        end
        
% PRE-TEST IN-SITU CALIBRATION CHECK --------------------------------------
    if doPreCalCheck == 1
        try
            if ~isempty(probeOutputL)
                labelL = probeOutputL.label;
                dL = dir([obj.savingGrace,'*',labelL,'_inSituCal_*']);
            else
                dL = [];
            end
            if ~isempty(probeOutputR)
                labelR = probeOutputR.label;
                dR = dir([obj.savingGrace,'*',labelR,'_inSituCal_*']);
            else
                dR = [];
            end
            iscCheckHandleL = [];
            if ~isempty(probeOutputL)
                dummy = load([obj.savingGrace,dL(1).name]); % first calibration chirp (loudspeaker 1)
                iscCheckHandleL = inSituCheck(probeOutputL.label,probeOutputL.ch(1),dummy.header,dummy.data,obj.fs,1,iscCheckHandleL);
                dummy = load([obj.savingGrace,dL(2).name]); % first calibration chirp (loudspeaker 2)
                iscCheckHandleL = inSituCheck(probeOutputL.label,probeOutputL.ch(2),dummy.header,dummy.data,obj.fs,2,iscCheckHandleL);
            end
            iscCheckHandleR = [];
            if ~isempty(probeOutputR)
                dummy = load([obj.savingGrace,dR(1).name]); % first calibration chirp (loudspeaker 1)
                iscCheckHandleR = inSituCheck(probeOutputR.label,probeOutputR.ch(1),dummy.header,dummy.data,obj.fs,1,iscCheckHandleR);
                dummy = load([obj.savingGrace,dR(2).name]); % first calibration chirp (loudspeaker 2)
                iscCheckHandleR = inSituCheck(probeOutputR.label,probeOutputR.ch(2),dummy.header,dummy.data,obj.fs,2,iscCheckHandleR);
            end
        catch ME
            disp('Start of measurement in-situ calibration check failed!')
        end
        answer = questdlg('Peaks at low and high frequencies should be between 1 and 2 Pa. Do NOT close figures!',...
            'In-situ Calibration Check',...
	        'Continue','Abort','Continue');
        switch answer
            case 'Continue'
            case 'Abort'
                return
        end
    end
    % END PRE-TEST IN-SITU CALIBRATION CHECK ---------------------------------------



        % get stimuli ---------------------------------------------------------
        fs = obj.fs; % get the system sampling rate
        % moved down into loop so that fRatio can change with level

        % COLLECT THE DATA --------------------------------------------------------
        disp('----- Collecting data -----')
        nLevels = length(targetL1); % number of stim
    end
        
    for ii=1:nLevels % loop over number of stimulus levels
        if runPostAnalysis == 0 % if actually recording data
            [Y1,Y2,F1,F2,Fdp,Time,Y1c,Y2c,Ydp,Ydpc] = ARLas_dpoaeStim_fixedRatioContinuous_vRT(fs,fmin,fmax,sweepRate,fRatio); % get the raw stimuli

            % Apply the in-situ calibration to the stimuli ----------------
            if strcmp(calType,'thev')
                if ~isempty(iscS1L)
                    [rows,cols] = size(Y1);
                    [S1L,clicksImpulse,errors1] = ARLas_applyISC(iscS1L,targetL1(ii)+adjust,targetCalType,F1(:),Y1(:));
                    S1L = reshape(S1L,rows,cols);
                    [S2L,clicksImpulse,errors1] = ARLas_applyISC(iscS2L,targetL2(ii)+adjust,targetCalType,F2(:),Y2(:));
                    S2L = reshape(S2L,rows,cols);
                else
                    S1L = [];
                    S2L = [];
                end
                if ~isempty(iscS1R)
                    [rows,cols] = size(Y1);
                    [S1R,clicksImpulse,errors1] = ARLas_applyISC(iscS1R,targetL1(ii)+adjust,targetCalType,F1(:),Y1(:));
                    S1R = reshape(S1R,rows,cols);
                    [S2R,clicksImpulse,errors1] = ARLas_applyISC(iscS2R,targetL2(ii)+adjust,targetCalType,F2(:),Y2(:));
                    S2R = reshape(S2R,rows,cols);
                else
                    S1R = [];
                    S2R = [];
                end
            elseif strcmp(calType,'long')         
                if ~isempty(probeInputL)
                    [rows,cols] = size(Y1);
                    [S1L,scaling] = ARLas_applyISC_LTglides(C1L.fo_spl,C1L.phi,C1L.freq,targetL1(ii),F1(:),Y1(:));
                    S1L = reshape(S1L,rows,cols);
                    S1L = S1L * multiplierL;
                    [S2L,scaling] = ARLas_applyISC_LTglides(C2L.fo_spl,C2L.phi,C2L.freq,targetL2(ii),F2(:),Y2(:));
                    S2L = reshape(S2L,rows,cols);
                    S2L = S2L * multiplierL;
                else
                    S1L = [];
                    S2L = [];
                end
                if ~isempty(probeInputR)
                    [rows,cols] = size(Y1);
                    [S1R,scaling] = ARLas_applyISC_LTglides(C1R.fo_spl,C1R.phi,C1R.freq,targetL1(ii),F1(:),Y1(:));
                    S1R = reshape(S1R,rows,cols);
                    S1R = S1R * multiplierR;
                    [S2R,scaling] = ARLas_applyISC_LTglides(C2R.fo_spl,C2R.phi,C2R.freq,targetL2(ii),F2(:),Y2(:));
                    S2R = reshape(S2R,rows,cols);
                    S2R = S2R * multiplierR;
                else
                    S1R = [];
                    S2R = [];
                end
            end

            % tell ARLas what to record and what to play
            obj.clearPlayList % clear out whatever was played previously
            if ~isempty(probeInputL) % load new signals to play
                obj.setPlayList(S1L,probeOutputL.ch(1));
                obj.setPlayList(S2L,probeOutputL.ch(2));
            end
            if ~isempty(probeInputR)
                obj.setPlayList(S1R,probeOutputR.ch(1));
                obj.setPlayList(S2R,probeOutputR.ch(2));
            end

            obj.clearRecList % clear out whatever was used previously 
            if ~isempty(probeInputL)
                if ii==1
                    probeInputL.label = [probeInputL.label,'_dpoae'];
                end
                obj.setRecList(probeInputL.ch,probeInputL.label,probeInputL.micSens,probeInputL.gain);    
            end
            if ~isempty(probeInputR)
                if ii==1
                    probeInputR.label = [probeInputR.label,'_dpoae'];
                end
                obj.setRecList(probeInputR.ch,probeInputR.label,probeInputR.micSens,probeInputR.gain);    
            end

            nReps = nSweeps(ii) * size(Y1,2); % number of reps (for ARLas, since stim folded into a matrix)
            obj.setNReps(nReps); % number of times to play stimulus
            obj.setFilter(1); % note: this is a highpass filter with a 75 Hz cutoff frequency.
            
            obj.objPlayrec.userInfo = []; % clear out previous, non-structure array info
            obj.objPlayrec.userInfo.testEar = testEar; % which ear is ACTUALLY being tested, regardless of what the probe info says
            obj.objPlayrec.userInfo.fmin = fmin;
            obj.objPlayrec.userInfo.fmax = fmax;
            obj.objPlayrec.userInfo.fs = fs;
            obj.objPlayrec.userInfo.dpoae_version = V;
            obj.objPlayrec.userInfo.sweepRate = sweepRate;
            obj.objPlayrec.userInfo.fRatio = fRatio;
            obj.objPlayrec.userInfo.Y1 = Y1;
            obj.objPlayrec.userInfo.Y2 = Y2;
            obj.objPlayrec.userInfo.Y1c = Y1c;
            obj.objPlayrec.userInfo.Y2c = Y2c;
            obj.objPlayrec.userInfo.Ydp = Ydp;
            obj.objPlayrec.userInfo.Ydpc = Ydpc;
            obj.objPlayrec.userInfo.Time = Time;
            obj.objPlayrec.userInfo.F1 = F1;
            obj.objPlayrec.userInfo.F2 = F2;
            obj.objPlayrec.userInfo.Fdp = Fdp;
            obj.objPlayrec.userInfo.nSweeps = nSweeps(ii);
            obj.objPlayrec.userInfo.nReps = nReps;
            obj.objPlayrec.userInfo.targetL1 = targetL1(ii);
            obj.objPlayrec.userInfo.targetL2 = targetL2(ii);
            obj.objPlayrec.userInfo.calType = calType;
            obj.objPlayrec.userInfo.targetCalType = targetCalType;
            obj.objPlayrec.userInfo.C1L = C1L;
            obj.objPlayrec.userInfo.C2L = C2L;
            obj.objPlayrec.userInfo.C1R = C1R;
            obj.objPlayrec.userInfo.C2R = C2R;
            if strcmp(calType,'thev')
                obj.objPlayrec.userInfo.iscS1L = iscS1L;
                obj.objPlayrec.userInfo.iscS2L = iscS2L;
                obj.objPlayrec.userInfo.iscS1R = iscS1R;
                obj.objPlayrec.userInfo.iscS2R = iscS2R;
            else
                iscS1L = [];
                iscS2L = [];
                iscS12L = [];
                iscS1R = [];
                iscS2R = [];
                iscS12R = [];
            end
            obj.objPlayrec.userInfo.probeR = probeR;
            obj.objPlayrec.userInfo.probeL = probeL;
            obj.objPlayrec.userInfo.micCorrectionL = micCorrectionL; 
            obj.objPlayrec.userInfo.micCorrectionR = micCorrectionR; 
            

            disp(['----- Recording data: ',num2str(ii),' of ',num2str(nLevels),' -----'])
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
                %[pl,Pl,phi,other,wf] = ARLas_convertPL(Data,iscS1L);
            end
            if ~isempty(probeInputR)
                [headerR,DataR] = obj.retrieveData(['Ch',num2str(probeInputR.ch)]); % retrive recorded data
                if doMicCorrection == 1
                    DataR = applyMicCorrection(DataR,micCorrectionR);
                end
                %[pl,Pl,phi,other,wf] = ARLas_convertPL(Data,iscS1R);
            end
        end
        
        % Analyze data ----------------------------------------------------
        disp('----- Analyzing data -----')
        % analysis prep -----
        if ii==1 % on the first iteration, initialize the analysis structure
            if ~isempty(probeInputL)
                DPOAE_L = getDPstruct(obj,fs,fmin,fmax,sweepRate,nSweeps(ii),targetL1(ii),targetL2(ii),calType,targetCalType,iscS1L,iscS2L,iscS12L,C1L,C2L,'Left');
                if ~strcmp(testEar,'Both')
                    DPOAE_L.ear = testEar;
                end
            end
            if ~isempty(probeInputR)
                DPOAE_R = getDPstruct(obj,fs,fmin,fmax,sweepRate,nSweeps(ii),targetL1(ii),targetL2(ii),calType,targetCalType,iscS1R,iscS2R,iscS12R,C1R,C2R,'Right');
                if ~strcmp(testEar,'Both')
                    DPOAE_R.ear = testEar;
                end
            end
        end 
        
        if runPostAnalysis == 1 % running post-hoc analyses on previously recorded data
            if ~isempty(dL)
                dummy = load([postPath,dL(ii).name]);
                headerL = dummy.header;
                DataL = dummy.data;
                clear dummy
                micCorrectionL = headerL.userInfo.micCorrectionL;
                DataL = applyMicCorrection(DataL,micCorrectionL);
                DPOAE_L.postHocTimeStamp = DPOAE_L.timeStamp; % note that this is a post-hoc analysis
                DPOAE_L.timeStamp = obj.timeStamp; % make sure that original time stamp shows up on the figues
            end
            if ~isempty(dR)
                dummy = load([postPath,dR(ii).name]);
                headerR = dummy.header;
                DataR = dummy.data;
                clear dummy
                micCorrectionR = headerR.userInfo.micCorrectionR;
                DataR = applyMicCorrection(DataR,micCorrectionR);
                DPOAE_R.postHocTimeStamp = DPOAE_R.timeStamp; % note that this is a post-hoc analysis
                DPOAE_R.timeStamp = obj.timeStamp; % make sure that original time stamp shows up on the figues
            end
        end    
        
        % ANALYSIS BEGINS HERE --------------------------------------------
        if ~isempty(probeInputL)
            DPOAE_L = dpoaeAnalysis_fixedRatioSweep_vRT(headerL,DataL,DPOAE_L,[]);
            %DPOAE_L = dpoae_sweepRatioANALYSIS_vRT3(headerL,DataL,DPOAE_L,[]);
        end
        if ~isempty(probeInputR)
            DPOAE_R = dpoaeAnalysis_fixedRatioSweep_vRT(headerR,DataR,DPOAE_R,[]);
            %DPOAE_R = dpoae_sweepRatioANALYSIS_vRT3(headerR,DataR,DPOAE_R,[]);
        end
        
        disp('----- Plotting data -----') % -------------------------------
        try
            if ~isempty(probeInputL)
                [h1L,h2L] = dpoaePlotting_fixedRatioSweep_vRT(DPOAE_L);
                %h1L.Position = [591   154  469  211];
                %h2L.Position = [591  442  469  314];
                pause(0.01)
            end
        catch
        end
        try
            if ~isempty(probeInputR)
                [h1R,h2R] = dpoaePlotting_fixedRatioSweep_vRT(DPOAE_R);
                %h1R.Position = [1060  154  469  211];
                %h2R.Position = [1060  442  469  314];
                pause(0.01);
            end
        catch
        end
    end
        
% PERFORM POST-TEST IN-SITU CALIBRATION -----------------------------------
        if doPostCalCheck == 1
            if ~isempty(probeInputL)
                if contains(probeInputL.label,'dpoae')
                    dummy = probeInputL.label;
                    dummy = dummy(1:end-6);
                    probeInputL.label = dummy;
                end
            end
            if ~isempty(probeInputR)
                if contains(probeInputR.label,'dpoae')
                    dummy = probeInputR.label;
                    dummy = dummy(1:end-6);
                    probeInputR.label = dummy;
                end
            end
            disp('----- Running post-test in-situ calibration -----')
            inSituReps = 6;
            doIndividual = 1;
            doSimultaneous = 0;
            [iscS1L,iscS2L,iscS12L] = ARLas_runISC(obj,probeInputL,probeOutputL,calPath1L,calPath2L,calType,fmin,fmax,micCorrectionL,inSituReps,doIndividual,doSimultaneous);
            [iscS1R,iscS2R,iscS12R] = ARLas_runISC(obj,probeInputR,probeOutputR,calPath1R,calPath2R,calType,fmin,fmax,micCorrectionR,inSituReps,doIndividual,doSimultaneous);       
            ISCL.(['Rec',num2str(ii+1)]).(['S1']) = iscS1L;
            ISCL.(['Rec',num2str(ii+1)]).(['S2']) = iscS2L;
            ISCR.(['Rec',num2str(ii+1)]).(['S1']) = iscS1R;
            ISCR.(['Rec',num2str(ii+1)]).(['S2']) = iscS2R;            
                
            try
                if ~isempty(probeOutputL)
                    labelL = probeOutputL.label;
                    dL = dir([obj.savingGrace,'*',labelL,'_inSituCal_*']);
                else
                    dL = [];
                end
                if ~isempty(probeOutputR)
                    labelR = probeOutputR.label;
                    dR = dir([obj.savingGrace,'*',labelR,'_inSituCal_*']);
                else
                    dR = [];
                end
                if ~isempty(probeOutputL)
                    figure(iscCheckHandleL)
                    dummy = load([obj.savingGrace,dL(3).name]); % first calibration chirp (loudspeaker 1)
                    iscCheckHandleL = inSituCheck(probeOutputL.label,probeOutputL.ch(1),dummy.header,dummy.data,obj.fs,1,iscCheckHandleL);
                    dummy = load([obj.savingGrace,dL(4).name]); % first calibration chirp (loudspeaker 2)
                    iscCheckHandleL = inSituCheck(probeOutputL.label,probeOutputL.ch(2),dummy.header,dummy.data,obj.fs,2,iscCheckHandleL);
                end
                if ~isempty(probeOutputR)
                    figure(iscCheckHandleR)
                    dummy = load([obj.savingGrace,dR(3).name]); % first calibration chirp (loudspeaker 1)
                    iscCheckHandleR = inSituCheck(probeOutputR.label,probeOutputR.ch(1),dummy.header,dummy.data,obj.fs,1,iscCheckHandleR);
                    dummy = load([obj.savingGrace,dR(4).name]); % first calibration chirp (loudspeaker 2)
                    iscCheckHandleR = inSituCheck(probeOutputR.label,probeOutputR.ch(2),dummy.header,dummy.data,obj.fs,2,iscCheckHandleR);
                end
        
            catch ME
                disp('End of measurement in-situ calibration check failed!')
            end
            pause(0.02)
        end
        % END POST-TEST IN-SITU CALIBRATION CHECK ---------------------------------------



    % Save figures ----------
    disp('----- Saving figures -----')
    try
        if ~isempty(probeInputL)
            if runPostAnalysis == 1 % saving post-hoc analyses on previously recorded data
                savePath = postPathSave;
            else
                savePath = [obj.savingGrace,obj.experimentID,'_analysis',obj.sep];
                %savePath = obj.objPlayrec.savedFilesPath;
                if exist(savePath,'dir') == 0 % if path does not exist
                    success = mkdir(savePath);
                    if ~success
                        warning('Unable to create new Experiment Directory: data save path')
                    end
                end 
            end
            %figureFileName = [obj.subjectID,'_dpoaePrimariesL.fig'];
            figureFileName = [obj.subjectID,'_dpoaePrimaries_',DPOAE_L.ear,'.fig'];
            figureFileName = ARLas_saveName(savePath,figureFileName);
            savefig(h1L,[savePath,figureFileName])
            %figureFileName = [obj.subjectID,'_dpoaePrimariesL.bmp'];
            figureFileName = [obj.subjectID,'_dpoaePrimaries_',DPOAE_L.ear,'.bmp'];
            figureFileName = ARLas_saveName(savePath,figureFileName);    
            saveas(h1L,[savePath,figureFileName])

            %figureFileName = [obj.subjectID,'_dpoaeL.fig'];
            figureFileName = [obj.subjectID,'_dpoae_',DPOAE_L.ear,'.fig'];
            figureFileName = ARLas_saveName(savePath,figureFileName);
            savefig(h2L,[savePath,figureFileName])
            figureFileName = [obj.subjectID,'_dpoae_',DPOAE_L.ear,'.bmp'];
            figureFileName = ARLas_saveName(savePath,figureFileName);    
            saveas(h2L,[savePath,figureFileName])
        end
    catch
        disp('Warning: One or more figures for DPOAE_L not saved!')        
    end

    try
        if ~isempty(probeInputR)
            if runPostAnalysis == 1 % saving post-hoc analyses on previously recorded data
                savePath = postPathSave;
            else % saving newly acquired data
                savePath = [obj.savingGrace,obj.experimentID,'_analysis',obj.sep];
                if exist(savePath,'dir') == 0 % if path does not exist
                    success = mkdir(savePath);
                    if ~success
                        warning('Unable to create new Experiment Directory: data save path')
                    end
                end 
            end
            figureFileName = [obj.subjectID,'_dpoaePrimaries_',DPOAE_R.ear,'.fig'];
            figureFileName = ARLas_saveName(savePath,figureFileName);
            savefig(h1R,[savePath,figureFileName])
            figureFileName = [obj.subjectID,'_dpoaePrimaries_',DPOAE_R.ear,'.bmp'];
            figureFileName = ARLas_saveName(savePath,figureFileName);    
            saveas(h1R,[savePath,figureFileName])

            figureFileName = [obj.subjectID,'_dpoae_',DPOAE_R.ear,'.fig'];
            figureFileName = ARLas_saveName(savePath,figureFileName);
            savefig(h2R,[savePath,figureFileName])
            figureFileName = [obj.subjectID,'_dpoae_',DPOAE_R.ear,'.bmp'];
            figureFileName = ARLas_saveName(savePath,figureFileName);    
            saveas(h2R,[savePath,figureFileName])
        end
    catch
       disp('Warning: One or more figures for DPOAE_R not saved!')
    end
    
    % Save analyses ------------------------
    disp('----- Saving analysis files -----')
    try 
        if ~isempty(probeInputL)
            if runPostAnalysis == 1 % saving post-hoc analyses on previously recorded data
                savePath = postPathSave;
            else % saving newly acquired data
                savePath = [obj.savingGrace,obj.experimentID,'_analysis',obj.sep];
                if exist(savePath,'dir') == 0 % if path does not exist
                    success = mkdir(savePath);
                    if ~success
                        warning('Unable to create new Experiment Directory: data save path')
                    end
                end 
            end
            %saveFileName = [obj.subjectID,'_analyzedDPOAE_L.mat'];
            saveFileName = [obj.subjectID,'_analyzedDPOAE_',DPOAE_L.ear,'.mat'];
            saveFileName = ARLas_saveName(savePath,saveFileName);
            save([savePath,saveFileName],'DPOAE_L','ISCL')
        end
    catch
        disp('Warning: DPOAE_L Analysis not saved!')
    end
    
    try 
        if ~isempty(probeInputR)
            if runPostAnalysis == 1 % saving post-hoc analyses on previously recorded data
                savePath = postPathSave;
            else % saving newly acquired data
                savePath = [obj.savingGrace,obj.experimentID,'_analysis',obj.sep];
                if exist(savePath,'dir') == 0 % if path does not exist
                    success = mkdir(savePath);
                    if ~success
                        warning('Unable to create new Experiment Directory: data save path')
                    end
                end 
            end
            %saveFileName = [obj.subjectID,'_analyzedDPOAE_R.mat'];
            saveFileName = [obj.subjectID,'_analyzedDPOAE_',DPOAE_R.ear,'.mat'];
            saveFileName = ARLas_saveName(savePath,saveFileName);
            save([savePath,saveFileName],'DPOAE_R','ISCR')
        end
    catch
        disp('Warning: DPOAE_R Analysis not saved!')
    end

    % Save in situ cal figures --------------------------------------------
    disp('----- Saving figures -----')
    try
        if ~isempty(probeInputL)
            savePath = [obj.savingGrace,obj.experimentID,'_analysis',obj.sep];
            if exist(savePath,'dir') == 0 % if path does not exist
                success = mkdir(savePath);
                if ~success
                    warning('Unable to create new Experiment Directory: data save path')
                end
            end 
            figureFileName = [obj.subjectID,'_inSituCal_',DPOAE_L.ear,'.fig'];
            figureFileName = ARLas_saveName(savePath,figureFileName);
            savefig(iscCheckHandleL,[savePath,figureFileName])
            figureFileName = [obj.subjectID,'_inSituCal_',DPOAE_L.ear,'.bmp'];
            figureFileName = ARLas_saveName(savePath,figureFileName);    
            saveas(iscCheckHandleL,[savePath,figureFileName])
        end
    catch
        if ~isempty(obj)
            disp('Warning: inSitu Cal L not saved!')        
        end
    end

    try
        if ~isempty(probeInputR)
            savePath = [obj.savingGrace,obj.experimentID,'_analysis',obj.sep];
            if exist(savePath,'dir') == 0 % if path does not exist
                success = mkdir(savePath);
                if ~success
                    warning('Unable to create new Experiment Directory: data save path')
                end
            end 
            figureFileName = [obj.subjectID,'_inSituCal_',DPOAE_R.ear,'.fig'];
            figureFileName = ARLas_saveName(savePath,figureFileName);
            savefig(iscCheckHandleR,[savePath,figureFileName])
            figureFileName = [obj.subjectID,'_inSituCal_',DPOAE_R.ear,'.bmp'];
            figureFileName = ARLas_saveName(savePath,figureFileName);    
            saveas(iscCheckHandleR,[savePath,figureFileName])
        end
    catch
       if ~isempty(obj)
           disp('Warning: inSitu Cal R not saved!')
       end
    end

    disp('----- Finished with DPOAE Fixed Ratio Sweep experiment -----')
    disp(' ')
    disp(' ')
toc/60
    
end

% INTERNAL FUNCTIONS ------------------------------------------------------
function [targetL1,nSweeps,targetL2,sweepRate,probeR,probeL,calType,targetCalType] = checkInputs(targetL1,nSweeps,targetL2,sweepRate,probeR,probeL,calType,targetCalType)
    if length(nSweeps) ~= length(targetL1) 
        if length(nSweeps) == 1
            nSweeps = ones(size(targetL1))*nSweeps;
        else
            error('Variables nSweeps and targetL1 must be same length.')
        end
    end
   
    % verify user modifiable parameters -----
    targetL1 = targetL1(:);
    targetL2 = targetL2(:);
    sweepRate = sweepRate(:);
    if length(targetL1) ~= length(targetL2)
        error('Variables targetL1 and targetL2 must be same length.')
    end
    if length(targetL1) ~= length(nSweeps)
        error('Variables targetL1 and nSweeps must be same length.')
    end
    if strcmp(probeR,probeL)
        error('Variables probeR and probeL cannot be the same letter.')
    end
    if ~(strcmp(probeR,'A') || strcmp(probeR,'B')) && ~isempty(probeR)
        error('Variable probeR contains unrecognized value. Must be string A, B, or empty set.')
    end
    if ~(strcmp(probeL,'A') || strcmp(probeL,'B')) && ~isempty(probeL)
        error('Variable probeL contains unrecognized value. Must be string A, B, or empty set.')
    end
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
function [DPOAE] = getDPstruct(obj,fs,fmin,fmax,sweepRate,nSweeps,targetL1,targetL2,calType,targetCalType,iscS1,iscS2,iscS12,C1,C2,ear)
    DPOAE.subjID = obj.subjectID;
    DPOAE.ear = ear;
    DPOAE.timeStamp = cellstr(datetime('now'));
    DPOAE.fs = fs;
    DPOAE.f2min = fmin;
    DPOAE.f2max = fmax;
    DPOAE.sweepRate = sweepRate;
    DPOAE.nSweeps = nSweeps;
    DPOAE.targetL1 = targetL1;
    DPOAE.targetL2 = targetL2;
    DPOAE.calType = calType;
    DPOAE.targetCalType = targetCalType;
    DPOAE.iscS1 = iscS1;
    DPOAE.iscS2 = iscS2;
    DPOAE.iscS12 = iscS12;
    DPOAE.C1 = C1;
    DPOAE.C2 = C2;

    DPOAE.L1 = [];
    DPOAE.N1 = [];
    DPOAE.P1 = [];
    DPOAE.L2 = [];
    DPOAE.N2 = [];
    DPOAE.P2 = [];
    DPOAE.Ldp = [];
    DPOAE.Ndp = [];
    DPOAE.Pdp = [];
    DPOAE.f1 = [];
    DPOAE.f2 = [];
    DPOAE.fdp = [];
    
    DPOAE.postHocTimeStamp = [];
    DPOAE.Ldp_3oct = [];
    DPOAE.Ndp_3oct = [];
    DPOAE.fc = [];
    
end
function [obj,nLevels,fs,dL,dR,iscS1L,iscS2L,iscS12L,iscS1R,iscS2R,iscS12R,postPath,postPathSave,testEar] = getSavedData(basePath,SUBJ,DATE,EXPRUN,probeInputL,probeInputR)      
    os = filesep; %'\';
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
    try
        dL = dir([postPath,'Ch',num2str(probeInputL.ch),'_',probeInputL.label,'_dpoae_*.mat']);
    catch
        dL = [];
    end
    try
        dR = dir([postPath,'Ch',num2str(probeInputR.ch),'_',probeInputR.label,'_dpoae_*.mat']);
    catch
        dR = [];
    end
    nLevels = max([size(dL,1),size(dR,1)]); %size(dL,1); % will be used to set iterations of the loop
    fs = 96000; % sampling rate (Hz)

    % get the saved in-situ calibrations
    if ~isempty(dL)
        dummy = load([postPath,dL(1).name]);
        q = dummy.header;
        %C1L = q.userInfo.C2L;
        %C2L = q.userInfo.C2L;
        iscS1L = q.userInfo.iscS1L;
        iscS2L = q.userInfo.iscS2L;        
        iscS12L = [];
    else
        iscS1L = [];
        iscS2L = [];        
        iscS12L = [];
    end

    if ~isempty(dR)
        dummy = load([postPath,dR(1).name]);
        q = dummy.header;
        %C1R = q.userInfo.C2R;
        %C2R = q.userInfo.C2R;
        iscS1R = q.userInfo.iscS1R;
        iscS2R = q.userInfo.iscS2R;        
        iscS12R = [];
    else
        iscS1R = [];
        iscS2R = [];        
        iscS12R = [];
    end

    obj.subjectID = SUBJ;
    obj.timeStamp = q.timeStamp;

    if ~isempty(dR) & ~isempty(dL)
        testEar = 'Both';
    elseif ~isempty(dR)
        testEar = 'Right';
    elseif ~isempty(dL)
        testEar = 'Left';
    end

end
function [Y1,Y2,F1,F2,Fdp,Time,Y1c,Y2c,Ydp,Ydpc] = ARLas_dpoaeStim_fixedRatioContinuous_vRT(fs,fmin,fmax,sweepRate,fRatio)
    % Create two vectors of DPOAE stimuli, f1 and f2.
    % The vectors are then folded into matrices for efficiency presenting through ARLas
    octPerSec = sweepRate;

    foldLength = 0.1; % fold into a matrix with columns this length (sec)
    foldN = round(foldLength * fs); % number of samples in each column

    nOctaves = log2(fmax)-log2(fmin); % numbmer of octaves to sweep over
    len = nOctaves / octPerSec; % length of total stimulus (sec)
    N = round(fs*len); % number of samples in total stimulus

    oct = linspace(0,nOctaves,N)'; % octave spacing
    f2 = 1000 * 2.^(oct); % put into linear space.
    k = fmin / 1000; % k is a multiplier to transform to the desired frequency range
    f2 = f2 * k; % f2 primary frequencies

    rampLen1 = 0.003; %5/fmin; % need some ramp on and off to avoid frequency splatter
    rampN1 = round(fs * rampLen1);
    rampLen2 = 0.003; %5/fmax;
    rampN2 = round(fs * rampLen2);
    pad1 = ones(rampN1,1)*fmin;
    pad2 = ones(rampN2,1)*fmax;
    f2 = [pad1;f2;pad2];

    % the following line adjsut for optimal (human) ratios
    %[f1,fdp,fRatio] = getRatios(targetL1,targetL2,f2);
    
    % the following lise are for truly fixed fRatio
    %fRatio = 1.22; % f2/fs ratio; now input argument
    f1 = f2 / fRatio; % f1 primary frequencies
    fdp = 2*f1 - f2; % cubic distortion dpoae frequencies

    N = length(f1); % recalculate new number of samples

    phi1 = 0; % starting phase
    phi2 = 0;
    phiDP = 0;
    y1 = zeros(N,1); % initialize output
    y2 = zeros(N,1); % initialize output
    ydp = zeros(N,1); % initialize output
    y1c = zeros(N,1); % initialize output
    y2c = zeros(N,1); % initialize output
    ydpc = zeros(N,1); % initialize output

    deltaT = 1/fs; % sampling period
    for ii=1:N
        phaseChange1 = 2*pi* (deltaT * f1(ii));
        phi1 = phi1 + phaseChange1;

        phaseChange2 = 2*pi* (deltaT * f2(ii));
        phi2 = phi2 + phaseChange2;

        phaseChangeDP = 2*pi* (deltaT * fdp(ii));
        phiDP = phiDP + phaseChangeDP;

        y1(ii,1) = sin(phi1);
        y2(ii,1) = sin(phi2);
        ydp(ii,1) = sin(phiDP);

        y1c(ii,1) = cos(phi1);
        y2c(ii,1) = cos(phi2);    
        ydpc(ii,1) = cos(phiDP);

    end

    h = hann(rampN1*2);
    h = h(1:rampN1);
    y1(1:rampN1) = y1(1:rampN1) .* h;
    y2(1:rampN1) = y2(1:rampN1) .* h;
    ydp(1:rampN1) = ydp(1:rampN1) .* h;
    y1c(1:rampN1) = y1c(1:rampN1) .* h;
    y2c(1:rampN1) = y2c(1:rampN1) .* h;
    ydpc(1:rampN1) = ydpc(1:rampN1) .* h;
    h = hann(rampN2*2);
    h = h(1:rampN2);
    h = flipud(h);
    y1(end-rampN2+1:end) = y1(end-rampN2+1:end) .* h;
    y2(end-rampN2+1:end) = y2(end-rampN2+1:end) .* h;
    ydp(end-rampN2+1:end) = ydp(end-rampN2+1:end) .* h;
    y1c(end-rampN2+1:end) = y1c(end-rampN2+1:end) .* h;
    y2c(end-rampN2+1:end) = y2c(end-rampN2+1:end) .* h;
    ydpc(end-rampN2+1:end) = ydpc(end-rampN2+1:end) .* h;

    NN = length(y1);
    nFolds = ceil(NN / foldN);
    extra = mod(NN,foldN);
    if extra ~= 0
        pad = zeros(foldN - extra,1);
        f1 = [f1;pad];
        f2 = [f2;pad];
        fdp = [fdp;pad];
        y1 = [y1;pad];
        y2 = [y2;pad];
        ydp = [ydp;pad];
        y1c = [y1c;pad];
        y2c = [y2c;pad];
        ydpc = [ydpc;pad];
    end
    Y1 = reshape(y1,foldN,nFolds);
    Y2 = reshape(y2,foldN,nFolds);
    Ydp = reshape(ydp,foldN,nFolds);
    Y1c = reshape(y1c,foldN,nFolds);
    Y2c = reshape(y2c,foldN,nFolds);
    Ydpc = reshape(ydpc,foldN,nFolds);
    F1 = reshape(f1,foldN,nFolds);
    F2 = reshape(f2,foldN,nFolds);
    Fdp = reshape(fdp,foldN,nFolds);
    NN = length(y1);
    time = (0:1:NN-1)'/fs;
    Time = reshape(time,foldN,nFolds);
end
function [f1,fdp,fRatio] = getRatios(L1,L2,f2)
% Values taken from Samantha Ginter dissertaion with Sumit Dhar, 2021
    %if L1==65 && L2==55
    if L1<75 && L2<75
        f = [0.5   0.75  1     1.5   2     3     4     6     8     10    11.2  12.5  14    15    16    17    18    19];
        r = [1.26  1.26  1.24  1.24  1.24  1.20  1.20  1.18  1.16  1.16  1.16  1.16  1.16  1.16  1.16  1.16  1.14  1.14];
    elseif L1>=75 && L2>=75
        f = [0.5   0.75  1     1.5   2     3     4     6     8     10    11.2  12.5  14    15    16    17    18    19];
        r = [1.28  1.28  1.28  1.28  1.28  1.24  1.22  1.20  1.18  1.18  1.18  1.18  1.18  1.18  1.18  1.18  1.18  1.18];
    else
        f = [0.5   0.75  1     1.5   2     3     4     6     8     10    11.2  12.5  14    15    16    17    18    19];
        r = [1.22  1.22  1.22  1.22  1.22  1.22  1.22  1.22  1.22  1.22  1.22  1.22  1.22  1.22  1.22  1.22  1.22  1.22];
    end
    f = f * 1000;
    fRatio = interp1(f,r,f2,'pchip');
    f1 = f2 ./ fRatio;
    fdp = 2*f1 - f2;
end

% OLD CODE ----------------------------------------------------------------
            % % look to see if mic correction exists as a field. If recording
            % % made before 3/23/2022, will not exist. 
            % if isfield(headerL,'micCorrectionL')
            %     micCorrectionL = headerL.micCorrectionL;
            % else % backwards compatability hack
            %     micCorrectionL = 1;
            % end
            % if isfield(headerL,'micCorrectionR')
            %     micCorrectionR = headerL.micCorrectionR;
            % else % backwards compatability hack
            %     micCorrectionR = 1;
            % end
            % 
            % DataL = applyMicCorrection(DataL,micCorrectionL);
            % DataR = applyMicCorrection(DataR,micCorrectionR);
            % 
            % DPOAE_L.postHocTimeStamp = DPOAE_L.timeStamp; % note that this is a post-hoc analysis
            % DPOAE_L.timeStamp = obj.timeStamp; % make sure that original time stamp shows up on the figues
            % DPOAE_R.postHocTimeStamp = DPOAE_R.timeStamp; % note that this is a post-hoc analysis
            % DPOAE_R.timeStamp = obj.timeStamp; % make sure that original time stamp shows up on the figues
