function [] = LRF(varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LRF(varargin)
%
% Measure DPOAEs using fixed F2 and sweeping F1 (i.e., sweep the frequency
% ratio).
%
% Author: Shawn Goodman, PhD
% Auditory Research Lab, the University of Iowa
% Original Date: January 3, 2025
% Last Updated: January 3, 2025 -- ssg -- 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    obj = varargin{1}; % get the arlas object
    V = '03JAN2025'; % this is the current version number of this program
    tic
%------ USER MODIFIABLE PARAMETERS ----------------------------------------
%--------------------------------------------------------------------------

    %  NOTE: f2 lists the frequencies to be tested (the higher frequency primary). The
    %        CURRENTLY YOU CAN ONLY ANALYZE 1 FREQUENCY!!
    f2 = [1925]; % f1 frequency (Hz) 
    Rmin = 1.01; % minimum f2/f1 ratio (must be >= 1), try 1.01
    Rmax = 1.3; % maximum f2/f1 ratio (must be >= 1), try 1.3
    nSweeps = 24; % number of test repeates 24 gives a good noise floor in humans when L2 = 65 and L1 = 55 dB FPL
    L1 = 65; % f1 level (dB FPL)
    L2 = 55; % f2 level (dB FPL)

    % specify which ER10X probes to use    
    probeL = []; % probe that is in LEFT ear; 'A', 'B', or []
    probeR = 'B'; % probe that is in RIGHT ear; 'A', 'B', or []
                  % Note: Either one or both ears can be tested at the same time.
                  % If either probe is empty set, will not be tested. 
                  
    % additional parameters that user usually doesn't need to mess with
    doMicCorrection = 1; % turn on and off microphone correction 
    doPreCalCheck = 1; % turn on and off the pre-test cal check plotting
    doPostCalCheck = 0; % turn on and off the post-test cal check check and plotting


%------ END USER MODIFIABLE PARAMETERS ------------------------------------
%--------------------------------------------------------------------------

    % edge effects mean that you need to extend the range a bit (sweepRate * 0.5)
    Rmin = Rmin - .01;
    Rmax = Rmax + .01;

    targetL1 = L1 * ones(size(f2)); % target levels for fixed L1 primaries (dB FPL);
    targetL2 = L2 * ones(size(f2)); % target levels for fixed L1 primaries (dB FPL);

    % specify calibration to use to set stimulus levels
    calType = 'thev'; % use 'thev' for Thevenin source
    targetCalType = 'fpl'; % 'fpl', 'spl', 'ipl'.

    disp(' ')
    disp('----- Starting dpoae experiment: fixed f2, swept f1 -----')
    ISCL = struct;
    ISCR = struct;
    
    if doPostCalCheck==1 && doPreCalCheck==0
        error('If doPostCalCheck=1, doPreCalcheck must also =1.')
    end

    if ~(isempty(probeL) | isempty(probeR)) % if both ears are being tested
        testProbe = 'Both';
    elseif ~(isempty(probeL))
        testProbe = 'Left';
    elseif ~(isempty(probeR))
        testProbe = 'Right';
    end

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
    if doMicCorrection==1
        if ~isempty(probeInputL)
            if micCorrectionL == 1
                warning('MIC CORRECTION IS NOT BEING APPLIED!')
            end
        end
        if  ~isempty(probeInputR)
            if micCorrectionR == 1
                warning('MIC CORRECTION IS NOT BEING APPLIED!')
            end
        end
    end
    [C1L,C2L,calPath1L,calPath2L] = getOutputCal(calType,probeOutputL,tempobj.map);
    [C1R,C2R,calPath1R,calPath2R] = getOutputCal(calType,probeOutputR,tempobj.map);
    

    % PERFORM IN-SITU CALIBRATION -----
    % Only perform in-situ calibration if using Thevenin source;
    % otherwise using long tube and constant voltage paradigm.
    if strcmp(calType,'thev')
        disp('----- Running in-situ calibration -----')
        inSituReps = 6;
        doIndividual = 1;
        doSimultaneous = 0;
        fmax = max(f2);
        fmin = min(f2./Rmax);
        [iscS1L,iscS2L,iscS12L] = ARLas_runISC(obj,probeInputL,probeOutputL,calPath1L,calPath2L,calType,fmin,fmax,micCorrectionL,inSituReps,doIndividual,doSimultaneous);
        [iscS1R,iscS2R,iscS12R] = ARLas_runISC(obj,probeInputR,probeOutputR,calPath1R,calPath2R,calType,fmin,fmax,micCorrectionR,inSituReps,doIndividual,doSimultaneous);
        ISCL.(['Rec',num2str(1)]).(['S1']) = iscS1L;
        ISCL.(['Rec',num2str(1)]).(['S2']) = iscS2L;
        ISCR.(['Rec',num2str(1)]).(['S1']) = iscS1R;
        ISCR.(['Rec',num2str(1)]).(['S2']) = iscS2R;            
    else
        error('Unrecognized calibration type. Aborting DPOAE program.')
    end

% PRE-TEST IN-SITU CALIBRATION CHECK -------------------------------------------
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
    % END PRE-TEST IN-SITU CALIBRATION CHECK ------------------------------

    % get stimuli ---------------------------------------------------------
    fs = obj.fs; % get the system sampling rate
    % moved down into loop so that fRatio can change with level

% COLLECT THE DATA --------------------------------------------------------
    disp('----- Collecting data -----')
    nLevels = length(targetL1); % number of stim
    
    % Apply the in-situ calibration to the stimuli ----------------
    nFreqs = length(f2);
    for kk=1:nFreqs

        fs = obj.fs;
        sweepLength = 8; % sweep length in seconds
        sweepN = round(fs*sweepLength);
        %[Y1,Y2,F1,F2,Fdp,Time,Y1c,Y2c,Ydp,Ydpc,fmin,fmax] = ARLas_dpoaeStim_sweepRatio(fs,f2(1),Rmin,Rmax,sweepN); % get the raw stimuli
        [y1,y2,F1,F2,Fdp,y1c,y2c,ydp,ydpc,fmin,fmax,fRatios] = ARLas_dpoaeStim_sweepRatio(fs,f2(1),Rmin,Rmax,sweepN); % get the raw stimuli
    
        sweepRate = (max(F1(:)) - min(F1(:))) / sweepLength; % f1 change in frequency per second. This will vary with f2 frequency 
        testLength = sweepLength * nSweeps; % total test length (sec)
        disp(['Test Length = ',num2str(testLength/60),' min.'])

        %f2(kk,1) = F2(1);
        %fdp(kk,1) = Fdp(1);
        if ~isempty(iscS1L)
            %[rows,cols] = size(Y1);
            [~,scaling,~] = ARLas_applyISC(iscS1L,targetL1(kk)+adjust,targetCalType,F1,y1);
            y1 = y1 .* abs(scaling);
            %S1L = Y1 .* abs(scaling);
            %Y2 = Y2 ./ max(Y2(:));
            [~,scaling,~] = ARLas_applyISC(iscS2L,targetL2(kk)+adjust,targetCalType,F2,y2);
            y2 = y2 .* abs(scaling);
            %S2L = Y2 .* abs(scaling);

            [F1,F2,Fdp,Y1,Y2,Ydp,Y1c,Y2c,Ydpc,foldLength,foldN] = foldme(F1,F2,Fdp,y1,y2,ydp,y1c,y2c,ydpc,fs);
            S1L = Y1;
            S2L = Y2;
            
            if kk==1
                S1Lfull = S1L;
                S2Lfull = S2L;
            else
                S1Lfull = S1L + S1Lfull;
                S2Lfull = S2L + S2Lfull;
            end                
        else
            S1L = [];
            S2L = [];
            S1Lfull = [];
            S2Lfull = [];
        end
        if ~isempty(iscS1R)
            [~,scaling,~] = ARLas_applyISC(iscS1R,targetL1(kk)+adjust,targetCalType,F1,y1);
            y1 = y1 .* abs(scaling); % multiply by magnitude only (scaling is complex; includes phase)
            [~,scaling,~] = ARLas_applyISC(iscS2R,targetL2(kk)+adjust,targetCalType,F2,y2);
            y2 = y2 .* abs(scaling);
            %S2R = Y2 .* abs(scaling);
            
            [F1,F2,Fdp,Y1,Y2,Ydp,Y1c,Y2c,Ydpc,foldLength,foldN] = foldme(F1,F2,Fdp,y1,y2,ydp,y1c,y2c,ydpc,fs);
            S1R = Y1;
            S2R = Y2;

            if kk==1
                S1Rfull = S1R;
                S2Rfull = S2R;
            else
                S1Rfull = S1R + S1Rfull;
                S2Rfull = S2R + S2Rfull;
            end
        else
            S1R = [];
            S2R = [];
            S1Rfull = [];
            S2Rfull = [];            
        end
    end
    S1L = S1Lfull;
    S2L = S2Lfull;
    S1R = S1Rfull;
    S2R = S2Rfull;
    %targetL2 = 20*log10(L2/.00002);
    
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
        probeInputL.label = [probeInputL.label,'_dpoae'];
        obj.setRecList(probeInputL.ch,probeInputL.label,probeInputL.micSens,probeInputL.gain);    
    end
    if ~isempty(probeInputR)
        probeInputR.label = [probeInputR.label,'_dpoae'];
        obj.setRecList(probeInputR.ch,probeInputR.label,probeInputR.micSens,probeInputR.gain);    
    end

    nSweeps = nSweeps(1); % all frequencies have the same number of sweeps
    nReps = nSweeps * size(Y1,2); % number of reps (for ARLas, since stim folded into a matrix)
    
    obj.setNReps(nReps); % number of times to play stimulus
    obj.setFilter(1); % note: this is a highpass filter with a 75 Hz cutoff frequency.
    
    % Populate userInfo fields --------------------------------------------
    obj.objPlayrec.userInfo = []; % clear out previous, non-structure array info
    obj.objPlayrec.userInfo.dpoae_version = V; % this version of the test
    obj.objPlayrec.userInfo.fs = fs; % sampling rate (Hz)
    obj.objPlayrec.userInfo.testProbe = testProbe; % which probe is being used (left,right, or both)
    obj.objPlayrec.userInfo.testEar = testEar; % which ear is ACTUALLY being tested (has probe insertion), regardless of which probe is being used
    

    obj.objPlayrec.userInfo.sweepLength = sweepLength; % length of each sweep (s)
    obj.objPlayrec.userInfo.sweepLength = sweepN; % number of samples in each sweep (samples)
    obj.objPlayrec.userInfo.sweepRate = sweepRate; % dB / sec
    obj.objPlayrec.userInfo.nSweeps = nSweeps; % number of sweeps to collect
    obj.objPlayrec.userInfo.nReps = nReps; % number of reps (for ARLas, since stim folded into a matrix)
    obj.objPlayrec.userInfo.fmin = fmin; % min frequency being tested
    obj.objPlayrec.userInfo.fmax = fmax; % max frequency being tested

    obj.objPlayrec.userInfo.targetL1 = targetL1; % target levels
    obj.objPlayrec.userInfo.targetL2 = targetL2; % target levels
    obj.objPlayrec.userInfo.L2min = Rmin; % 0 min L2 (dB FPL) 
    obj.objPlayrec.userInfo.L2max = Rmax; % 70 max L2 (dB FPL)
    obj.objPlayrec.userInfo.fRatios = fRatios; % f2 / f1 ratios
    obj.objPlayrec.userInfo.L1 = L1; % L1 vector of levels across sweep
    obj.objPlayrec.userInfo.L2 = L2; % L2 vector of levels across sweep
    obj.objPlayrec.userInfo.F1 = F1; % f1 vector of frequencies across sweep (these are constant in this program)
    obj.objPlayrec.userInfo.F2 = F2; % f2 vector of frequencies across sweep (these are constant in this program)
    obj.objPlayrec.userInfo.Fdp = Fdp; % fdp vector of frequencies across sweep (these are constant in this program)

    obj.objPlayrec.userInfo.Y1 = Y1;
    obj.objPlayrec.userInfo.Y2 = Y2;
    obj.objPlayrec.userInfo.Ydp = Ydp;
    obj.objPlayrec.userInfo.Y1c = Y1c;
    obj.objPlayrec.userInfo.Y2c = Y2c;
    obj.objPlayrec.userInfo.Ydpc = Ydpc;
    obj.objPlayrec.userInfo.foldLength = foldLength;
    obj.objPlayrec.userInfo.foldN = foldN;

    obj.objPlayrec.userInfo.calType = calType;
    obj.objPlayrec.userInfo.targetCalType = targetCalType;
    obj.objPlayrec.userInfo.C1L = C1L;
    obj.objPlayrec.userInfo.C2L = C2L;
    obj.objPlayrec.userInfo.iscS1L = iscS1L;
    obj.objPlayrec.userInfo.iscS2L = iscS2L;
    obj.objPlayrec.userInfo.probeL = probeL;
    obj.objPlayrec.userInfo.micCorrectionL = micCorrectionL; 
    obj.objPlayrec.userInfo.C1R = C1R;
    obj.objPlayrec.userInfo.C2R = C2R;
    obj.objPlayrec.userInfo.iscS1R = iscS1R;
    obj.objPlayrec.userInfo.iscS2R = iscS2R;
    obj.objPlayrec.userInfo.probeR = probeR;
    obj.objPlayrec.userInfo.micCorrectionR = micCorrectionR; 


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


% POST-TEST IN-SITU CALIBRATION CHECK -------------------------------------
    if doPostCalCheck == 1

        % PERFORM POST-TEST IN-SITU CALIBRATION -----
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
        ii=1;
        ISCL.(['Rec',num2str(ii+1)]).(['S1']) = iscS1L;
        ISCL.(['Rec',num2str(ii+1)]).(['S2']) = iscS2L;
        ISCR.(['Rec',num2str(ii+1)]).(['S1']) = iscS1R;
        ISCR.(['Rec',num2str(ii+1)]).(['S2']) = iscS2R;            
        toc/60
    
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
        
    % Analyze data ----------------------------------------------------
    disp('----- Analyzing data -----')
    if ~isempty(probeInputL)
        DPOAE_L = dpoae_sweepRatioANALYSIS_vRT3(headerL,DataL,obj);
    end
    if ~isempty(probeInputR)
        DPOAE_R = dpoae_sweepRatioANALYSIS_vRT3(headerR,DataR,obj);
    end

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
function [y1,y2,F1,F2,Fdp,y1c,y2c,ydp,ydpc,fmin,fmax,fRatios] = ARLas_dpoaeStim_sweepRatio(fs,f2,Rmin,Rmax,sweepN)
    % foldLength = 0.1; % fold into a matrix with columns this length (sec)
    % foldN = round(foldLength * fs); % number of samples in each column
    R = linspace(Rmin,Rmax,sweepN)'; % dB spacing, linear sweep

    % add padding for ramps on and off
    rampLen1 = 0.003; % 3 ms ramp on and off
    rampN1 = round(fs * rampLen1);
    rampLen2 = 0.003; 
    rampN2 = round(fs * rampLen2);
    pad1 = ones(rampN1,1)*min(R);
    pad2 = ones(rampN2,1)*max(R);
    fRatios = [pad1;R;pad2];

    F2 = ones(size(fRatios))*f2; % f2 is fixed
    F1 = F2 ./ fRatios; % get f1 vector of constant frequencies
    Fdp = 2*F1-F2;
   
    fmin = min(F1(:)); % minimum frequency being tested
    fmax = max(F2(:)); % maximum frequency being tested

    N = length(F2); % recalculate new number of samples (after padding added)
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
    % calculate stimuli
    for ii=1:N
        phaseChange1 = 2*pi* (deltaT * F1(ii));
        phi1 = phi1 + phaseChange1;

        phaseChange2 = 2*pi* (deltaT * F2(ii));
        phi2 = phi2 + phaseChange2;

        phaseChangeDP = 2*pi* (deltaT * Fdp(ii));
        phiDP = phiDP + phaseChangeDP;

        y1(ii,1) = sin(phi1);
        y2(ii,1) = sin(phi2);
        ydp(ii,1) = sin(phiDP);

        y1c(ii,1) = cos(phi1);
        y2c(ii,1) = cos(phi2);    
        ydpc(ii,1) = cos(phiDP);
    end
    % ramp on and off
    h = hann(rampN1*2);
    h = h(1:rampN1);
    y1(1:rampN1) = y1(1:rampN1) .* h; % stimulus in sine phase
    y2(1:rampN1) = y2(1:rampN1) .* h;
    ydp(1:rampN1) = ydp(1:rampN1) .* h;
    y1c(1:rampN1) = y1c(1:rampN1) .* h; % stimulus in cosine phase
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

    % % fold up the stimuli for efficient presentation
    % NN = length(y1);
    % nFolds = ceil(NN / foldN);
    % extra = mod(NN,foldN);
    % if extra ~= 0
    %     pad = zeros(foldN - extra,1);
    %     F1 = [F1;pad];
    %     F2 = [F2;pad];
    %     Fdp = [Fdp;pad];
    %     y1 = [y1;pad];
    %     y2 = [y2;pad];
    %     ydp = [ydp;pad];
    %     y1c = [y1c;pad];
    %     y2c = [y2c;pad];
    %     ydpc = [ydpc;pad];
    % end
    % Y1 = reshape(y1,foldN,nFolds);
    % Y2 = reshape(y2,foldN,nFolds);
    % Ydp = reshape(ydp,foldN,nFolds);
    % Y1c = reshape(y1c,foldN,nFolds);
    % Y2c = reshape(y2c,foldN,nFolds);
    % Ydpc = reshape(ydpc,foldN,nFolds);
    % F1 = reshape(F1,foldN,nFolds);
    % F2 = reshape(F2,foldN,nFolds);
    % Fdp = reshape(Fdp,foldN,nFolds);
    % NN = length(y1);
    % time = (0:1:NN-1)'/fs;
    % Time = reshape(time,foldN,nFolds);
end
function [F1,F2,Fdp,Y1,Y2,Ydp,Y1c,Y2c,Ydpc,foldLength,foldN] = foldme(F1,F2,Fdp,y1,y2,ydp,y1c,y2c,ydpc,fs)
    % fold up the stimuli for efficient presentation
    foldLength = 0.1; % fold into a matrix with columns this length (sec)
    foldN = round(foldLength * fs); % number of samples in each column    
    NN = length(y1);
    nFolds = ceil(NN / foldN);
    extra = mod(NN,foldN);
    if extra ~= 0
        pad = zeros(foldN - extra,1);
        F1 = [F1;pad];
        F2 = [F2;pad];
        Fdp = [Fdp;pad];
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
    F1 = reshape(F1,foldN,nFolds);
    F2 = reshape(F2,foldN,nFolds);
    Fdp = reshape(Fdp,foldN,nFolds);
    %NN = length(y1);
    %time = (0:1:NN-1)'/fs;
    %Time = reshape(time,foldN,nFolds);

end
function [fRatios] = getRatios(f2)
% NOTE: this code from work with Lori Dreisbach. Not for direct use here,
% but included for examples of reasonable ratios to use. This code was for
% use with newborns and may not work well for adults.
    n =5;
    if f2 == 2
        fRatios = linspace(1.19,1.27,n); %[1.1, 1.15, 1.2, 1.25, 1.3];
    elseif f2 == 4
        fRatios = linspace(1.16,1.24,n); %[1.1, 1.15, 1.2, 1.25, 1.3];
    elseif f2 == 6
        fRatios = linspace(1.23,1.15,n); %[1.1, 1.15, 1.2, 1.25, 1.3];
    elseif f2 == 8
        fRatios = linspace(1.14,1.22,n); %[1.1, 1.15, 1.2, 1.25, 1.3];
    elseif f2 == 10
        fRatios = linspace(1.13,1.21,n); %[1.05, 1.1, 1.15, 1.2, 1.25];
    elseif f2 == 12
        fRatios = linspace(1.15,1.23,n); %[1.1, 1.15, 1.2, 1.25, 1.3];
    elseif f2 == 14
        fRatios = linspace(1.13,1.21,n); %[1.1, 1.15, 1.2, 1.25, 1.3];
    elseif f2 == 16
        fRatios = linspace(1.14,1.22,n); %[1.1, 1.15, 1.2, 1.25, 1.3];
    else
        error('Unrecognized f2 value. Must be 2, 4, 6, 8, 10, 12, 14, or 16')
    end
end

% OLD CODE ----------------------------------------------------------------
