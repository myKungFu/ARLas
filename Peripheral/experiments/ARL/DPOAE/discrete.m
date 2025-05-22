function [] = discrete(varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% discrete(varargin)
%
% Measure DPOAEs using a fixed ratio and discrete levels and amplitudes.
% This version uses Thevenin source calibration and FPL.
%
% This version can be called from ARLas for new data collection, or it can
% be called stand-alone from the Matlab command prompt for post-hoc analysis, 
% like this: >> dpoae_fixedRatioDiscrete_vRT2([]);
% In that case, the program will use the information at the bottom of 
% USER MODIFIABLE PARAMETERS in order to do the analysis.
%
% Author: Shawn Goodman, PhD
% Auditory Research Lab, the University of Iowa
% Date: Original Date: July 19, 2024
% Last Updated: July 22, 2024 -- ssg -- Created a discrete experiment for JL.
% Last Updated: July 26, 2024 -- ssg -- fixed post cal check and saving
%                                       issues for figs and analyzed data.
% Last Updated: September 19, 2024 -- ssg -- added switch to do scissors or
%                                       fixed L1 paradigm
% Last Updated: December 5, 2024 -- ssg -- name change (_vRT2).
%               This version correctly handles multiple f2 and L2.
%               This version uses updated analysis/plotting codes.
%               This version adds option for repeating experiment k times.
% Last Updated: December 26, 2024 -- ssg -- new version number (vRT3).
%               Calls updated analysis and plotting functions. Removed
%               option for experimenting k times.
% Last Updated: January 3, 2025 -- ssg -- name changed from dpoae_fixedRatioDiscrete_vRT3
%               to dp_discrete_vRT3.m.
% Last Updated: Jan 28, 2025 -- ssg --Trying to get new IHS system to work
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    obj = varargin{1}; % get the arlas object
    V = '03JAN2025'; % this is the current version number of this program

%--------------------------------------------------------------------------
%------ USER MODIFIABLE PARAMETERS ----------------------------------------
    f2 = [1000,2000,4000,8000]; % frequencies to test (Hz)
    targetL2 = [55]; % L2 in dB FPL
    nSweeps = 12; %120; % Number of sweeps recorded at each level. 12 of these gives the same noise floor as 24 continuous sweeps

% additional parameters that user usually doesn't need to mess with
    doMicCorrection = 0; % turn on and off microphone correction
    doPreCalCheck = 0; % turn on and off the pre-test cal check
    doPostCalCheck = 0; % turn on and off the post-test cal check

% specify which ER10X probes to use    
    probeL = 'C'; %'A'; % probe that is in LEFT ear; 'A', 'B', or []
    probeR = []; %'B'; % probe that is in RIGHT ear; 'A', 'B', or []
                  % Note: Either one or both ears can be tested at the same time.
                  % If either probe is empty set, will not be tested. 
    
%------ END USER MODIFIABLE PARAMETERS ------------------------------------
%--------------------------------------------------------------------------

% GENERALLY, LEAVE THE FOLLOWING ALONE!
    fRatio = 1.22; % ratio of f2/f1. Ratio will remain fixed at this value
    stimParadigm = 'fixed'; % use 'fixed' or scissors'
    f1 = round(f2 ./ fRatio);
    fdp = 2*f1-f2;
    nLevels = length(targetL2);
    nFreqs = length(f2);
    fs = obj.fs; % get the system sampling rate

    if strcmp(stimParadigm,'fixed')
        % Fixed L1 paradigm
        targetL1 = ones(size(targetL2)) * 65; % target levels for primaries (dB FPL); 
    elseif strcmp(stimParadigm,'scissors')
        % Scissors paradigm
        targetL1 = 0.4*targetL2 + 39; % kummer scissors
        % targetL1 = 0.45.*targetL2 + 44; % neely scissors
    end

    if ~(isempty(probeL) | isempty(probeR)) % if both ears are being tested
        testProbe = 'Both';
    elseif ~(isempty(probeL))
        testProbe = 'Left';
    elseif ~(isempty(probeR))
        testProbe = 'Right';
    end

    disp(' ')
    disp('----- Starting DPOAE Discrete experiment -----')
    tic

    if doPostCalCheck==1 && doPreCalCheck==0
        error('If doPostCalCheck=1, doPreCalcheck must also =1.')
    end

    % specify calibration to use to set stimulus levels
    calType = 'thev'; % use 'thev' for Thevenin source
    targetCalType = 'fpl'; % 'fpl', 'spl', 'ipl'.
    ISCL = struct;
    ISCR = struct;

    % Create matrix of target conditions
    targetL1 = targetL1(:);
    targetL2 = targetL2(:);
    f2 = f2(:);
    f1 = f1(:);
    fdp = fdp(:);
    F1t = [];
    F2t = [];
    Fdpt = [];
    L1t = [];
    L2t = [];
    for jj=1:nFreqs
        dummy1 = ones(size(targetL1)) * f1(jj);
        dummy2 = ones(size(targetL1)) * f2(jj);
        dummy3 = ones(size(targetL1)) * fdp(jj);
        F1t = [F1t,dummy1];
        F2t = [F2t,dummy2];
        Fdpt = [Fdpt,dummy3];
        L1t = [L1t,targetL1];
        L2t = [L2t,targetL2];
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
    if strcmp(calType,'thev')
        disp('----- Running in-situ calibration -----')
        inSituReps = 6;
        doIndividual = 1;
        doSimultaneous = 0;
fmin = 100;
fmax = 8000; %18000;

        % get the new labels each time
        [probeInputL,probeOutputL,probeInputR,probeOutputR,~] = getProbeSettings(probeL,probeR,0);

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
        if kk==1
            answer = questdlg('Peaks at low and high frequencies should be between 1 and 2 Pa. Do NOT close figures!',...
                'In-situ Calibration Check',...
                'Continue','Abort','Continue');
            switch answer
                case 'Continue'
                case 'Abort'
                    return
            end
        end
    end
    % END PRE-TEST IN-SITU CALIBRATION CHECK ---------------------------------------

% COLLECT THE DATA --------------------------------------------------------
    disp('----- Collecting data -----')
    %nLevels = length(targetL1); % number of stim

    for JJ=1:nFreqs % loop over frequencies
        for II=1:nLevels % loop over stimulus levels
            disp(['Testing level ',num2str(II),' of ',num2str(nLevels),' at frequency ',num2str(JJ),' of ',num2str(nFreqs)])
            
            [Y1,Y2,F1,F2,Fdp,Time,Y1c,Y2c,Ydp,Ydpc] = ARLas_dpoaeStim_fixedRatioDiscrete_vRT3(fs,F1t(II,JJ),F2t(II,JJ),Fdpt(II,JJ)); % get the raw stimuli
            % Apply the in-situ calibration to the stimuli ----------------
            if strcmp(calType,'thev')
                if ~isempty(iscS1L)
                    [rows,cols] = size(Y1);
                    [S1L,clicksImpulse,errors1] = ARLas_applyISC(iscS1L,L1t(II,JJ)+adjust,targetCalType,F1(:),Y1(:));
                    S1L = reshape(S1L,rows,cols);
                    [S2L,clicksImpulse,errors1] = ARLas_applyISC(iscS2L,L2t(II,JJ)+adjust,targetCalType,F2(:),Y2(:));
                    S2L = reshape(S2L,rows,cols);
                else
                    S1L = [];
                    S2L = [];
                end
                if ~isempty(iscS1R)
                    [rows,cols] = size(Y1);
                    [S1R,clicksImpulse,errors1] = ARLas_applyISC(iscS1R,L1t(II,JJ)+adjust,targetCalType,F1(:),Y1(:));
                    S1R = reshape(S1R,rows,cols);
                    [S2R,clicksImpulse,errors1] = ARLas_applyISC(iscS2R,L2t(II,JJ)+adjust,targetCalType,F2(:),Y2(:));
                    S2R = reshape(S2R,rows,cols);
                else
                    S1R = [];
                    S2R = [];
                end
            else
                error('Thevenin source is the only calibration type that works with this experiment.')
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
                if II==1 & JJ==1
                    probeInputL.label = [probeInputL.label,'_fixedDP'];
                end
                obj.setRecList(probeInputL.ch,probeInputL.label,probeInputL.micSens,probeInputL.gain);    
            end
            if ~isempty(probeInputR)
                if II==1  & JJ==1
                    probeInputR.label = [probeInputR.label,'_fixedDP'];
                end
                obj.setRecList(probeInputR.ch,probeInputR.label,probeInputR.micSens,probeInputR.gain);    
            end

            nReps = nSweeps * size(Y1,2); % number of reps (for ARLas, since stim folded into a matrix)
            obj.setNReps(nReps); % number of times to play stimulus
            obj.setFilter(1); % note: this is a highpass filter with a 75 Hz cutoff frequency.
            
            obj.objPlayrec.userInfo = []; % clear out previous, non-structure array info
            obj.objPlayrec.userInfo.testEar = testEar; % which ear is ACTUALLY being tested, regardless of what the probe info says
            obj.objPlayrec.userInfo.testProbe = testProbe; % which probe is being used, regardless of ear in which it was placed
            obj.objPlayrec.userInfo.fmin = fmin;
            obj.objPlayrec.userInfo.fmax = fmax;
            obj.objPlayrec.userInfo.fs = fs;
            obj.objPlayrec.userInfo.dpoae_version = V;
            obj.objPlayrec.userInfo.L1t = L1t;
            obj.objPlayrec.userInfo.L2t = L2t;
            obj.objPlayrec.userInfo.F1t = F1t;
            obj.objPlayrec.userInfo.F2t = F2t;
            obj.objPlayrec.userInfo.Fdpt = Fdpt;
            obj.objPlayrec.userInfo.nFreqs = nFreqs;
            obj.objPlayrec.userInfo.nLevels = nLevels;
            obj.objPlayrec.userInfo.II = II; % the current level
            obj.objPlayrec.userInfo.JJ = JJ; % the current frequency
            obj.objPlayrec.userInfo.fRatio = fRatio;
            obj.objPlayrec.userInfo.nSweeps = nSweeps;
            obj.objPlayrec.userInfo.nReps = nReps;
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
        
% Analyze data ------------------------------------------------------------
            disp('----- Analyzing data -----')
            if II==1 & JJ==1 % on the first iteration, initialize the analysis structure
                if ~isempty(probeInputL)
                    if ~strcmp(testEar,'Both')
                        DPOAE_L.ear = testEar;
                    end
                    DPOAE_L = [];
                end
                if ~isempty(probeInputR)
                    if ~strcmp(testEar,'Both')
                        DPOAE_R.ear = testEar;
                    end
                    DPOAE_R = [];
                end
                h1 = [];
                h2 = [];
                h3 = [];
            end 
                    
            % ANALYSIS BEGINS HERE --------------------------------------------
            if ~isempty(probeInputL)
                if ~exist('h1','var')
                    h1 = [];
                    h2 = [];
                    h3 = [];
                end
                [DPOAE_L,h1,h2,h3] = dpoae_fixedRatioDiscreteANALYSIS_vRT3(headerL,DataL,DPOAE_L,h1,h2,h3,II,JJ);
            end
            if ~isempty(probeInputR)
                if ~exist('h1','var')
                    h1 = [];
                    h2 = [];
                    h3 = [];
                end
                [DPOAE_R,h1,h2,h3] = dpoae_fixedRatioDiscreteANALYSIS_vRT3(headerR,DataR,DPOAE_R,h1,h2,h3,II,JJ);
            end
        end
    end

% PERFORM POST-TEST IN-SITU CALIBRATION -----------------------------------
    if doPostCalCheck == 1  & runPostAnalysis == 0
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
    % END POST-TEST IN-SITU CALIBRATION CHECK -------------------------
    

% Save in situ cal figures ------------------------------------------------
    disp('----- Saving in-situ cal figures -----')
    try
        if ~isempty(probeInputL) & doPreCalCheck==1
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
        disp('Warning: inSitu Cal L not saved!')        
    end

    try
        if ~isempty(probeInputR)  & doPreCalCheck==1
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
       disp('Warning: inSitu Cal R not saved!')
    end

% %    disp(['Finished with repeat ',num2str(kk),' of ',num2str(experimentRepeats)])
%     % Have subject swallow between each run, press ok to continue
%     if experimentRepeats > 1
%         strng = 'Swallow. Press OK to continue.';
%         hh = msgbox(strng,'Alert');
%         th = findall(hh,'Type','Text'); 
%         th.FontSize = 14; 
%         hh.Position = [ 794.0000  674.0000  206.5000   62.8333];
%         uiwait(hh)
%     end
%end
    
disp('----- Finished with DPOAE Discrete experiment -----')
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
    %DPOAE.Ldp_3oct = [];
    %DPOAE.Ndp_3oct = [];
    %DPOAE.fc = [];
    
end
function [Y1,Y2,F1,F2,Fdp,Time,Y1c,Y2c,Ydp,Ydpc] = ARLas_dpoaeStim_fixedRatioDiscrete_vRT3(fs,f1,f2,fdp)
    % Create two vectors of DPOAE stimuli, f1 and f2.
    % The vectors are then folded into matrices for efficiency presenting through ARLas
    foldLength = 0.1; % fold into a matrix with columns this length (sec)
    foldN = round(foldLength * fs); % number of samples in each column
    len = 1; % length of total stimulus (sec)
    N = round(fs*len); % number of samples in total stimulus
    
    % the following lise are for truly fixed fRatio
    %fRatio = 1.22; % f2/fs ratio; now input argument
    %f1 = round(f2 / fRatio); % f1 primary frequencies
    %fdp = 2*f1 - f2; % cubic distortion dpoae frequencies

    time = (0:1:N-1)'/fs;
    
    y1 = sin(2*pi*f1*time);
    y2 = sin(2*pi*f2*time);
    ydp = sin(2*pi*fdp*time);
    y1c = cos(2*pi*f1*time);
    y2c = cos(2*pi*f2*time);
    ydpc = cos(2*pi*fdp*time);

    f1 = ones(size(y1))*f1;
    f2 = ones(size(y1))*f2;
    fdp = ones(size(y1))*fdp;

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
