function [] = goodmanMEMOC_v1(varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% goodmanMEMOC_v1(varargin)
%
% Measure MEMR (MOCR also possible) using a swept elicitor.
% This code modified from % Modified from goodmanMEMR_v2.m 
% This verison can test ipsilateral, contralalteral, or bilateral. Choose
% this by setting which probes to use. The noise has temporal gaps in it,
% during which the clicks are presented. The clicks are analyzed for MEMR,
% the CEOAEs are analyzed for MOCR. The levels of the clicks and noise need
% to be much lower to get an MOCR response (without MEMR), and thus, the
% averaging time needs to be much longer to obtain an adequate SNR.
%
% Authors: Shawn Goodman & Ehsan Khalili
% Auditory Research Lab, the University of Iowa
% Date: June 26, 2024
% Last Updated: June 26, 2024
% Last Updated: September 23, 2024 -- ssg. Code appears to be working fine
%           for MEMR. MOCR runs, but the levels and duration have not been optimized.
%           Need to implement a different analysis.
%
% TODO:
% Make noise scale appropriately. Now seems to stay at 120 dB SPL
% add post-hoc cal check
% new processing algorithm for MOCR
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    obj = varargin{1}; % get the arlas object
    V = '23SEP2024'; % this is the current version number of this program

%------ USER MODIFIABLE PARAMETERS ----------------------------------------
    testType = 'MEM' % set to 'MEM' or 'MOC'


%--------------------------------------------------------------------------
    if strcmp(testType,'MEM')
        targetNoise = 95;
        targetClick = 95;
        testLenMin = 3; % total test length in minutes
    elseif strcmp(testType,'MOC')
        targetNoise = 55;
        targetClick = 55;
        testLenMin = 8; % total test length in minutes
    else
        error('Unrecognized test type. Must be MOC or MEM')
    end


    multiplierNoise = 1/3; % this makes the levels come out right
    modLen = 4; % modulation duration in sec (from onset of noise to offset)

    % specify which ER10X probes to use    
    probeL = 'A'; % probe that is in LEFT ear; 'A', 'B', or []
    probeR = 'B'; % probe that is in RIGHT ear; 'A', 'B', or []
                  % Note: Either one or both ears can be tested at the same time.
                  % If either probe is empty set, will not be tested. 
                  
    % additional parameters that user usually doesn't need to mess with
    doMicCorrection = 1; % turn on and off microphone correction
    

    fmin =  100;
    fmax = 16000;
    
%------ END USER MODIFIABLE PARAMETERS ------------------------------------
%--------------------------------------------------------------------------
    disp(' ')
    disp('----- Starting memr experiment: memr -----')
    ISCL = struct;
    ISCR = struct;

    % specify calibration to use to set stimulus levels
    calType = 'thev'; % use 'thev' for Thevenin source
    targetCalType = 'fpl'; % 'fpl', 'spl', 'ipl'.
    % check to make sure all parameters are valid
    [probeR,probeL,calType,targetCalType] = checkInputs(probeR,probeL,calType,targetCalType);
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

    % PERFORM IN-SITU CALIBRATION -----
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

    % create stimuli ------------------------------------------------------
    
    % generate the clicks -----
    fs = obj.fs; % get the system sampling rate
    [nLoops,clickTrain,noiseSamples,nReps,H] = getStimuli(obj,modLen);
    % reshape the click stimuli into 100 ms chunks
    [rows,cols] = size(clickTrain);
    chunkLen = 0.1; % chunk size (sec)
    chunkSize = round(chunkLen*fs); % chunk size (samples)
    nChunks = floor(rows/chunkSize);
    time = (0:1:length(clickTrain)-1)'/fs;
    ClickTrain = reshape(clickTrain,chunkSize,nChunks);
    clickIndx = find(clickTrain>.01); % location of clicks
    %S1L = ClickTrain; % present clicks through left channel of probe 
    %S1R = S1L * 0;    % right channel is zeros

    nReps = ceil((testLenMin*60) / (length(clickTrain)/fs));
    nSweeps = nReps;
    nReps = nReps * nChunks;
    
    % generate noise -----
    [rows,cols] = size(ClickTrain);
    noise = randn(rows,cols*2);
    noise = noise(:)';
    [indx,nRejects] = AR(noise,'moderate',0);
    noise = AR_engine(noise,indx);
    noise = noise(1:rows*cols);
    noise = noise(:);
    h = H(:,1);
    Noise = noise .* h;
    Noise = Noise * multiplierNoise;
    %%%
    % in order to make ipsilateral noise, need to leave holes in the noise
    C = ClickTrain(:);
    N = Noise(:);
    %indx = find(C>.01); % location of each noise impulse
    holeDur1 = 0.003; % duration of hole (seconds)
    holeN1 = round(fs*holeDur1); % number of hole samplesii
    holeDur2 = 0.012; % duration of hole (seconds)
    holeN2 = round(fs*holeDur2); % number of hole samples
    Mask = ones(size(N));
    for ii=1:length(clickIndx)
        Mask(clickIndx(ii)-holeN1+1:clickIndx(ii)+holeN2) = 0;
    end
    Noise = Noise .* Mask;
    Noise = reshape(Noise,chunkSize,nChunks);
    %S2L = Noise;
    %S2R = Noise * 0;

%----
    % put noise on S2
    [rows,cols] = size(Noise);
    [Noise,~,~] = ARLas_applyISC(iscS1R,targetNoise,targetCalType,[fmin,fmax],Noise(:));
    Noise = reshape(Noise,rows,cols);
    %S2R = S2L * 0;

    % put clicks on S1
    %input = ClickTrain;
    %input = input * (max(abs(input(:))));
    [rows,cols] = size(ClickTrain);
    [Clicks,~,~] = ARLas_applyISC(iscS1R,targetClick,targetCalType,[fmin,fmax],ClickTrain(:));
    Clicks = reshape(Clicks,rows,cols);
    %S1R = S1L * 0;
    % if max(abs(S1L(:))) >= 1
    %     figure(111);
    %     plt(S1L(:))
    %     title('S1R')
    %     k = 1;
    % end
    % if max(abs(S1L(:))) >= 1
    %     figure(112);
    %     plt(S1L(:))
    %     title('S1L')
    %     k = 1;
    % end
    % if exist('k','var')
    %     if k==1
    %         keyboard
    %     end
    % end
    
    
% ----        
    
    % tell ARLas what to record and what to play % -----------------------------------
    obj.clearPlayList % clear out whatever was played previously
    if ~isempty(probeInputL) % Clicks are through inputL
        obj.setPlayList(Clicks,probeOutputL.ch(1));
        obj.setPlayList(Noise,probeOutputL.ch(2));
    end
    if ~isempty(probeInputR) % noise is thorugh inputR
        obj.setPlayList(Noise,probeOutputR.ch(1));
        obj.setPlayList(Noise*0,probeOutputR.ch(2));
    end

    obj.clearRecList % clear out whatever was used previously 
    if ~isempty(probeInputL)
        probeInputL.label = [probeInputL.label,'_memoc'];
        obj.setRecList(probeInputL.ch,probeInputL.label,probeInputL.micSens,probeInputL.gain);    
    end
    if ~isempty(probeInputR)
        probeInputR.label = [probeInputR.label,'_memoc'];
        obj.setRecList(probeInputR.ch,probeInputR.label,probeInputR.micSens,probeInputR.gain);    
    end

    obj.setNReps(nReps); % number of times to play stimulus
    obj.setFilter(1); % note: this is a highpass filter with a 75 Hz cutoff frequency.

    obj.objPlayrec.userInfo = []; % clear out previous, non-structure array info
    obj.objPlayrec.userInfo.fmin = fmin;
    obj.objPlayrec.userInfo.fmax = fmax;
    obj.objPlayrec.userInfo.fs = fs;
    obj.objPlayrec.userInfo.dpoae_version = V;
    obj.objPlayrec.userInfo.nSweeps = nSweeps;
    obj.objPlayrec.userInfo.nReps = nReps;
    obj.objPlayrec.userInfo.calType = calType;
    obj.objPlayrec.userInfo.targetCalType = targetCalType;
    % obj.objPlayrec.userInfo.target = target;
    obj.objPlayrec.userInfo.nChunks = nChunks;
    obj.objPlayrec.userInfo.chunkSize = chunkSize;
    obj.objPlayrec.userInfo.nSweeps = nSweeps;
    obj.objPlayrec.userInfo.clickIndx = clickIndx;
    obj.objPlayrec.userInfo.holeN1 = holeN1;
    obj.objPlayrec.userInfo.holeN2 = holeN2;
    obj.objPlayrec.userInfo.Mask = Mask;
    obj.objPlayrec.userInfo.targetNoise = targetNoise;
    obj.objPlayrec.userInfo.targetClick = targetClick;
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
    obj.objPlayrec.userInfo.micCorrectionL = micCorrectionL; % added 4/1/2022 -- ssg
    obj.objPlayrec.userInfo.micCorrectionR = micCorrectionR; % added 4/1/2022 -- ssg

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

    DataL = reshape(DataL,chunkSize*nChunks,nSweeps);
    DataR = reshape(DataR,chunkSize*nChunks,nSweeps);
    DataL = applyMicCorrection(DataL,micCorrectionL);
    DataR = applyMicCorrection(DataR,micCorrectionR);

    MEMOC.DataL = DataL;
    MEMOC.DataR = DataR;
    MEMOC.fs = fs;
    MEMOC.time = time;
    MEMOC.nSweeps = nSweeps;
    MEMOC.headerL = headerL;
    MEMOC.headerR = headerR;
    MEMOC.ModLen = modLen;
    MEMOC.clickIndx = clickIndx;
    MEMOC.iscS1L = iscS1L;
    MEMOC.iscS1R = iscS1R;
    
    % Data Analysis -------------------------------------------------------
    [MEMOC_inc,MEMOC_mem,MEMOC_moc,h1,h2] = analyzeMEMOC_v1(headerL,DataL,headerR,DataR,'Dummy',1);


    % Save analyses ------------------------
    disp('----- Saving analysis files -----')
    try 
        savePath = [obj.savingGrace,obj.experimentID,'_analysis',obj.sep];
        if exist(savePath,'dir') == 0 % if path does not exist
            success = mkdir(savePath);
            if ~success
                warning('Unable to create new Experiment Directory: data save path')
            end
        end 
        saveFileName = [obj.subjectID,'_analyzedMEMOC.mat'];
        saveFileName = ARLas_saveName(savePath,saveFileName);
        save([savePath,saveFileName],'MEMOC','MEMOC_inc','MEMOC_mem','MEMOC_moc')
    catch
        disp('Error saving analysis!')
    end

    % Save figures ----------
    disp('----- Saving figures -----')
    try
        savePath = [obj.savingGrace,obj.experimentID,'_analysis',obj.sep];
        if exist(savePath,'dir') == 0 % if path does not exist
            success = mkdir(savePath);
            if ~success
                warning('Unable to create new Experiment Directory: data save path')
            end
        end 
        
        figureFileName = [obj.subjectID,'_memr.fig'];
        figureFileName = ARLas_saveName(savePath,figureFileName);
        savefig(h1,[savePath,figureFileName])
        figureFileName = [obj.subjectID,'_memr.bmp'];
        figureFileName = ARLas_saveName(savePath,figureFileName);    
        saveas(h1,[savePath,figureFileName])
        
        figureFileName = [obj.subjectID,'_mocr.fig'];
        figureFileName = ARLas_saveName(savePath,figureFileName);
        savefig(h2,[savePath,figureFileName])
        figureFileName = [obj.subjectID,'_mocr.bmp'];
        figureFileName = ARLas_saveName(savePath,figureFileName);    
        saveas(h2,[savePath,figureFileName])
    catch
        disp('Error saving figures!')
    end
    
    %keyboard
    %return

    
end

% INTERNAL FUNCTIONS ------------------------------------------------------
function [nLoops,clickTrain,noiseSamples,nReps,H] = getStimuli(obj,modLen)
    % create click stimuli ---------------------------------------------------------
    fs = obj.fs; % get the system sampling rate
    clickN = 1; % number of samples in the click
    clickLen = 0.05; % desired stimulus length in seconds
    nSamples = round(clickLen * fs); % total number of samples
    if mod(nSamples,2) ~= 0 % make total number of samples even
        nSamples = nSamples + 1;
    end
    click = zeros(nSamples,1);
    startDelay = 0.003; % silence before the click onset (s)
    startSample = round(startDelay * fs);
    click(startSample:startSample + (clickN-1)) = 1;
     %dBatten = 23; % dB attenuation re: full out
     %multiplier = 10^(-dBatten/20);
     %click = click * multiplier;
    % Paige's threshold is 58 dB attenuation re: full output
    % Presentation level = 35 dB SL = 58 - 35 = 23 dB atten
    % IEC711 coupler results: 
    %   0.03 Pa peak
    %   0.03 + 0.024 = 0.054 peak to peak
    %   0.054 / 2 = 0.027 peak equivalent
    %   63.5 dB pSPL
    %   68.6 dB peak to peak SPL
    %   62.6 dB peSPL


    % create noise ------------------------------------------------------------
    clicksPerCondition = 1000; % total number of clicks recorded per condition (noise or no noise)
    %modLen = 8; % modulation duration in sec (from onset of noise to offset)
    loopDuration = modLen * 30; %240; %60; % seconds per loop
    clicksPerLoop = loopDuration; % because one second gives you 1 click at each level
    nLoops = ceil(clicksPerCondition / clicksPerLoop); % number of loops
    noiseSamples = round(modLen * fs); % total 
    nReps = loopDuration / (modLen * 2); % number of repeitions per loop
    %nReps = loopDuration / (modLen * 4); % number of repeitions per loop
    %h = linspace(0,70,noiseSamples/2)';
    h = linspace(0,70,noiseSamples)';
    h(1) = eps;
    h = [h;flipud(h)];
    h = 10.^(h/20); % noise amplitude in linear units
    h = h / max(abs(h)) * 1; % rescale to unit amplitude
    H = repmat(h,1,nReps);
    %dBattenN = 1; % dB attenuation re: full out
    %multiplierNoise = 10^(-dBattenN/20);
    % Paige's threshold is 65 dB attenuation re: full output
    % Presentation level = full output - 1 dB = 64 dB SL
    % IEC711 coupler results: 
    %   1.3 Pa peak
    %   96.3 dB SPL at peak intensity

    % finish making click train
    clicksPerTrain = (noiseSamples / nSamples)*2;
    clickTrain = repmat(click,1,clicksPerTrain);
    clickTrain = clickTrain(:);
end
function [probeR,probeL,calType,targetCalType] = checkInputs(probeR,probeL,calType,targetCalType)
    if strcmp(probeR,probeL)
        error('Variables probeR and probeL cannot be the same letter.')
    end
    if ~(strcmp(probeR,'A') || strcmp(probeR,'B') || ~isempty(probeR))
        error('Variable probeR contains unrecognized value. Must be string A, B, or empty set.')
    end
    if ~(strcmp(probeL,'A') || strcmp(probeL,'B') || ~isempty(probeL))
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
    dL = dir([postPath,'Ch',num2str(probeInputL.ch),'_',probeInputL.label,'_dpoae_*.mat']);
    dR = dir([postPath,'Ch',num2str(probeInputR.ch),'_',probeInputR.label,'_dpoae_*.mat']);
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
function [] = writeData2Excel(pathName,fileName,fc,Ldp_3oct_L,Ndp_3oct_L,Ldp_3oct_R,Ndp_3oct_R,L1,L2)
    warning off
    sheet = 'dpoae';
    nLevels = size(fc,2);
    nFreqs = size(fc,1);
    
    try
    % write LEFT ear -------------------------------
        range = 'A3';
        xlswrite([pathName,fileName],{'F2 Freq (Hz)'},sheet,range);
        range = 'A4';
        xlswrite([pathName,fileName],fc(:,1),sheet,range);

        N = nFreqs + 6;
        range = ['A',num2str(N+2)];
        xlswrite([pathName,fileName],{'F2 Freq (Hz)'},sheet,range);
        range = ['A',num2str(N+3)];
        xlswrite([pathName,fileName],fc(:,1),sheet,range);
        
        letters = {'B','C','D','E','F','G'};
        counter = 1;
        for ii=1:nLevels
            range = [char(letters(counter)),'1'];
            xlswrite([pathName,fileName],{['Stim Lvl:',num2str(L1(ii)),'\',num2str(L2(ii))]},sheet,range);
            range = [char(letters(counter)),'2'];
            status(1,1) = xlswrite([pathName,fileName],{'LEFT EAR'},sheet,range);
            range = [char(letters(counter)),'3'];
            status(3,1) = xlswrite([pathName,fileName],{'DPOAE Lvl (dB SPL)'},sheet,range);
            range = [char(letters(counter)),'4'];
            status(4,1) = xlswrite([pathName,fileName],Ldp_3oct_L(:,ii),sheet,range);
            counter = counter + 1;
            range = [char(letters(counter)),'1'];
            xlswrite([pathName,fileName],{['Stim Lvl:',num2str(L1(ii)),'\',num2str(L2(ii))]},sheet,range);
            range = [char(letters(counter)),'2'];
            status(1,1) = xlswrite([pathName,fileName],{'LEFT EAR'},sheet,range);
            range = [char(letters(counter)),'3'];
            status(3,1) = xlswrite([pathName,fileName],{'Noise Lvl (dB SPL)'},sheet,range);
            range = [char(letters(counter)),'4'];
            status(4,1) = xlswrite([pathName,fileName],Ndp_3oct_L(:,ii),sheet,range);
            counter = counter + 1;
        end
    catch
        disp('Warning: DPOAE_L Analysis not written to Excel file!')
    end
    
%     try
%     % write RIGHT ear -------------------------------
%         counter = 1;
%         N = nFreqs + 6;
%         for ii=1:nLevels
%             range = [char(letters(counter)),num2str(N)];
%             xlswrite([pathName,fileName],{['Stim Lvl:',num2str(L1(ii)),'\',num2str(L2(ii))]},sheet,range);
%             range = [char(letters(counter)),num2str(N+1)];
%             status(1,1) = xlswrite([pathName,fileName],{'Right EAR'},sheet,range);
%             range = [char(letters(counter)),num2str(N+2)];
%             status(3,1) = xlswrite([pathName,fileName],{'DPOAE Lvl (dB SPL)'},sheet,range);
%             range = [char(letters(counter)),num2str(N+3)];
%             status(4,1) = xlswrite([pathName,fileName],Ldp_3oct_R(:,ii),sheet,range);
%             counter = counter + 1;
%             range = [char(letters(counter)),num2str(N)];
%             xlswrite([pathName,fileName],{['Stim Lvl:',num2str(L1(ii)),'\',num2str(L2(ii))]},sheet,range);
%             range = [char(letters(counter)),num2str(N+1)];
%             status(1,1) = xlswrite([pathName,fileName],{'Right EAR'},sheet,range);
%             range = [char(letters(counter)),num2str(N+2)];
%             status(3,1) = xlswrite([pathName,fileName],{'Noise Lvl (dB SPL)'},sheet,range);
%             range = [char(letters(counter)),num2str(N+3)];
%             status(4,1) = xlswrite([pathName,fileName],Ndp_3oct_R(:,ii),sheet,range);
%             counter = counter + 1;
%         end
%     catch
%         disp('Warning: DPOAE_R Analysis not written to Excel file!')
%     end
    
    % do this last so cursor goes to top left of the sheet
    range = 'A1'; % center frequency ---
    xlswrite([pathName,fileName],{' '},sheet,range);
    
    
    % get rid of first, default tab
    excelFileName = fileName;
    excelFilePath = pathName;
    sheetName = 'Sheet';
    objExcel = actxserver('Excel.Application');
    objExcel.Workbooks.Open(fullfile(excelFilePath,excelFileName));
    try
        objExcel.ActiveWorkbook.Worksheets.Item([sheetName '1']).Delete;
    catch
    end
    
    % widen the column widths from default 8.26 to 12
    ws = objExcel.ActiveWorkbook.Worksheets.Item([sheet]);
    ws.Range('A1').ColumnWidth = 15;
    ws.Range('B1').ColumnWidth = 15;
    ws.Range('C1').ColumnWidth = 15;
    ws.Range('D1').ColumnWidth = 15;
    ws.Range('E1').ColumnWidth = 15;
    ws.Range('F1').ColumnWidth = 15;
    ws.Range('G1').ColumnWidth = 15;
    ws.Range('H1').ColumnWidth = 15;
    
    % close everything down
    objExcel.ActiveWorkbook.Save;
    objExcel.ActiveWorkbook.Close;
    objExcel.Quit;
    objExcel.delete;
    warning on
end

% OLD CODE ----------------------------------------------------------------

%     
%     [saveName] = ARLas_saveName(savePathName,saveName);
%     save([savePathName,saveName],'DataL','DataR','fs','time','nSweeps',...
%         'headerL','headerR','modLen','clickIndx','iscS1L','iscS1R')
%     disp('Done!')

%     save([savePathName,saveName],'DataL','DataR','fs','time','nSweeps',...
%         'headerL','headerR','modLen','clickIndx','iscS1L','iscS1R')


%         for jj=1:length(indx)
%             if indx(jj)-windowN/2 >= 1 & indx(jj)-windowN/2 < size(DataL,1)
%             chunk = DataL(indx(jj)-windowN/2:indx(jj)+windowN/2,:);
%             %chunk = mean(chunk,2);
%             timeChunk(1,jj) = time(indx(jj));
%             
%             %cc = chunk(238:300,:);
%             cc = chunk(477:540,:);
%             
%             [frequency,signal,noiseFloor] = ARLas_fda(cc(10:end,:),fs,0.00002,size(cc(10:end,:),1));
%             
%             
%             
%             %for kk=1:size(cc,2)
%             %    cc(:,kk) = cc(:,kk) - cc(1,kk);
%             %end
%             ccmean = mean(cc,2);
%             [~,ccMaxIndx] = min(ccmean);
%             %RMS(:,jj) = cc(ccMaxIndx,:)';
%             
%             RMS(:,jj) = sum(10.^(signal(2:14)/20)*0.00002);
%             
%             %RMS(:,jj) = max(cc(10:20));
%             %RMS(:,jj) = sqrt(mean(ccmean.^2,1));
%             %RMS(:,jj) = max(chunk);
%             end
%         end
%         
%         %[frequency,signal,noiseFloor] = ARLas_fda(chunk,fs,0.00002,size(chunk,1));
%  
%         
%         figure
%         subplot(2,1,1)
%         plot(time,Noise(:,1),'k')
%         xlim([time(1),time(end)])
%         subplot(2,1,2)
%         %plot(timeChunk(2:end)',RMS(:,2:end),'b')
%         hold on
%         %plot(timeChunk,mean(RMS,1),'b','LineWidth',2)
%         plot(timeChunk(2:end)',mean(RMS(:,2:end),1),'r')
%         xlim([time(1),time(end)])
%         

% October 10
        % Assume 76 microsecond travel time in ear canal.
        % Then there is a delay of 1.5 ms, round trip, before you see the
        % effects of the MEMR. Therefore, shift the window by this amount,
        % which is 146 samples.

        % first get template
%         for jj=1:6
%             if indx(jj)-windowN/2 >= 1 & indx(jj)-windowN/2 < size(DataL,1)
%                 chunk = DataL(indx(jj)-windowN/2:indx(jj)+windowN/2,:);
%                 cc = chunk(477:540,:);
%                 [Frequency,Signal(:,jj),NoiseFloor(:,jj)] = ARLas_fda(cc(10:end,:),fs,0.00002,size(cc(10:end,:),1));
%             end
%         end


% for kk=1:size(DataL,2)
%     [pl,Pl,phi,other,wf] = ARLas_convertPL(DataL,iscS1L);    
% end
% if ii==1
%     DATA_L = [];
%     DATA_R = [];
% end
% DATA_L = MOC_shawnIshan_analyzeV_1(headerL,DataL,'L',DATA_L);
% DATA_R = MOC_shawnIshan_analyzeV_1(headerR,DataR,'R',DATA_R);

%     Noise = Noise(:,1:nSweeps);
%     NNN = [];
%     for jj = 1:nSweeps
%         dummy = Noise(:,jj);
%         dummy = reshape(dummy,chunkSize,nChunks);
%         NNN = [NNN,dummy];
%     end
%     Noise = NNN;
%     clear NNN;
%     %NN = round(noiseSamples * nReps);
%     NN = round(noiseSamples*2 * nReps);
%     noise = randn(1,round(NN*1.1));  
%     [indx,nRejects] = AR(noise,'moderate',0);
%     noise = AR_engine(noise,indx);
%     noise = noise';
%     noise = noise(1:NN);
%     Noise = reshape(noise,noiseSamples*2,nReps);
%     Noise = Noise .* H;
%     Noise = Noise / max(abs(max(abs(Noise))));
%     %Noise = [Noise;Noise*0];
%     time = (0:1:size(Noise,1)-1)'/fs;
%     Noise = Noise * multiplierNoise;
%     %        windowN = round(fs*0.01);
%     %         for jj=2:length(indx)-1
%     %             Noise(indx(jj)-windowN:indx(jj)+windowN,:) = 0;
%     %         end
%     Noise = Noise(:,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %         SIGNAL = mean(Signal,2);
%          sigStart = 1;
%          sigFinish = 15;
% %         SIGNAL = SIGNAL(sigStart:sigFinish);
%          windowN = round(fs*0.01);
% 
%         
%         startOffset = 1; %13;
%         finishOffset = 60; % 60
%         whichOne = 1;
%         for jj=1:length(indx)
%             if indx(jj)-windowN/2 >= 1 & indx(jj)-windowN/2 < size(DataL,1)
%                 chunk = DataL(indx(jj)-windowN/2:indx(jj)+windowN/2,:);
%                 timeChunk(1,jj) = time(indx(jj));
%                 [~,hat] = max(abs(mean(chunk,2)));
%                 cc = chunk(hat+startOffset:hat+finishOffset,:); %chunk(477:540,:);
%                 [frequency,signal,noiseFloor] = ARLas_fda(cc,fs,0.00002,size(cc,1));
%                 signal = signal(sigStart:sigFinish);
%                 delta = sum(10.^(signal(whichOne)/20)*0.00002);
%                 RMS(:,jj) = 20*log10(delta/0.00002);
%             end
%         end
%         RMS = RMS(:,2:end);
%         RMS = RMS - RMS(1);
%         figure(10)
%         hold on
%         plot(timeChunk(2:end)',mean(RMS,1))
% 
% %---------------------------------------------------------------        
%  keyboard
%  
%         %figure
%         subplot(2,1,1)
%         plot(time,Noise(:,1),'k')
%         xlim([time(1),time(end)])
%         subplot(2,1,2)
%         hold on
%         plot(timeChunk(2:end)',mean(RMS(:,2:end),1))
%         xlim([time(1),time(end)])
%         
%         
% keyboard            
%         % PERFORM POST-TEST IN-SITU CALIBRATION -----
%         if ~isempty(probeInputL)
%             if contains(probeInputL.label,'dpoae')
%                 dummy = probeInputL.label;
%                 dummy = dummy(1:end-6);
%                 probeInputL.label = dummy;
%             end
%         end
%         if ~isempty(probeInputR)
%             if contains(probeInputR.label,'dpoae')
%                 dummy = probeInputR.label;
%                 dummy = dummy(1:end-6);
%                 probeInputR.label = dummy;
%             end
%         end
%         disp('----- Running post-test in-situ calibration -----')
%         inSituReps = 6;
%         doIndividual = 1;
%         doSimultaneous = 0;
%         [iscS1L,iscS2L,iscS12L] = ARLas_runISC(obj,probeInputL,probeOutputL,calPath1L,calPath2L,calType,fmin,fmax,micCorrectionL,inSituReps,doIndividual,doSimultaneous);
%         [iscS1R,iscS2R,iscS12R] = ARLas_runISC(obj,probeInputR,probeOutputR,calPath1R,calPath2R,calType,fmin,fmax,micCorrectionR,inSituReps,doIndividual,doSimultaneous);       
%         ISCL.(['Rec',num2str(ii+1)]).(['S1']) = iscS1L;
%         ISCL.(['Rec',num2str(ii+1)]).(['S2']) = iscS2L;
%         ISCR.(['Rec',num2str(ii+1)]).(['S1']) = iscS1R;
%         ISCR.(['Rec',num2str(ii+1)]).(['S2']) = iscS2R;            
%         
%         % Analyze data ----------------------------------------------------
%         disp('----- Analyzing data -----')
%         % analysis prep -----
%         if ii==1 % on the first iteration, initialize the analysis structure
%             if ~isempty(probeInputL)
%                 DPOAE_L = getDPstruct(obj,fs,fmin,fmax,sweepRate,nSweeps(ii),targetL1(ii),targetL2(ii),calType,targetCalType,iscS1L,iscS2L,iscS12L,C1L,C2L,'Left');
%             end
%             if ~isempty(probeInputR)
%                 DPOAE_R = getDPstruct(obj,fs,fmin,fmax,sweepRate,nSweeps(ii),targetL1(ii),targetL2(ii),calType,targetCalType,iscS1R,iscS2R,iscS12R,C1R,C2R,'Right');
%             end
%         end 
%         
%         if runPostAnalysis == 1 % running post-hoc analyses on previously recorded data
%             dummy = load([postPath,dL(ii).name]);
%             headerL = dummy.header;
%             DataL = dummy.data;
%             dummy = load([postPath,dR(ii).name]);
%             headerR = dummy.header;
%             DataR = dummy.data;
%             clear dummy
% 
%             % look to see if mic correction exists as a field. If recording
%             % made before 3/23/2022, will not exist. 
%             if isfield(headerL,'micCorrectionL')
%                 micCorrectionL = headerL.micCorrectionL;
%             else % backwards compatability hack
%                 micCorrectionL = 1;
%             end
%             if isfield(headerL,'micCorrectionR')
%                 micCorrectionR = headerL.micCorrectionR;
%             else % backwards compatability hack
%                 micCorrectionR = 1;
%             end
%             DataL = applyMicCorrection(DataL,micCorrectionL);
%             DataR = applyMicCorrection(DataR,micCorrectionR);
%             
%             DPOAE_L.postHocTimeStamp = DPOAE_L.timeStamp; % note that this is a post-hoc analysis
%             DPOAE_L.timeStamp = obj.timeStamp; % make sure that original time stamp shows up on the figues
%             DPOAE_R.postHocTimeStamp = DPOAE_R.timeStamp; % note that this is a post-hoc analysis
%             DPOAE_R.timeStamp = obj.timeStamp; % make sure that original time stamp shows up on the figues
%         end    
%         
%         % do the actual analysis here ---
%         if ~isempty(probeInputL)
%             DPOAE_L = ARLas_dpoaeAnalysis_fixedRatioContinuous_v1(headerL,DataL,DPOAE_L,[]);
%         end
%         if ~isempty(probeInputR)
%             DPOAE_R = ARLas_dpoaeAnalysis_fixedRatioContinuous_v1(headerR,DataR,DPOAE_R,[]);
%         end
%         
%         disp('----- Plotting data -----')
%         try
%             if ~isempty(probeInputL)
%                 [h1L,h2L] = ARLas_dpoaePlotting_fixedRatioContinuous_v1(DPOAE_L);
%                 %h1L.Position = [591   154  469  211];
%                 %h2L.Position = [591  442  469  314];
%                 pause(0.01)
%             end
%         catch
%         end
%         try
%             if ~isempty(probeInputR)
%                 [h1R,h2R] = ARLas_dpoaePlotting_fixedRatioContinuous_v1(DPOAE_R);
%                 %h1R.Position = [1060  154  469  211];
%                 %h2R.Position = [1060  442  469  314];
%                 pause(0.01);
%             end
%         catch
%         end
%         
% 
%     % write main results to excel file here ----------
%     disp('----- Writing Excel files -----')
%     fc = DPOAE_L.fc;
%     Ldp_3oct_L = DPOAE_L.Ldp_3oct;
%     Ndp_3oct_L = DPOAE_L.Ndp_3oct;
%     Ldp_3oct_R = DPOAE_R.Ldp_3oct;
%     Ndp_3oct_R = DPOAE_R.Ndp_3oct;
%     L1 = DPOAE_L.targetL1;
%     L2 = DPOAE_L.targetL2;
%     if runPostAnalysis == 1 % saving post-hoc analyses on previously recorded data
%         savePath = postPathSave;
%     else % saving newly acquired data
%         savePath = [obj.savingGrace,obj.experimentID,'_analysis',obj.sep];
%         if exist(savePath,'dir') == 0 % if path does not exist
%             success = mkdir(savePath);
%             if ~success
%                 warning('Unable to create new Experiment Directory: data save path')
%             end
%         end 
%     end
%     saveFileName = [obj.subjectID,'_analyzedDPOAE.xls'];
%     saveFileName = ARLas_saveName(savePath,saveFileName);
%     writeData2Excel(savePath,saveFileName,fc,Ldp_3oct_L,Ndp_3oct_L,Ldp_3oct_R,Ndp_3oct_R,L1,L2)
%     
% 
%     disp('----- Finished with memr experiment -----')
%     disp(' ')
%     disp(' ')

%     if modLen == 1
%         nSweeps = 80;
%     elseif modLen == 2
%         nSweeps = 40;
%     elseif modLen == 4
%         nSweeps = 20;
%     elseif modLen == 8
%         nSweeps = 10;
%     end

% specify nominal stimulus parameters:
%   Note: these will be modified slightly at low and high frequencies

% function [DPOAE] = getDPstruct(obj,fs,fmin,fmax,sweepRate,nSweeps,targetL1,targetL2,calType,targetCalType,iscS1,iscS2,iscS12,C1,C2,ear)
%     DPOAE.subjID = obj.subjectID;
%     DPOAE.ear = ear;
%     DPOAE.timeStamp = cellstr(datetime('now'));
%     DPOAE.fs = fs;
%     DPOAE.f2min = fmin;
%     DPOAE.f2max = fmax;
%     DPOAE.sweepRate = sweepRate;
%     DPOAE.nSweeps = nSweeps;
%     DPOAE.targetL1 = targetL1;
%     DPOAE.targetL2 = targetL2;
%     DPOAE.calType = calType;
%     DPOAE.targetCalType = targetCalType;
%     DPOAE.iscS1 = iscS1;
%     DPOAE.iscS2 = iscS2;
%     DPOAE.iscS12 = iscS12;
%     DPOAE.C1 = C1;
%     DPOAE.C2 = C2;
% 
%     DPOAE.L1 = [];
%     DPOAE.N1 = [];
%     DPOAE.P1 = [];
%     DPOAE.L2 = [];
%     DPOAE.N2 = [];
%     DPOAE.P2 = [];
%     DPOAE.Ldp = [];
%     DPOAE.Ndp = [];
%     DPOAE.Pdp = [];
%     DPOAE.f1 = [];
%     DPOAE.f2 = [];
%     DPOAE.fdp = [];
%     
%     DPOAE.postHocTimeStamp = [];
%     DPOAE.Ldp_3oct = [];
%     DPOAE.Ndp_3oct = [];
%     DPOAE.fc = [];
%     
% end
%         
%         if ~isempty(probeInputL)
%             if runPostAnalysis == 1 % saving post-hoc analyses on previously recorded data
%                 savePath = postPathSave;
%             else
%                 savePath = [obj.savingGrace,obj.experimentID,'_analysis',obj.sep];
%                 %savePath = obj.objPlayrec.savedFilesPath;
%                 if exist(savePath,'dir') == 0 % if path does not exist
%                     success = mkdir(savePath);
%                     if ~success
%                         warning('Unable to create new Experiment Directory: data save path')
%                     end
%                 end 
%             end
%             figureFileName = [obj.subjectID,'_dpoaePrimariesL.fig'];
%             figureFileName = ARLas_saveName(savePath,figureFileName);
%             savefig(h1L,[savePath,figureFileName])
%             figureFileName = [obj.subjectID,'_dpoaePrimariesL.bmp'];
%             figureFileName = ARLas_saveName(savePath,figureFileName);    
%             saveas(h1L,[savePath,figureFileName])
% 
%             figureFileName = [obj.subjectID,'_dpoaeL.fig'];
%             figureFileName = ARLas_saveName(savePath,figureFileName);
%             savefig(h2L,[savePath,figureFileName])
%             figureFileName = [obj.subjectID,'_dpoaeL.bmp'];
%             figureFileName = ARLas_saveName(savePath,figureFileName);    
%             saveas(h2L,[savePath,figureFileName])
%         end
%     catch
%         disp('Warning: One or more figures for DPOAE_L not saved!')        
%     end
% 
%     try
%         if ~isempty(probeInputR)
%             if runPostAnalysis == 1 % saving post-hoc analyses on previously recorded data
%                 savePath = postPathSave;
%             else % saving newly acquired data
%                 savePath = [obj.savingGrace,obj.experimentID,'_analysis',obj.sep];
%                 if exist(savePath,'dir') == 0 % if path does not exist
%                     success = mkdir(savePath);
%                     if ~success
%                         warning('Unable to create new Experiment Directory: data save path')
%                     end
%                 end 
%             end
%             figureFileName = [obj.subjectID,'_dpoaePrimariesR.fig'];
%             figureFileName = ARLas_saveName(savePath,figureFileName);
%             savefig(h1R,[savePath,figureFileName])
%             figureFileName = [obj.subjectID,'_dpoaePrimariesR.bmp'];
%             figureFileName = ARLas_saveName(savePath,figureFileName);    
%             saveas(h1R,[savePath,figureFileName])
% 
%             figureFileName = [obj.subjectID,'_dpoaeR.fig'];
%             figureFileName = ARLas_saveName(savePath,figureFileName);
%             savefig(h2R,[savePath,figureFileName])
%             figureFileName = [obj.subjectID,'_dpoaeR.bmp'];
%             figureFileName = ARLas_saveName(savePath,figureFileName);    
%             saveas(h2R,[savePath,figureFileName])
%         end
%     catch
%        disp('Warning: One or more figures for DPOAE_R not saved!')
%     end


%         
%         if ~isempty(probeInputL)
%             if runPostAnalysis == 1 % saving post-hoc analyses on previously recorded data
%                 savePath = postPathSave;
%             else % saving newly acquired data
%             end
%             saveFileName = [obj.subjectID,'_analyzedDPOAE_L.mat'];
%             saveFileName = ARLas_saveName(savePath,saveFileName);
%             save([savePath,saveFileName],'DPOAE_L','ISCL')
%         end
%     catch
%         disp('Warning: DPOAE_L Analysis not saved!')
%     end
%     
%     try 
%         if ~isempty(probeInputR)
%             if runPostAnalysis == 1 % saving post-hoc analyses on previously recorded data
%                 savePath = postPathSave;
%             else % saving newly acquired data
%                 savePath = [obj.savingGrace,obj.experimentID,'_analysis',obj.sep];
%                 if exist(savePath,'dir') == 0 % if path does not exist
%                     success = mkdir(savePath);
%                     if ~success
%                         warning('Unable to create new Experiment Directory: data save path')
%                     end
%                 end 
%             end
%             saveFileName = [obj.subjectID,'_analyzedDPOAE_R.mat'];
%             saveFileName = ARLas_saveName(savePath,saveFileName);
%             save([savePath,saveFileName],'DPOAE_R','ISCR')
%         end
%     catch
%         disp('Warning: DPOAE_R Analysis not saved!')
%     end
%     

