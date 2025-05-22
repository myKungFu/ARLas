function [] = LGF(varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LGF(varargin)
%
% Measure DPOAE leve growth functions (LGF) using fixed L1 and sweeping L2.
% See Gavin et al, 2023, JASA Express for details.
%
% Still to do:
% First Priority:
% 1) make this work when looping over multiple levels at once
% Second Priority:
% 5) Decide whether to use a "-3 dB offset"
% 6) Fix EPL issues 
%
%
% Author: Shawn Goodman, PhD
% Auditory Research Lab, the University of Iowa
% Original Date: March 14, 2022
% Last Updated: April 15, 2025 -- ssg
% Last Updated: May 21, 2025 -- ssg -- used new load calibration code
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    obj = varargin{1}; % get the arlas object
    V = '21MAY2025'; % this is the current version number of this program

%------ USER MODIFIABLE PARAMETERS ----------------------------------------
%--------------------------------------------------------------------------

    f2 = 2000; % f2 frequency (Hz) 
    nSweeps = 2; % number of test sweeps 24

    L2min = 0; % 0 min L2 (dB FPL) 
    L2max = 70; % 70 max L2 (dB FPL)

    doEPL = 0; % turn on and off EPL correction 

    % Note: You now specify which probes and mic cal to use in loadCalibration.m

%------ END USER MODIFIABLE PARAMETERS ------------------------------------
%--------------------------------------------------------------------------

    % perform load calibration
    doProbeFitCheck = 1;
    iscCheckHandleA = []; % figures not yet made so initialize handles to zero
    iscCheckHandleB = [];
    LC = []; % load calibration not yet done, so initialize to empty set
    [LC,iscCheckHandleA,iscCheckHandleB,OK] = loadCalibration(obj,doProbeFitCheck,LC,iscCheckHandleA,iscCheckHandleB);
    if OK == 0
        disp('Error: Probe fit not okay. Aborting program.')
        return
    end


    % edge effects mean that you need to extend the range a bit 
    L2min = L2min - 2;
    L2max = L2max + 3;

    targetL1 = 65 * ones(size(f2)); % target levels for fixed L1 primaries (dB FPL);
    fRatio = ones(size(f2))*1.22; % f2/f1 ratios 
    f2 = round(f2); % vector of f2 frequencies (force to be intergers)
    f1 = round(f2 ./ fRatio); % vector of f1 frequencies (force to be integers)
    fdp = 2*f1 - f2; % vector of dpoae frequencies
    fmin = min(f1); % minimum frequency being tested
    fmax = max(f2); % maximum frequency being tested

    sweepLength = 9; % number of seconds for each sweep
    sweepRate = sweepLength / (L2max - L2min); % dB per second
    sweepN = round(sweepLength * obj.fs); % number of samples in each sweep
    testLength = sweepLength * nSweeps; % total test length (sec)

    disp(' ')
    disp('----- Starting DPOAE experiment: LGF -----')
    disp(['Estimated Test Length = ',num2str(testLength/60),' min.'])
    
% If presenting X dB FPL, this is ~3-6 dB higher in dB SPL.
% In order to make more comparable with other studies that calibrate in
% dB SPL, can adjust the level accordingly:
% if strcmp(calType,'thev') 
%     if strcmp(targetCalType,'fpl')
%         adjust = -3; % if 55 dB fpl, this is more like 58 dB SPL. Therefore, subtract 3 to make the level 52 dB FPL ~= 55 dB SPL.
%     else
%         adjust = 0;
%     end
% end
adjust = 0; % hard coded 4/21/2025 ssg
        
    % COLLECT THE DATA --------------------------------------------------------
    disp('----- Collecting data -----')
    fs = obj.fs; % get the system sampling rate

    % Apply the in-situ calibration to the stimuli ------------------------
    iscSA1 = LC.iscSA1;
    iscSA2 = LC.iscSA2;
    iscSB1 = LC.iscSB1;
    iscSB2 = LC.iscSB2;
    probeInputA = LC.probeInputA;
    probeOutputA = LC.probeOutputA;
    probeInputB = LC.probeInputB;
    probeOutputB = LC.probeOutputB;

    nFreqs = length(f1);
    for kk=1:nFreqs
        [Y1,Y2,F1,F2,Fdp,Time,Y1c,Y2c,Ydp,Ydpc,L1,L2] = ARLas_dpoaeStim_sweepL2Continuous(fs,f1(kk),f2(kk),fdp(kk),L2min,L2max,sweepN); % get the raw stimuli
        if ~isempty(iscSA1)
            [rows,cols] = size(Y1);
            [SA1,scaling,errors1] = ARLas_applyISC(iscSA1,targetL1(kk)+adjust,'fpl',f1(kk),Y1(:));
            SA1 = Y1 .* abs(scaling);
            Y2 = Y2 ./ max(Y2(:));
            [SA2,scaling,errors1] = ARLas_applyISC(iscSA2,L2max+adjust,'fpl',f2(kk),Y2(:));
            SA2 = Y2 .* abs(scaling);
            if kk==1
                SA1full = SA1;
                SA2full = SA2;
            else
                SA1full = SA1 + SA1full;
                SA2full = SA2 + SA2full;
            end                
        else
            SA1 = [];
            SA2 = [];
            SA1full = [];
            SA2full = [];
        end
        if ~isempty(iscSB1)
            [SB1,scaling,errors1] = ARLas_applyISC(iscSB1,targetL1(kk)+adjust,'fpl',f1(kk),Y1(:));
            SB1 = Y1 .* abs(scaling);
            Y2 = Y2 ./ max(Y2(:));
            [SB2,scaling,errors1] = ARLas_applyISC(iscSB2,L2max+adjust,'fpl',f2(kk),Y2(:));
            SB2 = Y2 .* abs(scaling);
            if kk==1
                SB1full = SB1;
                SB2full = SB2;
            else
                SB1full = SB1 + SB1full;
                SB2full = SB2 + SB2full;
            end                
        else
            SB1 = [];
            SB2 = [];
            SB1full = [];
            SB2full = [];            
        end
    end
    SA1 = SA1full;
    SA2 = SA2full;
    SB1 = SB1full;
    SB2 = SB2full;
    targetL2 = 20*log10(L2/.00002);


    % tell ARLas what to record and what to play --------------------------
    obj.clearPlayList % clear out whatever was played previously
    if ~isempty(probeInputA) % load new signals to play
        obj.setPlayList(SA1,probeOutputA.ch(1));
        obj.setPlayList(SA2,probeOutputA.ch(2));
    end
    if ~isempty(probeInputB)
        obj.setPlayList(SB1,probeOutputB.ch(1));
        obj.setPlayList(SB2,probeOutputB.ch(2));
    end

    obj.clearRecList % clear out whatever was used previously 
    if ~isempty(probeInputA)
        probeInputA.label = [probeInputA.label,'_LGF'];
        obj.setRecList(probeInputA.ch,probeInputA.label,probeInputA.micSens,probeInputA.gain);    
    end
    if ~isempty(probeInputB)
        probeInputB.label = [probeInputB.label,'_LGF'];
        obj.setRecList(probeInputB.ch,probeInputB.label,probeInputB.micSens,probeInputB.gain);    
    end

    nSweeps = nSweeps(1); % all frequencies have the same number of sweeps
    nReps = nSweeps * size(Y1,2); % number of reps (for ARLas, since stim folded into a matrix)
    
    obj.setNReps(nReps); % number of times to play stimulus
    obj.setFilter(1); % note: this is a highpass filter with a 75 Hz cutoff frequency.
    
    % Populate userInfo fields --------------------------------------------
    obj.objPlayrec.userInfo = []; % clear out previous, non-structure array info
    obj.objPlayrec.userInfo.dpoae_version = V; % this version of the test
    obj.objPlayrec.userInfo.fs = fs; % sampling rate (Hz)
    obj.objPlayrec.userInfo.LC = LC; % load calibration structure
    obj.objPlayrec.userInfo.testEar = LC.testEar; % which ear is ACTUALLY being tested (has probe insertion), regardless of which probe is being used
    obj.objPlayrec.userInfo.sweepLength = sweepLength; % length of each sweep (s)
    obj.objPlayrec.userInfo.sweepLength = sweepN; % number of samples in each sweep (samples)
    obj.objPlayrec.userInfo.sweepRate = sweepRate; % dB / sec
    obj.objPlayrec.userInfo.nSweeps = nSweeps; % number of sweeps to collect
    obj.objPlayrec.userInfo.nReps = nReps; % number of reps (for ARLas, since stim folded into a matrix)
    obj.objPlayrec.userInfo.Time = Time; % time vector
    obj.objPlayrec.userInfo.fmin = fmin; % min frequency being tested
    obj.objPlayrec.userInfo.fmax = fmax; % max frequency being tested
    obj.objPlayrec.userInfo.targetL1 = targetL1; % target levels
    obj.objPlayrec.userInfo.targetL2 = targetL2; % target levels
    obj.objPlayrec.userInfo.L2min = L2min; % 0 min L2 (dB FPL) 
    obj.objPlayrec.userInfo.L2max = L2max; % 70 max L2 (dB FPL)
    obj.objPlayrec.userInfo.f1 = f1; % f1 primary frequencies
    obj.objPlayrec.userInfo.f2 = f2; % f2 primary frequencies
    obj.objPlayrec.userInfo.fdp = fdp; % dpoae frequencies
    obj.objPlayrec.userInfo.fRatio = fRatio; % f2 / f1 ratios
    obj.objPlayrec.userInfo.L1 = L1; % L1 vector of levels across sweep
    obj.objPlayrec.userInfo.L2 = L2; % L2 vector of levels across sweep
    obj.objPlayrec.userInfo.F1 = F1; % f1 vector of frequencies across sweep (these are constant in this program)
    obj.objPlayrec.userInfo.F2 = F2; % f2 vector of frequencies across sweep (these are constant in this program)
    obj.objPlayrec.userInfo.Fdp = Fdp; % fdp vector of frequencies across sweep (these are constant in this program)


    obj.objPlayrec.run % playback and record -----
    if obj.killRun
       return
    end

    % extract the recordings -----
    if ~isempty(probeInputA)
        [headerA,DataA] = obj.retrieveData(['Ch',num2str(probeInputA.ch)]); % retrive recorded data
        if LC.doMicCorrection == 1
            DataA = applyMicCorrection(DataA,LC.micCorrectionA);
        end
    end
    if ~isempty(probeInputB)
        [headerB,DataB] = obj.retrieveData(['Ch',num2str(probeInputB.ch)]); % retrive recorded data
        if LC.doMicCorrection == 1
            DataB = applyMicCorrection(DataB,LC.micCorrectionB);
        end
    end

    % perform post-hoc load calibration -----
    doProbeFitCheck = 1;
    if ~isgraphics(iscCheckHandleA) % check to make sure window hasn't been shut down
        iscCheckHandleA = []; % figures not yet made so initialize handles to zero
    end
    if ~isgraphics(iscCheckHandleB)
        iscCheckHandleB = []; 
    end
    [LC2,iscCheckHandleA,iscCheckHandleB] = loadCalibration(obj,doProbeFitCheck,LC,iscCheckHandleA,iscCheckHandleB);    
    % Where are calibration checks saved? They are saved in ...\Data\subjectID\loadCals\fileName
    % here we don't do anything with LC2 at this point, although it will
    % become the most recent load calibration if accessed using getMostRecentLoadCalibration.m

        
    % Analyze data ----------------------------------------------------
    disp('----- Analyzing data -----')

    % note: saving done by the analyzing program?
    % make analyze all levels
    % skip bootstrapping for now?
    
    try
    if ~isempty(probeInputA)
        DPOAE_A.subjID = obj.subjectID;
        DPOAE_A.testType = 'LGF';
        DPOAE_A.nSweeps = nSweeps;
        DPOAE_A.L2 = L2;
        DPOAE_A.f2 = f2;
        DPOAE_A.fdp = fdp;
        DPOAE_A.doEPL = doEPL;
        [DPOAE_A,h1A,h2A] = analyzeData(DataA,headerA,LC,DPOAE_A);
    else
        DPOAE_A = [];
    end
    catch ME
        keyboard
    end
    try
    if ~isempty(probeInputB)
        DPOAE_B.subjID = obj.subjectID;
        DPOAE_B.testType = 'LGF';
        DPOAE_B.nSweeps = nSweeps;
        DPOAE_B.L2 = L2;
        DPOAE_B.f2 = f2;
        DPOAE_B.fdp = fdp;
        DPOAE_B.doEPL = doEPL;
        [DPOAE_B,h1B,h2B] = analyzeData(DataB,headerB,LC,DPOAE_B);
    else
        DPOAE_B = [];
    end
    catch ME
        keyboard
    end


    % save analyses, including:
    % mat files
    % excel files
    % save to central location 

    saveProbeFitPlots(obj,LC,probeInputA,probeInputB,iscCheckHandleA,iscCheckHandleB)
    if ~isempty(probeInputA)
        saveAnalysisPlots(obj,LC,h1A,h2A)
    end
    if ~isempty(probeInputB)
        saveAnalysisPlots(obj,LC,h1B,h2B)
    end
    saveMatFiles(obj,DPOAE_A,DPOAE_B)

    % save excel sheets
    saveExcelFiles(obj,LC,DPOAE_A,DPOAE_B)

    disp('----- Finished DPOAE experiment: LGF -----')

end

% INTERNAL FUNCTIONS ------------------------------------------------------
function [Y1,Y2,F1,F2,Fdp,Time,Y1c,Y2c,Ydp,Ydpc,L1,L2] = ARLas_dpoaeStim_sweepL2Continuous(fs,f1,f2,fdp,L2min,L2max,sweepN)
    foldLength = 0.1; % fold into a matrix with columns this length (sec)
    foldN = round(foldLength * fs); % number of samples in each column
    
    %dbPerSec = sweepRate;
    %nDB = L2max - L2min; % number of dB to sweep over
    %len = nDB / dbPerSec; % length of total stimulus (sec)
    %N = round(fs*len); % number of samples in total stimulus
    
    dB = linspace(L2min,L2max,sweepN)'; % dB spacing, linear sweep
    L2 = 10.^(dB/20)*0.00002; % target L2 put into linear units.
    L1 = ones(sweepN,1); 

    % add padding for ramps on and off
    rampLen1 = 0.003; % 3 ms ramp on and off
    rampN1 = round(fs * rampLen1);
    rampLen2 = 0.003; 
    rampN2 = round(fs * rampLen2);
    pad1 = ones(rampN1,1)*min(L2);
    pad2 = ones(rampN2,1)*max(L2);
    L2 = [pad1;L2;pad2];
    pad3 = ones(rampN1,1);
    L1 = [pad3;L1;pad3];

    f1 = ones(size(L2))*f1; % get f1 vector of constant frequencies
    f2 = ones(size(L2))*f2;
    fdp = ones(size(L2))*fdp;
   
    N = length(L2); % recalculate new number of samples (after padding added)
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
        phaseChange1 = 2*pi* (deltaT * f1(ii));
        phi1 = phi1 + phaseChange1;

        phaseChange2 = 2*pi* (deltaT * f2(ii));
        phi2 = phi2 + phaseChange2;

        phaseChangeDP = 2*pi* (deltaT * fdp(ii));
        phiDP = phiDP + phaseChangeDP;

        y1(ii,1) = L1(ii).*sin(phi1);
        y2(ii,1) = L2(ii).*sin(phi2);
        ydp(ii,1) = sin(phiDP);

        y1c(ii,1) = L1(ii).*cos(phi1);
        y2c(ii,1) = L2(ii).*cos(phi2);    
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
    L1(1:rampN1) = L1(1:rampN1) .* h; % stimulus magnitude 
    L2(1:rampN1) = L2(1:rampN1) .* h; % stimulus magnitude
    
    h = hann(rampN2*2);
    h = h(1:rampN2);
    h = flipud(h);
    y1(end-rampN2+1:end) = y1(end-rampN2+1:end) .* h;
    y2(end-rampN2+1:end) = y2(end-rampN2+1:end) .* h;
    ydp(end-rampN2+1:end) = ydp(end-rampN2+1:end) .* h;
    y1c(end-rampN2+1:end) = y1c(end-rampN2+1:end) .* h;
    y2c(end-rampN2+1:end) = y2c(end-rampN2+1:end) .* h;
    ydpc(end-rampN2+1:end) = ydpc(end-rampN2+1:end) .* h;
    L1(end-rampN2+1:end) = L1(end-rampN2+1:end) .* h;
    L2(end-rampN2+1:end) = L2(end-rampN2+1:end) .* h;

    % fold up the stimuli for efficient presentation
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
        L1 = [L1;pad];
        L2 = [L2;pad];
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
function [DPOAE,h1,h2] = analyzeData(Data,header,LC,DPOAE)
    nSweeps = DPOAE.nSweeps;
    L2 = DPOAE.L2;
    fdp = DPOAE.fdp;
    fs = header.fs;

    % reshape and bandpass filter from 0.4-20 kHz
    Rows = length(Data(:))/nSweeps;
    Data = reshape(Data,Rows,nSweeps);
    b = getbpf; % fir filter coefficients                 
    [Rows,Cols] = size(Data);
    Data = Data(:);
    Data = fastFilter(b,Data);
    Data = reshape(Data,Rows,Cols);

    % perform swept LSF
    Lmin = header.userInfo.L2min;
    Lmax = header.userInfo.L2max;
    step = 1;
    x = 20*log10(L2/.00002);
    [signal,nf,snr,targets,zbar] = sweptLSF(x,Data,fdp,Lmin,Lmax,step,fs);
    
    % apply EPL correction ------------------------------------------------
    if DPOAE.doEPL == 1
        k1 = LC.iscSA1.eplCorrection;
        k2 = LC.iscSA2.eplCorrection;
        ff = LC.iscSA1.freq;
        k = abs((k1 + k2)/2); % take the mean from each channel as the correction
        k = meanSmoother(k,50); % smooth the epl correction
        [~,indx] = min(abs(ff - fdp)); % find the closest correction value
        k = k(indx);
        signal = 10.^(signal/20)*0.00002; % put signal magnitude back into Pascals
        signal = signal * abs(k); % multiply the signal magnitude by the epl magnitude correction
        signal = 20*log10(signal/.00002); % put signal back into dB (now EPL)
        % do the same thing to the noise floor
        nf = 10.^(nf/20)*0.00002; 
        nf = nf * abs(k); 
        nf = 20*log10(nf/.00002);
    end
    % end EPL correction --------------------------------------------------

    % fit only the low- to mid-level portion
    [~,indxMax] = max(signal);
    L2max = targets(indxMax);
    L2cut = L2max - 3; 
    if L2cut < 40
        L2cut = 40;
    end
    [~,indxCut] = min(abs(L2cut - targets));
    % perform Rician-based fit
    [jar,h1,h2] = ricianLGF_fit2d(signal(1:indxCut),nf(1:indxCut),targets(1:indxCut),nSweeps);
    
    DPOAE.jar = jar;
    DPOAE.L2max = L2max;
    DPOAE.L2cut = L2cut;
    DPOAE.indxCut = indxCut;
    DPOAE.signal = signal;
    DPOAE.targets = targets;
    DPOAE.nf = nf;
    DPOAE.Enf = mean(nf); % expected noise floor (mean value)
    DPOAE.LC = LC;
    DPOAE.nSweeps = nSweeps;
end
function [] = saveProbeFitPlots(obj,LC,probeInputA,probeInputB,iscCheckHandleA,iscCheckHandleB)
    disp('----- Saving Probe Fit Figures -----')
    tempPath = LC.tempResultsPath;
    try
        if ~isempty(probeInputA)
            savePath = [obj.savingGrace,obj.experimentID,'_analysis',obj.sep];
            if exist(savePath,'dir') == 0 % if path does not exist
                success = mkdir(savePath);
                addpath(savePath)
                if ~success
                    warning('Unable to create new Experiment Directory: data save path')
                end
            end 
            figureFileName = [obj.subjectID,'_inSituCal_','DPOAE_A_LGF','.fig'];
            %figureFileName = ARLas_saveName(savePath,figureFileName);
            savefig(iscCheckHandleA,[savePath,figureFileName])
            figureFileName = [obj.subjectID,'_inSituCal_','DPOAE_A_LGF','.bmp'];
            %figureFileName = ARLas_saveName(savePath,figureFileName);    
            saveas(iscCheckHandleA,[savePath,figureFileName])
            saveas(iscCheckHandleA,[tempPath,figureFileName])
        end
    catch
        disp('Warning: One or more figures for DPOAE_A not saved!')        
    end

    try
        if ~isempty(probeInputB)
            savePath = [obj.savingGrace,obj.experimentID,'_analysis',obj.sep];
            if exist(savePath,'dir') == 0 % if path does not exist
                success = mkdir(savePath);
                addpath(savePath)
                if ~success
                    warning('Unable to create new Experiment Directory: data save path')
                end
            end 
            figureFileName = [obj.subjectID,'_inSituCal_','DPOAE_B_LGF','.fig'];
            %figureFileName = ARLas_saveName(savePath,figureFileName);
            savefig(iscCheckHandleB,[savePath,figureFileName])
            figureFileName = [obj.subjectID,'_inSituCal_','DPOAE_B_LGF','.bmp'];
            %figureFileName = ARLas_saveName(savePath,figureFileName);    
            saveas(iscCheckHandleB,[savePath,figureFileName])
            saveas(iscCheckHandleB,[tempPath,figureFileName])
        end
    catch
       disp('Warning: One or more figures for DPOAE_B not saved!')
    end
end
function [] = saveAnalysisPlots(obj,LC,h1,h2)
    disp('----- Saving LGF Figures -----')
    tempPath = LC.tempResultsPath;
    try
        if ~isempty(h1)
            savePath = [obj.savingGrace,obj.experimentID,'_analysis',obj.sep];
            if exist(savePath,'dir') == 0 % if path does not exist
                success = mkdir(savePath);
                addpath(savePath)
                if ~success
                    warning('Unable to create new Experiment Directory: data save path')
                end
            end 
            figureFileName = [obj.subjectID,'DPOAE_LGF','.fig'];
            figureFileName = ARLas_saveName(savePath,figureFileName);
            savefig(h1,[savePath,figureFileName])
            figureFileName = [obj.subjectID,'DPOAE_LGF','.bmp'];
            figureFileName = ARLas_saveName(savePath,figureFileName);    
            saveas(h1,[savePath,figureFileName])
            saveas(h1,[tempPath,figureFileName])
        end
    catch
        disp('Warning: One or more figures for DPOAE_LGF not saved!')        
    end

    try
        if ~isempty(h2)
            savePath = [obj.savingGrace,obj.experimentID,'_analysis',obj.sep];
            if exist(savePath,'dir') == 0 % if path does not exist
                success = mkdir(savePath);
                addpath(savePath)
                if ~success
                    warning('Unable to create new Experiment Directory: data save path')
                end
            end 
            figureFileName = [obj.subjectID,'DPOAE_MLE','.fig'];
            figureFileName = ARLas_saveName(savePath,figureFileName);
            savefig(h2,[savePath,figureFileName])
            figureFileName = [obj.subjectID,'DPOAE_MLE','.bmp'];
            figureFileName = ARLas_saveName(savePath,figureFileName);    
            saveas(h2,[savePath,figureFileName])
            saveas(h2,[tempPath,figureFileName])
        end
    catch
       disp('Warning: One or more figures for DPOAE_MLE not saved!')
    end
end
function [] = saveMatFiles(obj,DPOAE_A,DPOAE_B)
    disp('----- Saving MAT files -----')
    % mat files don't get saved to temp location
    try
        if ~isempty(DPOAE_A)
            savePath = [obj.savingGrace,obj.experimentID,'_analysis',obj.sep];
            if exist(savePath,'dir') == 0 % if path does not exist
                success = mkdir(savePath);
                addpath(savePath)
                if ~success
                    warning('Unable to create new Experiment Directory: data save path')
                end
            end 
            fileName = [obj.subjectID,'DPOAE_A_LGF','.mat'];
            figureFileName = ARLas_saveName(savePath,fileName);
            save([savePath,fileName],'DPOAE_A')
        end
    catch
        disp('Warning: Mat file for DPOAE_A_LGF not saved!')        
    end
    try
        if ~isempty(DPOAE_B)
            savePath = [obj.savingGrace,obj.experimentID,'_analysis',obj.sep];
            if exist(savePath,'dir') == 0 % if path does not exist
                success = mkdir(savePath);
                addpath(savePath)
                if ~success
                    warning('Unable to create new Experiment Directory: data save path')
                end
            end 
            fileName = [obj.subjectID,'DPOAE_B_LGF','.mat'];
            figureFileName = ARLas_saveName(savePath,fileName);
            save([savePath,fileName],'DPOAE_B')
        end
    catch
        disp('Warning: Mat file for DPOAE_B_LGF not saved!')        
    end
end
function [] = saveExcelFiles(obj,LC,DPOAE_A,DPOAE_B)
    % Write analysis to spreadsheet
    disp('----- Writing Excel files -----')
    warning off
    tempPath = LC.tempResultsPath;
    try
        if ~isempty(DPOAE_A)
            savePath = [obj.savingGrace,obj.experimentID,'_analysis',obj.sep];
            if exist(savePath,'dir') == 0 % if path does not exist
                success = mkdir(savePath);
                addpath(savePath)
                if ~success
                    warning('Unable to create new Experiment Directory: data save path')
                end
            end 

            % ------
            fileName = [obj.subjectID,'DPOAE_A_LGF','.xls'];
            saveFileName = [fileName,'x'];
            saveFileNameCSV = [obj.subjectID,'DPOAE_A_LGF','.csv'];
            sheetName = [num2str(round(DPOAE_A.f2)),' Hz'];
            writematrix('F2_dB_FPL',[savePath,saveFileName],'Sheet',sheetName,'Range','A1');
            writematrix(DPOAE_A.targets,[savePath,saveFileName],'Sheet',sheetName,'Range','A2');
            writematrix('Ldp_dB_EPL',[savePath,saveFileName],'Sheet',sheetName,'Range','B1');
            writematrix(DPOAE_A.signal,[savePath,saveFileName],'Sheet',sheetName,'Range','B2');
            writematrix('Ln_dB_EPL',[savePath,saveFileName],'Sheet',sheetName,'Range','C1');
            writematrix(DPOAE_A.nf,[savePath,saveFileName],'Sheet',sheetName,'Range','C2');
            writematrix('X_dB_FPL',[savePath,saveFileName],'Sheet',sheetName,'Range','D1');
            writematrix(DPOAE_A.jar.X,[savePath,saveFileName],'Sheet',sheetName,'Range','D2');
            writematrix('Z_dB_EPL',[savePath,saveFileName],'Sheet',sheetName,'Range','E1');
            writematrix(DPOAE_A.jar.Z,[savePath,saveFileName],'Sheet',sheetName,'Range','E2');
            writematrix('Sz_dB_EPL',[savePath,saveFileName],'Sheet',sheetName,'Range','F1');
            writematrix(DPOAE_A.jar.Sz,[savePath,saveFileName],'Sheet',sheetName,'Range','F2');
            writematrix('Yhat_dB_EPL',[savePath,saveFileName],'Sheet',sheetName,'Range','G1');
            writematrix(DPOAE_A.jar.Yhat_mle,[savePath,saveFileName],'Sheet',sheetName,'Range','G2');
            writematrix('Slope',[savePath,saveFileName],'Sheet',sheetName,'Range','H1');
            writematrix(DPOAE_A.jar.Slope_mle,[savePath,saveFileName],'Sheet',sheetName,'Range','H2');
            writematrix('MLE',[savePath,saveFileName],'Sheet',sheetName,'Range','I1');
            writematrix(DPOAE_A.jar.mle(:),[savePath,saveFileName],'Sheet',sheetName,'Range','I2');
            % -----
            %Read the .xlsx file
            data = readtable([savePath,saveFileName]);
            % Save the data as a .csv file
            writetable(data,[tempPath,saveFileNameCSV]);

        end
    catch
        disp('Warning: Mat file for DPOAE_A XLSX not saved!')        
    end
    try
        if ~isempty(DPOAE_B)
            savePath = [obj.savingGrace,obj.experimentID,'_analysis',obj.sep];
            if exist(savePath,'dir') == 0 % if path does not exist
                success = mkdir(savePath);
                addpath(savePath)
                if ~success
                    warning('Unable to create new Experiment Directory: data save path')
                end
            end 

            % -----
            fileName = [obj.subjectID,'DPOAE_B_LGF','.xls'];
            saveFileName = [fileName,'x'];
            saveFileNameCSV = [obj.subjectID,'DPOAE_B_LGF','.csv'];
            sheetName = [num2str(round(DPOAE_B.f2)),' Hz'];
            writematrix('F2_dB_FPL',[savePath,saveFileName],'Sheet',sheetName,'Range','A1');
            writematrix(DPOAE_B.targets,[savePath,saveFileName],'Sheet',sheetName,'Range','A2');
            writematrix('Ldp_dB_EPL',[savePath,saveFileName],'Sheet',sheetName,'Range','B1');
            writematrix(DPOAE_B.signal,[savePath,saveFileName],'Sheet',sheetName,'Range','B2');
            writematrix('Ln_dB_EPL',[savePath,saveFileName],'Sheet',sheetName,'Range','C1');
            writematrix(DPOAE_B.nf,[savePath,saveFileName],'Sheet',sheetName,'Range','C2');
            writematrix('X_dB_FPL',[savePath,saveFileName],'Sheet',sheetName,'Range','D1');
            writematrix(DPOAE_B.jar.X,[savePath,saveFileName],'Sheet',sheetName,'Range','D2');
            writematrix('Z_dB_EPL',[savePath,saveFileName],'Sheet',sheetName,'Range','E1');
            writematrix(DPOAE_B.jar.Z,[savePath,saveFileName],'Sheet',sheetName,'Range','E2');
            writematrix('Sz_dB_EPL',[savePath,saveFileName],'Sheet',sheetName,'Range','F1');
            writematrix(DPOAE_B.jar.Sz,[savePath,saveFileName],'Sheet',sheetName,'Range','F2');
            writematrix('Yhat_dB_EPL',[savePath,saveFileName],'Sheet',sheetName,'Range','G1');
            writematrix(DPOAE_B.jar.Yhat_mle,[savePath,saveFileName],'Sheet',sheetName,'Range','G2');
            writematrix('Slope',[savePath,saveFileName],'Sheet',sheetName,'Range','H1');
            writematrix(DPOAE_B.jar.Slope_mle,[savePath,saveFileName],'Sheet',sheetName,'Range','H2');
            writematrix('MLE',[savePath,saveFileName],'Sheet',sheetName,'Range','I1');
            writematrix(DPOAE_B.jar.mle(:),[savePath,saveFileName],'Sheet',sheetName,'Range','I2');
            % -----
            %Read the .xlsx file
            data = readtable([savePath,saveFileName]);
            % Save the data as a .csv file
            writetable(data,[tempPath,saveFileNameCSV]);
            
        end
    catch
        disp('Warning: Mat file for DPOAE_B XLSX not saved!')        
    end
    warning on
end
function b = getbpf()
    Fs = 96000;
    Fstop1 = 100;             % First Stopband Frequency 
    Fpass1 = 400;            % First Passband Frequency -->for a low freq of 750, the dp is 450 Hz
    Fpass2 = 20000;           % Second Passband Frequency
    Fstop2 = 24000;           % Second Stopband Frequency
    Dstop1 = 0.001;           % First Stopband Attenuation
    Dpass  = 0.057501127785;  % Passband Ripple
    Dstop2 = 0.0001;          % Second Stopband Attenuation
    flag   = 'scale';         % Sampling Flag
    % Calculate the order from the parameters using KAISERORD.
    [N,Wn,BETA,TYPE] = kaiserord([Fstop1 Fpass1 Fpass2 Fstop2]/(Fs/2), [0 ...
                                 1 0], [Dstop1 Dpass Dstop2]);
    if mod(N,2)~= 0
        N = N + 1;
    end
    b  = fir1(N, Wn, TYPE, kaiser(N+1, BETA), flag)';
end


% OLD CODE ----------------------------------------------------------------
