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
% 
% 
% 
% load('Ch1_OP24A_DPG_0001.mat')
% 
% DPOAE.subjID = 'shawn';
% DPOAE.testType = 'GRM';
% DPOAE.f2min = header.userInfo.fmin;
% DPOAE.f2max = header.userInfo.fmax;
% DPOAE.nSweeps = header.userInfo.nSweeps;
% DPOAE.F1 = header.userInfo.F1;
% DPOAE.F2 = header.userInfo.F2;
% DPOAE.Fdp = header.userInfo.Fdp;
% DPOAE.L1 = header.userInfo.L1;
% DPOAE.L2 = header.userInfo.L2;
% DPOAE.fs = header.userInfo.fs;
% DPOAE.doEPL = 0;
% LC = header.userInfo.LC;
% 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    obj = varargin{1}; % get the arlas object
    V = '23MAY2025'; % this is the current version number of this program

%--------------------------------------------------------------------------
%------ USER MODIFIABLE PARAMETERS ----------------------------------------
    nSweeps = 24; %24; % 24 % Number of sweeps recorded at each level.
                    % Note: At sweep rate of 0.5 and 1-16 kHz, each sweep takes 8 sec.
                    %       Recommended nSweeps is 24, which should get noise floor to -10 dB SPL. This will take 24*8 = 192 sec = 3.2 min per
    
    % specify nominal stimulus parameters:
    %   Note: these will be modified slightly at low and high frequencies
    targetL1 = 65; %[65 65 65 65 65]; % target levels for primaries (dB FPL); Can be scalar or vector for more than one level
    targetL2 = 55; %[25 35 45 55 65]; % length of targetL2 must be same size as targetL1  

    fmin = 375; % min f2 value is 750 to ensure 1000 Hz
    fmax = 18000; % max f2 value is 18000 to ensure 16000 Hz (14000 is a 1/2 octave above 10 kHz)
    fRatio = 1.22; % nominal ratio of f2/f1. Sammi=1 will alter this
    sweepRate = 0.5; % sweep rate in octaves/s; should be 0.5.

    doEPL = 0; % turn on and off EPL correction 
    Sammi = 0; % turn on and off Sammi Ginter correction
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


    disp(' ')
    disp('----- Starting DPOAE experiment: DPG -----')

    % COLLECT THE DATA --------------------------------------------------------
    disp('----- Collecting data -----')
    fs = obj.fs; % get the system sampling rate
    nLevels = length(targetL1); % number of stim
    adjust = 0;
        
    for ii=1:nLevels % loop over number of stimulus levels
        if nLevels > 1
            disp(['Presenting stimulus level ',num2str(ii),' of ',num2sr(nLevels)])
        end

        %[f1,fdp,fRatio] = getRatios(L1,L2,f2);
        [Y1,Y2,F1,F2,Fdp,Time,Y1c,Y2c,Ydp,Ydpc,L1,L2,fRatio] = getStimuli(fs,fmin,fmax,sweepRate,fRatio,Sammi); % get the raw stimuli

        % Apply the in-situ calibration to the stimuli ----------------
        % Apply the in-situ calibration to the stimuli ------------------------
        iscSA1 = LC.iscSA1;
        iscSA2 = LC.iscSA2;
        iscSB1 = LC.iscSB1;
        iscSB2 = LC.iscSB2;
        probeInputA = LC.probeInputA;
        probeOutputA = LC.probeOutputA;
        probeInputB = LC.probeInputB;
        probeOutputB = LC.probeOutputB;

        if ~isempty(iscSA1)
            [rows,cols] = size(Y1);
            [SA1,clicksImpulse,errors1] = ARLas_applyISC(iscSA1,targetL1(ii)+adjust,'fpl',F1(:),Y1(:));
            SA1 = reshape(SA1,rows,cols);
            [SA2,clicksImpulse,errors1] = ARLas_applyISC(iscSA2,targetL2(ii)+adjust,'fpl',F2(:),Y2(:));
            SA2 = reshape(SA2,rows,cols);
        else
            SA1 = [];
            SA2 = [];
        end
        if ~isempty(iscSB1)
            [rows,cols] = size(Y1);
            [SB1,clicksImpulse,errors1] = ARLas_applyISC(iscSB1,targetL1(ii)+adjust,'fpl',F1(:),Y1(:));
            SB1 = reshape(SB1,rows,cols);
            [SB2,clicksImpulse,errors1] = ARLas_applyISC(iscSB2,targetL2(ii)+adjust,'fpl',F2(:),Y2(:));
            SB2 = reshape(SB2,rows,cols);
        else
            SB1 = [];
            SB2 = [];
        end

        if Sammi == 1
            D1 = L1(:) - targetL1(ii);
            D2 = L2(:) - targetL2(ii);
            zindx = find(Y1(:)==0);
            D1(zindx) = 0;
            D2(zindx) = 0;
            D1 = 10.^(D1/20);
            D2 = 10.^(D2/20);
            SA1 = SA1(:).*D1;
            SA2 = SA2(:).*D2;
            SA1 = reshape(SA1,rows,cols);
            SA2 = reshape(SA2,rows,cols);
        end

        % tell ARLas what to record and what to play
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
            if ii==1
                probeInputA.label = [probeInputA.label,'_DPG'];
            end
            obj.setRecList(probeInputA.ch,probeInputA.label,probeInputA.micSens,probeInputA.gain);    
        end
        if ~isempty(probeInputB)
            if ii==1
                probeInputB.label = [probeInputB.label,'_DPG'];
            end
            obj.setRecList(probeInputB.ch,probeInputB.label,probeInputB.micSens,probeInputB.gain);    
        end

        nReps = nSweeps(ii) * size(Y1,2); % number of reps (for ARLas, since stim folded into a matrix)
        obj.setNReps(nReps); % number of times to play stimulus
        obj.setFilter(1); % note: this is a highpass filter with a 75 Hz cutoff frequency.
        
        obj.objPlayrec.userInfo = []; % clear out previous, non-structure array info
        obj.objPlayrec.userInfo.dpoae_version = V;
        
        obj.objPlayrec.userInfo.LC = LC; % load calibration structure
        obj.objPlayrec.userInfo.testEar = LC.testEar; % which ear is ACTUALLY being tested (has probe insertion), regardless of which probe is being used

        obj.objPlayrec.userInfo.fmin = fmin;
        obj.objPlayrec.userInfo.fmax = fmax;
        obj.objPlayrec.userInfo.fs = fs;
        obj.objPlayrec.userInfo.sweepRate = sweepRate;
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
        obj.objPlayrec.userInfo.L1 = L1;
        obj.objPlayrec.userInfo.L2 = L2;
        obj.objPlayrec.userInfo.fRatio = fRatio;

        disp(['----- Recording data: ',num2str(ii),' of ',num2str(nLevels),' -----'])
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

        % Analyze data ----------------------------------------------------
        disp('----- Analyzing data -----')
        try
        if ~isempty(probeInputA)
            DPOAE_A.subjID = obj.subjectID;
            DPOAE_A.testType = 'DPG';
            DPOAE_A.f2min = fmin;
            DPOAE_A.f2max = fmax;
            DPOAE_A.nSweeps = nSweeps;
            DPOAE_A.F1 = F1;
            DPOAE_A.F2 = F2;
            DPOAE_A.Fdp = Fdp;
            DPOAE_A.L1 = headerA.userInfo.targetL1;
            DPOAE_A.L2 = headerA.userInfo.targetL2;
            DPOAE_A.fs = fs;
            DPOAE_A.doEPL = doEPL;
            [DPOAE_A,h1A] = analyzeData(DataA,headerA,LC,DPOAE_A);
        else
            DPOAE_A = [];
        end
        catch ME
            keyboard
        end
        try
        if ~isempty(probeInputB)
            DPOAE_B.subjID = obj.subjectID;
            DPOAE_B.testType = 'DPG';
            DPOAE_B.f2min = fmin;
            DPOAE_B.f2max = fmax;
            DPOAE_B.nSweeps = nSweeps;
            DPOAE_B.F1 = F1;
            DPOAE_B.F2 = F2;
            DPOAE_B.Fdp = Fdp;
            DPOAE_B.L1 = headerB.userInfo.targetL1;
            DPOAE_B.L2 = headerB.userInfo.targetL1;
            DPOAE_B.fs = fs;
            DPOAE_B.doEPL = doEPL;
            [DPOAE_B,h1B] = analyzeData(DataB,headerB,LC,DPOAE_B);
        else
            DPOAE_B = [];
        end
        catch ME
            keyboard
        end
    
        saveProbeFitPlots(obj,LC,probeInputA,probeInputB,iscCheckHandleA,iscCheckHandleB)
        if ~isempty(probeInputA)
            saveAnalysisPlots(obj,LC,h1A)
        end
        if ~isempty(probeInputB)
            saveAnalysisPlots(obj,LC,h1B)
        end
        saveMatFiles(obj,DPOAE_A,DPOAE_B)
    
        % save excel sheets
        saveExcelFiles(obj,LC,DPOAE_A,DPOAE_B)

    end
disp('----- Finished DPOAE experiment: DPG -----')
end

% INTERNAL FUNCTIONS ------------------------------------------------------

function [Y1,Y2,F1,F2,Fdp,Time,Y1c,Y2c,Ydp,Ydpc,L1,L2,fRatio] = getStimuli(fs,fmin,fmax,sweepRate,fRatio,Sammi)
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
    if Sammi==1
        %[f1,fdp,fRatio] = getRatios(targetL1,targetL2,f2);
        [f1,fdp,L1,L2,fRatio] = getSammiParameters(f2);
    else
        f1 = f2 ./ fRatio; % f1 primary frequencies
        fdp = 2*f1 - f2; % cubic distortion dpoae frequencies
    end

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

     if Sammi==0
         L1 = ones(size(y1));
         L2 = ones(size(y1));
         fRatio = ones(size(y1))*fRatio(1);
     end

    h = hann(rampN1*2); % ramp on and off the stimuli
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
    L1 = reshape(L1,foldN,nFolds);
    L2 = reshape(L2,foldN,nFolds);
    Fdp = reshape(fdp,foldN,nFolds);
    NN = length(y1);
    time = (0:1:NN-1)'/fs;
    Time = reshape(time,foldN,nFolds);
end
function [F1,Fdp,L1,L2,Ratio] = getSammiParameters(F2)
    % Values taken from Samantha Ginter dissertaion with Sumit Dhar, 2021
    Sf2 = [0.5,0.75,1,1.5,2,3,4,6,8,10,11.2,12.5,14,15,16]'*1000; % Sammi frequencies
    SL1 = [75,65,65,65,65,65,55,55,65,65,65,75,75,75,75]'; % Sammi L1
    SL2 = [75,55,55,55,55,55,40,40,55,55,55,75,75,75,75]'; % Sammi L2
    SRat= [1.24,1.22,1.22,1.18,1.2,1.2,1.18,1.18,1.18,1.16,1.16,1.16,1.14,1.14,1.14]'; % Sammi ratio
    L1 = interp1(Sf2,SL1,F2,'pchip');
    L2 = interp1(Sf2,SL2,F2,'pchip');
    Ratio = interp1(Sf2,SRat,F2,'pchip');
    F1 = F2 ./ Ratio;
    Fdp = 2*F1 - F2;
end
function [DPOAE,h1] = analyzeData(Data,header,LC,DPOAE)
    nSweeps = DPOAE.nSweeps;
    %L2 = DPOAE.L2;
    fdp = DPOAE.Fdp;
    fs = header.fs;

    % reshape and bandpass filter from 0.4-20 kHz
    Rows = length(Data(:))/nSweeps;
    Data = reshape(Data,Rows,nSweeps);
    b = getbpf; % fir filter coefficients                 
    [Rows,Cols] = size(Data);
    Data = Data(:);
    Data = fastFilter(b,Data);
    Data = reshape(Data,Rows,Cols);
    
    Y1 = header.userInfo.Y1;
    Y1c = header.userInfo.Y1c;
    Y2 = header.userInfo.Y2;
    Y2c = header.userInfo.Y2c;
    Ydp = header.userInfo.Ydp;
    Ydpc = header.userInfo.Ydpc;
    Y1 = reshape(Y1,Rows,1);
    Y1c = reshape(Y1c,Rows,1);
    Y2 = reshape(Y2,Rows,1);
    Y2c = reshape(Y2c,Rows,1);
    Ydp = reshape(Ydp,Rows,1);
    Ydpc = reshape(Ydpc,Rows,1);
    F1 = DPOAE.F1;
    F2 = DPOAE.F2;
    Fdp = DPOAE.Fdp;
    F1 = reshape(F1,Rows,1);
    F2 = reshape(F2,Rows,1);
    Fdp = reshape(Fdp,Rows,1);


    % -------------------------------------------
    frameLen = 0.5; % analysis frame length (s)
    frameN = round(fs*frameLen); % number of samples in each frame
    downsampledN = 100; % downsample the results to this number of samples
    
    fmin = DPOAE.f2min;
    fmax = DPOAE.f2max;
    step = (fmax - fmin) / (downsampledN-1); % step size for downsampled result
    
    targets = (fmin:step:fmax)';
    for ii=1:length(targets) % find sample location of each target
        [~,indices(ii,1)] = min(abs(targets(ii)-F2));
    end
    nIterations = length(indices);
    DPOAE.pcOverlap = median(diff(indices)) / frameN; % percent overlap between each successive analysis frame

    % ----------------------------------------------
    % create solution matrix for LSF
    W = blackman(frameN); %ones(frameN,1); %hann(frameN);
    frameN2 = ceil(frameN/2); % half-frame size (number of samples)
    nSamples = size(Data,1); % total number of samples in each sweep
    % bad results happen when the frames are not full. Remove partially
    % filled frames
    [~,lossIndx] = min(abs(indices - frameN2));
    indices = indices(lossIndx:end-lossIndx);
    
    DPOAE.targetF2 = F2(indices);
    DPOAE.targetF1 = F1(indices);
    nIterations = length(indices);

    hitCounter = 0;
    for ii=1:nIterations % loop across frames
        if mod(ii,10)==0
            disp(['Analyzing frame ',num2str(ii),' of ',num2str(nIterations)])
        end
        start = indices(ii)-frameN2+1; % starting sample of current frame
        finish = indices(ii)+frameN2; % ending sample of current frame
        if start < 0 % handle frame run on and run off of the data matrix
            finish = indices(ii)+frameN2;
            Chunk = Data(1:finish,:);
            Pad1 = zeros(abs(start)+1,nSweeps);
            Chunk = [Pad1;Chunk];
            w = W;
            
            pad1 = zeros(abs(start)+1,1);
            chunkY1 = Y1(1:finish,:);
            chunkY1 = [pad1;chunkY1];            
            chunkY1c = Y1c(1:finish,:);
            chunkY1c = [pad1;chunkY1c]; 
            chunkY2 = Y2(1:finish,:);
            chunkY2 = [pad1;chunkY2];            
            chunkY2c = Y2c(1:finish,:);
            chunkY2c = [pad1;chunkY2c];            
            chunkYdp = Ydp(1:finish,:);
            chunkYdp = [pad1;chunkYdp];            
            chunkYdpc = Ydpc(1:finish,:);
            chunkYdpc = [pad1;chunkYdpc]; 
            %w(1:abs(start)) = 0;
        elseif finish > nSamples
            hitCounter = hitCounter + 1;
            disp('hit')
            extra = finish - nSamples;
            Chunk = Data(start:end,:);
            Pad1 = zeros(extra,nSweeps);
            Chunk = [Chunk;Pad1];
            w = W;

            %pad1 = zeros(abs(start)+1,1);
            pad1 = Pad1(:,1);
            chunkY1 = Y1(start:end,:);
            chunkY1 = [chunkY1;pad1]; 
            chunkY1c = Y1c(start:end,:);
            chunkY1c = [chunkY1c;pad1]; 
            chunkY2 = Y2(start:end,:);
            chunkY2 = [chunkY2;pad1];            
            chunkY2c = Y2c(start:end,:);
            chunkY2c = [chunkY2c;pad1];            
            chunkYdp = Ydp(start:end,:);
            chunkYdp = [chunkYdp;pad1];            
            chunkYdpc = Ydpc(start:end,:);
            chunkYdpc = [chunkYdpc;pad1]; 
            %w(nSamples:end) = 0;
        else
            Chunk = Data(start:finish,:);
            w = W;

            chunkY1 = Y1(start:finish);
            chunkY1c = Y1c(start:finish);
            chunkY2 = Y2(start:finish);
            chunkY2c = Y2c(start:finish);
            chunkYdp = Ydp(start:finish);
            chunkYdpc = Ydpc(start:finish);
        end
        X = [chunkY1,chunkY1c,chunkY2,chunkY2c,chunkYdp,chunkYdpc];
        
        % perform wLSF for each frame in each sweep separately.
        try
        for jj=1:nSweeps
            [B(:,jj),yhat(:,jj),residuals(:,jj)] = OLSfit_internal(X,Chunk(:,jj),w);
        end
        catch ME
            keyboard
        end
        % compute weighted coherent mean, noise floor (standard error) and snr
        [signal(1,ii),nf(1,ii),snr(1,ii)] = bswSNR(B(1,:)); % for primary f1
        [signal(2,ii),nf(2,ii),snr(2,ii)] = bswSNR(B(2,:)); % for primary f2
        [signal(3,ii),nf(3,ii),snr(3,ii),ns] = bswSNR(B(3,:)); % for dpoae fdp
        %bigB = [bigB;ns(:)];
    end
    targets = targets(1:size(signal,2));
    
    % apply EPL correction ------------------------------------------------
    if DPOAE.doEPL == 1
        keyboard
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

    DPOAE.Ldp = signal(3,:)';
    DPOAE.Ndp = nf(3,:)';
    DPOAE.f2 = targets;
    DPOAE.L1_measured = signal(1,:)';
    DPOAE.N1_measured = nf(1,:)';
    DPOAE.L2_measured = signal(2,:)';
    DPOAE.N2_measured = nf(2,:)';
    DPOAE.LC = LC;
    DPOAE.nSweeps = nSweeps;
    DPOAE.hitCounter = hitCounter;

    % PLOTTING ------------------------------------------------------------
    h1 = figure;
    hold on
    plot(DPOAE.f2/1000,DPOAE.Ldp,'b')
    hold on
    plot(DPOAE.f2/1000,DPOAE.Ndp,'k')
    xmin = min(DPOAE.f2/1000);
    xmax = max(DPOAE.f2/1000);
    xlim([xmin,xmax])
    xlabel('Target F2 (kHz)')
    ylabel('DPOAE Level (dB SPL)')
    grid on
    title([DPOAE.subjID,'DPGRAM'])

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
            figureFileName = [obj.subjectID,'_inSituCal_','DPOAE_A_DPG','.fig'];
            %figureFileName = ARLas_saveName(savePath,figureFileName);
            savefig(iscCheckHandleA,[savePath,figureFileName])
            figureFileName = [obj.subjectID,'_inSituCal_','DPOAE_A_DPG','.bmp'];
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
            figureFileName = [obj.subjectID,'_inSituCal_','DPOAE_B_DPG','.fig'];
            %figureFileName = ARLas_saveName(savePath,figureFileName);
            savefig(iscCheckHandleB,[savePath,figureFileName])
            figureFileName = [obj.subjectID,'_inSituCal_','DPOAE_B_DPG','.bmp'];
            %figureFileName = ARLas_saveName(savePath,figureFileName);    
            saveas(iscCheckHandleB,[savePath,figureFileName])
            saveas(iscCheckHandleB,[tempPath,figureFileName])
        end
    catch
       disp('Warning: One or more figures for DPOAE_B not saved!')
    end
end
function [] = saveAnalysisPlots(obj,LC,h1)
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
            figureFileName = [obj.subjectID,'DPOAE_DPG','.fig'];
            figureFileName = ARLas_saveName(savePath,figureFileName);
            savefig(h1,[savePath,figureFileName])
            figureFileName = [obj.subjectID,'DPOAE_DPG','.bmp'];
            figureFileName = ARLas_saveName(savePath,figureFileName);    
            saveas(h1,[savePath,figureFileName])
            saveas(h1,[tempPath,figureFileName])
        end
    catch
        disp('Warning: One or more figures for DPOAE_LGF not saved!')        
    end

    % try
    %     if ~isempty(h2)
    %         savePath = [obj.savingGrace,obj.experimentID,'_analysis',obj.sep];
    %         if exist(savePath,'dir') == 0 % if path does not exist
    %             success = mkdir(savePath);
    %             addpath(savePath)
    %             if ~success
    %                 warning('Unable to create new Experiment Directory: data save path')
    %             end
    %         end 
    %         figureFileName = [obj.subjectID,'DPOAE_MLE','.fig'];
    %         figureFileName = ARLas_saveName(savePath,figureFileName);
    %         savefig(h2,[savePath,figureFileName])
    %         figureFileName = [obj.subjectID,'DPOAE_MLE','.bmp'];
    %         figureFileName = ARLas_saveName(savePath,figureFileName);    
    %         saveas(h2,[savePath,figureFileName])
    %         saveas(h2,[tempPath,figureFileName])
    %     end
    % catch
    %    disp('Warning: One or more figures for DPOAE_MLE not saved!')
    % end
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
            fileName = [obj.subjectID,'DPOAE_B_DPG','.mat'];
            figureFileName = ARLas_saveName(savePath,fileName);
            save([savePath,fileName],'DPOAE_B')
        end
    catch
        disp('Warning: Mat file for DPOAE_B_DPG not saved!')        
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
            fileName = [obj.subjectID,'DPOAE_A_DPG','.xls'];
            saveFileName = [fileName,'x'];
            saveFileNameCSV = [obj.subjectID,'DPOAE_A_DPG','.csv'];
            sheetName = ['DPGRAM'];
            writematrix('F2_Hz',[savePath,saveFileName],'Sheet',sheetName,'Range','A1');
            writematrix(DPOAE_A.f2,[savePath,saveFileName],'Sheet',sheetName,'Range','A2');
            writematrix('Ldp_dB_EPL',[savePath,saveFileName],'Sheet',sheetName,'Range','B1');
            writematrix(DPOAE_A.Ldp,[savePath,saveFileName],'Sheet',sheetName,'Range','B2');
            writematrix('Ndp_dB_EPL',[savePath,saveFileName],'Sheet',sheetName,'Range','C1');
            writematrix(DPOAE_A.Ndp,[savePath,saveFileName],'Sheet',sheetName,'Range','C2');
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
            fileName = [obj.subjectID,'DPOAE_B_DPG','.xls'];
            saveFileName = [fileName,'x'];
            saveFileNameCSV = [obj.subjectID,'DPOAE_B_DPG','.csv'];
            sheetName = ['DPGRAM'];
            writematrix('F2_Hz',[savePath,saveFileName],'Sheet',sheetName,'Range','A1');
            writematrix(DPOAE_B.f2,[savePath,saveFileName],'Sheet',sheetName,'Range','A2');
            writematrix('Ldp_dB_EPL',[savePath,saveFileName],'Sheet',sheetName,'Range','B1');
            writematrix(DPOAE_B.Ldp,[savePath,saveFileName],'Sheet',sheetName,'Range','B2');
            writematrix('Ndp_dB_EPL',[savePath,saveFileName],'Sheet',sheetName,'Range','C1');
            writematrix(DPOAE_B.Ndp,[savePath,saveFileName],'Sheet',sheetName,'Range','C2');
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
function [B,yhat,residuals] = OLSfit_internal(X,data,w)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [B,yhat,residuals] = OLSfit_internal(X,data,w,t)
%
% Ordinary Least-Squares Regression for the model
% y = a + b1*X1 + b2*X2 + ... + bk*Xk.
% X = matrix of independent variables.
% y = vector of dependent variables.
% w = vector of weighting values (optional input)
%
% Author: Shawn Goodman, PhD
% Auditory Research Lab, the University of Iowa
% Date: March 7, 2022
% Last Updated: July 20, 2023 -- ssg
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    y = data;
    X = [X,ones(size(X,1),1)]; % solution matrix
    [rows,columns] = size(X); % size of the solution matrix
    %k = columns-1; % k is the number of IV (the column of ones for y-intercept/dc offset is not counted)
    %n = rows; % number of observations
    
    if ~isempty(w) % if a weighting vector exists...
        w = w(:); % vectorize as a column vector...
        sqrtw = sqrt(w); % ...and take the square root
        X = repmat(sqrtw,1,columns).*X; % replicate sqrtw to matrix the size of X; J is the weighted version of X
        y = sqrtw.*y; % Dy is the weighted version of y
    end
    
    % Solve the least squares problem
    [Q,R] = qr(X,0); % orthogonal-triangular decomposition; R is the Cholesky factor of the X matrix
                     % the input argument 0 makes an "economy sized" decomposition, so that
                     % [nSamples,nIVs] = size(Q), and 
                     % [nIvs,nIVs] = size(R).
    b = full(R\(Q'*y)); % Same as p = D*X\(D*y); b is a vector of coefficients; also same as b = (X'*X)\(X'*Y).
    yhat = X*b;      % Predicted responses at each data point.
    residuals = y - yhat;

    b = b(1:end-1); % throw away dc offset
    % loop to put cosine and sine parts into complex form
    counter = 1; 
    for ii=1:2:length(b)
        B(counter,1) = b(ii) + 1i*b(ii+1); 
        counter = counter + 1; 
    end 

end

% OLD CODE ----------------------------------------------------------------

        % if runPostAnalysis == 1 % running post-hoc analyses on previously recorded data
        %     if ~isempty(dL)
        %         dummy = load([postPath,dL(ii).name]);
        %         headerL = dummy.header;
        %         DataL = dummy.data;
        %         clear dummy
        %         micCorrectionL = headerL.userInfo.micCorrectionL;
        %         DataL = applyMicCorrection(DataL,micCorrectionL);
        %         DPOAE_L.postHocTimeStamp = DPOAE_L.timeStamp; % note that this is a post-hoc analysis
        %         DPOAE_L.timeStamp = obj.timeStamp; % make sure that original time stamp shows up on the figues
        %     end
        %     if ~isempty(dR)
        %         dummy = load([postPath,dR(ii).name]);
        %         headerR = dummy.header;
        %         DataR = dummy.data;
        %         clear dummy
        %         micCorrectionR = headerR.userInfo.micCorrectionR;
        %         DataR = applyMicCorrection(DataR,micCorrectionR);
        %         DPOAE_R.postHocTimeStamp = DPOAE_R.timeStamp; % note that this is a post-hoc analysis
        %         DPOAE_R.timeStamp = obj.timeStamp; % make sure that original time stamp shows up on the figues
        %     end
        % end    

                % % Analyze data ----------------------------------------------------
        % disp('----- Analyzing data -----')
        % % analysis prep -----
        % if ii==1 % on the first iteration, initialize the analysis structure
        %     if ~isempty(probeInputL)
        %         DPOAE_L = getDPstruct(obj,fs,fmin,fmax,sweepRate,nSweeps(ii),targetL1(ii),targetL2(ii),calType,targetCalType,iscS1L,iscS2L,iscS12L,C1L,C2L,'Left');
        %         if ~strcmp(testEar,'Both')
        %             DPOAE_L.ear = testEar;
        %         end
        %     end
        %     if ~isempty(probeInputR)
        %         DPOAE_R = getDPstruct(obj,fs,fmin,fmax,sweepRate,nSweeps(ii),targetL1(ii),targetL2(ii),calType,targetCalType,iscS1R,iscS2R,iscS12R,C1R,C2R,'Right');
        %         if ~strcmp(testEar,'Both')
        %             DPOAE_R.ear = testEar;
        %         end
        %     end
        % end 
        
        % % ANALYSIS BEGINS HERE --------------------------------------------
        % if ~isempty(probeInputL)
        %     DPOAE_L = dpoaeAnalysis_fixedRatioSweep_vRT(headerL,DataL,DPOAE_L,[]);
        %     %DPOAE_L = dpoae_sweepRatioANALYSIS_vRT3(headerL,DataL,DPOAE_L,[]);
        % end
        % if ~isempty(probeInputR)
        %     DPOAE_R = dpoaeAnalysis_fixedRatioSweep_vRT(headerR,DataR,DPOAE_R,[]);
        %     %DPOAE_R = dpoae_sweepRatioANALYSIS_vRT3(headerR,DataR,DPOAE_R,[]);
        % end
        % 
        % disp('----- Plotting data -----') % -------------------------------
        % try
        %     if ~isempty(probeInputL)
        %         [h1L,h2L] = dpoaePlotting_fixedRatioSweep_vRT(DPOAE_L);
        %         %h1L.Position = [591   154  469  211];
        %         %h2L.Position = [591  442  469  314];
        %         pause(0.01)
        %     end
        % catch
        % end
        % try
        %     if ~isempty(probeInputR)
        %         [h1R,h2R] = dpoaePlotting_fixedRatioSweep_vRT(DPOAE_R);
        %         %h1R.Position = [1060  154  469  211];
        %         %h2R.Position = [1060  442  469  314];
        %         pause(0.01);
        %     end
        % catch
        % end


    % % Save figures ----------
    % disp('----- Saving figures -----')
    % try
    %     if ~isempty(probeInputL)
    %         if runPostAnalysis == 1 % saving post-hoc analyses on previously recorded data
    %             savePath = postPathSave;
    %         else
    %             savePath = [obj.savingGrace,obj.experimentID,'_analysis',obj.sep];
    %             %savePath = obj.objPlayrec.savedFilesPath;
    %             if exist(savePath,'dir') == 0 % if path does not exist
    %                 success = mkdir(savePath);
    %                 if ~success
    %                     warning('Unable to create new Experiment Directory: data save path')
    %                 end
    %             end 
    %         end
    %         %figureFileName = [obj.subjectID,'_dpoaePrimariesL.fig'];
    %         figureFileName = [obj.subjectID,'_dpoaePrimaries_',DPOAE_L.ear,'.fig'];
    %         figureFileName = ARLas_saveName(savePath,figureFileName);
    %         savefig(h1L,[savePath,figureFileName])
    %         %figureFileName = [obj.subjectID,'_dpoaePrimariesL.bmp'];
    %         figureFileName = [obj.subjectID,'_dpoaePrimaries_',DPOAE_L.ear,'.bmp'];
    %         figureFileName = ARLas_saveName(savePath,figureFileName);    
    %         saveas(h1L,[savePath,figureFileName])
    % 
    %         %figureFileName = [obj.subjectID,'_dpoaeL.fig'];
    %         figureFileName = [obj.subjectID,'_dpoae_',DPOAE_L.ear,'.fig'];
    %         figureFileName = ARLas_saveName(savePath,figureFileName);
    %         savefig(h2L,[savePath,figureFileName])
    %         figureFileName = [obj.subjectID,'_dpoae_',DPOAE_L.ear,'.bmp'];
    %         figureFileName = ARLas_saveName(savePath,figureFileName);    
    %         saveas(h2L,[savePath,figureFileName])
    %     end
    % catch
    %     disp('Warning: One or more figures for DPOAE_L not saved!')        
    % end

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
%             figureFileName = [obj.subjectID,'_dpoaePrimaries_',DPOAE_R.ear,'.fig'];
%             figureFileName = ARLas_saveName(savePath,figureFileName);
%             savefig(h1R,[savePath,figureFileName])
%             figureFileName = [obj.subjectID,'_dpoaePrimaries_',DPOAE_R.ear,'.bmp'];
%             figureFileName = ARLas_saveName(savePath,figureFileName);    
%             saveas(h1R,[savePath,figureFileName])
% 
%             figureFileName = [obj.subjectID,'_dpoae_',DPOAE_R.ear,'.fig'];
%             figureFileName = ARLas_saveName(savePath,figureFileName);
%             savefig(h2R,[savePath,figureFileName])
%             figureFileName = [obj.subjectID,'_dpoae_',DPOAE_R.ear,'.bmp'];
%             figureFileName = ARLas_saveName(savePath,figureFileName);    
%             saveas(h2R,[savePath,figureFileName])
%         end
%     catch
%        disp('Warning: One or more figures for DPOAE_R not saved!')
%     end
% 
%     % Save analyses ------------------------
%     disp('----- Saving analysis files -----')
%     try 
%         if ~isempty(probeInputL)
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
%             %saveFileName = [obj.subjectID,'_analyzedDPOAE_L.mat'];
%             saveFileName = [obj.subjectID,'_analyzedDPOAE_',DPOAE_L.ear,'.mat'];
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
%             %saveFileName = [obj.subjectID,'_analyzedDPOAE_R.mat'];
%             saveFileName = [obj.subjectID,'_analyzedDPOAE_',DPOAE_R.ear,'.mat'];
%             saveFileName = ARLas_saveName(savePath,saveFileName);
%             save([savePath,saveFileName],'DPOAE_R','ISCR')
%         end
%     catch
%         disp('Warning: DPOAE_R Analysis not saved!')
%     end
% 
%     % Save in situ cal figures --------------------------------------------
%     disp('----- Saving figures -----')
%     try
%         if ~isempty(probeInputL)
%             savePath = [obj.savingGrace,obj.experimentID,'_analysis',obj.sep];
%             if exist(savePath,'dir') == 0 % if path does not exist
%                 success = mkdir(savePath);
%                 if ~success
%                     warning('Unable to create new Experiment Directory: data save path')
%                 end
%             end 
%             figureFileName = [obj.subjectID,'_inSituCal_',DPOAE_L.ear,'.fig'];
%             figureFileName = ARLas_saveName(savePath,figureFileName);
%             savefig(iscCheckHandleL,[savePath,figureFileName])
%             figureFileName = [obj.subjectID,'_inSituCal_',DPOAE_L.ear,'.bmp'];
%             figureFileName = ARLas_saveName(savePath,figureFileName);    
%             saveas(iscCheckHandleL,[savePath,figureFileName])
%         end
%     catch
%         if ~isempty(obj)
%             disp('Warning: inSitu Cal L not saved!')        
%         end
%     end
% 
%     try
%         if ~isempty(probeInputR)
%             savePath = [obj.savingGrace,obj.experimentID,'_analysis',obj.sep];
%             if exist(savePath,'dir') == 0 % if path does not exist
%                 success = mkdir(savePath);
%                 if ~success
%                     warning('Unable to create new Experiment Directory: data save path')
%                 end
%             end 
%             figureFileName = [obj.subjectID,'_inSituCal_',DPOAE_R.ear,'.fig'];
%             figureFileName = ARLas_saveName(savePath,figureFileName);
%             savefig(iscCheckHandleR,[savePath,figureFileName])
%             figureFileName = [obj.subjectID,'_inSituCal_',DPOAE_R.ear,'.bmp'];
%             figureFileName = ARLas_saveName(savePath,figureFileName);    
%             saveas(iscCheckHandleR,[savePath,figureFileName])
%         end
%     catch
%        if ~isempty(obj)
%            disp('Warning: inSitu Cal R not saved!')
%        end
%     end
% 
%     disp('----- Finished with DPOAE Fixed Ratio Sweep experiment -----')
%     disp(' ')
%     disp(' ')
% toc/60
