function [peakLocations1,peakLocations2,calParams] = DW_capLevelSeries_Thev(varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [peakLocations1,peakLocations2,calParams] = DW_capLevelSeries_Thev(varargin)
%
% Measure CAPs at mulitple frequencies and levels.
% Input Arguments: obj,preInjectionMin,pumpRateSeries,counter,peakLocations1,peakLocations2
% Required functions: DW_analyzeSeries.m
%
% Authors: Shawn S Goodman, PhD & Jeffery T Lichtenhan PhD
% Auditory Research Lab, Dept. of Comm. Sciences & Disorders, The University of Iowa
% Auditory Physiology Lab, Dept. of Otolaryngology, Washington University School of Medicine, St. Louis
% Date of Original Function: January 8-10, 2019
% Major Update: September 22-24, 2019 -- ssg
%                                     Includes numerous updates. Note this
%                                     no longer does MOC measurements.
% Updated: September 25-26, 2019 -- ssg
%                                Upgraded to save version numbers to header file.
% Updated: November 7-12, 2019 -- ssg
%                                Upgraded to Thevenin calibration
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    experimentVersion = '11.12.2019'; % this is the current version number of this program    
    obj = varargin{1}; % get the arlas object
    fs = obj.fs; % sampling rate in Hz
    preInjectionMin = varargin{2}; % number of requested pre-injection minutes
    pumpRateSeries = varargin{3}; % vector of pump rates
    pumpMasterCounter = varargin{4}; % current location in the pump rates vector
    tau = (pumpMasterCounter - preInjectionMin) -1; % assumed experiment time (re: start of injection)
    peakLocations1 = varargin{5}; % location of the n1 CAP peak
    peakLocations2 = varargin{6}; % location of the p1 CAP peak
    pumpMasterVersion = varargin{7}; % the version of pumpMaster that is being used
    calParams = varargin{8}; % how many times pump master has called this program
    
    % User controlled parameters -----------------------------------------------
    pipProbe = 'B';            % which 10x probe to present the pips through
    doVerify = 0;              % verify the output levels in coupler (1) or regular test mode (0)
    
    % IMPORTANT NOTE: Specify low- and high-cut frequencies (Hz) for ALL stimuli.
    %                 If low- and high-cut values are the same, will create
    %                 a pip. If different, will create a click.
    % IMPORTANT NOTE: Make sure that the requested values match pip and
    %                 click audiogram values (from DW_fastCapAudio_Thev.m and from
    %                 DW_fastCapAudioClicks_Thev.m.
    lowCut  = [ 200,  200, 4000, 8000, 16000, 20000]; % low-frequency cutoff value (Hz)
    highCut = [4000,32000, 4000, 8000, 16000, 20000]; % high-frequency cutoff values (Hz)
    levels = [10,20,30,40,50]; % stimulus levels to present (dB SL)

    % for doVerify:
    %lowCut  = [  200,  200]; % low-frequency cutoff value (Hz)
    %highCut = [32000,32000]; % high-frequency cutoff values (Hz)
    %levels = [50,60]; % stimulus levels to present (dB SL)
    
    minReps = 12;               % minimum number of stimulus repetitions allowed; must be multiple of 6. 12 is the suggested number
    targetSNR = 15;             % desired SNR to achieve (move on once this level is reached or max reps is reached). At threshold this program is expecting 6 dB. So, targetSNR can't be above 6 if you want 0 dB SL
    maxTotalLength = 50;        % maximum total length to complete series (seconds).

    pathName =  ['C:\myWork\ARLas\Data\',obj.experimentID,'\'];  % location of the saved cap audiogram 
    fileNamePIP= 'CAP_Audiogram_CAPaudio1_1.mat';  % name of the saved PIP cap audiogram file
    fileNameCLK= 'CAP_AudiogramClicks_CAPaudioClicks2_1.mat';  % name of the saved CLICK cap audiogram file
    %--------------------------------------------------------------------------

    freqs = [lowCut;highCut]; % the frequency content of each stimulus (Hz)
    nFreqs = size(freqs,2);   % the number of stimuli to test (including pips and clicks)
    doRandomize = 1;  % turn on (1) or off (0) stimulus order randomization
    maxStimSPL = 100; % PIP maximum stimulus level in dB SPL--will not test SL levels that exceed this value
    maxStimSPLreFO = -15; % CLICK maximum output level to present (dB re: full out)
    
    [inputs,outputs] = hardwareSetup;
    inputG = inputs{3};               % GRASS Bioamp
    if strcmp(pipProbe,'A')
        input = inputs{1};         % ER10xA microphone
        output = outputs{1};       % ER10xA loudspeakers
    elseif strcmp(pipProbe,'B')
        input = inputs{2};         % ER10xB microphone
        output = outputs{2};       % ER10xB loudspeakers
    end
    if doVerify == 1
        inputRef = inputs{4};         % GRAS 1/8" reference microphone
    end
    
    [pathName_mic,folderName_mic,fileName_mic] = mostRecentCalibration('mic',input.label,[]);
    if ~isempty(fileName_mic)
        dummy = load([pathName_mic,folderName_mic,fileName_mic]);
        micCorrection = dummy.micCorrection; % microphone correction
    else
        micCorrection = 1; % convolving by this changes nothing
    end
    
    if isempty(calParams) % this will be empty on the first function call, but not subsequently
        calParams.fmin = 200;   
        calParams.fmax = 32000; 
        calParams.testProbe = pipProbe;
        calParams.usingChannels = 2; % use this channel (1 or 2; only one channel needed for this experiment)
        calParams = ARLas_earCanalRecordings_DW10x(obj,calParams); % run ear canal recordings to get in-situ measurement
    end
    
    dummyPIP = load([pathName,fileNamePIP]);
    dummyCLK = load([pathName,fileNameCLK]);
    for ii=1:nFreqs % calculate threshold for each requested stimulus
        if freqs(1,ii) ~= freqs(2,ii) % if stimulus is a click
            for kk=1:size(dummyCLK.FRQ,2)
                if (freqs(1,ii)==(dummyCLK.FRQ(kk,1)*1000) & freqs(2,ii)==(dummyCLK.FRQ(kk,2))*1000) == 1
                    try
                        indx = kk;
                    catch
                        keyboard
                        end
                    break
                end
            end
            thds(ii,1) = dummyCLK.THD(indx); % click threshold re full out
            stimType(ii,1) = 1; % designate stim type as click (1; because shape of a click)
            originalReps(ii,1) = dummyCLK.originalReps; % number of recordings used to get original thresholds            
        elseif freqs(1,ii) == freqs(2,ii) % if stimulus is a pip
            thds(ii,1) = interp1(dummyPIP.FRQ *1000,dummyPIP.THD,freqs(1,ii),'pchip');
            stimType(ii,1) = 0; % designate stim type as pip (0)
            originalReps(ii,1) = dummyPIP.originalReps; % number of recordings used to get original thresholds            
        end
    end        
    levels = levels(:);      % force to be a column vector
    thds = thds(:)';         % force to a row vector
    stimType = stimType(:)';
    originalReps = originalReps(:)';

    [stim,cuts] = getStim(1000,fs);            % get sample cap stimulus (tone burst)
    stimLength = length(stim)/fs;              % stimulus length (seconds). This includes both polarities and zero padding. Typically ~13 ms
    maxTotalReps = floor(maxTotalLength / stimLength); % maximum possible number of stimulus presentations in allotted time
    
    stimSPL = repmat(levels,1,size(freqs,2)) + repmat(thds,size(levels,1),1); % SPL = SL + thd; requested level of each stimulus (dB SPL)
    stimSL = stimSPL - repmat(thds,size(levels,1),1);                         % requested level of each stimulus (dB SL)
    expectedSNR = stimSL + 6;                                                 % expected SNR that would be achieved using originalReps;
                                                                              %  here, we are assuming that threshold is close to 6 dB SNR
                                                                              % amount of change in SNR that is desired to make everything meet targetSNR
                                                                              % code for direction of change (increase or decrease)
    dBchange = targetSNR - expectedSNR;                                       
    pos = dBchange >= 0;                                                      
    neg = (dBchange < 0) *-1;
    direction = pos + neg;
    doublings = abs(dBchange/3);                                              % do this many doublings/halvings in order to get desired snr
    nReps = floor(2.^(log2(originalReps) + (direction.*doublings)));           % requested number of repetitions
    nReps(find(nReps < minReps)) = minReps;                                   % make sure no reps are less than minReps
    nReps = nReps .* (stimSPL < maxStimSPL);                                  % do not test levels that are too high
    NReps = sum(sum(nReps));                                                  % total number of reps to be tested
    if NReps > maxTotalReps                                                   % check to make sure the total number of possible reps is not exceeded
       %keyboard
       warning('Requested number of repetitions exceeds the number possible. Choose fewer levels or frequencies, or max dB SPL.')
    end
    nReps = floor(nReps * (maxTotalReps/NReps));                              % re-allocate the number of reps that are left over, so as to fill up the time allotted time

    
    for ii=1:nFreqs % calculate threshold for each requested stimulus
        if stimType(ii) == 1 % if stimulus is a click
            % requested values with attenuation > 0 cannot be created.
            % therefore, simply turn off
            stimSPL(find(stimSPL(:,ii) > 0),ii) = -300;
            %disp('WARNING: a requested click SL results in a click amplitude >1.')
        end
    end

    % randomization needs to be done in blocks of 6 repetitions.
    % Sub-averaged in blocks of 6, the artifact (multiples of 60 Hz) will
    % maximally cancel. This means that when performing randomization, blocks
    % of 6 must stay together. (Note: actual randomization done further down in
    % code.)
    nReps = nReps - mod(nReps,6); % constrain reps to be muliples of 6

    NReps = sum(sum(nReps));                                                  % re-calculate the total number of reps to be tested
    if NReps > maxTotalReps                                                   % check to make sure the total number of possible reps is not exceeded
        error('Requested number of repetitions exceeds the number possible.')
    end
    stimSPL = stimSPL .* (stimSPL<maxStimSPL);                                % zero out levels that will not be presented due to overexposure


    nLevels = length(levels);
    nFreqs = length(freqs);
    Stim = zeros(length(stim),NReps);                                         % initialize the stimulus matrix
    counter = 1;
    for ii=1:nLevels
        for jj=1:nFreqs
            if stimType(jj) == 0 % if stimulus is a pip
                f = freqs(1,jj);             % current frequency
                [stim,~] = getStim(f,fs);    % get sample cap stimulus (tone burst)
                desiredLvl = stimSPL(ii,jj); % stimulus level to present
                desiredLvl = desiredLvl + 3; % the values shown above refer rms levels; however, they are
                     % implemented using stimulus peak levels. 
                     % Therefore, add 3 dB to make them come out correctly.
                if doVerify == 1
                    type = 'ipl';
                else
                    type = 'fpl';
                end
                try
                [stimScaled,a,errors] = ARLas_applyISC(calParams.isc,desiredLvl,type,f,stim);
                catch 
                    keyboard
                end
                % stim = stim * abs(a);
                % The stimulus is in peak amplitude, but we want the target to
                % be in rms. Therefore, add 3 dB:
                %stim = stimScaled * (1./sqrt(.5));
                stim = stimScaled;
                nn = nReps(ii,jj);
                if nn > 0
                    for kk=1:nn
                        if max(abs(stim)) > 1
                            keyboard
                        elseif max(abs(stim)) == 1
                            stim = stim * 0.99;
                        end
                        Stim(:,counter) = stim;
                        counter = counter + 1;
                    end
                end
            elseif stimType(jj) == 1 % if stimulus is a click
                f = freqs(:,jj)';               % current frequency cutoff values
                [stim,~] = getStimClick(fs);    % get sample cap stimulus (click)
                desiredLvl = stimSPL(ii,jj); % stimulus level to present
                if doVerify == 1
                    type = 'ipl';
                else
                    type = 'fpl';
                end
                DL = 60; % desired level. Scaling will take place after obtained,
                         % so simply use a constant level that should always be
                         % available (60 dB here)
                try
                [~,a,errors] = ARLas_applyISC(calParams.isc,DL,type,f,stim);
                catch
                    keybard
                end
                a = a ./ max(abs(a)); % scale to 100%
                multiplier = 10^(desiredLvl/20);
                a = a .* multiplier;
                stim = fastFilter(a,stim);

                if any(isnan(stim))
                    keyboard
                end

                nn = nReps(ii,jj);
                if nn > 0
                    for kk=1:nn
                        if max(abs(stim)) > 1.05
                            keyboard
                        elseif max(abs(stim)) >= 1
                            stim = stim * 0.99;
                        end
                        Stim(:,counter) = stim;
                        counter = counter + 1;
                    end
                end
            end
        end
    end

    if doRandomize == 1                                                       % turn on and off randomization
        r = rand(1,NReps/6);
        [~,ri] = sort(r);
        ri = (ri-1)*6+1;
        Ri = [ri',ri'+1,ri'+2,ri'+3,ri'+4,ri'+5];
        randIndx = reshape(Ri',1,NReps);    
        Stim = Stim(:,randIndx);                                              % randomize the stimuli
        [~,sortIndx] = sort(randIndx(:));                                     % index to re-sort the randomized data
    else
       randIndx = (1:1:NReps);                                                % if no randomization is desired, leave in order
       sortIndx = randIndx;
    end

    obj.objPlayrec.userInfo = [];                                             % clear old header info  
    obj.objPlayrec.userInfo.nReps = nReps;                                    % load new header info
    obj.objPlayrec.userInfo.stimSPL = stimSL; 
    obj.objPlayrec.userInfo.stimSPL = stimSPL; 
    obj.objPlayrec.userInfo.freqs = freqs;
    obj.objPlayrec.userInfo.thds = thds;
    obj.objPlayrec.userInfo.levels = levels;
    obj.objPlayrec.userInfo.randIndx = randIndx;
    obj.objPlayrec.userInfo.sortIndx = sortIndx;
    obj.objPlayrec.userInfo.maxStimSPL = maxStimSPL;
    obj.objPlayrec.userInfo.minReps = minReps;
    obj.objPlayrec.userInfo.targetSNR = targetSNR;
    obj.objPlayrec.userInfo.maxTotalLength = maxTotalLength;
    obj.objPlayrec.userInfo.maxTotalLength = maxTotalReps;
    obj.objPlayrec.userInfo.NReps = NReps;
    obj.objPlayrec.userInfo.cuts = cuts;
    obj.objPlayrec.userInfo.preInjectionMin = preInjectionMin;
    obj.objPlayrec.userInfo.pumpRateSeries = pumpRateSeries;
    obj.objPlayrec.userInfo.pumpMasterCounter = pumpMasterCounter;
    obj.objPlayrec.userInfo.pumpMasterVersion = pumpMasterVersion;
    obj.objPlayrec.userInfo.experimentVersion = experimentVersion;
    obj.objPlayrec.userInfo.tau = tau;
    obj.objPlayrec.userInfo.calParams = calParams;
    %--------------------------------------------------------------------------

    obj.clearRecList % clear out whatever was used previously for recordings
    obj.setRecList(inputG.ch,inputG.label,inputG.micSens,inputG.gain); % load the Gras bioamp first
    obj.setRecList(input.ch,input.label,input.micSens,input.gain); % load the pip stimulus probe second
    if doVerify == 1
        obj.setRecList(inputRef.ch,inputRef.label,inputRef.micSens,inputRef.gain); % load the reference mic third
    end
    
    obj.clearPlayList; % clear out whatever was used previously for playback
    obj.setPlayList(Stim,output.ch(1));  % load output for pip stimuli. This is always used
    obj.setNReps(NReps); % number of stimulus reps to play and record. This will run through the stimulus matrix ONCE!

    obj.objPlayrec.run   % playback and record
    if obj.killRun
        return  
    end
    
    if doVerify == 1
        [headerAudio,DataAudio] = obj.retrieveData(input.label); % get acoustic pip recordings (audio from ER10X)
        DataAudio = fastFilter(micCorrection,DataAudio);
        [headerRef,DataRef] = obj.retrieveData(inputRef.label);     % get reference data
        
        keyboard        
        [magA,nfA,snrA,freqA,capA,timeA,stimulusSPLA,phiA,signalA,noiseA,frequencyA] = DW_fastCapAnalyzeInternal(headerAudio,DataAudio,1,stimType);
        [magR,nfR,snrR,freqR,capR,timeR,stimulusSPLR,phiR,signalR,noiseR,frequencyR] = DW_fastCapAnalyzeInternal(headerRef,DataRef,1,stimType);
        disp(' ')
        disp('Difference between reference and ipl target:')
        [freqs'/1000,(max(magR,[],1)' - stimSPL)]
    else
        [header,Data] = obj.retrieveData(inputG.label);   % get electrical data (cap from gras amplifier)
    end
    
    % analyze the data ----------------------------------------------------
% take the following out!!!
%[header,Data] = obj.retrieveData(input.label); % get acoustic pip recordings (audio from ER10X)
%Data = fastFilter(micCorrection,Data);
    
    [mag,nf,snr,freq,cap,time,stimulusSPL,phi,signal,noise,frequency] = DW_fastCapAnalyzeInternal(header,Data,0,stimType);
    %[FF,SPL,SL,MAG,NF,nF,nMaxLevels,peakLocations1,peakLocations2] = DW_analyzeSeriesInternal(header,Data,peakLocations1,peakLocations2);
    
    % plot data -----
    try 
        delete(h1)
    end
    h1 = figure(2000);
    plt(time,cap)
    xlabel('Time (ms)')
    ylabel('Amplitude (uV)')
    title('CAP Waveforms (All)')
    
    %offset = 5;
    %h1 = plotData(h1,mag,freqs(1,:),nf,nFreqs,nMaxLevels,offset);

end % EOF: end of experiment file


% internal functions ------------------------------------------------------
function [stim,cuts] = getStimClick(fs) % get the cap stimulus (tone burst)
    % This stim designed October 25, 2019
    badHarmonic = 180; % what frequency to cancel (Hz). Here, 60 and 120 are already filtered out by the GRASS filter settings
    nCycles = 1.5; % number of cycles. must be n+.5, where n is any positive integer.
    nSamples = round(((1/badHarmonic)*nCycles)*fs); % number of samples in the pip
    %t = (0:1:nSamples-1)'/fs; % time vector
    %center = round(length(t)/2); % center point of the stimulus
    stim = zeros(nSamples,1);
    %stim(center) = 1;
    stim(144) = 1;
    zpad = zeros(nSamples,1); % make zero padding
    stim = [zpad;stim;zpad;-stim;zpad]; % concatenate to make two pips, one out of phase with the other
    cuts = [nSamples*1+1,nSamples*2+1;  nSamples*3+1,nSamples*4+1]; % samples for where to cut the pips
end
function [stim,cuts] = getStim(f,fs) % get the cap stimulus (tone burst)
    % Phis stim designed June 8, 2017
    % Previously, we were using something 50 ms in length, which is an integer
    % of 60 Hz noise and harmonics. This stimulus is designed to cancel out
    % such electrical noise
    badHarmonic = 180; % what frequency to cancel (Hz). Here, 60 and 120 are already filtered out by the GRASS filter settings
    nCycles = 1.5; % number of cycles. must be n+.5, where n is any positive integer.
    nSamples = round(((1/badHarmonic)*nCycles)*fs); % number of samples in the pip
    pipLen = nSamples/fs; % pip length (seconds)
    t = (0:1:nSamples-1)'/fs; % time vector
    stim = sin(2*pi*f*t); % sinusoid
    rampLen = 0.001; % onset and offset ramp lengths (s) <-----------------------------------------------------------
    rampSamples = round(rampLen * fs); % number of samples in stimulus
        if mod(rampSamples,2) ~= 0 % force to be even number
            rampSamples = rampSamples + 1;
        end
    h = hann(rampSamples*2); % apply hanning window ramps
    h = h(1:rampSamples);
    stim(1:rampSamples,:) = stim(1:rampSamples,:) .* h;
    stim = flipud(stim);
    stim(1:rampSamples,:) = stim(1:rampSamples,:) .* h;
    stim = flipud(stim);
    zpad = zeros(nSamples,1); % make zero padding
    stim = [zpad;stim;zpad;-stim;zpad]; % concatenate to make two pips, one out of phase with the other
    cuts = [nSamples*1+1,nSamples*2+1;  nSamples*3+1,nSamples*4+1]; % samples for where to cut the pips
end
function [h1] = plotData(h1,P2P,FF,NF,nF,nL,offset)
    figure(h1)
        for ii=1:nF
            plot((FF(:,ii)+offset)/1000,20*log10(P2P(:,ii)/20),'*-','LineWidth',1)
            hold on
        end
        for ii=1:nL
            plot(FF(ii,:)/1000,20*log10(NF(ii,:)/20),'--.k')
        end
        xlabel('Frequency (kHz)')
        ylabel('Cap Amplitude (dB re: 20uV)')
        hold off
end
function [mag,nf,snr,freq,cap,time,stimulusSPL,phi,signal,noise,frequency] = DW_fastCapAnalyzeInternal(header,Data,audioSwitch,stimType)
    % ---------------------------------------------------------------------
    % set to 1 to analyze ER10X audio, rather than CAP recordings
    % ---------------------------------------------------------------------

    % de-interleave the data
    sortIndx = header.userInfo.sortIndx;
    Data = Data(:,sortIndx);  % put data back into sorted order

    % unpack needed information from the header
    levels = header.userInfo.levels;
    freqs = header.userInfo.freqs;
    nReps = header.userInfo.nReps;
    stimSPL = header.userInfo.stimSPL;
    fs = header.fs;
    cuts = header.userInfo.cuts;
    %stimSL = stimSPL;

    % Analyze Data for each level and frequency.
    nLevels = length(levels); % number of levels tested
    nFreqs = size(freqs,2);   % number of frequencies tested
    nSamples = size(Data,1);  % number of samples in each recorded buffer
    counter = 1; 
    start = 1;
    snrCutP1 = 8;
    snrCutN1 = 8;
    
    for ii=1:nLevels
        for jj=1:nFreqs
            f = freqs(:,jj);              % current frequency
            nn = nReps(ii,jj);          % number of recordings at this combintaion
            finish = start + nn-1;
            chunk = Data(:,start:finish);
            if nn > 0                     % if data were actually recorded at this level and frequency combination
                nSubs = size(chunk,2)/6;  % number of sub-averages to perform. Always in sets of 6 to cancel electrical artifact
                step2 = 6;
                start2 = 1;
                finish2 = step2;
                counter2 = 1;
                chunk2 = zeros(nSamples,nSubs);
                for kk=1:nSubs           % loop to create sub-averages
                    chunk2(:,counter2) = mean(chunk(:,start2:finish2),2);
                    start2 = start2 + step2;
                    finish2 = finish2 + step2;
                    counter2 = counter2 + 1;
                end
                % this is the new way; it uses Fourier-based RMS magnitudes
                if audioSwitch == 0 % if analyzing cap waveforms
                    [Cap,time] = getCAP(chunk2,cuts,fs);
                    peakLocations1 = NaN;
                    peakLocations2 = NaN;
                    [mag(1,counter),nf(1,counter),snr(1,counter),cap(:,counter),time,peakLocations1,peakLocations2] = DW_capAnalysisPeaks(Cap,time,fs,peakLocations1,peakLocations2,snrCutP1,snrCutN1);
                    %[mag(1,counter),nf(1,counter),snr(1,counter),phi(1,counter),cap(:,counter),signal(:,counter),noise(:,counter),frequency] = DW_capAnalysis(Cap,time,fs);
                    freq(:,counter) = f;
                    %stimulusSL(1,counter) = stimSL(ii,jj);   % stimulus sensation level
                    stimulusSPL(1,counter) = stimSPL(ii,jj); % stimulus sound pressure level
                    phi = [];
                    signal = [];
                    noise = [];
                    frequency = [];                    
                else % if analyzing audio waveforms
                    q = chunk2;
                    q1 = q(cuts(1,1):cuts(1,2),:);
                    q2 = q(cuts(2,1):cuts(2,2),:);
                    q = [q1,-q2];
                    len = 0.001;
                    q = ARLas_ramp(q,fs,len);
                    if stimType(jj) == 0 % for pip audio
                        [fff,signal,noiseFloor,phase] = ARLas_fda(q,fs,0.00002,fs);
                        [~,indx] = min(abs(fff-f(1)));
                        mag(ii,jj) = signal(indx);
                        phi(ii,jj) = phase(indx);
                        nf(ii,jj) = noiseFloor(indx);
                        freq(1,counter) = f(1);
                        stimulusSPL(1,counter) = stimSPL(ii,jj);
                        snr(ii,jj) = mag(ii,jj) - nf(ii,jj);
                    elseif stimType(jj) == 1 % for click audio
                        mag(ii,jj) = 20*log10(max(abs((mean(q,2))))/.00002);
                        phi(ii,jj) = nan;
                        nf(ii,jj) = nan;
                        freq(:,counter) = f;
                        stimulusSPL(1,counter) = stimSPL(ii,jj);
                        snr(ii,jj) = nan;

                        if stimSPL(ii,jj) == -15
                            [fff,signal,noiseFloor,phase] = ARLas_fda(q,fs,0.00002);
                            %figure
                            plt(fff,signal)
                            hold on
                            plt(fff,noiseFloor)
                            %keyboard
                        end

                    end

                    cap = [];
                    time = [];
                    signal = [];
                    noise = [];
                    frequency = [];
                end
                counter = counter + 1;
                start = finish + 1;
            end
        end
    end
end
function [Cap,time] = getCAP(chunk2,cuts,fs)
    D1 = chunk2(cuts(1,1):cuts(2,1)-1,:); % pip
    D2 = chunk2(cuts(2,1):end,:);         % invertd pip
    Cap = (D1 + D2)/2;                  % compound action potential
    Cap = Cap * 1000000; % cap in uV
    nSamples = size(Cap,1);           % number of samples in each response
    time = (0:1:nSamples-1)'/fs;      % time vector associated with each CAP waveform
    time = time *1000;                % time in ms                
end
function b = lpf(fs,Fpass) % lowpass filter coefficients
    Fstop = 1.7 * Fpass;    % Stopband Frequency
    Dpass = 0.057501127785;  % Passband Ripple
    Dstop = 0.0001;          % Stopband Attenuation
    flag  = 'scale';         % Sampling Flag
    [N,Wn,BETA,TYPE] = kaiserord([Fpass Fstop]/(fs/2), [1 0], [Dstop Dpass]);
    if mod(N,2)~= 0
        N = N + 1;
    end
    b  = fir1(N, Wn, TYPE, kaiser(N+1, BETA), flag)';
end

% OLD CODE
%--------------------------------------------------------------------------------------------------------
