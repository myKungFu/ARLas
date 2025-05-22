function [] = DW_fastCapAudio_Thev(varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DW_fastCapAudio_Thev(varargin)
%
% Measure CAPs at mulitple frequencies and levels. Get fast audiogram.
% Required functions: DW_analyzeCapLevelSeries.m
%                     DW_plotCapLevelSeries.m
%
% Authors: Shawn S Goodman, PhD & Jeffery T Lichtenhan PhD
% Auditory Research Lab, Dept. of Comm. Sciences & Disorders, The University of Iowa
% Auditory Physiology Lab, Dept. of Otolaryngology, Washington University School of Medicine, St. Louis
% Date: August 12, 2019
% Updated: August 12-15, 2019 -- ssg
% Updated: November 5-6, 2019 -- ssg; Changed to work with Thevenin source calibration.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    doSimulation = 0; % use previously-collected data instead of recording
    obj = varargin{1}; % get the arlas object
    fs = obj.fs;
    tic

    % User controlled parameters -----------------------------------------------
    pipProbe = 'B';             % which 10x probe to present the pips through
    targetReps = 66; % number of stimulus repetitions per condition (make this a multiple of 6!)
            % Suggested Values:
            %  66 for a full set; Matches old method; Will take about 5.7 minutes to run.
            %  12 for a "quick and dirty" look to see if everything is working okay. Takes about 70 seconds to run.

    freqs = [1000,2000,3000,4000,8000,10000,12000,14000,18000,24000]; % stimulus frequencies to test (Hz)
    retl =  [  35,  30,  30,  30,  20,   20,   15,   20,   15,   30   ]; % reference equivalent level for guinea pigs; lowest level expected to hear (our data)

    doIndividualPlots = 1; % turn on (1) and off (0) plots of inididual threshold searches
    plotAudiogram = 1;     % turn on (1) and off (0) plot of final audiogram
    doVerify = 0; % verify the output levels using a reference mic (1) or run in normal mode (0)
    %--------------------------------------------------------------------------

    doRandomize = 1; % randomize the presentation order
    stepSize = 5;    % stepsize in dB
    maxStimSPL = 70; % maximum level to present I moved this from 80 to 75 to compensate for lowering retl by 5 dB, above.

    [inputs,outputs] = hardwareSetup; % an updated version of DWio
    inputG = inputs{3};               % GRASS Bioamp
    if strcmp(pipProbe,'A')
        pipInput = inputs{1};         % ER10xA microphone
        pipOutput = outputs{1};       % ER10xA loudspeakers
    elseif strcmp(pipProbe,'B')
        pipInput = inputs{2};         % ER10xB microphone
        pipOutput = outputs{2};       % ER10xB loudspeakers
    end
    if doVerify == 1
        inputRef = inputs{4};         % GRAS 1/8" reference microphone
    end
    
    [pathName_mic,folderName_mic,fileName_mic] = mostRecentCalibration('mic',pipInput.label,[]);
    if ~isempty(fileName_mic)
        dummy = load([pathName_mic,folderName_mic,fileName_mic]);
        micCorrection = dummy.micCorrection; % microphone correction
    else
        micCorrection = 1; % convolving by this changes nothing
    end
    
    calParams.fmin = 200;   
    calParams.fmax = 32000; 
    calParams.testProbe = pipProbe;
    calParams.usingChannels = 2; % use this channel (1 or 2; only one channel needed for this experiment)

    calParams = ARLas_earCanalRecordings_DW10x(obj,calParams); % run ear canal recordings to get in-situ measurement

    % CREATE STIMULUS MATRIX --------------------------------------------------
    nFreqs = length(freqs);
    dummy = repmat((0:stepSize:100)',1,nFreqs); % vector of levels to test
    nLevels = size(dummy,1);
    Retl = repmat(retl,nLevels,1);
    levels = Retl + dummy;
    levels(levels>maxStimSPL) = 0;
    indx = sum(levels,2) == 0;
    levels(indx,:) = [];
    nLevels = size(levels,1);

    [stim,cuts] = getStim(1000,fs);  % get sample cap stimulus (tone burst)
    stimLength = length(stim)/fs;    % stimulus length (seconds). This includes both polarities and zero padding. Typically ~13 ms
    maxTotalLength = ceil(stimLength * targetReps * nFreqs * nLevels);
    stimSPL = levels;
    [Stim,sortIndx,randIndx,nReps,NReps] = getStimMatrix(calParams,levels,freqs,stim,stimSPL,maxTotalLength,stimLength,fs,doRandomize,doVerify);

    obj.objPlayrec.userInfo = [];                                             % clear old header info  
    obj.objPlayrec.userInfo.nReps = nReps;                                    % load new header info
    obj.objPlayrec.userInfo.stimSPL = stimSPL; 
    obj.objPlayrec.userInfo.freqs = freqs;
    obj.objPlayrec.userInfo.levels = levels;
    obj.objPlayrec.userInfo.randIndx = randIndx;
    obj.objPlayrec.userInfo.sortIndx = sortIndx;
    obj.objPlayrec.userInfo.maxStimSPL = maxStimSPL;
    obj.objPlayrec.userInfo.maxTotalLength = maxTotalLength;
    obj.objPlayrec.userInfo.NReps = NReps;
    obj.objPlayrec.userInfo.cuts = cuts;
    obj.objPlayrec.userInfo.calParams = calParams;

    obj.clearRecList % clear out whatever was used previously for recordings
    obj.setRecList(inputG.ch,inputG.label,inputG.micSens,inputG.gain); % load the Gras bioamp first
    obj.setRecList(pipInput.ch,pipInput.label,pipInput.micSens,pipInput.gain); % load the pip stimulus probe second
    if doVerify == 1
        obj.setRecList(inputRef.ch,inputRef.label,inputRef.micSens,inputRef.gain); % load the reference mic third
    end

    obj.clearPlayList; % clear out whatever was used previously for playback
    obj.setPlayList(Stim,pipOutput.ch(1)); % load output for pip stimuli. This is always used
    obj.setNReps(NReps); % number of stimulus reps to play and record. This will run through the stimulus matrix ONCE!

    tic
    obj.objPlayrec.run  % playback and record
    if obj.killRun
        return  
    end
    toc

    if doSimulation == 1
        pathName = 'C:\myWork\ARLas\Data\DW226\CAPaudio1\';
        stub = 'Ch5_GRASSbio_';
        d = dir([pathName,stub,'*.mat']);
        N = size(d,1);
        counter = 1;
        for ii=5:N
            disp(num2str(ii))
            pause(0.001)
            dummy = load([d(ii).folder,'\',d(ii).name]);
            [Cap,time,fs,timeStamp] = DW_capPreAnalysis(dummy.data,dummy.header);
            [mag(counter,1),nf(counter,1),snr(counter,1),phi(counter,1),cap(:,counter),signal(:,counter),noise(:,counter),frequency(:,counter)] = DW_capAnalysis(Cap,time,fs);
            freq(counter,1) = dummy.header.userInfo.frequency;
            stimulusSPL(counter,1) = dummy.header.userInfo.level;
            counter = counter + 1;
        end
    else
        if doVerify == 1
            [headerAudio,DataAudio] = obj.retrieveData(pipInput.label); % get acoustic pip recordings (audio from ER10X)
            [headerRef,DataRef] = obj.retrieveData(inputRef.label);     % get reference data
            [magA,nfA,snrA,freqA,capA,timeA,stimulusSPLA,phiA,signalA,noiseA,frequencyA] = DW_fastCapAnalyze(headerAudio,DataAudio,1);
            [magR,nfR,snrR,freqR,capR,timeR,stimulusSPLR,phiR,signalR,noiseR,frequencyR] = DW_fastCapAnalyze(headerRef,DataRef,1);
            disp(' ')
            disp('Difference between reference and ipl target:')
            %[freqs'/1000,(max(magR,[],1)' - maxStimSPL +3)]
            [freqs'/1000,(max(magR,[],1)' - maxStimSPL)]
            keyboard
            %[mag,nf,snr,freq,cap,time,stimulusSPL,phi,signal,noise,frequency] = DW_fastCapAnalyze(headerRef,DataRef); % extract the cap waveforms
        else
            [header,Data] = obj.retrieveData(inputG.label);   % get electrical data (cap from gras amplifier)
            [mag,nf,snr,freq,cap,time,stimulusSPL,phi,signal,noise,frequency] = DW_fastCapAnalyze(header,Data,0); % extract the cap waveforms
            % NOTE: the "fast analyze" version above is used (instead of DW_capPreAnalysis) because of the need to de-interleave everything.
        end
    end

    % Pack results into structures and analyze --------------------------------
    f = 1000;
    cap1000 = struct;
    cap1000 = packStructure(cap1000,f,mag,nf,snr,freq,cap,time,stimulusSPL,phi,fs,signal,noise,frequency);
    f = 2000;
    cap2000 = struct;
    cap2000 = packStructure(cap2000,f,mag,nf,snr,freq,cap,time,stimulusSPL,phi,fs,signal,noise,frequency);
    f = 3000;
    cap3000 = struct;
    cap3000 = packStructure(cap3000,f,mag,nf,snr,freq,cap,time,stimulusSPL,phi,fs,signal,noise,frequency);
    f = 4000;
    cap4000 = struct;
    cap4000 = packStructure(cap4000,f,mag,nf,snr,freq,cap,time,stimulusSPL,phi,fs,signal,noise,frequency);
    f = 8000;
    cap8000 = struct;
    cap8000 = packStructure(cap8000,f,mag,nf,snr,freq,cap,time,stimulusSPL,phi,fs,signal,noise,frequency);
    f = 10000;
    cap10000 = struct;
    cap10000 = packStructure(cap10000,f,mag,nf,snr,freq,cap,time,stimulusSPL,phi,fs,signal,noise,frequency);
    f = 12000;
    cap12000 = struct;
    cap12000 = packStructure(cap12000,f,mag,nf,snr,freq,cap,time,stimulusSPL,phi,fs,signal,noise,frequency);
    f = 14000;
    cap14000 = struct;
    cap14000 = packStructure(cap14000,f,mag,nf,snr,freq,cap,time,stimulusSPL,phi,fs,signal,noise,frequency);
    f = 18000;
    cap18000 = struct;
    cap18000 = packStructure(cap18000,f,mag,nf,snr,freq,cap,time,stimulusSPL,phi,fs,signal,noise,frequency);
    f = 24000;
    cap24000 = struct;
    cap24000 = packStructure(cap24000,f,mag,nf,snr,freq,cap,time,stimulusSPL,phi,fs,signal,noise,frequency);
%     f = 32000;
%     cap32000 = struct;
%     cap32000 = packStructure(cap32000,f,mag,nf,snr,freq,cap,time,stimulusSPL,phi,fs,signal,noise,frequency);
    aGram.cap1000 = cap1000;
    aGram.cap2000 = cap2000;
    aGram.cap3000 = cap3000;
    aGram.cap4000 = cap4000;
    aGram.cap8000 = cap8000;
    aGram.cap10000 = cap10000;
    aGram.cap12000 = cap12000;
    aGram.cap14000 = cap14000;
    aGram.cap18000 = cap18000;
    aGram.cap24000 = cap24000;
%     aGram.cap32000 = cap32000;
    aGram.freqs = freqs';
    aGram.thd = [cap1000.THD,cap2000.THD,cap3000.THD,cap4000.THD,cap8000.THD,cap10000.THD,cap12000.THD,cap14000.THD,cap18000.THD,cap24000.THD]';
    aGram.retl = retl';

    % save data ---------------------------------------------------------------
    disp('Saving data...')
    % save full data set
    try
        pathName = obj.objPlayrec.savedFilesPath;
        fileName = ['fastCapAudio_',obj.subjectID,'.mat'];
        fileName = ARLas_saveName(pathName,fileName);
        save([pathName,fileName],'aGram')
    catch
        disp('Error saving full data set!')
    end
    % save just thresholds, for use with subsequent level series
    originalReps = nReps(1,1);
    FRQ = aGram.freqs / 1000;
    THD = aGram.thd;
    try
        fileName = ['CAP_Audiogram_',obj.subjectID,'.mat'];
        fileName = ARLas_saveName(pathName,fileName);
        pathName =  ['C:\myWork\ARLas\Data\',obj.experimentID,'\'];
        save([pathName,fileName],'FRQ','THD','originalReps')
    catch
        disp('Error saving partial data set!')
    end
    % save the audiogram in Excel
    try
        pathNameForXLS = obj.objPlayrec.savedFilesPath; % TEST  TEST TEST TEST
        fileName4XLS = ['CAP_Audiogram_',obj.subjectID,'.xls'];
        fileName4XLS = ARLas_saveName(pathNameForXLS,fileName4XLS);
        writeAudiogramData(FRQ,THD,pathNameForXLS,fileName4XLS)
    catch
        disp('Error saving Excel data!')
    end
    disp('...Finished saving.')
    toc

    % plotting
    if doIndividualPlots == 1 
        plotFastFit12(aGram.cap1000.THD,aGram.cap1000.cap,aGram.cap1000.time,1000,fs,aGram.cap1000.lvl',aGram.cap1000.plotInfo)
        plotFastFit12(aGram.cap2000.THD,aGram.cap2000.cap,aGram.cap2000.time,2000,fs,aGram.cap2000.lvl',aGram.cap2000.plotInfo)
        plotFastFit12(aGram.cap3000.THD,aGram.cap3000.cap,aGram.cap3000.time,3000,fs,aGram.cap3000.lvl',aGram.cap3000.plotInfo)
        plotFastFit12(aGram.cap4000.THD,aGram.cap4000.cap,aGram.cap4000.time,4000,fs,aGram.cap4000.lvl',aGram.cap4000.plotInfo)
        plotFastFit12(aGram.cap8000.THD,aGram.cap8000.cap,aGram.cap8000.time,8000,fs,aGram.cap8000.lvl',aGram.cap8000.plotInfo)
        plotFastFit12(aGram.cap10000.THD,aGram.cap10000.cap,aGram.cap10000.time,10000,fs,aGram.cap10000.lvl',aGram.cap10000.plotInfo)
        plotFastFit12(aGram.cap12000.THD,aGram.cap12000.cap,aGram.cap12000.time,12000,fs,aGram.cap12000.lvl',aGram.cap12000.plotInfo)
        plotFastFit12(aGram.cap14000.THD,aGram.cap14000.cap,aGram.cap14000.time,14000,fs,aGram.cap14000.lvl',aGram.cap14000.plotInfo)
        plotFastFit12(aGram.cap18000.THD,aGram.cap18000.cap,aGram.cap18000.time,18000,fs,aGram.cap18000.lvl',aGram.cap18000.plotInfo)
        plotFastFit12(aGram.cap24000.THD,aGram.cap24000.cap,aGram.cap24000.time,24000,fs,aGram.cap24000.lvl',aGram.cap24000.plotInfo)
       % plotFastFit12(aGram.cap32000.THD,aGram.cap32000.cap,aGram.cap32000.time,32000,fs,aGram.cap32000.lvl',aGram.cap32000.plotInfo)
    end

    if plotAudiogram == 1
        h = figure;
        h = DW_plotAudiogram12(aGram.freqs/1000,aGram.thd,h);
        title({'CAP Audiogram','Norms from N = 23 using A Salt Rig'})
    end
end % EOF: end of experiment file

% internal functions ------------------------------------------------------
function [st] = packStructure(st,f,mag,nf,snr,freq,cap,time,stimulusSPL,phi,fs,signal,noise,frequency)
    indx = find(freq == f);
    st.lvl = stimulusSPL(indx);
    st.mag = mag(indx);
    st.nf = nf(indx);
    st.snr = snr(indx);
    st.phi = phi(indx);
    dummy = cap(:,indx);
    st.wave = dummy;
    st.time = time;
    [THD,cap,time,plotInfo] = fastFit12(st.wave,time,st.lvl,fs,signal,noise,frequency);
    st.THD = THD;
    st.cap = cap;
    st.time = time;
    st.plotInfo = plotInfo;
end
function [THD,cap,time,plotInfo] = fastFit12(CAP,time,LVL,fs,signal,noise,frequency)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [cap,THD,h1,h2] = fastFit12(CAP,time,F,LVL,fs);
%
% Performs threshold analysis of CAP level series
%
% CAP = matrix of mean cap waveforms
%       Should include the cap proper, as well as time afterwards to
%       estimate the noise floor
% time = corresponding matrix of time vectors (ms)
% f = stimulus frequency (Hz)
% LVL = vector of stimulus levels (dB SPL)
% fs = sampling rate (Hz)
%
% Auditory Research Lab, The University of Iowa
% Deptartment of Communication Sciences & Disorders
% The University of Iowa
% Author: Shawn S. Goodman, PhD
% Date: June 16, 2017
% Updated: August 14, 2019 ssg
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    [LVL,indx] = sort(LVL,'descend'); % sort levels from highest to lowest
    CAP = CAP(:,indx);
    level = LVL(:); % stimulus levels as a column vector
    nLevels = length(level); % number of stimulus levels used
    n = size(CAP,1); % number of samples in cap
    if mod(n,2)~=0 % force to be an even number of samples
        CAP = CAP(1:end-1,:);
        n = size(CAP,1);
    end

    startTime = 0;  % earliest time over which to look for cap (ms)
    finishTime = 5; % latest time over which to look for cap (ms)
    startTimeNoise = 7; % time at which to start looking for noise (ms)
    [~,startIndx] = min(abs(time - startTime));   % define starting index
    [~,finishIndx] = min(abs(time - finishTime)); % define ending index
    [~,startIndxNoise] = min(abs(time - startTimeNoise));   % define starting index
    cap = CAP(startIndx:finishIndx,:); % cap is the first part of the recording
    noise = CAP(startIndxNoise:end,:); % noise estimate occurs in time after the cap
    timeNoise = time(startIndxNoise:end); % cut time vector to match the noise
    time = time(startIndx:finishIndx); % cut time vector to match the cap
%     if size(noise,1)< 100
%         nn = floor(size(cap,2)/2);
%         n1 = cap(:,1:nn,:);
%         n2 = cap(:,nn+1:nn*2,:);
%         noise = (n1 + n2)/2;
%         timeNoise = time;
%     end
    rampLen = 0.0005; % ramp on and off the cap and noise measurements (0.5 ms ramp)
    cap = ARLas_ramp(cap,fs,rampLen);
    noise = ARLas_ramp(noise,fs,rampLen);

    b = bpf(fs); % get filter coefficients for a bandpass filter (100-1250 Hz)
    cap = fastFilter(b,cap); % apply the filter to the cap
    noise = fastFilter(b,noise); % apply the filter to the noise

    % METHOD 1 -----
    % Look for time domain peak shifts (peaks should be delayed as level drops)
    for ii=1:nLevels
        Q(:,ii) = xcorr(cap(:,ii),cap(:,1)); % cross correlate all caps with the respone the highest stimulus level
    end
    [pAmp,pIndx] = max(Q); % find the max of the cross correlation functions
    pIndx = pIndx - pIndx(1); % make the delays relative to the high-level comparision
    pAmp = pAmp ./ pAmp(1);   % make the amplitudes relative to the high-level comparison
    pIndx = pIndx'; % transpose to column vectors
    pAmp = pAmp';
    % Range over which the delays are expected to vary. Can tune these parameters
    pUpper = pIndx + 35; % suggest 35 samples
    pLower = pIndx - 2.5; % suggest 5 samples
    pSig = zeros(size(pIndx)); % initialize a significance vector
    pSig(1) = 1; % the first (high-level comparison) is significant
    for ii=2:length(pIndx) % look over the others.
        if pIndx(ii) < pUpper(ii-1) && pIndx(ii) > pLower(ii-1)
            pSig(ii,1) = 1; % peak delay is significant falls within expected values
        else
            break; % once the delay falls outside the significant range, stop looking.
        end
    end
    thdIndxP = max(find(pSig==1)); % find the lowest level that was significant
    thdP = level(thdIndxP); % the stimulus level at cap threshold obtained by this measure

    % METHOD 2 -----
    % In the complex frequency domain, determine whether SNR falls within 95% tolerance region
    lowCut = 100; % lowest frequency to consider (Hz); determined emperically
    highCut = 1250; % highest frequency to consider (Hz); determined emperically
    nfft = fs; % number of points in the fft
    freq = (0:1:nfft-1)'*(fs/nfft); % frequency vector for the fft
    [~,indxLow] = min(abs(freq-lowCut));
    [~,indxHigh] = min(abs(freq-highCut));
    Noise = fft(noise,fs); % fft of noise
    Noise = Noise(indxLow:indxHigh,:); % cut to include only desired frequencies
    Cap = fft(cap,fs); % fft of the cap
    Cap = Cap(indxLow:indxHigh,:); % cut to include only desired frequencies
    R_alpha = toleranceRegion(Noise); % calculate a 95% tolerance region using the noise
    sig = zeros(size(Cap)); % initialize a significance vector
    for ii=1:size(cap,2) % loop over levels
        sig(:,ii) = calculateSignificance(Cap(:,ii),R_alpha); % whether the cap falls inside the tolerance region
    end

    pSig = sum(sig,1)'/size(Cap,1); % proportion of significant points
    cutoff = 0.2; % cutoff for saying cap is present; determined empierically
    Sig = zeros(size(cap,2),1); % initialize to not siginficant
    for ii=1:size(cap,2)
        if pSig(ii) > cutoff
            Sig(ii,1) = 1;
        end
    end
    % loop to find the lowest stimulus used that yielded a signficant result
    ok = 1;
    counter = 1;
    while ok==1 && counter <= length(Sig)
        if Sig(counter) == 1
            counter = counter + 1;
        else
            ok = 0;
        end
    end
    thdIndx = counter - 1;
    if thdIndx == 0 % if there were no significant levels
        thd = nan;
    else 
        thd = level(thdIndx); % the stimulus level at cap threshold obtained by this measure
    end

    % Combine Methods -----
    % No response is coded as NaN.
    THD = mean([thdP,thd]); % the final threshold is the mean of the two methods 
    THDindx = mean([thdIndx,thdIndxP]);
    
    plotInfo.noise = noise;
    plotInfo.timeNoise = timeNoise;
    plotInfo.THDindx = THDindx;
    plotInfo.thdIndxP = thdIndxP;
    plotInfo.thdIndx = thdIndx;
    plotInfo.pIndx = pIndx;
    plotInfo.pAmp = pAmp;
    plotInfo.nLevels = nLevels;
    plotInfo.pSig = pSig;
    plotInfo.Sig = Sig;
    plotInfo.cutoff = cutoff;
end
function [Stim,sortIndx,randIndx,nReps,NReps] = getStimMatrix(calParams,levels,freqs,stim,stimSPL,maxTotalLength,stimLength,fs,doRandomize,doVerify)
    % MAKE STIMULUS MATRIX
    maxTotalReps = floor(maxTotalLength / stimLength); % maximum possible number of stimulus presentations in allotted time
    % divide up the reps equally among all frequencies:
    nLevels = size(levels,1);
    nFreqs = length(freqs);
    nReps = floor(maxTotalReps / (nFreqs*nLevels));
    nReps = nReps - mod(nReps,6); % constrain reps to be muliples of 6
    nReps = repmat(nReps,nLevels,nFreqs);
    nReps = nReps .* (stimSPL ~= 0);                                          % do not test levels that are too high
    %nReps = floor(nReps * (maxTotalReps/NReps));                              % re-allocate the number of reps that are left over, so as to fill up the time allotted time
    NReps = sum(sum(nReps));                                                  % re-calculate the total number of reps to be tested
    if NReps > maxTotalReps                                                   % check to make sure the total number of possible reps is not exceeded
        warning('Requested number of repetitions exceeds the number possible.')
        keyboard
    end
    Stim = zeros(length(stim),NReps);                                         % initialize the stimulus matrix
    counter = 1;
    for ii=1:nLevels
        for jj=1:nFreqs
            f = freqs(jj);               % current frequency
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
            [stimScaled,a,errors] = ARLas_applyISC(calParams.isc,desiredLvl,type,f,stim);
            % a = applyCal(pathName_speakerCal,fileName_speakerCal,fileName_earBarCal,f,1,desiredLvl); % level multiplier 
            % stim = stim * abs(a);
            % The stimulus is in peak amplitude, but we want the target to
            % be in rms. Therefore, add 3 dB:
            stim = stimScaled * (1./sqrt(.5));
            nn = nReps(ii,jj);
            if nn > 0
                for kk=1:nn
                    Stim(:,counter) = stim;
                    counter = counter + 1;
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
end
function [stim,cuts] = getStim(f,fs) % get the cap stimulus (tone burst)
    % This stim designed June 8, 2017
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
function [R_alpha] = toleranceRegion(noiseSamples)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% R_alpha = toleranceRegion(noiseSamples);
%
% Calculate the ellipse that contains 95% of the data points.
%
% noiseSamples = Samples of noise in the frequency domain. Values are in
%            complex rectangular form (real and imaginary parts).
%            noiseSamples is a matrix of size n by m, where each n is a
%            freqeuency between 100 and 3000 Hz, and each m is a different
%            stimulus level
% R_alpha = the ellipse that contains 100*(1-alpha) percent of the data,
%            alpha = 0.05.
%
% Authors: Shawn S. Goodman, PhD & Ian B. Mertes, PhD
%          Dept. of Communication Sciences & Disorders
%          University of Iowa
% Date: July 13, 2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
deltaHat = noiseSamples(:); % force to column vector
barDeltaHat = mean(deltaHat,1); % mean of the bootstrap estimates
deltaHat = deltaHat - barDeltaHat; % remove the mean to center around zero
deltaHat = [real(deltaHat),imag(deltaHat)]; % convert to two-column matrix
K = size(deltaHat,1);
S = (1/(K-1))*(deltaHat'*deltaHat); % Eq. 3 in manuscript
[lambda,upsilon] = eig(S); 
m = 2000; % number of samples in the ellipse
theta = linspace(0,2*pi,m); % Linearly-spaced vector 
C = [cos(theta);sin(theta)]; % Create unit circle
R = (lambda*sqrt(upsilon)*C)'; % Eq. 4 in manuscript
alpha = 0.05; % significance level
df = 2; % degrees of freedom
invChiSq = 5.9915; % Inverse of chi-squared cumulative distribution function when alpha = 0.05 and df = 2.
                   % Note that invChiSq can be computed for other values of alpha and df using 
                   % the chi2inv.m function of the MATLAB Statistics Toolbox, 
                   % where invChiSq = chi2inv(1-alpha,df);
Ra = R * sqrt(invChiSq); % Eq. 5 in manuscript
Ra = Ra(:,1) + 1i*Ra(:,2); % convert back to complex vector
R_alpha = Ra + barDeltaHat; % put back to the original location
end
function [sig] = calculateSignificance(signal,R_alpha)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [sig,deltaNorm] = calculateSignificance(signal,R_alpha);
%
% Determine whether delta is significantly different from zero.
%
% signal = fft of the cap signal. Each sample is a complex frequency value
%           from 100 to 1250 Hz.
% R_alpha = the elliptical tolerance interval that contains 100*(1-alpha) 
%            percent of the data
%
% Authors: Shawn S. Goodman, PhD & Ian B. Mertes, PhD
%          Dept. of Communication Sciences & Disorders
%          University of Iowa
% Date: July 13, 2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
signal = signal(:); % force to be a column vector
ns = length(signal); % number of samples
sig = zeros(1,ns);
for ii=1:ns
    phiR = angle(R_alpha); % compute phase values 
    phiS = angle(signal(ii));
    magR = abs(R_alpha); % compute magnitude values
    magS = abs(signal(ii));
    [~,indx] = min(abs(phiR - phiS)); % find the phase on the ellipse closest to 
        % the phase of delta
    if magS > magR(indx)
        sig(1,ii) = 1; % if the magnitude of the data point is larger than the ellipse
    else
        sig(1,ii) = 0; % if the magnitude of the data point falls within the ellipse
    end
end
end
function b = bpf(fs)
%bandpass filter
Fstop1 = 50;              % First Stopband Frequency
Fpass1 = 100;             % First Passband Frequency
Fpass2 = 1250;            % Second Passband Frequency
Fstop2 = 1500;            % Second Stopband Frequency
Dstop1 = 0.001;           % First Stopband Attenuation
Dpass  = 0.057501127785;  % Passband Ripple
Dstop2 = 0.0001;          % Second Stopband Attenuation
flag   = 'scale';         % Sampling Flag
[N,Wn,BETA,TYPE] = kaiserord([Fstop1 Fpass1 Fpass2 Fstop2]/(fs/2), [0 ...
                             1 0], [Dstop1 Dpass Dstop2]);
if mod(N,2)~= 0
    N = N + 1;
end
b  = fir1(N, Wn, TYPE, kaiser(N+1, BETA), flag)';
end
function [] = plotFastFit12(THD,cap,time,f,fs,level,plotInfo)
    noise = plotInfo.noise;
    timeNoise = plotInfo.timeNoise;
    THDindx = plotInfo.THDindx;
    thdIndxP = plotInfo.thdIndxP;
    thdIndx = plotInfo.thdIndx;
    pIndx = plotInfo.pIndx;
    pAmp = plotInfo.pAmp;
    nLevels = plotInfo.nLevels;
    pSig = plotInfo.pSig;
    Sig = plotInfo.Sig;
    cutoff = plotInfo.cutoff;

    h1 = figure; %(555); % plot raw cap waveforms
    hold off
    offset = (max(cap(:,1)) - min(cap(:,1))) * 0.1;
    for ii=1:size(cap,2)
        ytick(ii,1) = -((ii-1)*offset);
        if ii-THDindx <= 0.5
            if ii-THDindx <= 0.5 && ii-THDindx >= -0.5
                LW = 2;
            else
                LW = 0.5;
            end
            plot(time,cap(:,ii)-((ii-1)*offset),'Color',[1 0 0],'LineWidth',LW)
            hold on
        else
            plot(time,cap(:,ii)-((ii-1)*offset),'Color',[0 0 0])
            hold on
        end
    end
    plot(timeNoise,noise-((ii)*offset),'Color',[.7 .7 .7])
    hold on
    xlabel('Time (ms)','FontSize',12)
    ylabel('Stimulus SPL (dB)','FontSize',12)
    xlim([time(1),time(end)])
    title(['CAP ',num2str(f/1000),' kHz: Threshold = ',num2str(THD),' dB SPL'])
    ytick = [ytick;-(ii)*offset];
    lvl = [level;0];
    for ii=1:length(ytick)
        decimate = 2; % restrict to this number of digits
        ytickstr(ii,1) = cellstr(num2str(lvl(ii),decimate)); % make tick labels
    end
    set(gca,'YTick',flipud(ytick),'YTickLabel',flipud(ytickstr));

    [amp,indx] = max(cap(:,1));
    for ii=1:size(cap,2)
        tP1(ii,1) = time(indx)+(pIndx(ii)/(fs/1000));
        aP1(ii,1) = amp*pAmp(ii)-((ii-1)*offset);
    end
    [amp,indx] = min(cap(:,1));
    for ii=1:size(cap,2)
        tP2(ii,1) = time(indx)+(pIndx(ii)/(fs/1000));
        aP2(ii,1) = amp*pAmp(ii)-((ii-1)*offset);
    end

    if thdIndxP > 1
        plot(tP1(1:thdIndxP),aP1(1:thdIndxP),'*-','Color',[0 0 1]) % P1 of CAP
        plot(tP2(1:thdIndxP),aP2(1:thdIndxP),'*-','Color',[0 1 0]) % N1 of CAP
    end

    if thdIndxP <= nLevels
        for ii=thdIndxP+1:nLevels
            plot(tP1(ii),aP1(ii),'o','Color',[0 0 1])
            plot(tP2(ii),aP2(ii),'o','Color',[0 1 0])
        end
    end

    h2 = figure; %(556); % plot percentage of significant points for each cap waveform
    hold off
    for ii=1:length(pSig)
        if Sig(ii,1) == 1
            plot(level(ii),pSig(ii),'*','Color',[0 0 1])
            hold on
        else
            plot(level(ii),pSig(ii),'*','Color',[0 0 0])
            hold on
        end
    end
    plot(level(1:round(THDindx)),pSig(1:round(THDindx)),'-','Color',[0 0 1])
    hold on
    line([min(level)-5,max(level)+5],[cutoff,cutoff],'LineWidth',1,'LineStyle','--','Color',[.7 .7 .7])
    line([THD,THD],[0,1],'LineWidth',2,'LineStyle','-','Color',[1 0 0])
    xlim([min(level)-5,max(level)+5])
    xlabel('Stimulus SPL (dB)','FontSize',12)
    ylabel('Proportion of CAP > Noise Floor','FontSize',12)
    title(['CAP ',num2str(f/1000),' kHz: Threshold = ',num2str(THD),' dB SPL'])
    set(h2,'Position',[403 883 560 420]);
end
function [] = writeAudiogramData(FRQ,THD,pathName,fileName4XLS)
    FRQ = FRQ(:); % force to column vectors
    THD = THD(:);
    warning off
    sheet = 'audiogram'; % sheet name is current frequency
    range = 'A1';
    status(1,1) = xlswrite([pathName,fileName4XLS],{'Frequency (kHz)'},sheet,range);
    range = 'A2';
    status(2,1) = xlswrite([pathName,fileName4XLS],FRQ,sheet,range);
    range = 'B1';
    status(3,1) = xlswrite([pathName,fileName4XLS],{'CAP thd (dB SPL)'},sheet,range);
    range = 'B2';
    status(4,1) = xlswrite([pathName,fileName4XLS],THD,sheet,range);
    warning on
end

% old code ----------------------------------------------------------------
% function [multiplier] = applyCal(pathName_speakerCal,fileName_speakerCal,fileName_earBarCal,f,ch,targetAmp)
%     % f = stimulus frequency (Hz)
%     % ch = probe output (1 or 2; probe is 1, suppressor is 2)
%     % targetAmp = desired output level in dB SPL rms
%     dummy = load([pathName_speakerCal,fileName_speakerCal]);
%     fo = interp1(dummy.d.freq,dummy.d.FO(:,ch),f,'pchip');
%     dOut = targetAmp - fo;
%     dummy2 = load([pathName_speakerCal,fileName_earBarCal]);
%     c2 = interp1(dummy2.d.freq,dummy2.d.correction2(:,ch),f,'pchip');
%     dOut = dOut - c2;
%     multiplier = 10^(dOut/20);
%     if multiplier >= 1 % if cannot achieve target output
%         hit = 1;
%         multiplier = .99;
%         disp('hit! hit! hit! hit! hit! hit! hit! hit! hit! hit!')
%     else
%         hit = 0;
%     end
% end

% ------------------- start hack
% % this is a bit of a hack, but by reducing the number of averages at higher
% % levels (still must be multiples of 6) we can save ~20% recording time.
% mask = nReps ~= 0;
% strip = ones(1,nFreqs);
% nReps = [strip*66;strip*66;strip*60;strip*60;strip*54;strip*54;strip*48;strip*48;strip*42;strip*42;strip*36;strip*36;strip*36] .* mask;
% NReps = sum(nReps(:));
% ------------------- end hack

% NEW CALIBRATION re. SHAWN & LAURA VISIT ON 8/9/17 --------
% pathName_speakerCal = 'C:\myWork\ARLas\Peripheral\calibrations\speakerCals\'; % location of loudspeaker calibrations
% fileName_speakerCal = 'speakerCal_10xB_1.mat';
% fileName_earBarCal = 'speakerCalEarBar_10xB_1.mat';            
% pathName_micCal = 'C:\myWork\ARLas\Peripheral\calibrations\micCals\'; % location of microphone calibrations
% fileName_micCal = 'micCal_10xB_avg12-15PLASTIC.mat';

%[Stim,sortIndx,randIndx,nReps,NReps] = getStimMatrix(calParams,levels,freqs,stim,stimSPL,maxTotalLength,stimLength,fs,doRandomize,pathName_speakerCal,fileName_speakerCal,fileName_earBarCal);
%obj.objPlayrec.userInfo.calFile_speaker = fileName_speakerCal;
%obj.objPlayrec.userInfo.calFile_earBar = fileName_earBarCal;
%obj.objPlayrec.userInfo.calFile_mic = fileName_micCal;
