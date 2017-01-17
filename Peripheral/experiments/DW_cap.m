function [] = DW_cap(varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DW_cap(varargin)
%
% Measure CAP
%
% Auditory Research Lab, The University of Iowa
% Deptartment of Communication Sciences & Disorders
% The University of Iowa
% Author: Shawn S. Goodman, PhD
% Date: October 27, 2016
% Last Updated: November 29, 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

obj = varargin{1}; % get the arlas object


% USER CAN ADUST PARAMETERS HERE: -----------------------------------------
freqs = [1000,2000,3000,4000,6000,8000,12000,16000]; % stimulus frequencies to test

% ----- Load Electrical Input:
label = 'GRASS'; % GRASS CP511 
micSens = 1; % sensitivity in V/Pa
gain = 80; % amplifier gain in dB (10,000x = 80 dB)
ch = 2; % CHANGE THIS, IF NEEDED
obj.setRecList(ch,label,micSens,gain);
% make a structure to save
input1.label = label;
input1.micSens = micSens;
input1.gain = gain;
input1.ch = ch;
% label = 'GRAS8'; % GRAS 1/8-inch pressure mic, SN:266921
% micSens = 0.00092; % sensitivity in V/Pa
% gain = -0.25; % quarter inch GRAS pre-amp gain
% ch = 5; % CHANGE THIS, IF NEEDED
% obj.setRecList(ch,label,micSens,gain);
% % make a structure to save
% input1.label = label;
% input1.micSens = micSens;
% input1.gain = gain;
% input1.ch = ch;

% ----- Load Secondary Input:
label = 'ER10C_AOS';
micSens = 0.05;
gain = 20;
ch = 1; % CHANGE THIS, IF NEEDED!
obj.setRecList(ch,label,micSens,gain);
% make a structure to save
input2.label = label;
input2.micSens = micSens;
input2.gain = gain;
input2.ch = ch;

% Specify output:
% make a structure to save
output1.label = 'ER10C_AOS';
output1.gain = 20;
output1.ch = 1; % CHANGE THIS< IF NEEDED (channel on the sound card)
output1.receiver = 1; % which receiver on the ER10C is being used
% -------------------------------------------------------------------------

audiogram = figure; 
title('CAP Audiogram')

fs = obj.fs; % sampling rate being used
nFreqs = length(freqs); % number of frequencies to measure
fullOutSPL = getCal(freqs,fs); % get the calibration values

obj.objPlayrec.nReps = 24; % number of times to play stimulus
snrCriterion = 9; % snr in dB; set critierion for present CAP
startLvl = 80; % desired starting amplitude
thd = 0;
fs = obj.fs;
for ii=1:nFreqs % loop over frequencies
    f = freqs(ii); % get the current frequency
    stim = getStim(f,fs); % get full out stimulus
    foSPL = fullOutSPL(ii); % get the full out value
    
%     if ii>1
%         startLvl = thd + 30;
%         if startLvl > 90
%             startLvl = 90;
%         end
%     end
    
    disp(' ')
    disp(' ')
    disp(['Testing CAP: ',num2str(f),' Hz.'])
    [levels,p2p,snr,cap,time,thd] = findThreshold(obj,f,stim,foSPL,startLvl,output1,snrCriterion); % get the threshold
    THD(1,ii) = thd;
    FRQ(1,ii) = f/1000;
    
    % write data ------------------------
    disp(' Writing data to XLS file...')
    fileName = ['CAP_thresholds_',obj.subjectID,'.xls'];
    pathName = obj.objPlayrec.savedFilesPath;    
    warning off
    if isempty(pathName)
        disp('Empty path name!')
        keyboard
    end
    sheet = num2str(f); % sheet name is current frequency
        range = 'A1';
        status(1,1) = xlswrite([pathName,fileName],{'CAP thd (dBSPL)'},sheet,range);
        range = 'A2';
        status(2,1) = xlswrite([pathName,fileName],thd,sheet,range);          
        range = 'B1';
        status(3,1) = xlswrite([pathName,fileName],{'stimLvl (dBSPL)'},sheet,range);
        range = 'B2';
        status(4,1) = xlswrite([pathName,fileName],levels,sheet,range);
        range = 'C1';
        status(5,1) = xlswrite([pathName,fileName],{'CAP p2p amp (mV)'},sheet,range);
        range = 'C2';
        status(6,1) = xlswrite([pathName,fileName],p2p*1000,sheet,range);
        range = 'D1';
        status(7,1) = xlswrite([pathName,fileName],{'CAP snr (dB)'},sheet,range);
        range = 'D2';
        status(8,1) = xlswrite([pathName,fileName],snr,sheet,range);
        range = 'E1';
        warning off
        status(9,1) = xlswrite([pathName,fileName],{'time (ms)'},sheet,range);
        range = 'E2';
        status(10,1) = xlswrite([pathName,fileName],time,sheet,range);
        range = 'F1';
        status(11,1) = xlswrite([pathName,fileName],{'CAP waveform (V)'},sheet,range);
        range = 'F2';
        status(12,1) = xlswrite([pathName,fileName],cap,sheet,range);
    if min(status) == 1
        disp(' ...Successful data write.')
    else
        disp(' ...Error writing data!')
        status
    end
    warning on
    
    % plot audiogram
    figure(audiogram)
    plot(FRQ,THD,'*-b')
    xlim([.5 40])
    grid on
    xlabel('Frequency (kHz)')
    ylabel('Threshold (dB SPL)')
    title('CAP Audiogram')
    
end

% save audiogram as a jpeg
figFileName = ['CAPaudiogram_',obj.subjectID,'.tif'];
saveas(audiogram,[pathName,figFileName])

% save the audiogram in Excel
FRQ = FRQ';
THD = THD';
warning off
sheet = 'audiogram'; % sheet name is current frequency
range = 'A1';
status(13,1) = xlswrite([pathName,fileName],{'Frequency (kHz)'},sheet,range);
range = 'A2';
status(14,1) = xlswrite([pathName,fileName],FRQ,sheet,range);
range = 'B1';
status(15,1) = xlswrite([pathName,fileName],{'CAP thd (dB SPL)'},sheet,range);
range = 'B2';
status(16,1) = xlswrite([pathName,fileName],THD,sheet,range);
warning on
end % end of experiment file



% internal functions ------------------------------------------------------
function [levels,p2p,snr,cap,timeCut,thd] = findThreshold(obj,f,stim,foSPL,startLvl,output,snrCriterion)
    outputCh = output.ch; % get the output channel
    fs = obj.fs;
%     a = getAmp(startLvl,foSPL); % get the multiplier to give desired output level
%     obj.setPlayList(stim*a,outputCh);% Load Output
    t0 = 480; % time zero, the onset of the tone burst, is time zero
    time = (((0:1:length(stim)-1)'-t0)/fs)*1000; % time in ms
    doPlot = 110; % plot to this figure

    % loop to find threshold
    lvl = startLvl;
    stepSizeDn = 10; % if response is present, go down this many dB
    stepSizeUp = 5; % if response is absent, go up this many dB
    maxLvl = 96; % dB SPL
    minLvl = 10; % dB SPL
    maxReversals = 8; % maximum number of reversals (up to down or vice versa) allowed
    maxEdges = 4; % maximum number of presentations at the max or min levels
    maxRecordings = 100; % maximum numbe rof recordings allowed
    
    done = 0; % when finished, this will change to 1.
    direction = -1; % code negative for decreasing
    nReversals = 0; % initilze
    nEdges = 0;
    counter = 0; % count the recordings
    while done == 0
        
        a = getAmp(lvl,foSPL); % get the multiplier to give desired output level
        obj.setPlayList(stim*a,outputCh);% Load Output
        
        disp(['   testing CAP: ',num2str(lvl),' dB SPL.'])
        counter = counter + 1;
        levels(counter,1) = lvl;
        
        
        
        obj.objPlayrec.run % playback and record
        ok = obj.checkForErrors;
        if ~ok
           return
        end    
        %[header,data] = obj.retrieveData(['Ch',num2str(outputCh)]); % get raw data
        [header,data] = obj.retrieveData('Ch2'); % get raw data
        [p2p(counter,1),snr(counter,1),cap(:,counter),timeCut] = analyzeCAP(data,time,f,startLvl,fs,doPlot,snrCriterion); % analyze cap data        
        
        % adjust next level based on whether a response was present
        if snr(counter,1) >= snrCriterion % if response present
            lvl = lvl - stepSizeDn; % decrease level ----------
            if direction > 0 % if current direction is up
                nReversals = nReversals + 1; % code a reversal
                direction = -1; % set current direction as down
                %disp(['Reversal!',num2str(snr(counter,1))])
            else % current direction is down
                % no change is necessary
            end
            if lvl < minLvl
                lvl = minLvl;
                nEdges = nEdges + 1;
            end
        else % if response absent
            lvl = lvl + stepSizeUp; % increase lvl -----------
            if direction < 0 % if current direction is down
                nReversals = nReversals + 1; % code a reversal
                direction = +1; % set current direction as up
                %disp(['Reversal!',num2str(snr(counter,1))])
            else % current direction is up
                % no change is necessary
            end
            if lvl > maxLvl
                lvl = maxLvl;
                nEdges = nEdges + 1;
            end            
        end
        % check to see if end condition has been reached
        if nReversals > maxReversals
            done = 1;
            disp('     Max reversals reached.')
        end
        if nEdges > maxEdges
            done = 1;
            disp('     Max edge presentations reached.')
        end
        if counter > maxRecordings
            done = 1;
            disp('     Max recordings reached.')
        end            
    end % while loop finishes here
    disp(['   Number of reversals = ',num2str(nReversals)])
    thd = mean(levels(end-maxReversals+1:end)); % CAP threshold
end

function [p2p,snr,cap,time] = analyzeCAP(data,time,f,lvl,fs,doPlot,snrCriterion)
    totalSamples = size(data,1); % number of total samples in buffer 
    nSamples = totalSamples / 2; % number in each half
    d1 = data(1:nSamples,:); % first half
    d2 = data(nSamples+1:end,:); % second half
    Cap = d1+d2; % compound action potential
    %cm = d1-d2; % cochlear microphonic
    % cut down to look only at the cap itself
    start = 500; % start sample
    finish = 1100; % finish sample
    Cap = Cap(start:finish,:);
    time = time(start:finish);
    rampSamples = 100; % apply onset/offset ramps
    h = hann(rampSamples * 2);
    h = h(1:rampSamples);
    H = repmat(h,1,size(Cap,2));
    Cap(1:100,:) = Cap(1:100,:) .* H;
    Cap = flipud(Cap);
    Cap(1:100,:) = Cap(1:100,:) .* H;
    Cap = flipud(Cap);
    Cap = ARLas_artifactReject(Cap); % two-pass artifact rejection
    Cap = ARLas_artifactReject(Cap);
    cap = mean(Cap,2);
    % frequency-domain analysis
    ref = 1; % reference is 1 Volt
    [frequency,signal,noiseFloor] = ARLas_fda(Cap,fs,ref);
    % find root mean square SNR
    snr = (10.^(signal/20) ./ 10.^(noiseFloor/20));
    snr = snr(2:8);
    snr = sqrt(mean(snr.^2));
    snr = 20*log10(snr);
    % time-domain analysis, peak-to-peak amplitude
    if snr < 6
        p2p = 0;
    else
        p2p = max(cap) - min(cap);
    end
    if doPlot > 0 % if a plot is requested
        figure(doPlot)
        subplot(2,1,1) % time subplot
        plot(time,Cap*1000,'Color',[.7 .7 .7]) % plot all traces
        hold on
        plot(time,mean(cap*1000,2),'r') % plot mean trace
        xlim([time(1) time(end)])
        hold off
        xlabel('Time (ms)')
        ylabel('Amplitude (mV)')
         title(['CAP:   ',num2str(f),' Hz, ',num2str(lvl),' dB SPL     p2pAmp = ',num2str(p2p*1000),' mV'])
        subplot(2,1,2) % frequency subplot
        plot(frequency,signal,'r')
        hold on
        plot(frequency,noiseFloor,'k')
        plot(frequency,noiseFloor+snrCriterion,'k:')
        xlim([0 2000])
        hold off
        xlabel('Frequency (Hz)')
        ylabel('Magnitude (dB re:1V)')
    end
end

function [a] = getAmp(desiredLvl,fullOut) % return the desired amplitude as a multiplier
    d = desiredLvl - fullOut; % difference in dB
    a = 10^(d/20); % muliplier
    if a > 1
        disp('ERROR: desiredLvl cannot exceed full output.')
        a = 1;
    end
end

function [fullOutSPL] = getCal(freqs,fs) % get calibration values
    pathName = 'C:\myWork\ARLas\Peripheral\calibrations\DW_AAcal\'; % location of calibration files
    d = dir(pathName); % get the names in the directory
    calFile = d(end).name; % get the name of the most recent chirp calibration
    dummy = load(calFile); % load the calibration
    fullOutSPL = dummy.fullOutSPL; % the dB SPL produced using the maximum stimulus amplitude (scaled +/- 1)
    frequency = dummy.frequency; % associated frequency vector
    % smoothe the calibration by making a relatively low-order FIR filter
    M = round(.02 * fs); % filter order is 2% of the sampling rate
    if mod(M,2)~= 0 % filter order must be an even number
        M = M + 1;
    end
    n = length(fullOutSPL);
    X = [fullOutSPL;flipud(fullOutSPL(2:end-1))]; % add aliased portion of the DFT
    X = 10.^(X/20); % convert to linear units
    x = ifft(X); % create in impulse response
    N = length(x);
    x = [x(N/2+2:end);x(1:N/2+1)]; % make filter causal
    x = x((N/2)-(M/2):(N/2)+(M/2)); % take only the desired center portion
    w = hann(M+1); % create hann window
    x = x .* w; % window the filter to smooth edges
    X = fft(x,N); % put back into the frequency domain
    X = 20*log10(abs(X)); % put back into decibels
    X = X(1:n);
    %figure
    %plot(fullOutSPL,'b')
    %hold on
    %plot(X,'r')
    %fullOutSPL = meanSmoother(fullOutSPL,10); this is another way to smooth
    nFreqs = length(freqs);
    foSPL = zeros(nFreqs,1); % initialize full out SPL values for the freqs being tested
    for jj=1:nFreqs
        [~,indx] = min(abs(freqs(jj)-frequency));
        foSPL(jj,1) = fullOutSPL(indx);
    end
    fullOutSPL = foSPL;
end

function [stim] = getStim(f,fs) % get the cap stimulus (tone burst)
    len = 0.01; % desired stimulus length (s)
        nSamples = round(len * fs); % number of samples in stimulus
        if mod(nSamples,2) ~= 0 % force to be even number
            nSamples = nSamples + 1;
        end
    t = (0:1:nSamples-1)'/fs;
    stim = sin(2*pi*f*t);
    rampLen = 0.001;
    rampSamples = round(rampLen * fs); % number of samples in stimulus
        if mod(rampSamples,2) ~= 0 % force to be even number
            rampSamples = rampSamples + 1;
        end
    h = hann(rampSamples*2);
    h = h(1:rampSamples);
    stim(1:rampSamples,:) = stim(1:rampSamples,:) .* h;
    stim = flipud(stim);
    stim(1:rampSamples,:) = stim(1:rampSamples,:) .* h;
    stim = flipud(stim);
    padLen = 0.005;
    padSamples = round(padLen * fs); % number of samples in stimulus
        if mod(padSamples,2) ~= 0 % force to be even number
            padSamples = padSamples + 1;
        end
    zpad = zeros(padSamples,1);
    stim = [zpad;stim;zpad;zpad];
    stim = [stim;-stim];
end

