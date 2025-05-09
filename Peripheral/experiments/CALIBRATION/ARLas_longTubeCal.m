function [] = ARLas_longTubeCal(varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ARLas_longTubeCal(varargin)
%
% Perform incident pressure calibration in a long, "reflectionless" tube.
% Saves results into a structure called LTC (Long Tube Calibration),
% located in 'C:\myWork\ARLas\Peripheral\calibrations\LTCals\' in a
% sub-folder under the date of the calibration.
%
% Author: Shawn Goodman, PhD
% Date: June 7, 2021
% Last Updated: June 7, 2021 -- ssg
% Last Updated: June 9, 2021 -- ssg
% Last Updated: June 10, 2021 -- ssg -- changed the saved waveform to the filter impulse
% Last Updated: September 8, 2021 -- ssg 
% Last Updated: September 13, 2021 -- ssg -- cleaned up code; calculate
%                                            surge impedance and save value
% Last Updated: March 15, 2022 -- ssg -- added minor comments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    obj = varargin{1}; % get the arlas object
    V = '13SEP2021'; % this is the current version number of this program
    
    %------ USER MODIFIABLE PARAMETERS ----------------------------------------
    probe = 'B'; % which ER10X probe to use (A = left; B = right)
    cavityDiameter = 0.9; % long tube cavity diameter (cm) 
                     % IMPORTANT! Make sure this is the correct value for
                     % the tube that you are using!!!
    cavityTemperature = 30; % degrees Celcius--read this off of the 10X screen
    
    fmin = 100; % minimum frequency to include
    fmax = 20000; % maximum frequency to include
    applyMicCorrection = 1; % turn on (1) and off (0) mic correction
    nReps = 32; % Number of stimulus repetitions to record
    level = 0; % stimulus level in dB re: full output (must be zero or negative)
    doPlot = 1; % turn on (1) and off (0) plotting
    %--------------------------------------------------------------------------

    [inputs,outputs] = hardwareSetup;
    if strcmp(probe,'A')
        probeInput = inputs{1};         % ER10xA microphone
        probeOutput = outputs{1};       % ER10xA loudspeakers
    elseif strcmp(probe,'B')
        probeInput = inputs{2};         % ER10xB microphone
        probeOutput = outputs{2};       % ER10xB loudspeakers
    end
    
    disp(' '); disp(' '); disp(' '); disp(' ')
    alertTxt = {'Starting ARLas Long Tube Calibration Routine.'
         ['  Calibrating probe ',probe,'.']
        };
    nn = size(alertTxt,1);
    for ii=1:nn
        cprintf([0,0,.4],[alertTxt{ii},'\n']);
    end

    % -------------------------------------------------------------------------
    % get most recent calibration files (mic only)
    [pathName_mic,folderName_mic,fileName_mic] = mostRecentCalibration('mic',probeInput.label,[]);
    if applyMicCorrection == 1
        if ~isempty(fileName_mic)
            dummy = load([pathName_mic,folderName_mic,fileName_mic]);
            micCorrection = dummy.micCorrection; % microphone correction
        else
            micCorrection = 1; % convolving by this changes nothing
        end
    else
        micCorrection = 1; % convolving by this changes nothing
    end
    %--------------------------------------------------------------------------
    
    % calculate tube surge impedance
    d = cavityTemperature - 26.85; % from Keefe (1984)
    C = 3.4723e4 * (1 + 0.00166 * d); % speed of sound, adjusted for temperature, from Keefe (1984)
    rho = 1.1769e-3 * (1 - 0.00335 * d);
    r = cavityDiameter ./ 2; % radius in cm
    z0 = (rho * C) ./ (pi * r.^2); % z0 = Ro; characteristic impedance (z0) is constant & real; 
    
    % create stimuli ----------------------------------------------------------
    fs = obj.fs; % get the system sampling rate
    len = 0.1; % stimulus length (sec) -->use 100 ms
    stimN = round(len * fs); % number of samples in stimulus
    stim = zeros(stimN,1); % initialze to all zeros
    stimOffset_sec = 0.003; % stimulus offset from time zero (sec)
    stimOffset_samples = round(stimOffset_sec * fs);
    stim(stimOffset_samples,1) = 1; % impulse click, full output
    if level > 0
        error('variable level must be negative or zero.')
    end
    if level == 0
        multiplier = 0.99;
    else
        multiplier = 10^(level/20); % convert from dB to linear
    end
    stim = stim * multiplier; % scale output to desired level
    
    % include instructions here... -------------------
    alertTxt = {['  Place ER10X probe ',probeInput.label,' in the end of a long lossy tube.']};
    nn = size(alertTxt,1);
    for ii=1:nn
        cprintf([0,0,.4],[alertTxt{ii},'\n']);
    end
    disp(' ')        

    nChannelsOut = length(probeOutput.ch); % number of output channnels being calibrated (typically 2)
    for kk=1:nChannelsOut % loop over output channels; test output channels one at a time
        alertTxt = {['     Calibrating loudspeaker ',num2str(kk)]};
        nn = size(alertTxt,1);
        for ii=1:nn
            cprintf([0,0,.4],[alertTxt{ii},'\n']);
        end
    
        % COLLECT THE DATA ------------------------------------------------
        obj.clearRecList % clear out whatever was used previously 
        obj.setRecList(probeInput.ch,probeInput.label,probeInput.micSens,probeInput.gain);
        obj.setNReps(nReps); % number of times to play stimulus
        obj.setFilter(0); % note: this is a highpass filter with a 75 Hz cutoff frequency.

        obj.clearPlayList % clear out whatever was used previously
        obj.setPlayList(stim,probeOutput.ch(kk)); % load the currently tested output channel

        obj.objPlayrec.userInfo = []; % clear out previous, non-structure array info
        obj.objPlayrec.userInfo.Lvl = level;
        obj.objPlayrec.userInfo.fs = fs;
        obj.objPlayrec.userInfo.stim = stim;
        obj.objPlayrec.userInfo.caps_version = V;

        obj.objPlayrec.run % playback and record
        if obj.killRun
           return
        end

        % RETRIEVE DATA ---------------------------------------------------
        [header,Data] = obj.retrieveData(probeInput.label); % get raw data
        if applyMicCorrection == 1
            Data = fastFilter(micCorrection,Data);
        end
        
        [LTRecording,stimulusFiltered,Data,time] = LTCleanData(Data,fs,stim,header);
        if kk == 1
            [fo_spl,phi,freq,Amax,pRef,nf] = XferFunction(Data,stimulusFiltered,fs,fmin,fmax);
            LTC1 = buildStruct(stim,time,stimulusFiltered,LTRecording,Data,probeInput,probeOutput,...
                Amax,pRef,nReps,len,fs,level,fmin,fmax,micCorrection,cavityDiameter,...
                cavityTemperature,z0,fo_spl,phi,freq,nf,V);
        elseif kk == 2
            [fo_spl,phi,freq,Amax,pRef,nf] = XferFunction(Data,stimulusFiltered,fs,fmin,fmax);
            LTC2 = buildStruct(stim,time,stimulusFiltered,LTRecording,Data,probeInput,probeOutput,...
                Amax,pRef,nReps,len,fs,level,fmin,fmax,micCorrection,cavityDiameter,...
                cavityTemperature,z0,fo_spl,phi,freq,nf,V);
        end
    end
    
% now do calibration for fractional click bands ---------------------------
    plotBands = 0; % but don't plot these fractional bands
    freqList = [1	16
                1	8
                1	4
                1	2
                2	16
                2	8
                2	4
                4	16
                4	8
                8	16] * 1000;
    nSubBands = size(freqList,1);

    for kk=1:nChannelsOut % loop over output channels; test output channels one at a time
        alertTxt = {['     Calibrating loudspeaker ',num2str(kk)]};
        nn = size(alertTxt,1);
        for ii=1:nn
            cprintf([0,0,.4],[alertTxt{ii},'\n']);
        end
    
        for jj=1:nSubBands %-----------------------------------------------
            fmin = freqList(jj,1);
            fmax = freqList(jj,2);

            alertTxt = {['         calibrating sub-band ',num2str(fmin),' to ',num2str(fmax),' Hz']};
            nn = size(alertTxt,1);
            for ii=1:nn
                cprintf([0,0,.4],[alertTxt{ii},'\n']);
            end

            % COLLECT THE DATA ---------------------------
            obj.clearRecList % clear out whatever was used previously 
            obj.setRecList(probeInput.ch,probeInput.label,probeInput.micSens,probeInput.gain);
            obj.setNReps(nReps); % number of times to play stimulus
            obj.setFilter(0); % note: this is a highpass filter with a 75 Hz cutoff frequency.

            obj.objPlayrec.userInfo = []; % clear out previous, non-structure array info
            obj.objPlayrec.userInfo.Lvl = level;
            obj.objPlayrec.userInfo.fs = fs;
            obj.objPlayrec.userInfo.stim = stim;
            obj.objPlayrec.userInfo.caps_version = V;

            targetLvl = 10; % this is a junk value
            stimFreq = [fmin*.75,fmax*1.1];
            if kk==1
                [stimScaled,impulse,errors] = ARLas_applyLongTubeCal(LTC1,targetLvl,stimFreq,stim);
            elseif kk==2
                [stimScaled,impulse,errors] = ARLas_applyLongTubeCal(LTC2,targetLvl,stimFreq,stim);
            end
            obj.clearPlayList % clear out whatever was used previously
            obj.setPlayList(stimScaled,probeOutput.ch(kk)); % load the currently tested output channel   
            obj.objPlayrec.run % playback and record
                if obj.killRun
                   return
                end
            [header,Data] = obj.retrieveData(probeInput.label); % get raw data
            if applyMicCorrection == 1
                Data = fastFilter(micCorrection,Data);
            end
            [flatRecording,stimulusFiltered,DataFlt,time] = LTCleanData(Data,obj.fs,stimScaled,header);

            pSPL = 20*log10(max(flatRecording)/.00002); % peak SPL
            ppSPL = 20*log10((max(flatRecording)-min(flatRecording))/.00002); % peak to peak SPL
            
            origN = 400;
            nfft = 1024;
            fr = flatRecording(1:origN);
            mag = abs(fft(fr,nfft));
            mag = mag / (origN/2);
            mag = 20*log10(mag/.00002);
            mag = mag + 20*log10(sqrt(.5)); % since we want rms instead of peak values (given this is broadband)
            frq = (0:1:nfft-1)'*(fs/nfft);
            
            [~,fminIndx] = min(abs(frq-fmin));
            [~,fmaxIndx] = min(abs(frq-fmax));
            averageOut = mean(mag(fminIndx:fmaxIndx));
            stdOut = std(mag(fminIndx:fmaxIndx));
            maxDev = max(abs(mag(fminIndx:fmaxIndx)-averageOut));
            
            if plotBands == 1
                h = figure;
                subplot(2,1,1)
                plot(fr)
                subplot(2,1,2)
                plt(frq,mag)
                hold on
                line([fmin fmax],[averageOut,averageOut],'Color',[1 0 0])
                title(['maxOut = ',num2str(averageOut),' dB SPL    std = ',num2str(stdOut),' dB    ','maxAbsDev = ',num2str(maxDev),' dB'])
                xlim([fmin,fmax])
                ymin = min(mag(fminIndx:fmaxIndx));
                ymax = max(mag(fminIndx:fmaxIndx));
                ylim([ymin,ymax])
                xlabel('Frequency (Hz)')
                ylabel('Magnitude (dB SPL)')
              %keyboard
                pause(0.5)
                close(h)
            end
            
            if kk==1
                LPC1(jj,1) = averageOut;
                PP1(jj,1) = ppSPL;
                
                LTC1.(['BW_',num2str(fmin),'_',num2str(fmax)]).('foStimulus') = impulse;
                LTC1.(['BW_',num2str(fmin),'_',num2str(fmax)]).('maxOutLPC') = averageOut;
                LTC1.(['BW_',num2str(fmin),'_',num2str(fmax)]).('stdOutLPC') = stdOut;
                LTC1.(['BW_',num2str(fmin),'_',num2str(fmax)]).('maxDevLPC') = maxDev;
                LTC1.(['BW_',num2str(fmin),'_',num2str(fmax)]).('pSPL') = pSPL;
                LTC1.(['BW_',num2str(fmin),'_',num2str(fmax)]).('ppSPL') = ppSPL;
            elseif kk==2
                try
                LPC2(jj,1) = averageOut;
                PP2(jj,1) = ppSPL;
                catch ME
                    keyboard
                end
                
                LTC2.(['BW_',num2str(fmin),'_',num2str(fmax)]).('foStimulus') = impulse;
                LTC2.(['BW_',num2str(fmin),'_',num2str(fmax)]).('maxOutLPC') = averageOut;
                LTC2.(['BW_',num2str(fmin),'_',num2str(fmax)]).('stdOutLPC') = stdOut;
                LTC2.(['BW_',num2str(fmin),'_',num2str(fmax)]).('maxDevLPC') = maxDev;
                LTC2.(['BW_',num2str(fmin),'_',num2str(fmax)]).('pSPL') = pSPL;
                LTC2.(['BW_',num2str(fmin),'_',num2str(fmax)]).('ppSPL') = ppSPL;
            end
        end
    end    
    LTC1.LPC_summary = LPC1;
    LTC1.PP_summary = PP1;
    LTC2.LPC_summary = LPC2;
    LTC2.PP_summary = PP2;
    
    % Save LTC here -----------------------------
    alertTxt = {['  Saving calibrations']};
    nn = size(alertTxt,1);
    for ii=1:nn
        cprintf([0,0,.4],[alertTxt{ii},'\n']);
    end

    pathName_LTC = 'C:\myWork\ARLas\Peripheral\calibrations\LTCals\';
    if exist(pathName_LTC,'dir') ~= 7 % check to make sure that the save path exists
        success = mkdir(pathName_LTC); % if not, try to create it
        if success ~= 1
            warning('Specified folder for saving data does not exist. Entering debug mode.')
            keyboard
        else
            addpath(genpath(pathName_LTC)) % add all directories and subdirectories
        end
    end
    
    t = datetime('now'); % create the folder name for saving calibrations
    folderName = datestr(t,'mm_dd_yyyy');
    folderName = [folderName,'\']; 
    if exist([pathName_LTC,folderName],'dir') ~= 7 % check to make sure that the save path exists
        success = mkdir([pathName_LTC,folderName]); % if not, try to create it
        if success ~= 1
            warning('Specified folder for saving data does not exist. Entering debug mode.')
            keyboard
        else
            addpath(genpath([pathName_LTC,folderName])) % add all directories and subdirectories
        end
    end
    pathName = [pathName_LTC,folderName];
    
    channel = 1;
    fileName = [LTC1.probeOutput.label,'_Ch',num2str(channel),'.mat'];
    fileName = ARLas_saveName(pathName,fileName);
    LTC = LTC1;
    save([pathName,fileName],'LTC');
    
    if nChannelsOut == 2
        channel = 2;
        fileName = [LTC2.probeOutput.label,'_Ch',num2str(channel),'.mat'];
        fileName = ARLas_saveName(pathName,fileName);
        LTC = LTC2;
        save([pathName,fileName],'LTC');
    end
    
    % do plotting here ------------------------------------------
    if doPlot == 1
        ch = 1;
        h1 = makePlots(LTC1,ch);
        if nChannelsOut == 2
            ch = 2;
            h2 = makePlots(LTC2,ch);
        end
    end
    
    alertTxt = {'Long Tube calibration finished.'};
    nn = size(alertTxt,1);
    for ii=1:nn
        cprintf([0,0,.4],[alertTxt{ii},'\n']);
    end
    disp(' ')    
    
end

% INTERNAL FUNCTIONS ------------------------------------------------------
function [LTC] = buildStruct(stim,time,stimulusFiltered,LTRecording,Data,probeInput,probeOutput,...
    Amax,pRef,nReps,len,fs,level,fmin,fmax,micCorrection,cavityDiameter,cavityTemperature,z0,fo_spl,phi,freq,nf,V)
    LTC.stimOrig = stim;
    LTC.time = time;
    LTC.stimFiltered = stimulusFiltered;
    LTC.recording = LTRecording;
    LTC.Recording = Data;
    LTC.probeInput = probeInput;
    LTC.probeOutput = probeOutput;
    LTC.Amax = Amax;
    LTC.pRef = pRef;
    LTC.nReps = nReps;
    LTC.origStimLen = len;
    LTC.fs = fs;
    LTC.presentLvlReFullOut = level;
    LTC.fmin = fmin;
    LTC.fmax = fmax;
    LTC.micCorrection = micCorrection;
    LTC.cavityDiameter = cavityDiameter;
    LTC.cavityTemperature = cavityTemperature;
    LTC.z0 = z0;
    LTC.fo_spl = fo_spl;
    LTC.phi = phi;
    LTC.freq = freq;
    LTC.nf = nf;
    LTC.V = V;
end
            
function [fo_spl,phi,freq,Amax,pRef,nf] = XferFunction(Recording,stimulus,fs,fmin,fmax)
    recording = mean(Recording,2);
    stimulus = stimulus(:);
    N = length(stimulus);
    nfft = 1024;
    if N > nfft
        error('Calibration stimulus is longer than nfft. Cannot compute values.')
    end    
    S = fft(stimulus,nfft); % stimulus vector
    R = fft(recording,nfft); % recordings vector
    PL = R ./ S; % load pressure
    
    freq = (0:1:nfft-1)'*(fs/nfft);
    [~,fminIndx] = min(abs(fmin - freq));
    [~,fmaxIndx] = min(abs(fmax - freq));
    PL = PL(fminIndx:fmaxIndx); % cut to highest calibrated frequency
    freq = freq(fminIndx:fmaxIndx);
    
    Amax = 1; % system maximum output (usually 1)
    pRef = 0.00002; % pressure reference (20 micropascals)
    fo_spl = 20*log10(Amax*abs(PL)/pRef);  % load sound pressure magnitude (dB spl)
    phi = angle(PL); % phase of the response
    
    originalN = N;
    ref = 0.00002;
    [frequency,signal,noiseFloor] = ARLas_fda(Recording,fs,ref,nfft,originalN);
    snr = signal(fminIndx:fmaxIndx) - noiseFloor(fminIndx:fmaxIndx);
    nf = fo_spl - snr;

end

function [data,stimulus,Data,time] = LTCleanData(X,fs,stim,header)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [data,stimulus,Data,time] = LTcleanData(X,fs,stim,header);
%
% For use with ARLas (Auditory Research Laboratory auditory software)
% Used in applying long tube (LT) calibration to cavity and/or ear canal recordings.
%
% X (input argument) = matrix of input waveforms
% fs = sampling rate (Hz)
% stim (input argument) = waveform vector; the original electrical stimulus that was used to obtain X
% cut1 = cutoff value to remove excess zero padding at the start of the signal
% cut2 = cutoff value to remove excess zero padding at the end of the signal
%
% X (output argument) = the lowpass Butterworth filtered (22 kHz cutoff), 
%       truncated waveform w/ artifact rejection and 3 ms onset/offset windowing
% stim (output argument) the stimulus with the same filtering and windowing applied.
%
% Auditory Research Lab, The University of Iowa
% Deptartment of Communication Sciences & Disorders
% The University of Iowa
% Author: Shawn Goodman
% Original Date: October 16, 2018; cleanData.m
% Adapted June 7, 2021 -- ssg
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    [rows,cols] = size(X); % orignal size of the recording
    % it is important to do the same processing to the stimulus as to the
    % recorded waveforms.
    stim = repmat(stim,1,cols+2); % replicate stimulus to matrix the same size as recording
    stim = stim(:); % reshape into a long vector
    X = X(:);
    stim = filtfilt(header.SOS,header.G,stim); % apply the playrec highpass filter to stim (same as was for recording)
    [SOS,G] = getLPiir(fs); % get coefficients for a lowpass filter
    stim = filtfilt(SOS,G,stim); % apply the lowpass filter
    X = filtfilt(SOS,G,X);
    X = reshape(X,rows,cols); % reshape into matrices
    stim = reshape(stim,rows,cols+2);
    stimulus = stim(:,2:end-1);
    stimulus = mean(stimulus,2);

    Data = ARLas_artifactReject(X); % perform basic artifact reject
    data = mean(Data,2); % averaged data
    [~,indx] = max(abs(data)); % location of time zero
    preLen = 0.001; % amount of time to include before time zero
    postLen = 0.004; % amount of time to include after time zero
    preSamples = round(preLen * fs); % number of samples before indx zero
    postSamples = round(postLen * fs); % number of samples after indx zero
    
    data = data(indx-preSamples:indx+postSamples); % truncate data; averaged vector
    Data = Data(indx-preSamples:indx+postSamples,:); % truncate data; full data set
    stimulus = stimulus(indx-preSamples:indx+postSamples);
    N = length(data); % number of samples in truncated vector
    if mod(N,2) ~= 0 % force length to be even
        data = data(1:end-1);
        Data = Data(1:end-1,:);
        stimulus = stimulus(1:end-1);
        N = length(data);
    end
    
    data = ARLas_ramp(data,fs,0.001); % ramp onset and offset
    Data = ARLas_ramp(Data,fs,0.001);
    n = (0:1:N-1)';
    time = n/fs*1000;
    time = time - (preLen*1000);
end

function [h] = makePlots(LTC,ch)
    time = LTC.time;
    Recording = LTC.Recording;
    recording = LTC.recording;
    nf = LTC.nf;
    fo_spl = LTC.fo_spl;
    freq = LTC.freq;

    [~,stimPeakIndx] = max(abs(LTC.stimFiltered));

    h = figure;
    subplot(2,1,1)
    plt(time,Recording,'Color',[.7 .7 .7],'LineWidth',0.5)
    hold on
    plt(time,recording,'b','LineWidth',0.5)
    xlim([time(1),time(end)])
    xlabel('Time (ms)')
    ylabel('Amplitude (Pa)')
    title(['LongTube Click   ',LTC.probeOutput.label,' Ch ',num2str(ch)])
    ymax = max(abs(recording))*1.05;
    ylim([-ymax ymax])
    line([time(stimPeakIndx),time(stimPeakIndx)],[-ymax,ymax],'LineStyle','-','Color',[1 0 0],'LineWidth',0.5)

    subplot(2,1,2)
    plt(freq/1000,fo_spl,'r')
    hold on
    plt(freq/1000,nf,'k--')
    xlim([freq(1)/1000,freq(end)/1000])
    xlabel('Frequency (kHz)')
    ylabel('Full Output Level (dB SPL)')
    legend('signal','noise floor')

end

% OLD CODE ----------------------------------------------------------------
