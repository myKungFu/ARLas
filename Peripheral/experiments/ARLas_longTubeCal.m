function [] = ARLas_longTubeCal(varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ARLas_longTubeCal(varargin)
%
% Perform a "long tube" calibration for stimulus levels. 
% Time windowing used to avoid any standing waves due to reflections from
% the end of the tube.
%
% Insrt the probe tip into a long, lossy tube. 
% The incident pressure wave is measured in response to two stimuli:
% a click and a tone pip. The click is used to get the broadband transfer
% function. The pip (at 1000 Hz) is used to scale the overall output
% correctly. 
%
% Auditory Research Lab, The University of Iowa
% Deptartment of Communication Sciences & Disorders
% The University of Iowa
% Author: Shawn S. Goodman, PhD
% Date: April 7, 2017
% Revised: May 4, 2017
% Revised June 6, 2017
% Updated: August 15, 2017  Updated a few lines to make compatible with
%                           current ARLas version
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

obj = varargin{1}; % get the arlas object

% ------------------- USER ADJUSTABLE PARAMETERS --------------------------
tubeLength = 91.44; % long bass tube length (cm)
tubeDiameter =  0.322; % long tube inner diameter (cm)
fmin = 100; % minimum frequency over which to calibrate (Hz)
fmax = 20000; % maximum frequency over which to calibrate (Hz)
nReps = 24; % number of times to play the each stimulus
nTrips = 10; % number of round trip delays to wait before repeat presentation (4-5 should be sufficient)
testing = 'A'; % which probe we are testing

% INPUTS -----
    % ER10x probe A
    inputA.label = '0141'; % change this to be whatever label you want
    inputA.micSens = 0.05; % sensitivity in V/Pa
    inputA.gain = 20; % amplifier gain in dB
    inputA.ch = 1; % input channel on the sound card

    % ER10x probe B
    inputB.label = '0142'; % change this to be whatever label you want
    inputB.micSens = 0.05; % sensitivity in V/Pa
    inputB.gain = 20; % amplifier gain in dB
    inputB.ch =  2; % input channel on the sound card
    
    % G.R.A.S 1/8" Pressure Mic, Type 40DP (SHAWN'S MIC)
    inputC.label = '109338'; % change this to be whatever label you want
    inputC.micSens = 0.0012; % sensitivity in V/Pa
    inputC.gain = 0; % amplifier gain in dB
    inputC.ch = 5; % input channel on the sound card

% OUTPUTS -----    
    % ER10x probe A
    outputA.label = '0141';
    outputA.ch = [1,2]; % output channels on the sound cardthe ER10C is being used

    % ER10x probe B
    outputB.label = '0142';
    outputB.ch = [3,4]; % output channels on the sound cardthe ER10C is being used
    
if strcmp(testing,'A')
    driverProbe = outputA;      % probe being used for playing OUTPUT
    measurementProbe = inputA;  % probe whose mic is being calibrated INPUT
    fileName_inputCal = []; % microphone correction file; if empty set, will skip correction
    clickAmp = 0.75; 
    pipAmp = 0.5; 
elseif strcmp(testing,'B')
    driverProbe = outputB;      % probe being used for playing OUTPUT
    measurementProbe = inputB;  % probe whose mic is being calibrated INPUT
    fileName_inputCal = [];
    clickAmp = 0.75; 
    pipAmp = 0.5; 
else
    error('Unrecognized probe!')
end
%--------------------------------------------------------------------------

% Calibration files are typically locate here and shouldn't change
pathName_inputCal = 'C:\myWork\ARLas\Peripheral\calibrations\micCals\'; % location of mic calibration files
pathName_outputCal = 'C:\myWork\ARLas\Peripheral\calibrations\speakerCals\'; % where to save the long-tube cal files

fs = obj.fs; % get the system sampling rate
roundTripTime = 2*(tubeLength/3.4723e4); % approx round trip delay for long tube (seconds)
roundTripN = round(fs * roundTripTime); % samples

nSamples = roundTripN;
if mod(nSamples,2)~= 0
    nSamples = nSamples +1;
end
click = zeros(nSamples,1); % make the click
click(100,1) = clickAmp; % set the click amplitude to a value >0 and <1
padN = roundTripN * nTrips; % number of samples of zero padding
pad = zeros(padN,1); % zero padding
clickStim = [click;pad]; % zero pad the end of the stimulus 

pipFreq = 1000; % pip frequency in Hz
time = (0:1:nSamples-1)'/fs; % time vector
a = pipAmp; % pip amplitude (re: full out = 1). make this < 1
pip = a * sin(2*pi*pipFreq*time); % create the sinusoid
w = hann(length(time)); % make a hann window
correction = 1+(1-mean(w.^2)); % hann window reduces overall amplitude; correct by this much
w = w * correction; % boost the level of the window
pip = pip .* w; % apply the hann window to the sinusoid
pipStim = [pip;pad]; % zero pad the end of the pip stimulus 

disp(' ')
disp(' ')
disp('----- Starting Long Tube Calibration -----')

obj.clearRecList % clear out whatever was used previously 
obj.setRecList(measurementProbe.ch,measurementProbe.label,measurementProbe.micSens,measurementProbe.gain);
obj.clearPlayList % clear out whatever was used previously
obj.objPlayrec.nReps = nReps; % number of times to play stimulus

nOuts = length(driverProbe.ch); % number of outputs to calibrate
for ii=1:nOuts % loop across output channels
    
    disp(['Calibrating channel ',num2str(ii),' of ',num2str(nOuts),': Ch=',num2str(driverProbe.ch(ii))])
    if ii>1
        obj.clearPlayList % clear out whatever was used previously
    end
    obj.setPlayList(clickStim,driverProbe.ch(ii)); % load the click stimulus
    obj.objPlayrec.userInfo = []; % clear out previous, non-structure array info
    obj.objPlayrec.userInfo.channel = driverProbe.ch(ii);
    obj.objPlayrec.userInfo.fs = fs;
    
    obj.setNReps(nReps); % number of times to play stimulus
    obj.objPlayrec.run % playback and record
    if obj.killRun
       return  
    end
    [header,DataClick] = obj.retrieveData(['Ch',num2str(measurementProbe.ch)]); % retrive recorded data
    
    obj.clearPlayList 
    obj.setPlayList(pipStim,driverProbe.ch(ii));
    obj.objPlayrec.userInfo = []; % clear out previous, non-structure array info
    obj.objPlayrec.userInfo.channel = driverProbe.ch(ii);
    obj.objPlayrec.userInfo.fs = fs;
    obj.objPlayrec.run % playback and record
    if obj.killRun
       return
    end
    [~,DataPip] = obj.retrieveData(['Ch',num2str(measurementProbe.ch)]); % retrieve recorded data
            
    cutoff = fmin-50; % highpass filter cutoff
    if cutoff < 50
        cutoff = 50;
    end
    DataClick = ARLas_hpFilter(DataClick,fs,cutoff); % highpass filter recording
    click = ARLas_hpFilter(clickStim,fs,cutoff); % highpass filter stimulus
    DataPip = ARLas_hpFilter(DataPip,fs,cutoff);
    pip = ARLas_hpFilter(pipStim,fs,cutoff);

    % Apply microphone correction:
    if ~isempty(fileName_inputCal)
        try
            dummy = load([pathName_inputCal,fileName_inputCal]);
            DataClick = fastFilter(dummy.d.h,DataClick);
            DataPip = fastFilter(dummy.d.h,DataPip);
        catch
            disp('ERROR: Problem applying microphone correction.')
        end
    end
    
    % cut down to the original length (get rid of zero padding)
    start = 1; % start sample
    finish = nSamples; % finish sample
    DataClick = DataClick(start:finish,:); 
    click = click(start:finish);
    DataPip = DataPip(start:finish,:); 
    pip = pip(start:finish);
    
    DataClick = ARLas_artifactReject(DataClick); % perform post-hoc artifact reject
    DataPip = ARLas_artifactReject(DataPip);
    if mod(length(clickStim),2)~=0 % force to be an even number of samples
        DataClick = DataClick(2:end,:);
        click = click(2:end);
        DataPip = DataPip(2:end,:);
        pip = pip(2:end);
    end
    
    figure(532) % plot raw click and pip
    subplot(2,1,1)
    if ii==1
        plot(DataClick,'b')
    elseif ii == 2
        plot(DataClick,'r')
    end
    hold on
    xlabel('Time (samples)','FontSize',12)
    ylabel('Amplitude (Pa)','FontSize',12)
    xlim([1,size(DataClick,1)])
    title(['Probe: ',measurementProbe.label])
    subplot(2,1,2)
    if ii==1
        plot(DataPip,'b')
    elseif ii==2
        plot(DataPip,'r')
    end
    hold on
    xlabel('Time (samples)','FontSize',12)
    ylabel('Amplitude (Pa)','FontSize',12)
    xlim([1,size(DataPip,1)])
    
    nfft = size(click,1); % fft size
    dataClick = mean(DataClick,2); % take mean
    dataPip = mean(DataPip,2);
    d.dataClick = dataClick; % save original waveforms
    d.dataPip = dataPip;
    
    % Scale output--The size that would have been obtained if we had
    % presented the data at full output amplitude (1), assuming system
    % linearity:
    dataPip = dataPip ./ a; 
    P = abs(fft(dataPip,nfft)); % dft magnitude output of the scaled pip
    P = P / (nfft/2); % scale magnitude appropriately
    P = 20*log10(P/.00002); % put into dB SPL
    Pmax = max(P); % the maximum output possible at this frequency
    
    H = abs(fft(dataClick,nfft) ./ fft(click,nfft)); % magnitude transfer function
    freq = (0:1:nfft-1)'*(fs/nfft); % frequency (Hz)
    [~,indxMin] = min(abs(freq-fmin)); % location of frequency min and max
    [~,indxMax] = min(abs(freq-fmax));
    [~,indxNyquist] = min(abs(freq-(fs/2))); % location of Nyquist frequency
    H = H(1:indxNyquist); % cut to include only unaliased parts
    freq = freq(1:indxNyquist);
    H(1:indxMin) = 0; % zero out frequencies outside of fmin and fmax
    H(indxMax:end) = 0;
    
    % this is just a sanity check
    %Hdb = 20*log10((H+eps)/.00002);
    %[~,indx] = min(abs(freq - pipFreq)); % location of the pip frequency
    %q = interp1(freq,Hdb,pipFreq,'pchip');
    %correction = Pmax - q; % amount of difference between pip and click
    %Hdb = Hdb + correction; % apply dB correction    
    
    % smooth the calibration by making a relatively low-order FIR filter
    pRef = 0.0002; % 20 uPa is pressure reference
    n = length(H);
    X = [H;flipud(H(2:end-1))]; % add aliased portion of the DFT
    x = ifft(X); % create in impulse response
    N = length(x);
    x = [x(N/2+2:end);x(1:N/2+1)]; % make filter causal
    
    M = round(.02 * fs); % filter order is 2% of the sampling rate
    if mod(M,2)~= 0 % filter order must be an even number
        M = M + 1;
    end
    if M > length(x)
        M = round(length(x)*.75);
        if mod(M,2)~= 0 % filter order must be an even number
            M = M + 1;
        end
    end        
    
    x = x((N/2)-(M/2):(N/2)+(M/2)); % take only the desired center portion
    w = hann(M+1); % create hann window
    h = x .* w; % window the filter to smooth edges

    H2 = fft(h,fs); % put back into the frequency domain, but with 1-Hz spacing
    H2 = 20*log10(abs(H2)/pRef); % put into dB SPL
    H2 = H2(1:fs/2); % take only the un-aliased portion
    freq = (0:1:fs-1)'*(fs/fs);
    freq = freq(1:fs/2);
    [~,indx] = min(abs(freq - pipFreq)); % location of the pip frequency
    correction = Pmax - H2(indx); % amount of difference between pip and click
    H2 = H2 + correction; % apply dB correction
    FO(:,ii) = H2; % save output as FO (stands for "Full Output" matrix)
                   % FO is a spectrum showing the full output possible at
                   % each freuqency
end

% additional information to save:
d.FO = FO; % full out amplitude (dB SPL)
d.freq = freq; % frequency (Hz);
d.fmin = fmin; % minimum frequency over which to calibrate (Hz)
d.fmax = fmax; % maximum frequency over which to calibrate (Hz)
d.tubeLength = tubeLength; % long tube length (cm)
d.tubeDiameter = tubeDiameter; % long tube diameter (cm)
d.timeStamp = header.timeStamp; % time at which data were collected
d.micCorr = fileName_inputCal; % mic correction file used

% save data
disp('Saving data. Please wait...')
if exist(pathName_outputCal,'dir') ~= 7
    try 
        mkdir(pathName_outputCal)
        addpath(genpath(pathName_outputCal)) 
    catch
    end
end
fileName = ['speakerCal_', driverProbe.label,'.mat'];
fileName = ARLas_saveName(pathName_outputCal,fileName);
save([pathName_outputCal,fileName],'d') % save as a binary file    

figure % plot full output
plot(d.freq/1000,d.FO(:,1),'b')
hold on
if nOuts == 2
    plot(d.freq/1000,d.FO(:,2),'r')
end
xlim([d.fmin/1000,d.fmax/1000])
xlabel('Frequency (kHz)','FontSize',12)
ylabel('Full Output (dB SPL)','FontSize',12)
title(['Probe: ',measurementProbe.label])
if nOuts == 1
    legend('out1')
elseif nOuts == 2
    legend('out1','out2')
end

disp('----- Finished Long Tube Calibration -----')
disp(' ')
disp(' ')
