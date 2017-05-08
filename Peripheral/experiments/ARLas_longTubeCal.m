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
% function. The pip (at 1000 Hz) is used to scale the overal output
% correctly. 
%
% Auditory Research Lab, The University of Iowa
% Deptartment of Communication Sciences & Disorders
% The University of Iowa
% Author: Shawn S. Goodman, PhD
% Date: April 7, 2017
% Revised: May 4, 2017

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

obj = varargin{1}; % get the arlas object

% ------------------- USER ADJUSTABLE PARAMETERS --------------------------
tubeLength = 475; % long tube length (cm) NOTE: Must be > 347 cm (3.47 m, 136.6 in, 11.33 ft)
tubeDiameter = 0.8; % long tube diameter (cm) NOTE: make this similar to human ear canal diameter, 0.7 or 0.8

nReps = 24; % number of times to play the each stimulus
nTrips = 4; % number of round trip delays to wait before repeat presentation (4-5 should be sufficient)


    % Specify Inputs and Outputs:
    %  NOTE: It is assumed that a 1-channel input configuration is being used--
    %   1 channel in and N channels out.
    input.label = 'ER10xA'; % change this to be whatever label you want
    input.micSens = 0.05; % sensitivity in V/Pa
    input.gain = 20; % amplifier gain in dB
    input.ch = 1; % input channel on the sound card --can only be one channel

    output.label = 'ER10xA';
    output.ch = [1,2]; % output channels on the sound card --can be one or more
%--------------------------------------------------------------------------


fmin = 100; % minimum frequency over which to calibrate (Hz)
fmax = 25000; % maximum frequency over which to calibrate (Hz)
if tubeLength < 347
    disp('ERROR: tube length must be >= 347 cm (3.47 m, 136.6 in, 11.33 ft).')
    return
end

fs = obj.fs; % get the system sampling rate
roundTripTime = 2*(tubeLength/3.4723e4); % approx round trip delay for long tube (seconds)
roundTripN = round(fs * roundTripTime); % samples

if roundTripTime > 0.03 % if round trip time exceeds 30 ms
    nSamples = round(0.03 * fs); % set stim length to 30 ms
else
    nSamples = roundTripN; % othrwise, make as long as possible;
end
if mod(nSamples,2)~= 0
    nSamples = nSamples +1;
end
click = zeros(nSamples,1); % make the click
click(200,1) = .5; % set the click amplitude to a value >0 and <1
padN = roundTripN * nTrips; % number of samples of zero padding
pad = zeros(padN,1); % zero padding
clickStim = [click;pad]; % zero pad the end of the stimulus 

pipFreq = 1000; % pip frequency in Hz
time = (0:1:nSamples-1)'/fs; % time vector
a = .5; % pip amplitude (re: full out = 1). make this < 1
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
obj.setRecList(input.ch,input.label,input.micSens,input.gain);
obj.clearPlayList % clear out whatever was used previously
obj.objPlayrec.nReps = nReps; % number of times to play stimulus

nOuts = length(output.ch); % number of outputs to calibrate
for ii=1:nOuts % loop across output channels
    
    disp(['Calibrating channel ',num2str(ii),' of ',num2str(nOuts),': Ch=',num2str(output.ch(ii))])
    obj.clearPlayList % clear out whatever was used previously
    obj.setPlayList(clickStim,output.ch(ii)); % load the click stimulus
    obj.objPlayrec.userInfo = []; % clear out previous, non-structure array info
    obj.objPlayrec.userInfo.channel = output.ch(ii);
    obj.objPlayrec.userInfo.fs = fs;
    obj.objPlayrec.run % playback and record
    if obj.killRun
       return  
    end
    [header,DataClick] = obj.retrieveData(['Ch',num2str(input.ch)]); % retrive recorded data
    
    obj.setPlayList(pipStim,output.ch(ii));
    obj.objPlayrec.userInfo = []; % clear out previous, non-structure array info
    obj.objPlayrec.userInfo.channel = output.ch(ii);
    obj.objPlayrec.userInfo.fs = fs;
    obj.objPlayrec.run % playback and record
    if obj.killRun
       return
    end
    [header2,DataPip] = obj.retrieveData(['Ch',num2str(input.ch)]); % retrive recorded data
    
    
    cutoff = fmin-50; % highpass filter cutoff
    if cutoff < 50
        cutoff = 50;
    end
    DataClick = ARLas_hpFilter(DataClick,fs,cutoff); % highpass filter recording
    click = ARLas_hpFilter(clickStim,fs,cutoff); % highpass filter stimulus
    DataPip = ARLas_hpFilter(DataPip,fs,cutoff);
    pip = ARLas_hpFilter(pipStim,fs,cutoff);
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
    
    nfft = size(click,1); % fft size
    dataClick = mean(DataClick,2); % take mean to increase SNR
    dataPip = mean(DataPip,2);
    
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
    % smooth the calibration by making a relatively low-order FIR filter
    M = round(.02 * fs); % filter order is 2% of the sampling rate
    if mod(M,2)~= 0 % filter order must be an even number
        M = M + 1;
    end
    pRef = 0.0002; % 20 uPa is pressure reference
    n = length(H);
    X = [H;flipud(H(2:end-1))]; % add aliased portion of the DFT
    x = ifft(X); % create in impulse response
    N = length(x);
    x = [x(N/2+2:end);x(1:N/2+1)]; % make filter causal
    x = x((N/2)-(M/2):(N/2)+(M/2)); % take only the desired center portion
    w = hann(M+1); % create hann window
    h = x .* w; % window the filter to smooth edges
    H2 = fft(h,N); % put back into the frequency domain
    H2 = 20*log10(abs(H2)/pRef); % put into dB SPL
    H2 = H2(1:n); % take only the un-aliased portion
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

% save data
disp('Saving data. Please wait...')
pathName = obj.objPlayrec.savedFilesPath;      
fileName = ['LongTubeCal_',obj.subjectID,'.mat'];
fileName = ARLas_saveName(pathName,fileName);
disp(['   Saving binary data (',fileName,')'])
save([pathName,fileName],'d') % save as a binary file    

disp('----- Finished Long Tube Calibration -----')
disp(' ')
disp(' ')



