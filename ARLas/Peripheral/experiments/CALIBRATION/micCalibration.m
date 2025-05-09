function micCalibration(varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% micCalibration;
%
% Microphone calibration routine for 
% Follows Siegel and Neely calibrations.
%
% All of these calibrations are saved in 'C:\myWork\ARLas\Peripheral\calibrations\micCals'
% Each calibration is further saved in a sub-directory named by date in the form 'mm_dd_yyyy'.
%
% Author: Shawn S Goodman, PhD
% Auditory Research Lab, Dept. of Comm. Sciences & Disorders, The University of Iowa
% Original Date: April 28, 2025
% Last Updated: April 28, 2025 -- ssg
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    obj = varargin{1};
    V = '28APR2025'; % version 

% USER MODIFIABLE PARAMTERS -----------------------------------------------

    doCalCheck = 1; % set = 1 to check new calibration.

    [inputs,outputs] = hardwareSetup;
    inputIHS = inputs{1}; % OAE probe IHS24 mic
    inputRef = inputs{2}; % reference 1/4 inch condensor mic
    inputPT = inputs{3}; % ER7 probe tube mic
    output = outputs{1}; % IHS24 loudspeaker
                         % NOTE: Here I am using the left channel (1) to
                         % play the chirp and the right channel shold be
                         % playing nothing.

    nReps = 128; % number of times to play the chirp

% END USER MODIFIABLE PARAMTERS -------------------------------------------
    
    fs = obj.fs; % sampling rate (Hz)
    [stimulus,F,cut1,cut2] = getStimulus(fs); % get logarithmic chirp
    fminStim = min(F); % minimum frequency in the stimulus
    fmaxStim = max(F); % maximum frequency in the stimulus
    stimulus = adjustOutput(F,stimulus); % adjust for better SNR
    stimGain_dB = 0; % reduce so you don't overdrive the system; 
    if stimGain_dB >0
        error('stimGain_dB must be <= 0')
    end
    stimulus = stimulus * 10^(stimGain_dB/20);

% -------------------------------------------------------------------------
    txt = ({'Place the sound source at one end of a coupler';'';...
        'Place the condenser mice in the other end with the mic screen very close to the probe tube';'';...
        ''});
        choice = questdlg(txt,'Mic Cal','Continue','Cancel','Continue');
    if strcmp(choice,'Continue')
        % do nothing; continue
    elseif strcmp(choice,'Cancel')
        return
    else % if user shut down box using the x
        return
    end
% -------------------------------------------------------------------------

    % LOAD SETUP FOR RECORDING 1-------------------------------------------
    obj.clearRecList % clear the previously used recording list
    %obj.setRecList(inputIHS.ch,inputIHS.label,inputIHS.micSens,inputIHS.gain); % load the recording info
    obj.setRecList(inputRef.ch,inputRef.label,inputRef.micSens,inputRef.gain); % load the recording info
    obj.setRecList(inputPT.ch,inputPT.label,inputPT.micSens,inputPT.gain); % load the recording info
 
    obj.setNReps(nReps); % number of times to play stimulus
    obj.setFilter(1); % turn on default highpass filter

    obj.clearPlayList % clear the previously used playback list
    obj.setPlayList(stimulus,output.ch(1)); % playback is through the first (left) channel of the IHS
    
    % PLAYBACK & RECORD ---------------------------------------------------
    obj.objPlayrec.run % run the stimulus
    if obj.killRun
       return
    end    

    % RETRIEVE DATA ----------------------------------------------------
    %[headerIHS,DataIHS] = obj.retrieveData(inputIHS.label); % get chirp recorded by IHS microphone
    [headerRef,DataRef] = obj.retrieveData(inputRef.label); % get chirp recorded by reference microphone
    [headerPT1,DataPT1] = obj.retrieveData(inputPT.label); % get chirp recorded by probe tube microphone

% -------------------------------------------------------------------------
    txt = ({'Be very careful not to disturbe the sound source or probe tube microphone!';'';...
        'Remove the reference condenser mic and replace it with the OAE probe mic.';'';...
        'Be sure to put the probe mic opening in the same plane as where the reference mic was.'
        ''});
        choice = questdlg(txt,'Mic Cal','Continue','Cancel','Continue');
    if strcmp(choice,'Continue')
        % do nothing; continue
    elseif strcmp(choice,'Cancel')
        return
    else % if user shut down box using the x
        return
    end
% -------------------------------------------------------------------------

    % LOAD SETUP FOR RECORDING 2-------------------------------------------
    obj.clearRecList % clear the previously used recording list
    obj.setRecList(inputIHS.ch,inputIHS.label,inputIHS.micSens,inputIHS.gain); % load the recording info
    %obj.setRecList(inputRef.ch,inputRef.label,inputRef.micSens,inputRef.gain); % load the recording info
    obj.setRecList(inputPT.ch,inputPT.label,inputPT.micSens,inputPT.gain); % load the recording info
 
    obj.setNReps(nReps); % number of times to play stimulus
    obj.setFilter(1); % turn on default highpass filter

    obj.clearPlayList % clear the previously used playback list
    obj.setPlayList(stimulus,output.ch(1)); % playback is through the first (left) channel of the IHS
    
    % PLAYBACK & RECORD ---------------------------------------------------
    obj.objPlayrec.run % run the stimulus
    if obj.killRun
       return
    end    

    % RETRIEVE DATA ----------------------------------------------------
    [headerIHS,DataIHS] = obj.retrieveData(inputIHS.label); % get chirp recorded by IHS microphone
    %[headerRef,DataRef] = obj.retrieveData(inputRef.label); % get chirp recorded by reference microphone
    [headerPT2,DataPT2] = obj.retrieveData(inputPT.label); % get chirp recorded by probe tube microphone

% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
    % CALCULATE TRANSFER FUNCTION
    pRef = mean(DataRef,2);
    pPT1 = mean(DataPT1,2);
    pPT2 = mean(DataPT2,2);
    pIHS = mean(DataIHS,2);

    nfft = length(pRef)*10;
    PRef = fft(pRef,nfft);
    PPT1 = fft(pPT1,nfft);
    PPT2 = fft(pPT2,nfft);
    PIHS = fft(pIHS,nfft);

    H = (PIHS ./ PPT2) .* (PPT1 ./ PRef); % transfer function

    % Apply smoothing splines o real and imaginary parts of H.
    % H is a complex vector in the frequency domain.
    N = length(pRef); % number of samples
    t = (0:1:N-1)'/fs * 1000; % time vector (ms)
    f = (0:1:nfft-1)'*(fs/nfft) / 1000; % frequency vector (kHz)
    smoothing = 0.5; % smoothing factor--0 is a straight line, 1 is cubic spline
    w = ones(size(f));
    mr = real(H);
    mi = imag(H);
    ppr = csaps(f,mr,smoothing,[],w); % piecewise polynomial real coefficients
    ppi = csaps(f,mi,smoothing,[],w); % piecewise polynomial imaginary coefficients
    mr_sm = ppval(ppr,f);
    mi_sm = ppval(ppi,f);
    Hsm = mr_sm + 1i*mi_sm;

    % store results in a microphone calibration structure called MCAL
    MCAL.pRef = pRef;
    MCAL.pPT1 = pPT1;
    MCAL.pPT2 = pPT2;
    MCAL.pIHS = pIHS;
    MCAL.PRef = PRef;
    MCAL.PPT1 = PPT1;
    MCAL.PPT2 = PPT2;
    MCAL.PIHS = PIHS;
    MCAL.H = H;
    MCAL.Hsm = Hsm;
    MCAL.smoothing = smoothing;
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------    

    % plot results --------------------
    figure % plot the time waveforms recorded by the microphones
    subplot(2,2,1) % pressure recorded by the reference microphone 
     plot(t,pRef,'b')
     xlabel('Time (ms)')
     ylabel('Pressure (Pa)')
     title('Reference Mic')
     grid on
    subplot(2,2,2) % pressure recorded by the probe tube microphone 
     plot(t,pPT1,'g')
     xlabel('Time (ms)')
     ylabel('Pressure (Pa)')
     title('Probe Tube Mic 1')
     grid on
    subplot(2,2,3) % pressure recorded by the IHS microphone 
     plot(t,pIHS,'r')
     xlabel('Time (ms)')
     ylabel('Pressure (Pa)')
     title('OAE Probe Mic')
     grid on
    subplot(2,2,4) % pressure recorded by the probe tube microphone  
     plot(t,pPT2,'g')
     xlabel('Time (ms)')
     ylabel('Pressure (Pa)')
     title('Probe Tube Mic 2')
     grid on

    figure % plot of the transfer function, raw and smoothed
    subplot(2,1,1) % magnitude
     plt(f,20*log10(abs(H)),'Color',[.7 .7 .7])
     hold on
     plt(f,20*log10(abs(Hsm)),'Color',[0 0 1])
     xlabel('Frequency (kHz)')
     ylabel('Magnitude (dB)')
     title('Transfer Function')
     grid on
     xlim([.1 20])
    subplot(2,1,2) % phase
     plt(f,unwrap(angle(H)),'Color',[.7 .7 .7])
     hold on
     plt(f,unwrap(angle(Hsm)),'Color',[0 0 1])
     xlabel('Frequency (kHz)')
     ylabel('Phase (rad)')
     xlim([.1 20])
     grid on

% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
    % CREATE INVERSE IMPULSE RESPONSE
    tau = 0.01; %0.0015; % filter group delay (sec)
    L = round(tau*2 * fs); % calculate filter order
    if mod(L,2)~=1 % ensure filter order is even
        L = L + 1;
    end
    M = L-1; % filter order
    M2 = M/2; % filter group delay in samples
    N = fs; % number of samples in design (using FIR from the "window" method)
    F = (0:1:N-1)'*(fs/N); % make a frequency axis corresponding to X
    f1 = 100; % minimum frequency to include (Hz)
    f2 = 22000; % maximum frequency to include (Hz)
    [~,indx1] = min(abs(F - f1)); % index of minimum frequency
    [~,indx2] = min(abs(F - f2)); % index of maximum frequency
    [~,indx4] = min(abs(F-(fs - F(indx1)))); % the locations of the alised portions of the spectrum
    [~,indx3] = min(abs(F-(fs - F(indx2))));
    Fcut = F(indx1:indx2); % cut frequency vector down to size    
    
    % normally here you use fo (full output). But in this case, you just
    % want what you measured (in dB, not dB SPL!)
    fo = interp1(f*1000,20*log10(abs(Hsm)),Fcut,'pchip');
    phi = interp1(f*1000,angle(Hsm),Fcut,'pchip');

    % calculate the maximum possible output as an impulse
    Mag = zeros(fs,1); % initialize magnitude spectrum
    Phase = zeros(fs,1); % initialize phase spectrum
    bar = 0; % level of constant output (dB SPL)
    foInv = bar - fo; % inverse full output
    foInv = foInv + bar; % put back around min value
    foInv = 10.^(foInv/20); % put into linear units
        
    % create an impulse response that will flatten the mic response
    Mag(indx1:indx2) = foInv;
    Mag(indx3:indx4) = flipud(foInv);
    Phase(indx1:indx2) = -phi;
    Phase(indx3:indx4) = flipud(phi);
    X = Mag .* cos(Phase) + 1i*Mag.*sin(Phase); % spectrum, complex rectangular form
    x = real(ifft(X)); % time domain impulse
    x = [x(end-M2:end);x(1:M2)]; % make filter causal and cut down to size
    h = hann(length(x)); % window
    b = x .* h; % apply window to the impulse
    
    tt = (0:1:length(b)-1)'/fs*1000;
    tt = tt - (length(b)/fs*1000/2);

    MCAL.b = b; % mic impulse response
    MCAL.t = tt; % time vector
    MCAL.f = f; % frequency vector
    MCAL.nfft = nfft; % fft size (samples)
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------     

    % plot results --------------------
    figure(10) % plot microphone filter
    subplot(2,1,1) % time domain impulse
     hold on
     plt(tt,b,'k')
     xlabel('Time(ms)')
     ylabel('Amplitude')
     title('Transfer Function')
     grid on
    subplot(2,1,2) % frequency domain magnitude spectrum
     hold on
     plt(f,20*log10(abs(fft(b,nfft))),'k')
     xlabel('Frequency (kHz)')
     ylabel('Magnitude (dB)')
     xlim([.1 20])
     grid on

% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
    % SAVE NEW CALIBRATION
    fsep = filesep;
    try
        [pathName_mic,~,~] = mostRecentCalibration('mic',inputIHS.label,[],obj.map);
    catch
        keyboard
    end
    
    tt = datetime('now'); % create the folder name for saving microphone calibrations
    folderName = datestr(tt,'mm_dd_yyyy');
    folderName = [folderName,fsep]; 
    if exist([pathName_mic,folderName],'dir') ~= 7 % check to make sure that the save path exists
        success = mkdir([pathName_mic,folderName]); % if not, try to create it
        if success ~= 1
            warning('Specified folder for saving mic cal does not exist. Entering debug mode.')
            keyboard
        else
            addpath(genpath(pathName_mic)) % add all directories and subdirectories
        end
    end
    fileName_mic = [inputIHS.label,'.mat'];
    saveName = ARLas_saveName([pathName_mic,folderName],fileName_mic);

    fmin = fminStim;
    fmax = fmaxStim;
    micCorrection = b;
    probeInput = inputIHS;
    probeSN = 101;
    
    MCAL.fmin = fmin;
    MCAL.fmax = fmax;
    MCAL.fs = fs;
    MCAL.probeInput = probeInput;

    save([pathName_mic,folderName,saveName],'micCorrection','fmin','fmax','fs','probeInput','folderName','probeSN','MCAL')
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------



% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
if doCalCheck == 1
    % Now run again to see if the mic correction worked:

    % LOAD SETUP FOR RECORDING 2b-------------------------------------------
    obj.clearRecList % clear the previously used recording list
    obj.setRecList(inputIHS.ch,inputIHS.label,inputIHS.micSens,inputIHS.gain); % load the recording info
    %obj.setRecList(inputRef.ch,inputRef.label,inputRef.micSens,inputRef.gain); % load the recording info
    obj.setRecList(inputPT.ch,inputPT.label,inputPT.micSens,inputPT.gain); % load the recording info
 
    obj.setNReps(nReps); % number of times to play stimulus
    obj.setFilter(1); % turn on default highpass filter

    obj.clearPlayList % clear the previously used playback list
    obj.setPlayList(stimulus,output.ch(1)); % playback is through the second (right) channel of the IHS
    
    % PLAYBACK & RECORD ---------------------------------------------------
    obj.objPlayrec.run % run the stimulus
    if obj.killRun
       return
    end    

    % RETRIEVE DATA ----------------------------------------------------
    [headerIHS,DataIHS] = obj.retrieveData(inputIHS.label); % get chirp recorded by IHS microphone
   DataIHS = fastFilter(b,DataIHS);
    %[headerRef,DataRef] = obj.retrieveData(inputRef.label); % get chirp recorded by IHS microphone
    [headerPT2,DataPT2] = obj.retrieveData(inputPT.label); % get chirp recorded by probe tube microphone

% -------------------------------------------------------------------------
    txt = ({'Place the condenser mice in the other end with the mic screen very close to the probe tube';'';...
        ''});
        choice = questdlg(txt,'Mic Cal','Continue','Cancel','Continue');
    if strcmp(choice,'Continue')
        % do nothing; continue
    elseif strcmp(choice,'Cancel')
        return
    else % if user shut down box using the x
        return
    end
% -------------------------------------------------------------------------

    % LOAD SETUP FOR RECORDING 1b-------------------------------------------
    obj.clearRecList % clear the previously used recording list
    %obj.setRecList(inputIHS.ch,inputIHS.label,inputIHS.micSens,inputIHS.gain); % load the recording info
    obj.setRecList(inputRef.ch,inputRef.label,inputRef.micSens,inputRef.gain); % load the recording info
    obj.setRecList(inputPT.ch,inputPT.label,inputPT.micSens,inputPT.gain); % load the recording info
 
    obj.setNReps(nReps); % number of times to play stimulus
    obj.setFilter(1); % turn on default highpass filter

    obj.clearPlayList % clear the previously used playback list
    obj.setPlayList(stimulus,output.ch(1)); % playback is through the second (right) channel of the IHS
    
    % PLAYBACK & RECORD ---------------------------------------------------
    obj.objPlayrec.run % run the stimulus
    if obj.killRun
       return
    end    

    % RETRIEVE DATA ----------------------------------------------------
    %[headerIHS,DataIHS] = obj.retrieveData(inputIHS.label); % get chirp recorded by IHS microphone
    [headerRef,DataRef] = obj.retrieveData(inputRef.label); % get chirp recorded by IHS microphone
    [headerPT1,DataPT1] = obj.retrieveData(inputPT.label); % get chirp recorded by IHS microphone
% -------------------------------------------------------------------------

    % CALCULATE TRANSFER FUNCTION a second time
    pRef = mean(DataRef,2);
    pPT1 = mean(DataPT1,2);
    pPT2 = mean(DataPT2,2);
    pIHS = mean(DataIHS,2);

    nfft = length(pRef)*10;
    PRef = fft(pRef,nfft);
    PPT1 = fft(pPT1,nfft);
    PPT2 = fft(pPT2,nfft);
    PIHS = fft(pIHS,nfft);

    H = (PIHS ./ PPT2) .* (PPT1 ./ PRef);

    % Create the transfer function impulse response -----------------------
    w = ones(size(f));
    mr = real(H);
    mi = imag(H);
    ppr = csaps(f,mr,smoothing,[],w); % piecewise polynomial real coefficients
    ppi = csaps(f,mi,smoothing,[],w); % piecewise polynomial imaginary coefficients
    mr_sm = ppval(ppr,f);
    mi_sm = ppval(ppi,f);
    Hsm = mr_sm + 1i*mi_sm;
    h = real(ifft(Hsm));

    NN = length(h);
    nyquist = NN/2 + 1;
    hh = [h(nyquist:end);h(1:nyquist-1)]; % make filter causal
    
    filterGroupDelay = 0.01; % desired filter group delay (s)
    m = round(filterGroupDelay * fs); % filter order
    if mod(m,2)~= 0 % force group delay to be even so that can remove it after filtering
        m = m + 1;
    end
    m2 = m/2; % half group delay
    hh = hh(nyquist-m2:nyquist+m2); % cut down filter to desired size
    h = hann(length(hh)); % window to smooth edges
    b = hh .* h; % filter coefficients in the time domain.    
    tt = (0:1:length(b)-1)'/fs*1000;
    tt = tt - (length(b)/fs*1000/2);


    figure(10)
    subplot(2,1,1)
     plt(tt,b,'r')
     xlabel('Time(ms)')
     ylabel('Amplitude')
     title('Transfer Function')
     grid on
    subplot(2,1,2)
     plt(f,20*log10(abs(fft(b,nfft))),'r')
     xlabel('Frequency (kHz)')
     ylabel('Magnitude (dB)')
     xlim([.1 20])
     grid on
end



    disp('----- Finished microphone calibration -----')
end

% INTERNAL FUNCTIONS ------------------------------------------------------
function [stim,F,cut1,cut2] = getStimulus(fs)
    fmin = 100; % starting frequency (hz)
    fmax = 24000; % ending frequency (Hz); this is 1.2 times the desired ending frequency. 
                  % fs = 96000 Hz; so this is still below nyquist (48000 Hz)
%     fmin = 100; % starting frequency (hz)
%         % find ending frequency:
%         %   2 samples = (half a cycle / f) * fs;
%         %   2 = (0.5/x)*fs
%         %   0.5*fs/2 = x
%         %   2400 = x;
%     fmax = 0.5*fs/2;
    dur = 0.125; % desired chirp duration is 125 ms
    N = round(fs * dur); % number of samples in recording
    X = (0:1:N-1)';
    X = X / max(X);
    % formula = fmin*exp(-a*X)
    % X(end) = 1, so the value of a is as follows:
    % fmax = 100*exp(-aX)
    % fmax / 100 = exp(-a)
    % -log(fmax / 100) = a
    a = -log(fmax / fmin);
    F = 100*exp(-a*X);
    lead = round(0.5/F(1)*fs); % lead samples for half a cycle
    lag = round(0.25/F(end)*fs); % lead samples for a quarter cycle
    F = [ones(lead,1)*F(1);F;ones(lag,1)*F(end)]; % instantateous frequencies
    % Rotate the phase:
    % At each time sample, time has advanced by the sampling period (1/fs).
    % At each time sample, the phase has rotated by the frequency (f) times
    % the sampling period: phase advance (in cycles) = f * (1/fs), or f/fs.
    % To make this a radian change, multiply by 2pi: phi = (2*pi*f)/fs;
    % Finally, add the newly accumulated phase to the previous phase value.
    stim = zeros(length(F),1);
    phi = 0; % starting phase is zero
    for jj = 1:length(F)
        stim(jj,1) = sin(phi);
        phi = phi + ((2*pi*F(jj))/fs);
    end

    rampN = round(lead*2); % ramp stimulus on
    h = hann(rampN*2);
    h = h(1:rampN);
    stim(1:rampN) = stim(1:rampN) .* h;
    rampN = round(lead*.5); % ramp stimulus off
    h = hann(rampN*2);
    h = h(1:rampN);
    stim = flipud(stim);
    stim(1:rampN) = stim(1:rampN) .* h;
    stim = flipud(stim);
    % zero pad stimulus with 5 ms on both sides. This is to allow for
    % filtering later. Define cut1 and cut2 so you can remove padding after
    % recording
    extra = round(fs*0.002);
    cut2 = length(stim) + extra;
    padLen = 0.050; % zero-pad the stimulus with 5 ms on both sides
    padN = round(padLen * fs);
    cut1 = padN - extra;
    cut2 = cut2 + padN;
    pad = zeros(padN,1);
    stim = [pad;stim;pad];
    F = [pad;F;pad];
    stim = stim / max(abs(stim));  % rescale so 1 is max out
    stim = stim * .99;
end
function [stim] = adjustOutput(F,stim)
   % adjust the output amplitudes to give approximately flat output across
    % the frequency range. These values based on "average" 10x driver
    % response, re: the supplied operating instructions manual from ER.
    % (plus a few emperical tweaks)
    fff =    [100,200,400,600,1000,2000,4000,6000,10000,12000,16000,20000,30000]'; % frequency (Hz)
    %output = [ 72, 55, 50, 50,  40,  40,  40,  40,   50,   60,   65,   65,   65]';
    %output = [ 50, 50, 50, 50,  50,  50,  50,  50,   50,   50,   50,   50,   50]';
    % adjusted during NU visit 3/5/20205 Shawn and Mary
    output = [ 60, 55, 50, 50,  40,  40,  40,  40,   40,   30,   60,   60,   50]'; 
    relativeOut = (output - max(output));
    a = zeros(length(stim),1);
    for jj=1:length(stim)
        a(jj,1) = interp1(fff,relativeOut,F(jj),'pchip'); % relative output amplitude (dB)
    end
    a = 10.^(a/20); % put a into linear units
    stim = a .* stim; % amplitude adjusted stimulus
    stim = stim / max(abs(stim));
    stim = stim * .99;
end

% OLD CODE ----------------------------------------------------------------
    % -----------------------------------------------
    % h = real(ifft(Hsm));
    % 
    % NN = length(h);
    % nyquist = NN/2 + 1;
    % hh = [h(nyquist:end);h(1:nyquist-1)]; % make filter causal
    % 
    % filterGroupDelay = 0.0015; % desired filter group delay (s)
    % m = round(filterGroupDelay * fs); % filter order
    % if mod(m,2)~= 0 % force group delay to be even so that can remove it after filtering
    %     m = m + 1;
    % end
    % m2 = m/2; % half group delay
    % hh = hh(nyquist-m2:nyquist+m2); % cut down filter to desired size
    % h = hann(length(hh)); % window to smooth edges
    % b = hh .* h; % filter coefficients in the time domain.    
