function [] = ARLas_sfoae(varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ARLas_sfoae(varargin)
%
% Generic SFOAE program
%
% Auditory Research Lab, The University of Iowa
% Deptartment of Communication Sciences & Disorders
% The University of Iowa
% Author: Shawn S. Goodman, PhD
% Date: February 20, 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
obj = varargin{1}; % get the arlas object

% IMPORTANT! UPDATE THE FOLLOWING LINES OF CODE BEFORE RUNNING!! -----------------------------------------
fProbe = [1000,1250,1430]; % list probe frequencies (in Hz)
suppFreq = 50; % suppressor frequencies will be this many Hz higher than probe
suppLvl = 20; % suppressor level will be this many dB higher than probe

% ----- Specify Input:
input.label = 'ER10C'; % label that will be used for naming saved files 
input.micSens = 0.05; % microphone sensitivity in V/Pa
input.gain = 20; % amplifier gain in dB 
input.ch = 1; % channel on the sound card through which the recording is coming

% ----- Specify Output:
output.label = {'ER10C_Ch1','ER10C_Ch2'}; % label to identify the transducer being calibrated
output.gain = [20,20]; % gain in dB; this value not used, but just documents the physical configuration of the probe, 
output.ch = [1,2]; % channel on the sound card through which the playback is sent
output.receiver = [1,2]; % which receiver on the ER10C is being used
%--------------------------------------------------------------------------

nFreqs = length(fProbe); % number of frequencies to test
fSupp = fProbe + suppFreq; % suppressor frequencies in Hz
fs = obj.fs; % get the system sampling rate

% initialize output variables
sfoaeMag = zeros(nFreqs,1); % sfoae magnitudes (dB SPL)
sfoaeNf = zeros(nFreqs,1); % sfoae noise floors (dB SPL)
sfoaeSnr = zeros(nFreqs,1); % sfoae signal to noise ratio (dB)
probeMag = zeros(nFreqs,1); % probe level (dB SPL)
suppMag = zeros(nFreqs,1); % suppressor level (dB SPL)

for ii=1:nFreqs % loop across probe frequencies
    disp(' ')
    disp(['Testing ',num2str(fProbe(ii)),' Hz.'])
    disp(['Frequency ',num2str(ii),' of ',num2str(nFreqs)])
    disp(' ')

    % 1) CREATE THE STIMULUS --------------------------------------------------
    len = .25; % stimulus length in seconds
    N = round(len * fs); % total number of samples
    if mod(N,2)~= 0 % force to be an even numbe of samples
        N = N + 1;
    end
    t = (0:1:N-1)'/fs; % time vector
    probe = sin(2*pi*fProbe(ii)*t); % first primary tone
    supp = sin(2*pi*fSupp(ii)*t); % second primary tone
    rampLen = 0.03; % ramp length in seconds
    probe = ARLas_ramp(probe,fs,rampLen);
    supp = ARLas_ramp(supp,fs,rampLen);
    probe = probe * 10^(-suppLvl/20); % make probe smaller than suppressor
    padLen = 0.02; % ramp length in seconds
    N = round(padLen * fs); % total number of samples
    zpad = zeros(N,1);
    probe = [zpad;probe;zpad];
    supp = [zpad;supp;zpad];
    a = 10^(-40/20); % turn down stimuli 40 db from full output on the card
    supp = supp * a;
    probe = probe * a;
    zpad = zeros(size(supp));
    stimTrainCh1 = [probe;zpad;probe];
    stimTrainCh2 = [zpad;supp;supp];
    
  % LOAD THE CHANNELS ----------------------------------------------------
    obj.clearRecList % clear any previously used recording list
    obj.setRecList(input.ch,input.label,input.micSens,input.gain) % load the recording info for ARLas to use
    obj.clearPlayList % clear any previously used playback list
    obj.setPlayList(stimTrainCh1,output.ch(1)) % load the currently tested output channel
    obj.setPlayList(stimTrainCh2,output.ch(2))

  % SET USER DATA FOR HEADER FILE
    obj.objPlayrec.userInfo = []; % clear out previous, non-structure array info
    obj.objPlayrec.userInfo.fProbe = fProbe(ii);
    obj.objPlayrec.userInfo.fSupp = fSupp(ii);
    obj.objPlayrec.userInfo.lvlProbe = 20*log10(a);
    obj.objPlayrec.userInfo.lvlSuppressor = 20*log10(a) + suppLvl;
    
    % PLAYBACK & RECORD ----------------------------------------------------
    obj.objPlayrec.nReps = 24; % number of times to play stimulus
    obj.objPlayrec.run % run the stimulus
    ok = obj.checkForErrors;
    if ~ok
       return
    end   

    % RETRIEVE & ANALYZE DATA ----------------------------------------------------
    [header,data] = obj.retrieveData(input.label); % get raw data
    cutoff = 125;
    data = ARLas_hpFilter(data,fs,cutoff);
    
    %   Separate into the three parts
    nSamples = size(data,1)/3;
    start = 1;
    finish = nSamples;
    P1 = data(start:finish,:);
    Probe = P1;
    probe = mean(Probe,2);
    start = start + nSamples;
    finish = finish + nSamples;
    P2 = data(start:finish,:);
    start = start + nSamples;
    finish = finish + nSamples;
    P12 = data(start:finish,:);
    Residual = P1 + P2 - P12;

    %   Ramp Edges
    len = 0.001; % use 1 ms ramps
    Residual = ARLas_ramp(Residual,fs,len);
    %   Artifact Rejection
    Residual = ARLas_artifactReject(Residual);
    ref = 0.00002; % reference is 20 uPa
    [frequency,signal,noiseFloor] = ARLas_fda(Residual,fs,ref);
    [~,signalP1,noiseFloorP1] = ARLas_fda(P1,fs,ref);
    [~,signalP2,noiseFloorP2] = ARLas_fda(P2,fs,ref);
    
    expected = fProbe(ii); % expected frequency (Hz)
    [~,indx] = min(abs(expected - frequency));
    sfoaeMag(ii) = signal(indx);
    sfoaeNf(ii) = noiseFloor(indx);
    sfoaeSnr(ii) = sfoaeMag(ii) - sfoaeNf(ii);
    probeMag(ii) = signalP1(indx);
    expected = fSupp(ii);
    [~,indx2] = min(abs(expected - frequency));
    suppMag(ii) = signalP2(indx2);
    
    figure % analysis figure
    plot(frequency/1000,signal)
    hold on
    plot(frequency/1000,noiseFloor,'k')
    plot(frequency(indx)/1000,signal(indx),'r*')
    plot(frequency(indx)/1000,noiseFloor(indx),'k*')
    xlim([0.5,4])
    xlabel('Frequency (kHz)','FontSize',12)
    ylabel('Magnitude (dB SPL)','FontSize',12)
    title('SFOAE','FontSize',12)    
    
    disp('NOTE: The analyzed data are not currently saved. This still needs to be written...')
    
end
audiogram = figure; % sfoae-gram
plot(fProbe,sfoaeMag,'r*-')
hold on
plot(fProbe,sfoaeNf,'k:')
xlim([min(fProbe)*.95,max(fProbe)*1.05])
xlabel('Frequency (kHz)','FontSize',12)
ylabel('Magnitude (dB SPL)','FontSize',12)
title('SFOAE-gram','FontSize',12)    
    
save([obj.objPlayrec. savedFilesPath,'sfoaeData.mat'],'sfoaeMag','sfoaeNf',...
    'sfoaeSnr','probeMag','suppMag','fs')
figFileName = ['SFOAEaudiogram_',obj.subjectID,'.tif'];
saveas(audiogram,[obj.objPlayrec. savedFilesPath,figFileName])


end % end of experiment file

