function [] = ARLas_dpoae(varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ARLas_dpoae(varargin)
%
% Measure distortion product otoacoustic emissions; standard paradigm.
%
% Auditory Research Lab, The University of Iowa
% Deptartment of Communication Sciences & Disorders
% The University of Iowa
% Author: Shawn S. Goodman, PhD
% Date: October 27, 2016
% Last Updated: November 29, 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

obj = varargin{1}; % get the arlas object

fs = obj.fs; % get the system sampling rate
len = 1; % desired stimulus length (s)
    nSamples = round(len * fs); % number of samples in stimulus
    if mod(nSamples,2) ~= 0 % force to be even number
        nSamples = nSamples + 1;
    end
    
f2 = [2000;4000;8000]; % list f2 frequencies (in Hz) to test here
nFreqs = length(f2); % number of frequencies to test
fRatio = 1.22; % f2/f1 ratio
f1 = round(f2 ./ fRatio); % calculate f1 frequencies (Hz)
fdp = 2*f1 - f2; % expected frequency of cubic distortion tone (Hz)
L1 = .025; 
L2 = .025;
time = (0:1:nSamples-1)'/fs; % time in s

DPOAE = zeros(nSamples,nFreqs); % initialize a container for the data
for ii=1:nFreqs % loop across f2 frequencies

    % 1) CREATE THE STIMULUS --------------------------------------------------
    % Each column is one output channel. The stimulus must be a matrix the same
    % size as the number of initialized output channels
    S1 = L1 * sin(2*pi*f1(ii)*time);
    S2 = L2 * sin(2*pi*f2(ii)*time);

    % 2) LOAD THE STIMULUS ----------------------------------------------------
    obj.objPlayrec.stimTrain.Ch1 = S1; % load the stimulus
    obj.objPlayrec.stimTrain.Ch2 = S2; % load the stimulus
    obj.objPlayrec.nReps = 5; % number of times to play stimulus
    
    % 3) PLAYBACK & RECORD ----------------------------------------------------
    obj.objPlayrec.run % run the stimulus
    if obj.killRun
       return
    end    
    % 4) RETRIEVE DATA ----------------------------------------------------
    channel = 1; % which channel to retrieve data from
    [header,data] = obj.retrieveData(channel);
    DPOAE(:,ii) = mean(data,2);
end

% 5) ANALYZE & PLOT DATA ----------------------------------------------------
% Common Data Dperations: (code is available in ARLas/Core/general/)
for ii=1:nFreqs
    % Frequency Domain Analysis
    ref = 0.00002; % reference is 20 uPa
    [frequency,signal,noiseFloor] = ARLas_fda(DPOAE(:,ii),fs,ref);
    
    [~,indx] = min(abs(frequency - fdp(ii)));
    Ldp(ii,1) = signal(indx);
    Ndp(ii,1) = mean([signal(indx-10),signal(indx+10)]);
    
    figure(101)
    hold off
    plot(frequency/1000,signal)
    hold on
    plot(frequency(indx)/1000,Ldp(ii,1),'*g')
    plot(frequency(indx)/1000,Ndp(ii,1),'*r')
    xlim([0.8 10])
    xlabel('Frequency (kHz)','FontSize',12)
    ylabel('Magnitude (dB SPL)','FontSize',12)
    title('DPOAE SPECTRUM','FontSize',12)    
    pause(1)
end

figure
plot(f2,Ldp,'b-*')
hold on
plot(f2,Ndp,'k:*')
xlim([f2(1)-100 f2(end)+100])
xlabel('Frequency (Hz)','FontSize',12)
ylabel('Magnitude (dB SPL)','FontSize',12)
title('DP-GRAM','FontSize',12)    

end % end of experiment file