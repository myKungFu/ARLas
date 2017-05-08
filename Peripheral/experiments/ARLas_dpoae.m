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
% Date: May 4, 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ------------------- USER ADJUSTABLE PARAMETERS --------------------------
fmin = 500; % minimum f2 frequency to test (Hz)
fmax = 20000; % maximum f2 frequency to test (Hz)
stepSize = 1/3; % f2 frequency step size (octaves)
fRatio = 1.22; % f2/f1 ratio
L1 = 75; % level of f1 (dB SPL)
L2 = 65; % level of f2 (dB SPL)
nReps = 20; % number of times to play the 1-second stimulus

    % Specify Inputs and Outputs:
    %   It is assumed that a 1-channel configuration is being used--
    %   1 channel in and 2 channels out.
    input.label = 'ER10xA'; % change this to be whatever label you want
    input.micSens = 0.05; % sensitivity in V/Pa
    input.gain = 20; % amplifier gain in dB
    input.ch = 1; % input channel on the sound card

    output.label = 'ER10xA';
    output.ch = [1,2]; % output channels on the sound card
    
calFile = 'LongTubeCal_q_2.mat'; % set to a specific file, using ARLas_longTubeCal.m works well.
%--------------------------------------------------------------------------

disp(' ')
disp(' ')
disp('----- Starting DPOAE experiment -----')
obj = varargin{1}; % get the arlas object
fs = obj.fs; % get the system sampling rate
len = 1; % desired stimulus length (s). Should be 1 second
    nSamples = round(len * fs); % number of samples in stimulus
    if mod(nSamples,2) ~= 0 % force to be even number
        nSamples = nSamples + 1;
    end
time = (0:1:nSamples-1)'/fs; % time in s
f2 = 2.^((log2(fmin):stepSize:log2(fmax))'); % vector of f2 frequencie to test
f2 = round(f2); % force to be integers (so they fall precisely on an fft bin)
nFreqs = length(f2); % number of frequencies to test
f1 = f2 ./ fRatio; % calculate f1 frequencies (Hz)
f1 = round(f1); % force to be integers
fdp_cubic = 2*f1 - f2; % frequencies of cubic distortion tone (Hz)
fdp_diff = f2-f1; % frequencies of difference tones

% CREATE THE STIMULI ------------------------------------------------------
S1 = zeros(nSamples,nFreqs); % initialize stimulus matrix
S2 = zeros(nSamples,nFreqs);
for ii=1:nFreqs
    dummy = load(calFile);
    fo = dummy.d.FO; % full card output (dB SPL)
    freq = dummy.d.freq; % frequency (Hz)
    
    [~,indx] = min(abs(freq-f1(ii))); % indx of f1
    a1 = 10^((L1 - fo(indx,1)) /20); % multiplier for first primary
    [~,indx] = min(abs(freq-f2(ii))); % indx for f2
    a2 = 10^((L2 - fo(indx,2)) /20); % multiplier for second primary
    if a1 > 1
        disp('ERROR: desired output esceeds max output for f1!')
        return
    end
    if a2 > 1
        disp('ERROR: desired output esceeds max output for f2!')
        return
    end
    S1(:,ii) = a1 * sin(2*pi*f1(ii)*time);
    S2(:,ii) = a2 * sin(2*pi*f2(ii)*time);
end
DPOAE = zeros(nSamples,nFreqs); % initialize a container for the saved data

% COLLECT THE DATA --------------------------------------------------------
obj.clearRecList % clear out whatever was used previously 
obj.setRecList(input.ch,input.label,input.micSens,input.gain);
obj.clearPlayList % clear out whatever was used previously
obj.objPlayrec.nReps = nReps; % number of times to play stimulus

h1 = figure(111); 
    subplot(2,1,1) % analysis figure
    hold off
    plot(1,1,'.w')
    xlabel('Frequency (kHz)','FontSize',12)
    ylabel('Magnitude (dB SPL)','FontSize',12)
    title('DPOAE: ANALYSIS','FontSize',12)    
    subplot(2,1,2) % dp-gram figure
    hold off
    plot(1,1,'.w')
    xlabel('Frequency (Hz)','FontSize',12)
    ylabel('Magnitude (dB SPL)','FontSize',12)
    title('DP-GRAM','FontSize',12)    

d = [];
t0 = clock;
for ii=1:nFreqs % loop across f2 frequencies
    t1 = clock;
    disp(['Testing frequency ',num2str(ii),' of ',num2str(nFreqs),': f2=',num2str(f2(ii)),' Hz'])
    if ii > 1
        disp(['  ',num2str(estRemain),' min remaining (est.)'])
    end
    obj.setPlayList(S1(:,ii),output.ch(1));
    obj.setPlayList(S2(:,ii),output.ch(2));
    
    obj.objPlayrec.userInfo = []; % clear out previous, non-structure array info
    obj.objPlayrec.userInfo.f1 = f1(ii);
    obj.objPlayrec.userInfo.f2 = f2(ii);
    obj.objPlayrec.userInfo.fdp_cubic = fdp_cubic(ii);
    obj.objPlayrec.userInfo.fdp_diff = fdp_diff(ii);
    obj.objPlayrec.userInfo.L1 = L1;
    obj.objPlayrec.userInfo.L2 = L2;
    obj.objPlayrec.userInfo.fs = fs;
    
    obj.objPlayrec.run % playback and record
    if obj.killRun
       return
    end    
    [header,Data] = obj.retrieveData(['Ch',num2str(input.ch)]); % retrive recorded data
    d = analyzeDPOAE_HEARD(header,Data,d,h1);
    plotDPOAE_HEARD(d,h1,0) % update plot for each frequency tested
    t2 = clock;
    g = etime(t2,t1)/60;
    estRemain = g*(nFreqs-ii); % estimated time remaining
end
% additional information to save:
d.L1_target = L1; % the target level
d.L2_target = L2; % the target level
d.nReps = nReps; % the number of stimulus repetitions
d.fRatio = fRatio; % primary frequency ratio
d.stepSize = stepSize; % frequency step size (octaves)
d.calFile = calFile; % which calibration file was used for this data set
d.subjName = obj.subjectID; % the subject ID
d.timeStamp = header.timeStamp; % time at which data were collected

t2 = clock;
g = etime(t2,t0)/60;
disp(['Data collection finished. Test time was ',num2str(g),' minutes.'])

figure(h1)
subplot(1,1,1)
plotDPOAE_HEARD(d,h1,1) % final dp-gram

% save data
disp('Saving data. Please wait...')
pathName = obj.objPlayrec.savedFilesPath;      
fileName = ['DPOAE_',obj.subjectID,'.mat'];
fileName = ARLas_saveName(pathName,fileName);
disp(['   Saving binary data (',fileName,')'])
save([pathName,fileName],'d') % save as a binary file    
% save figure
figFileName = fileName;
figFileName = figFileName(1:end-3);
figFileName = [figFileName,'tif'];
disp(['   Saving figure (',figFileName,')'])
saveas(h1,[pathName,figFileName])
% write date to excell spreadsheet
xlFileName = fileName;
xlFileName = xlFileName(1:end-3);
xlFileName = [xlFileName,'xls'];
disp(['   Saving Excel data (',xlFileName,'('])
writeDPgram(pathName,xlFileName,d)

disp('----- Finished DPOAE experiment -----')
disp(' ')
disp(' ')
end % end of experiment file

% Internal Functions ------------------------------------------------------
function [] = writeDPgram(pathName,fileName,d)
    warning off
    disp('    .')
    sheet = 'dpgram';
    range = 'A1'; % F2 data
    xlswrite([pathName,fileName],{'F2 (Hz)'},sheet,range);
    range = 'A2';
    xlswrite([pathName,fileName],d.f2',sheet,range);
    range = 'B1';
    xlswrite([pathName,fileName],{'L2 (Hz)'},sheet,range);
    range = 'B2';
    xlswrite([pathName,fileName],d.L2',sheet,range);
    range = 'C1';
    xlswrite([pathName,fileName],{'Phi_F2 (rad)'},sheet,range);
    range = 'C2';
    xlswrite([pathName,fileName],d.P2',sheet,range);
    disp('    .')
    range = 'D1'; % F1 data
    xlswrite([pathName,fileName],{'F1 (Hz)'},sheet,range);
    range = 'D2';
    xlswrite([pathName,fileName],d.f1',sheet,range);
    range = 'E1';
    xlswrite([pathName,fileName],{'L1 (Hz)'},sheet,range);
    range = 'E2';
    xlswrite([pathName,fileName],d.L1',sheet,range);
    range = 'F1';
    xlswrite([pathName,fileName],{'Phi_F1 (rad)'},sheet,range);
    range = 'F2';
    xlswrite([pathName,fileName],d.P1',sheet,range);
    disp('    .')
    range = 'G1'; % 2f1-f2 data
    xlswrite([pathName,fileName],{'Fdp_cubic (Hz)'},sheet,range);
    range = 'G2';
    xlswrite([pathName,fileName],d.fdp_cubic',sheet,range);
    range = 'H1';
    xlswrite([pathName,fileName],{'Ldp_cubic (Hz)'},sheet,range);
    range = 'H2';
    xlswrite([pathName,fileName],d.Ldp_cubic',sheet,range);
    range = 'I1';
    xlswrite([pathName,fileName],{'Phi_dp_cubic (rad)'},sheet,range);
    range = 'I2';
    xlswrite([pathName,fileName],d.Pdp_cubic',sheet,range);    
    range = 'J1';
    xlswrite([pathName,fileName],{'Noise_dp_cubic (rad)'},sheet,range);
    range = 'J2';
    xlswrite([pathName,fileName],d.Ndp_cubic',sheet,range);    
    disp('    .')
    range = 'K1'; % f2-f1 data
    xlswrite([pathName,fileName],{'Fdp_diff (Hz)'},sheet,range);
    range = 'K2';
    xlswrite([pathName,fileName],d.fdp_diff',sheet,range);
    range = 'L1';
    xlswrite([pathName,fileName],{'Ldp_diff (Hz)'},sheet,range);
    range = 'L2';
    xlswrite([pathName,fileName],d.Ldp_diff',sheet,range);
    range = 'M1';
    xlswrite([pathName,fileName],{'Phi_dp_diff (rad)'},sheet,range);
    range = 'M2';
    xlswrite([pathName,fileName],d.Pdp_cubic',sheet,range);    
    range = 'N1';
    xlswrite([pathName,fileName],{'Noise_dp_diff (rad)'},sheet,range);
    range = 'N2';
    xlswrite([pathName,fileName],d.Ndp_diff',sheet,range);    
    disp('    .')
    sheet = 'Exp Parameters';
    range = 'A1'; 
    xlswrite([pathName,fileName],{'Subj Name'},sheet,range);
    range = 'B1';
    xlswrite([pathName,fileName],d.subjName,sheet,range);    
    range = 'A2'; 
    xlswrite([pathName,fileName],{'Time Stamp'},sheet,range);
    range = 'B2';
    xlswrite([pathName,fileName],d.timeStamp,sheet,range);
    range = 'A3'; 
    xlswrite([pathName,fileName],{'Sample Rate'},sheet,range);
    range = 'B3';
    xlswrite([pathName,fileName],d.fs,sheet,range);
    range = 'A4'; 
    xlswrite([pathName,fileName],{'Stim Reps'},sheet,range);
    range = 'B4';
    xlswrite([pathName,fileName],d.nReps,sheet,range);       
    range = 'A5'; 
    xlswrite([pathName,fileName],{'L1 Target'},sheet,range);
    range = 'B5';
    xlswrite([pathName,fileName],d.L1_target,sheet,range);
    range = 'A6';
    xlswrite([pathName,fileName],{'L2 Target'},sheet,range);
    range = 'B6';
    xlswrite([pathName,fileName],d.L2_target,sheet,range);
    range = 'A7';
    xlswrite([pathName,fileName],{'fRatio'},sheet,range);
    range = 'B7';
    xlswrite([pathName,fileName],d.fRatio,sheet,range);       
    range = 'A8';
    xlswrite([pathName,fileName],{'Calibration File'},sheet,range);
    if ~isempty(d.calFile)
         range = 'B8';
        xlswrite([pathName,fileName],d.calFile,sheet,range);
    end
    
    warning on
end

