function [] = fPTCinput(varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Overview: 
    % Creates input data files (.txt) that can be loaded into the fPTC program 
    % of the Psychoacoustics toolbox. These input data files define fPTC 
    % parameters. fPTCinput has two main functions:
    %
    % 1) Converts a participant's behavioral tracking threshold measured at
    % any frequency/frequencies and from FPL to SPL based on
    % the participant's unique FPL/SPL transfer function obtained during
    % in situ calibration (isc) IF probeFPL = 1 (otherwise tracking thresholds are kept in 
    % % dB FPL). Then adds 10 dB (+10 SL). This determines the 
    % probe SPL value to be used for fast
    % psychophysical tuning curves (fPTC) at each test frequency. This can
    % information can be used separately if you want to define your own fPTC input
    % parameters.
    %
    % 2) Creates a .in file that defines input data for the fPTC program
    % for each possible probe frequency for a given individual participant 
    % (the possible probe frequencies are based on the absolute thresholds
    % measured using behavioral tracking). 
%% Notes/Caveats:
    % 1) SPL to FPL transfer functions are saved in the participant's main
    % directory during isc. *Currently if you run multiple iscs, the
    % transfer functions will be overwritten with the newest data. Need to
    % think about how to handle this*
    %
    % 2) Calculates FPL to SPL conversion for both channels. However, only
    % data from the first channel are saved right now (to reduce confusion
    % for the tester). It is assumed that fPTCs will be measured in channel 1. This requires 
    % loading settings into the Psychoacoustics toolbox; otherwise the toolbox will default to
    % play through channel 2. Note that behavioral tracking thresholds
    % are measured through channel 1 so it is most appropriate to use
    % channel 1 for all measures.
    %
    % 3) The psychoacoustics toolbox requires a .in file (rather than .txt)
    % when loading input data. 
%% Dependencies:
    %
    % 1) Behavioral tracking thresholds must be measured and stored in
    % folder titled 'Thresholds' within the participant's main directory
    %
    % 2) isc must have been measured using inSituCalibration_v2.m and a
    % .txt file containing the participant's SPL/FPL transfer function(s) 
    % saved to the participant's directory
    %
    % 3) Participant info file must be saved in the participant's main
    % directory 

% Author: Courtney Glavin 
% Date: May 13, 2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

inputfPTC = 1;      % set = 1 to save fPTC input files (for psychoacoustics toolbox) for each probe frequency
probeFPL = 1;       % set = 1 to save probe level in FPL; set = 0 to save probe level in SPL 

% ----------------------------------------------------------------------
% Calculate absolute thresholds in SPL -----------

% Select participant's directory 
h = waitbar(0,'Select the participant main directory');
wax = findobj(h, 'type','axes');
tax = get(wax,'title');
set(tax,'fontsize',10);
steps = 2500;
for step = 1:steps
    % computations take place here
    waitbar(step / steps)
end
close(h)

folderName = uigetdir('C:\myWork2\ARLas\Data',['Select participant folder']);
cd(folderName)

% Get list of participant's .mat files (from all folders and subfolders in directory)
files = dir('**/*.xlsx');

% % Find participant's information file(s)
% for n=1:length(files)
%     if ~isempty(strfind(files(n).name,'_participantInfo'))
%         infoFiles{n} = files(n);
%     else
%     end
% end

M = readmatrix(files(1).name);
L = M(:,2);
R = M(:,3);
frequency = M(:,1);
if ~all(isnan(L)) & ~all(isnan(R)) % if both ears have data 
    disp('Pick me!')
    keyboard
elseif all(isnan(L)) % if right ear only has data
    threshFPL = R; 
elseif all(isnan(R)) % if left ear only has data
    threshFPL = L;
else all(isnan(L)) % if neither ear has data
    error('All data are NaN.')
end

% clean: only keep rows with data
indx = find(~isnan(threshFPL));
frequency = frequency(indx);
threshFPL = threshFPL(indx);


% % Get rid of empty cells
% infoFiles_new = infoFiles(~cellfun(@isempty, infoFiles));
% i=1;

% % Find the newest participantInfo file if there are multiple
% if length(infoFiles_new) > 1;
%     for i=1:length(infoFiles_new)
%         f(i) = infoFiles_new{i}.datenum;
%     end 
% 
%     % Find index of newest file based on datenum
%     [~,indx] = max(f);
%     newestFile = infoFiles_new{indx};
% else length(infoFiles_new) == 1;
%     newestFile = infoFiles_new{i};
% end

% % Load participant's data from most recent file
% load([newestFile.folder,filesep,newestFile.name]);
% 
% clear i;
% clear n;
% clear indx;
% 
% % Select participant's tracking thresholds
% h = waitbar(0,'Select the tracking thresholds (.mat file)');
% wax = findobj(h, 'type','axes');
% tax = get(wax,'title');
% set(tax,'fontsize',10);
% steps = 2500;
% for step = 1:steps
%     % computations take place here
%     waitbar(step / steps)
% end
% close(h)
% 
% thresholdfile = uigetfile(fullfile([folderName,filesep, 'Thresholds',filesep,'*trackingThresholdsFPL*.mat']));
% trackingthresholds = load(deblank([folderName,filesep,'Thresholds',filesep,thresholdfile]));
% 
% % Load in threshold frequencies and FPL values
% threshFPL = [trackingthresholds.frequency, trackingthresholds.all(isnan(L))];

% if probeFPL == 0
    % Load participant's SPL -> FPL transfer functions for each channel
    % t1 = dir([folderName,filesep,'*t1.mat']);   % channel 1
    % t2 = dir([folderName,filesep,'*t2.mat']);   % channel 2
    % load(t1.name)
    % load(t2.name)

    % Interpolate transfer function over range of threshold test frequencies
    % Channel 1
    % spl1 = interp1(t1(:,1),t1(:,2),(threshFPL(:,1)'*1000)');    % tracking threshold frequency in kHz; change to Hz
    % % Channel 2
    % spl2 = interp1(t2(:,2),t2(:,2),(threshFPL(:,1)'*1000)');
% 

    path = pwd;
    k = strfind(path,'Data');
    path = path(k+5:end);
    k = strfind(path,'\');
    subjectID = path(1:k(1)-1);


    LCpath = ['C:\myWork2\ARLas\Data\',subjectID,'\loadCals\'];
    fPTCpath = ['C:\myWork2\ARLas\Data\',subjectID,'\fPTC\'];
    cd(LCpath)
    d = dir([LCpath,'*.mat']);
    nFiles = length(d);
    for ii=1:nFiles
        Q{ii,:} = d(ii).date; 
    end
    if isempty(Q)
        disp(['No valid load calibration found for this subject.'])
        disp(['Run the program "loadCalibration.m" first.'])
        return
    end
    [pick,iPick] = max(datetime(Q));
    recall = load([LCpath,d(iPick).name]);
    LC = recall.LC;



    spl1 = LC.iscSA1.spl;
    spl2 = LC.iscSA2.spl;
    ff = LC.iscSA1.freq;
    fpl1 = LC.iscSA1.fpl;
    fpl2 = LC.iscSA2.fpl;

    t1 = spl1 - fpl1;
    t2 = spl2 - fpl2;
    for ii=1:length(threshFPL)
        [~,Indx(ii,1)] = min(abs(frequency(ii) - ff));
    end
    probePTCSPL1 = [frequency,threshFPL+t1(Indx)];
    probePTCSPL2 = [frequency,threshFPL+t2(Indx)];


    % % Apply transfer function to thresholds (go from FPL -> SPL by adding
    % % transfer function to the threshold)
    % % Channel 1
    % probePTCSPL1 = [threshFPL(:,1)*1000,(threshFPL(:,2)+spl1)];
    % % Channel 2
    % probePTCSPL2 = [threshFPL(:,1)*1000,(threshFPL(:,2)+spl2)];

    % Save data for tester to reference
    %writematrix(probePTCSPL1,[folderName,filesep,'Thresholds',filesep,'probePTCSPL_Ch1','.txt']); 
    writematrix(probePTCSPL1,[fPTCpath,'probePTCSPL_Ch1','.txt']); 
    
    %writematrix(probePTCSPL2,[folderName,filesep,'Thresholds',filesep,'probePTCSPL_Ch2','.txt']); 
% 
% end 

% ----------------------------------------------------------------------
% Create input files for each possible PTC probe frequency -----------

if inputfPTC == 1

    % Create inputData cell array to store data 
    % note: parameters outside of for loop should not change across probe frequencies

    inputData{1,1} = 'Input data to SWPTC software';
    inputData{2,1} = 95.5;                  % Calibration data (dB SPL value for xx V RMS; value of xx defined below; currently using 95.5 based on ER10xA0130 data from Steve Viryani;
                                            % these data are stored under C:\myWork\ARLas\Peripheral\experiments\CCG\Dissertation\fPTC\711data_original_5.17.2022.xlsx. Note that these data are the 
                                            % output of our probe mic's driver measured by the probe mic and in the coupler. The sound card used to measure the data was not able to drive to 1V, 28.5 dB should
                                            % be added to all data to approximate a 1V drive. 95.5 came from the data (+28.5 dB) at 1 kHz.
    inputData{3,1} = 1 ;                    % Calibration data (V RMS) - hardcode here based on system measurements; the V RMS needed to produce the output in the line above
    inputData{4,1} = 96000;                 % Sample rate 
    channel = 1;                            % Define channel to play out; default to channel 1
    inputData{5,1} = channel;               % Channel
    
    
    % subjectID = subjectData.subjectID(1,:);      % Participant ID; based on participant info. file loaded
    %     if ~isempty(strfind(subjectID,'_')) % Get rid of '_' character in subjectID if necessary. fPTC will not accept '_' character
    %         subjectID = erase(subjectID,'_');
    %     end 
    inputData{6,1} = subjectID;         
    %testear = subjectData.ear(1,:);          % Participant test ear; based on participant info. file loaded
    testear = LC.testEar;
    inputData{7,1} = testear; 
    
    % Make new directory if needed
    if ~exist(fPTCpath, 'dir')
       mkdir(fPTCpath)
       addpath(fPTCpath)
    end

    % Create a new input data .txt file (for the Psychoacoustics software)
    % for each possible fPTC probe frequency (based on absolute threshold
    % frequencies tested) 
    for i=1:length(threshFPL(:,1)) 
        
        %fProbe = threshFPL(i,1);   % probe frequency
        fProbe = frequency(i,1);   % probe frequency
        durProbe = 0.2;                 % probe signal duration (200ms per Sek & Moore 2011)
        interval = 0.2;                 % interval between signal pulses (200ms per Sek & Moore 2011)
        
        % Probe signal level (absolute threshold in FPL or SPL + 10 dB)
        if probeFPL == 0
            if channel == 1                 % define by channel 
                levelProbe = probePTCSPL1(i)+10;
            elseif channel == 2
                levelProbe = probePTCSPL2(i)+10; 
            end 
        elseif probeFPL == 1
            if channel == 1                 % define by channel
                levelProbe = threshFPL(i)+10;
            elseif channel == 2
                levelProbe = threshFPL(i)+10; 
            end 
        end 
        
        % Set start and stop frequency of masker (based on probe frequency)
        if fProbe < 10000
            fStartMask = 1.5 * fProbe;  % Charaziak et al (2012)
            fStopMask = 0.5 * fProbe;   % Charaziak et al (2012); Sek et al (2005)
        elseif fProbe >= 10000
            fStartMask = 1.3 * fProbe;  % Buus et al (1986)
            fStopMask = 0.7 * fProbe; 
        end 
    
        % Noise duration 
        durMask = 240;                  % 240 seconds = 4 minutes; from KC's code & from Sek et al. (2011) 
    
        % Select bandwidth of masker (based on probe frequency)
        % Balance to minimize beat detection and ability to define tip frequency. Should not be consequential for normal hearing participants (CCG)
        if fProbe <= 1500            
            bw = 0.2 * fProbe;
        elseif fProbe > 1500
            bw = 320;                   % (Sek et al 2005; Baiduc et al 2012)
        end 
        
        sweepDir = 'Downward';          % Upward; Masker sweep direction; default downward (high to low) (Charaziak et al 2012; Sek & Moore 2005; Sek & Moore 2015)
        
        % Switch definition of fStart/Stop based on direction of masker
        % First number printed in text file should be lowest frequency of
        % masker; second number should be highest frequency of masker 
        if strcmp(sweepDir, 'Downward') == 1
            fStart = fStopMask;
            fStop = fStartMask; 
        end 
    
        delta = 2;                      % rate of change of masker level (in dB/s) 
        levelStartMask = 50;            % initial noise level (Sek et al 2005)
        
        % Combine into matrix
        inputData{8,1} = fProbe; 
        inputData{9,1} = durProbe;
        inputData{10,1} = interval;
        inputData{11,1} = levelProbe; 
        inputData{12,1} = fStart;
        inputData{13,1} = fStop;
        inputData{14,1} = durMask;
        inputData{15,1} = bw; 
        inputData{16,1} = sweepDir;
        inputData{17,1} = delta;
        inputData{18,1} = levelStartMask;
        % Lines 19-22 can be set to add extra noise (LPN; TEN)
        inputData{19,1} = 'Noise Section: No additional noise';
        inputData{20,1} = 'False';
        inputData{21,1} = 0; 
        inputData{22,1} = 0;
        inputData{23,1} = 'End of Input file';
    
        % Save in psychoacoustics toolbox folder
        savePsychPath1 = 'C:\Users\Public\PSYCHOACOUSTICS\SWPTC\Input_data\nOAExt\';
        writecell(inputData,[savePsychPath1,subjectID,'_',num2str(fProbe),'_','Ch1','.in'],'FileType','text');
        
        % Save in participant's folder 
        %writecell(inputData,[folderName,filesep,'fPTC',filesep,subjectID,'_',testear,'_','inputData','_',num2str(fProbe),'kHz','.txt']);

    end 
end
       
% Sample text file format (for reference only)
%Input data to SWPTC software 
%95 - calibration data (driver output for 1V RMS input)
%1 - calibration (V RMS) 
%48000 - fs
%1 - channel number
%Example - subject ID
%L - test ear
%500 - signal frequency
%0.3 - signal duration
%0.2 - interval between signal pulses
%67 - signal level
%250 - lowest frequency of masker
%750 - highest frequency of masker
%40 - noise duration  of masker
%100 - bandwidth of masker
%Upward - sweep direction
%4 - rate of change of masker level
%50 - initial noise level
%Noise Section: No additional noise
%False
%0
%0
%End of Input file