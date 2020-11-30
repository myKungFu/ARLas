function [] = ARLas_micCorrection_DW10x(varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ARLas_micCorrection_DW10x(varargin);
%
% This code is used to make microphone corrections.
% This version (DW10x) is for the ER10X in Jeff Lichtenhan's lab. 
%
% Calibration tube is 36 inches in length (36 * 2.54 = 91.44 cm)
% Travel time in the tube is 91.44 / 34400 = 0.0027 s = 2.7 ms (one way)
% Round-trip travel time = 2.7*2 = 5.4 ms.
% Assume 10 round trip travel times are needed to allow the stimulus to
% completely decay. The stimulus duration should be 5.4*10 = 54 ms.
%
%
% OPTIONAL INPUT ARGUMENT:
% inputRef = struture specifying reference input
%
% OUTPUT ARGUMENT:
% isc = in-situ calibration structure.
%
% Authors: Shawn S Goodman, PhD & Jeffery T Lichtenhan PhD
% Auditory Research Lab, Dept. of Comm. Sciences & Disorders, The University of Iowa
% Auditory Physiology Lab, Dept. of Otolaryngology, Washington University School of Medicine, St. Louis
% Date: October 31, 2019
% Updated: October 31, 2019 -- ssg
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    obj = varargin{1};
    params = varargin{2};

    % unpack the paramter structure
    %applyMicCorrection = params.applyMicCorrection;
    fmin = params.fmin;
    fmax = params.fmax;
    testProbe = params.testProbe;
    %cavityTemperature = params.cavityTemperature;
    fs = obj.fs; % get the system sampling rate


% ---  IMPORTANT! UPDATE THE FOLLOWING LINES OF CODE BEFORE RUNNING!!   ---
    nRepsClick = 1024; % number of times to play the stimulus
% -------------------------------------------------------------------------

    [inputs,outputs] = hardwareSetup; %hardwareSetup; % read in the saved hardware setup
    inputA = inputs{1};        % ER10xA microphone
    outputA = outputs{1};      % ER10xA loudspeakers
    inputB = inputs{2};        % ER10xB microphone
    outputB = outputs{2};      % ER10xB loudspeakers
    inputRef = inputs{4};         % GRAS 1/8" reference microphone
   
    if strcmp(testProbe,'A') % make mic correction for probe A (use probe B to present)
        input = inputA;
        output = outputB;
        [pathName_mic,folderName_mic,fileName_mic] = mostRecentCalibration('mic',inputA.label,[]);
        %[pathName_thev1,folderName_thev1,fileName_thev1] = mostRecentCalibration('thev',outputA.label,1);
        %[pathName_thev2,folderName_thev2,fileName_thev2] = mostRecentCalibration('thev',outputA.label,2);
    elseif strcmp(testProbe,'B') % make mic correction for probe B (use probe A to present)
        input = inputB;
        output = outputA;
        [pathName_mic,folderName_mic,fileName_mic] = mostRecentCalibration('mic',inputB.label,[]);
        %[pathName_thev1,folderName_thev1,fileName_thev1] = mostRecentCalibration('thev',outputB.label,1);
        %[pathName_thev2,folderName_thev2,fileName_thev2] = mostRecentCalibration('thev',outputB.label,2);
    end    

    %----------------------------------------------------------------------
    % create a click stimulus with a high-frequency emphasis 
    % in order to compensate for high-frequency rolloff in the recievers.
    stimLen = 0.054; % length of stimulus
    stimN = round(stimLen * fs); % number of samples in stimulus
    stimulus = zeros(stimN,1); % initialize stimulus
    stimulus(200,1) = 1; % impulse is here
    % create the high-frequency emphasis
    N = fs/2;
    mag = zeros(N,1);
    mag(200:10000) = 1;
    dBGain = 15; % gain above 10 kHz
    gain = 10^(dBGain/20);
    mag(10000:34000) = gain;
    mag = [mag;flipud(mag)];
    phi = zeros(size(mag));
    X = mag.*cos(phi) + 1i*mag.*sin(phi);
    x = real(ifft(X));
    x = [x(N+1:end);x(1:N)];
    order = round(fs*0.001);
    if mod(order,2)~=0
        order = order + 1;
    end
    M2 = order/2;
    x = x(N-M2:N+M2);
    stimulus = fastFilter(x,stimulus); % filter to give the desired gain
    stimulus = stimulus / max(abs(stimulus)); % re-scale to full out
    stimulus = stimulus * .99; % keep from overdriving the system
    
    
    % -----------------------------------------------------------------
    obj.clearRecList % clear the previously used recording list
    obj.setRecList(input.ch,input.label,input.micSens,input.gain); % load the recording info for the test probe
    obj.setNReps(nRepsClick); % number of times to play stimulus
    obj.setFilter(1); % turn on default highpass filter
    obj.clearPlayList % clear the previously used playback list
    obj.setPlayList(stimulus,output.ch(1)); % load the playback probe
    obj.objPlayrec.run % playback and record
    if obj.killRun
       return
    end    
    [header,data] = obj.retrieveData(input.label); % get test probe recording
    
    disp('Insert reference microphone')
    pause
        
    obj.clearRecList % clear the previously used recording list
    obj.setRecList(inputRef.ch,inputRef.label,inputRef.micSens,inputRef.gain); % load the recording info for reference mic
    obj.objPlayrec.run % playback and record
    if obj.killRun
       return
    end    
    [headerRef,dataRef] = obj.retrieveData(inputRef.label); % get reference recording


    % -----------------------------------------------------------------
    % Calibration tube is 36 inches in length (36 * 2.54 = 91.44 cm)
    % Travel time in the tube is 91.44 / 34400 = 0.0027 s = 2.7 ms (one way)
    % Round-trip travel time = 2.7*2 = 5.4 ms.
    % Assume 10 round trip travel times are needed to allow the stimulus to
    % completely decay. The stimulus duration should be 5.4*10 = 54 ms.
    %
    % The recording should be time windowed to allow only the incident,
    % recorded waveform in. One travel time (1) is required for the
    % incident wave to reach the microphone. The first reflection travels
    % back to the source for a second travel time (2). The wave is then
    % re-reflected and travels back to the recording micrphone for a third
    % travel time (3). Therefore, the recording window should be a bit less
    % than 3 travel times. We will sue 2.5 travel times, or 
    % 2.7 * 2.5 = 6.75 ms.
    tubeLength = 91.44; % cm
    speedOfSound = 34400; % cm/s
    travelTime = tubeLength / speedOfSound; % s
    nTravelTimes = 2.5;
    winLen = nTravelTimes * travelTime; % length of time window (s)
    winN = round(winLen * fs); % length of time window (samples)
    [~,indxStart] = max(abs(stimulus)); % window indices
    indxFinish = indxStart + winN;
    probe = mean(data,2);
    probe = probe(indxStart:indxFinish);
    ref = mean(dataRef,2);
    ref = ref(indxStart:indxFinish);
    probe = ARLas_ramp(probe,fs,0.001);
    ref = ARLas_ramp(ref,fs,0.001);
    b = ARLas_makeMicCorrection_DW10x(probe,ref,fmin,fmax,fs); % microphone transfer function
    micCorrection = b;

    % -----------------------------------------------------------------
    tt = datetime('now'); % create the folder name for saving microphone calibrations
    folderName = datestr(tt,'mm_dd_yyyy');
    folderName = [folderName,'\']; 
    if exist([pathName_mic,folderName],'dir') ~= 7 % check to make sure that the save path exists
        success = mkdir([pathName_mic,folderName]); % if not, try to create it
        if success ~= 1
            warning('Specified folder for saving mic cal does not exist. Entering debug mode.')
            keyboard
        else
            addpath(genpath(pathName_mic)) % add all directories and subdirectories
        end
    end
    fileName_mic = [input.label,'.mat'];
    saveName = ARLas_saveName([pathName_mic,folderName],fileName_mic);
    save([pathName_mic,folderName,saveName],'micCorrection','fmin','fmax','fs','testProbe','folderName')

end

% Internal Functions ------------------------------------------------------

% OLD CODE ----------------------------------------------------------------
