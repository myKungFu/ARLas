function [] = ARLas_cavityRecordings_IHS(varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ARLas_cavityRecordings_IHS(varargin);
%
% Make recordings from ER10x calibrator for calculation of Thevenin source
% parameters. This version (ER10X) is optimized for calibrations in humans.
% It is assumed that you are using the built-in ER10X calibrator tube OR
% a set of custom stainless tubes 0.8 cm in diameter, of the lengths [2.9; 3.6; 4.15; 5.2; 6.9] cm
%
% For these calibrations, data in a folder named
% 'CalFiles'. Each calibration is saved in a sub-directory named by the
% date in the form 'mm_dd_yyyy'.
%
% Authors: Shawn S Goodman, PhD
% Auditory Research Lab, Dept. of Comm. Sciences & Disorders, The University of Iowa
% Original Date: January 8-10, 2019
% Updated: October 10-22, 2019 -- ssg
% Updated: May 28, 2021 -- ssg
% Updated: July 1, 2021 -- ssg -- fixed bad function call on channel 2 (EB)
%                                 fixed the way that mic cal is applied
% Updated: January 29, 2025 -- ssg -- updated for IHS stuff
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    obj = varargin{1};
    params = varargin{2};
    V = '20JAN2025'; % this is the current version number of this program


    % unpack the paramter structure
    applyMicCorrection = params.applyMicCorrection;
    fmin = params.fmin;
    fmax = params.fmax;
    testProbe = params.testProbe;
    cavityTemperature = params.cavityTemperature;
    
    nReps = 64; %128; % number of stimulus repetitions; recommend ~48
    cavityLengths = [2.9; 3.6; 4.15; 5.2; 6.9]; % cavity lengths (cm) -- the custom set is fixed at these lengths
    % extra value is added to the desired length to account for insertion of the tip 
    adjust = -0.4; % (cm); reduction due to insertion of the probe tip
    cavityLengths = cavityLengths + adjust;

    cavityDiameter = 0.8; % ER10x cavity diameter (cm)
    stimGain_dB = 0; % reduce so you don't overdrive the system; 
                      % -5 dB works for ER10X with 20 dB input gain and no
                      % output gain. If you use amplifiers, may need to reduce
                      % the output by more. If you disable the output
                      % limiter, set this value to -20.
    if params.outputLimiterEnabled == 0 % if limiter is disabled, 
        stimGain_dB = stimGain_dB; % -15; % reduce output by another 15 dB -- NOT FOR IHS
    end
    
    [inputs,outputs] = hardwareSetup; %hardwareSetup; % read in the saved hardware setup
    if strcmp(testProbe,'B')
        input = inputs{2};        % ER10xA microphone
        output = outputs{2};      % ER10xA loudspeakers
    % elseif strcmp(testProbe,'B')
    %     input = inputs{2};         % ER10xB microphone
    %     output = outputs{2};       % ER10xB loudspeakers
    end

    [pathName_mic,folderName_mic,fileName_mic] = mostRecentCalibration('mic',input.label,[],obj.map);
    [pathName_thev1,folderName_thev1,fileName_thev1] = mostRecentCalibration('thev',output.label,1,obj.map);
    [pathName_thev2,folderName_thev2,fileName_thev2] = mostRecentCalibration('thev',output.label,2,obj.map);
    
    % -------------------------------------------------------------------------
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

    fs = obj.fs; % sampling rate (Hz)
    [stimulus,F,cut1,cut2] = getStimulus(fs); % get logarithmic chirp
    fminStim = min(F); % minimum frequency in the stimulus
    fmaxStim = max(F); % maximum frequency in the stimulus
    stimulus = adjustOutput(F,stimulus); % adjust for better SNR
    stimulus = stimulus * 10^(stimGain_dB/20);


    cavityLengths = cavityLengths(:)'; % force these to be a column vector
    nCavities = length(cavityLengths); % loop over number of cavity lengths specified
    nChannelsOut = length(output.ch); % number of output channnels being calibrated (typically 2)
    for jj=1:nCavities % loop over number of cavities
        extra = -adjust;
        showInstructions(cavityLengths,extra,jj);
        for kk=1:nChannelsOut % loop over output channels; test output channels one at a time
            % LOAD THE STIMULUS ----------------------------------------------------
            obj.clearRecList % clear the previously used recording list
            obj.setRecList(input.ch,input.label,input.micSens,input.gain); % load the recording info for ARLas to use
            obj.setNReps(nReps); % number of times to play stimulus
            obj.setFilter(1); % turn on default highpass filter
            obj.objPlayrec.userInfo.caps_version = V;


            % PLAYBACK & RECORD ----------------------------------------------------
            obj.clearPlayList % clear the previously used playback list
            obj.setPlayList(stimulus,output.ch(kk)); % load the currently tested output channel

            obj.objPlayrec.run % run the stimulus
            if obj.killRun
               return
            end    
            
            % RETRIEVE DATA ----------------------------------------------------
            [header,data] = obj.retrieveData(input.label(kk)); % get raw data
            if applyMicCorrection == 1
                [Rows,Cols] = size(data);
                data = data(:);
                data = fastFilter(micCorrection,data);
                data = reshape(data,Rows,Cols);                
            end
            if kk == 1
                [thevRecording1(:,jj),stimf] = cleanData(data,obj.fs,stimulus,header,cut1,cut2);
            elseif kk == 2
                [thevRecording2(:,jj),stimf] = cleanData(data,obj.fs,stimulus,header,cut1,cut2);
            end
        end
    end

    % save recording data -------- --------------------------------------------
    recordingParams.fs = obj.fs; % sampling rate (Hz)
    recordingParams.cavityTemperature = cavityTemperature; % degrees C
    recordingParams.cavityDiameter = cavityDiameter; % (cm)
    recordingParams.cavityLengths_nominal = cavityLengths;
    recordingParams.fmin = fmin; % calibrate over this range: fmin-fmax
    recordingParams.fmax = fmax;
    recordingParams.timeStamp = datestr(clock); % time stamp for when this calibration was done
    recordingParams.delay_now = obj.objInit.delay_now;
    recordingParams.card2volts_now = obj.objInit.card2volts_now;
    recordingParams.hostAPI_now = obj.objInit.hostAPI_now;
    recordingParams.input = input;
    recordingParams.output = output;
    recordingParams.thevCalPathName = obj.objPlayrec.savedFilesPath;
    recordingParams.cavityRecordingsPathName = obj.objPlayrec.savedFilesPath;
    recordingParams.fminStim = fminStim; % minimum frequency in calibration stimulus
    recordingParams.fmaxStim = fmaxStim; % maximum frequency in calibration stimulus
    recordingParams.stimOrig = stimulus; % stimulus presented to the sound card DAC
    recordingParams.cut1 = cut1; % cut sample for start of stimulus
    recordingParams.cut2 = cut2; % cut sample for end of stimulus
    recordingParams.stimFilt = stimf; % stimulus filtered and cut to size
    recordingParams.stimGain_dB = stimGain_dB; % gain used (dB re: full out)
    recordingParams.micCorrection = micCorrection; % the microphone correction filter
    recordingParams.params = params;
    recordingParams.F = F; % vector of instantaneous frequencies (Hz)
    stimulus = stimf; % stimulus is the filtered version of the stimulus
    
    channel = 1;
    recordings = thevRecording1;
    pathName = recordingParams.cavityRecordingsPathName;
    fileName1 = [recordingParams.output.label,'_Ch',num2str(channel),'.mat'];
    fileName1 = ARLas_saveName(pathName,fileName1);
    save([pathName,fileName1],'recordingParams','stimulus','recordings');
    if nChannelsOut == 2
        channel = 2;
        recordings = thevRecording2;
        pathName = recordingParams.cavityRecordingsPathName;
        fileName2 = [recordingParams.output.label,'_Ch',num2str(channel),'.mat'];
        fileName2 = ARLas_saveName(pathName,fileName2);
        save([pathName,fileName2],'recordingParams','stimulus','recordings');
    end
    
    % perform Thevenin calculations ---------------------------------------
    alertTxt = {'Making Thevenin Source calculations.'
        '  This process may take a few minutes.'
        '  Please wait patiently (or wait impatiently, whatever :) ...'
        };
    nn = size(alertTxt,1);
    for ii=1:nn
        cprintf([0,0,.4],[alertTxt{ii},'\n']);
    end
    disp(' ')    
    
    
    t = datetime('now'); % create the folder name for saving calibrations
    folderName = datestr(t,'mm_dd_yyyy');
    fsep = filesep;
    folderName = [folderName,fsep]; 
    if exist([pathName_thev1,folderName],'dir') ~= 7 % check to make sure that the save path exists
        success = mkdir([pathName_thev1,folderName]); % if not, try to create it
        if success ~= 1
            error('Specified folder for saving data does not exist. Entering debug mode.')
            keyboard
        else
            addpath(genpath(pathName_thev1)) % add all directories and subdirectories
        end
    end
    
    dummy = load([pathName,fileName1]);
    t = thevCal_new(dummy.recordingParams,dummy.stimulus,dummy.recordings);
    t.fmin = fmin;
    t.fmax = fmax;
    t.calculate
    fileName1 = ARLas_saveName([pathName_thev1,folderName],fileName1);
    save([pathName_thev1,folderName,fileName1],'t')
    
    if nChannelsOut == 2
        dummy = load([pathName,fileName2]);
        t = thevCal_new(dummy.recordingParams,dummy.stimulus,dummy.recordings);
        t.fmin = fmin;
        t.fmax = fmax;
        t.calculate
        fileName2 = ARLas_saveName([pathName_thev2,folderName],fileName2);
        save([pathName_thev2,folderName,fileName2],'t')
    end
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
    output = [ 60, 55, 50, 50,  40,  40,  40,  40,   40,   40,   60,   60,   60]'; 
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

function [] = showInstructions(cavityLengths,extra,jj)
    gui.height = 200;
    gui.width = 497;    

    scrsz = get(0,'ScreenSize'); % get the current screen size
    gui.left = round(scrsz(4) * .05);
    gui.bottom = scrsz(4) * .1;    
    %[left, bottom,width, height]
    rect = [gui.left,gui.bottom+125,gui.width,gui.height]; % give instructions
        h = figure('Position',rect);
        set(h,'Name',' ','NumberTitle','off','MenuBar','none','Resize','off','Color',[1 1 1])
        h1 = uicontrol('Parent',h,'Style','togglebutton','BackgroundColor',[1 1 1],...
            'Position',[1 1 gui.width gui.height],'Visible','on');
        set(h1,'CData',imread(['thevCavities_',num2str(jj),'.jpg']))
    txt = ({['1. Place probe into tube of length = ',(num2str((cavityLengths(jj)+extra)*10)),' mm.'];'';...
        '2. Press OK when ready.'});
    uiwait(msgbox(txt,'Thevenin Calibration'));
    try
        delete(h)
    catch
    end % end instructions
end

% OLD CODE ----------------------------------------------------------------

% function [stim,F,cut1,cut2] = getStimulus(fs)
%     fmin = 100; % starting frequency (hz)
%     fmax = 24000; % ending frequency (Hz); this is 1.2 times the desired ending frequency. 
%                   % fs = 96000 Hz; so this is still below nyquist (48000 Hz)
%     dur = 0.150; % desired chirp duration is 125 ms
%     N = round(fs * dur); % number of samples in recording
%     X = (0:1:N-1)';
%     X = X / max(X);
%     % formula = fmin*exp(-a*X)
%     % X(end) = 1, so the value of a is as follows:
%     % fmax = 100*exp(-aX)
%     % fmax / 100 = exp(-a)
%     % -log(fmax / 100) = a
%     a = -log(fmax / fmin);
%     F = 100*exp(-a*X);
%     lead = round(0.5/F(1)*fs); % lead samples for half a cycle
%     lag = round(0.25/F(end)*fs); % lead samples for a quarter cycle
%     F = [ones(lead,1)*F(1);F;ones(lag,1)*F(end)]; % instantateous frequencies
%     % Rotate the phase:
%     % At each time sample, time has advanced by the sampling period (1/fs).
%     % At each time sample, the phase has rotated by the frequency (f) times
%     % the sampling period: phase advance (in cycles) = f * (1/fs), or f/fs.
%     % To make this a radian change, multiply by 2pi: phi = (2*pi*f)/fs;
%     % Finally, add the newly accumulated phase to the previous phase value.
%     stim = zeros(length(F),1);
%     phi = 0; % starting phase is zero
%     for jj = 1:length(F)
%         stim(jj,1) = sin(phi);
%         phi = phi + ((2*pi*F(jj))/fs);
%     end
% 
%     rampN = round(lead*2); % ramp stimulus on
%     h = hann(rampN*2);
%     h = h(1:rampN);
%     stim(1:rampN) = stim(1:rampN) .* h;
%     rampN = round(lead*.5); % ramp stimulus off
%     h = hann(rampN*2);
%     h = h(1:rampN);
%     stim = flipud(stim);
%     stim(1:rampN) = stim(1:rampN) .* h;
%     stim = flipud(stim);
%     % zero pad stimulus with 5 ms on both sides. This is to allow for
%     % filtering later. Define cut1 and cut2 so you can remove padding after
%     % recording
%     extra = round(fs*0.002);
%     cut2 = length(stim) + extra;
%     padLen = 0.050; % zero-pad the stimulus with 5 ms on both sides
%     padN = round(padLen * fs);
%     cut1 = padN - extra;
%     cut2 = cut2 + padN;
%     pad = zeros(padN,1);
%     stim = [pad;stim;pad];
%     F = [pad;F;pad];
%     stim = stim / max(abs(stim));  % rescale so 1 is max out
%     stim = stim * .99;
%     % There are extra frequencies added to the stimulus to ensure that the
%     % desired range is met. The actual minimum and maximum frequencies that
%     % can be tested accurately are recoded here:
%     %fmin = 200;   % minimum testable frequency (Hz)
%     %fmax = 32000; % maximum testable frequency (Hz)
% end
% function [stim] = adjustOutput(F,stim)
%    % adjust the output amplitudes to give approximately flat output across
%     % the frequency range. These values based on "average" 10x driver
%     % response, re: the supplied operating instructions manual from ER.
%     % (plus a few emperical tweaks)
%     fff =    [100,200,400,600,1000,2000,4000,6000,10000,12000,16000,20000,30000,45000]'; % frequency (Hz)
%     output = [ 72, 55, 50, 50,  40,  40,  40,  40,   50,   65,   70,   75,   80,   80]';
%     %output = [ 72, 55, 50, 50,  40,  40,  40,  40,   50,   60,   65,   65,   65,   70]';
%     relativeOut = (output - max(output));
%     a = zeros(length(stim),1);
%     for jj=1:length(stim)
%         a(jj,1) = interp1(fff,relativeOut,F(jj),'pchip'); % relative output amplitude (dB)
%     end
%     a = 10.^(a/20); % put a into linear units
%     stim = a .* stim; % amplitude adjusted stimulus
%     stim = stim / max(abs(stim));
%     stim = stim * .99;
% end





%     if exist(pathName_thev1,'dir') ~= 7 % check to make sure that the save path exists
%         success = mkdir(pathName_thev1); % if not, try to create it
%         if success ~= 1
%             error('Specified path for saving data does not exist. Entering debug mode.')
%             keyboard
%         else
%             addpath(genpath(savePath)) % add all directories and subdirectories
%         end
%     end
%     % ensure the path names end with \
%     if strcmp(savePath(end),'\') == 0
%        savePath = [savePath,'\']; 
%     end
%     savedFile1 = saveData(recordingParams,stimulus,thevRecording1,channel,recordingParams.cavityRecordingsPathName);
%     if nChannelsOut == 2
%         channel = 2;
%         savedFile2 = saveData(recordingParams,stimulus,thevRecording2,channel,recordingParams.cavityRecordingsPathName);
%     end
    %saveThev(recordingParams.thevCalPathName,savedFile1,t)
    %saveThev(recordingParams.thevCalPathName,savedFile2,t)
    %disp('Finished runing Thevenin calibration.')
    %disp('.....................................')
    % location to save Thevenin calibration to
    %savePath = 'C:\myWork\ARLas\Peripheral\calibrations\thevCals\'; 
    % location of mic calibration files
    %pathName_micCorr = 'C:\myWork\ARLas\Peripheral\calibrations\micCals\'; 
% function [] = saveThev(pathName,fileName,obj)
%     % assume fileName ends with file extension
%     fileName = fileName(1:end-4); % remove the extension
%     counter = 1;
%     while exist([pathName,fileName,'_thevCal_',num2str(counter),'.mat'],'file') == 2
%         counter = counter + 1;
%     end
%     t = obj;
%     try
%         save([pathName,fileName,'_thevCal_',num2str(counter),'.mat'],'t')
%         OK = 1;
%     catch ME
%         OK = 0;
%          errorTxt = {'  Issue: Error saving data.'
%              '  Action: Save aborted.'
%              '  Location: in thevCal.saveThev.'
%             };
%         errorMsgARLas(errorTxt);
%     end
%     if OK
%          alertTxt = {'  Issue: Calibration.'
%              '  Action: Calibration successfully saved.'
%              '  Location: in thevCal.saveThev.'
%             };
%         alertMsgARLas(alertTxt);
%     end
% end
% function [savedFile] = saveData(recordingParams,stimulus,recordings,channel,pathName)
%     if nargin < 5
%         pathName = recordingParams.cavityRecordingsPathName;
%     end
%     if ~exist(pathName,'dir')
%         mkdir(pathName)
%         addpath(genpath(pathName))
%     end
%     try
%         [fileName,pathName] = uiputfile([pathName,recordingParams.output.label,'_Ch',num2str(channel),'.mat'],'Save FPL Calibration.');
%         if fileName == 0 % user chose not to create a file
%             dlgtxt = {'You chose NOT to save the current calibration recordings.' ' ' 'Are you sure you want to do this?'};
%             choice = questdlg(dlgtxt,'Save Initialization','Yes','No','No');
%             if strcmp(choice,'Yes') == 1
%                 return
%             else % give one more chance
%                 [fileName,pathName] = uiputfile([pathName,recordingParams.output.label,'_Ch',num2str(channel),'.mat'],'Save FPL Calibration.');
%             end
%         end
%         savedFile = fileName;
%         if isempty(fileName) % user chose not to create a file
%             return
%         end
% 
%         if exist(pathName,'dir') == 0 % if correct folder does not exist
%             success = mkdir(pathName); % try to create it
%             if success ~= 1
%                 errorTxt = {'  Issue: Error creating the directory for saving cavity recordings.'
%                      '  Action: Cavity recordings NOT saved.'
%                      '  Location: cavityRecordings.'
%                     };
%                 errorMsgARLas(errorTxt);
%                 return
%             end
%         end
%         try % save file here ----------------------------------------------
%             recordingParams.fileName = fileName;
%             save([pathName,fileName],'recordingParams','stimulus','recordings');
%         catch
%             alertTxt = {'  Issue: Error saving cavity recordings.'
%                  '  Action: Cavity recordings NOT saved.'
%                  '  Location: cavityRecordings.'
%                 };
%             alertMsgARLas(alertTxt);
%         end
%     catch ME
%         errorTxt = {'  Issue: Error saving cavity recordings.'
%              '  Action: Cavity recordings NOT saved..'
%              '  Location: cavityRecordings.'
%             };
%         errorMsgARLas(errorTxt);
%         objInit.printError(ME)
%     end
% end

   
    

