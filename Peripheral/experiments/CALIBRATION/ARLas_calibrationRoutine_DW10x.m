function [] = ARLas_calibrationRoutine_DW10x(varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ARLas_calibrationRoutine_DW10x(varargin);
%
% Calibration routine for Lichtenhan lab. This code to be loaded into ARLas
% as an experiment and run. It does the following:
%   1. Performs Thevenin source calibration on the 10X with no mic correction.
%   2. Calculates microphone correction.
%   3. Performs Thevenin source calibration on the 10X with new mic correction.
%   4. Verifies final Thevenin and microphone calibrations.
%
% This version (DW10x) is optimized for calibrations in the Lichtenhan lab:
% Calibration is done in a hollow ear bar from 200 - 32000 Hz.
% Custom stainles tube set used for Thevenin source calibration.
% A custom brass coupler is used to get mic correction and to verify.
% A 1/8" condensor mic is used as the reference.
%
% For these calibrations, data in a folder named
% 'CalFiles'. Each calibration is saved in a sub-directory named by the
% date in the form 'mm_dd_yyyy'.
%
% Authors: Shawn S Goodman, PhD & Jeffery T Lichtenhan PhD
% Auditory Research Lab, Dept. of Comm. Sciences & Disorders, The University of Iowa
% Auditory Physiology Lab, Dept. of Otolaryngology, Washington University School of Medicine, St. Louis
% Date: October 12, 2019
% Updated: October 12-22, 2019 -- ssg
% Updated: February 4, 2020 -- ssg
%                              Clarified in comments  that doNew calibration type does
%                              not work well and should be avoided.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    obj = varargin{1};

% ---  IMPORTANT! UPDATE THE FOLLOWING LINES OF CODE BEFORE RUNNING!!   ---
    testProbe = 'A';        % set to 'A' or 'B';
    cavityTemperature = 26; % 30.5 (C); read this off of the 10X and update as needed    
    fmin = 200;     % minimum frequency (Hz) over which to calculate fit
    fmax = 32000;   % maximum frequency (Hz) over which to calculate fit
                            % NOTE: These are different from the stimulus chirp range.
                            % The above values allow you to perform a calibration over a sub-range of the
                            % frequencies that are in the calibration chirp, which always goes from 200-32000 Hz.
                            % In the 10x internal cavity, the calibration falls apart above 22000 Hz.
                            % In the small stainless tubes, the calibration is good up to 32000 Hz.
    doBackup = 1;   % turn on (1) and off (0) backup of calibration files.
    doMicCal = 1;   % turn on (1) and off (0) microphone calibration.
    doThevCal = 1;  % turn on (1) and off (0) Thevenin calibration.
    doVerify = 1;   % turn on (1) and off (0) verification of calibration.
                    % NOTE: Thevenin calibrations must exist in order to do
                    % mic calibration. 
    doNew = 0; % new mic correction type. Does not work; use old (set to 0)
% -------------------------------------------------------------------------

    params.fmin = fmin;
    params.fmax = fmax;
    params.testProbe = testProbe;
    params.cavityTemperature = cavityTemperature;
    params.folderName_mic = []; % no file names or folders typically given to start
    params.folderName_thev1 = [];
    params.folderName_thev2 = [];
    params.fileName_mic = []; % no file names or folders typically given to start
    params.fileName_thev1 = [];
    params.fileName_thev2 = [];
    
    disp(' '); disp(' '); disp(' '); disp(' ')
    alertTxt = {'Starting DW Calibration Routine.'
         ['  Calibrating ER10X probe ',testProbe,'.']
         '..................................'
        };
    nn = size(alertTxt,1);
    for ii=1:nn
        cprintf([0,0,.4],[alertTxt{ii},'\n']);
    end
    disp(' ')    


% 0. Back up existing calibration files. ----------------------------------
if doBackup == 1
    ARLas_calibrationBackup % backup the old calibration files
end

if doNew == 1
    ARLas_earBarCorrection_DW10x(obj,params)
end

% 2. Perform Thevenin source calibration on the 10X WITH mic correction, if it exists.
if doThevCal == 1
    alertTxt = {'Performing Thevenin Source calibration.'
        };
    nn = size(alertTxt,1);
    for ii=1:nn
        cprintf([0,0,.4],[alertTxt{ii},'\n']);
    end
    disp(' ')
    params.applyMicCorrection = 0;
    ARLas_cavityRecordings_DW10x(obj,params) % perform Thevenin calibration
end

% 1. Calculate microphone correction.--------------------------------------
if doMicCal == 1
    if strcmp(testProbe,'A')
        playProbe = 'B';
    elseif strcmp(testProbe,'B')
        playProbe = 'A';
    else
        error('Unrecognized value for testProbe: must be A or B.')
    end
    alertTxt = {'Performing calibration to find microphone correction.'
        ['  Place ER10X probe ',testProbe,' in the end of the long brass tube with a clear connector.']
        ['  Place ER10X probe ',playProbe,' in the end of the long brass tube with a colored connector.']
        '  Ensure that the Reference microphone is turned on.'
        };
    nn = size(alertTxt,1);
    for ii=1:nn
        cprintf([0,0,.4],[alertTxt{ii},'\n']);
    end
    disp(' ')    
    
    showInstructions
    params.applyMicCorrection = 0;
    params.createMicCorrection = 1;
    %removeMicCorrection(testProbe)
    ARLas_couplerRecordings_DW10x(obj,params); % make recordings in coupler
    %ARLas_micCorrection_DW10x(obj,params)
    
    alertTxt = {'Applying microphone correction to Thevenin Source calculations.'
        '  This process may take a few minutes.'
        '  Please wait patiently (or wait impatiently, if you prefer!)...'
        };
    nn = size(alertTxt,1);
    for ii=1:nn
        cprintf([0,0,.4],[alertTxt{ii},'\n']);
    end
    disp(' ')    
    insertMicCorrection(testProbe) % put in the mic correction, recalculate, and save
end


% 3. Verify Thevenin and microphone calibrations. -------------------------   
if doVerify == 1
    alertTxt = {'Performing Coupler verification of Thevenin Source and microphone calibrations.'
        '  Place ER10X probe in the brass coupler.'
        '  Ensure that the Reference microphone is also attached and turned on.'
        };
    nn = size(alertTxt,1);
    for ii=1:nn
        cprintf([0,0,.4],[alertTxt{ii},'\n']);
    end
    disp(' ')    
    showInstructions
    params.applyMicCorrection = 1;
    params.createMicCorrection = 0;
    ARLas_couplerRecordings_DW10x(obj,params);
end

% 4. Double check system distortion and noise floor. ----------------------
    % NOT ACTUALLY IMPLEMENTED HERE
    alertTxt = {'..................................'
         'Finished DW Calibration Routine'
        };
    nn = size(alertTxt,1);
    disp(' ')
    for ii=1:nn
        cprintf([0,0,.4],[alertTxt{ii},'\n']);
    end
    disp(' '); disp(' ');
    
    
end

% internal functions ------------------------------------------------------
function [] = insertMicCorrection(testProbe)
    % put the microphone correction into the Thevenin calibration, re-calculate, and re-save.
    [inputs,outputs] = hardwareSetup; %hardwareSetup; % read in the saved hardware setup
    if strcmp(testProbe,'A')
        input = inputs{1};        % ER10xA microphone
        output = outputs{1};      % ER10xA loudspeakers
    elseif strcmp(testProbe,'B')
        input = inputs{2};        % ER10xB microphone
        output = outputs{2};      % ER10xB loudspeakers
    end
    inputRef = inputs{4};         % GRAS 1/8" reference microphone
   
    [pathName_mic,folderName_mic,fileName_mic] = mostRecentCalibration('mic',input.label,[]);
    [pathName_thev1,folderName_thev1,fileName_thev1] = mostRecentCalibration('thev',output.label,1);
    [pathName_thev2,folderName_thev2,fileName_thev2] = mostRecentCalibration('thev',output.label,2);
    
    dummy = load([pathName_mic,folderName_mic,fileName_mic]); % load the microphone recording
    micCorrection = dummy.micCorrection;

    dummy = load([pathName_thev1,folderName_thev1,fileName_thev1]); % load the thev cal
    t = dummy.t;
    t.applyMicCorrection = 1;
    t.micCorrection = micCorrection; % write in the new mic correction value
    t.calculate; % recalculate thevenin source
    save([pathName_thev1,folderName_thev1,fileName_thev1],'t') % save file (overwrites old value)

    dummy = load([pathName_thev2,folderName_thev2,fileName_thev2]); % load the thev cal
    t = dummy.t;
    t.applyMicCorrection = 1;
    t.micCorrection = micCorrection;
    t.calculate;
    save([pathName_thev2,folderName_thev2,fileName_thev2],'t')

end
function [] = removeMicCorrection(testProbe)
    % put the microphone correction into the Thevenin calibration, re-calculate, and re-save.
    [inputs,outputs] = hardwareSetup; %hardwareSetup; % read in the saved hardware setup
    if strcmp(testProbe,'A')
        input = inputs{1};        % ER10xA microphone
        output = outputs{1};      % ER10xA loudspeakers
    elseif strcmp(testProbe,'B')
        input = inputs{2};        % ER10xB microphone
        output = outputs{2};      % ER10xB loudspeakers
    end
    inputRef = inputs{4};         % GRAS 1/8" reference microphone
   
    [pathName_mic,folderName_mic,fileName_mic] = mostRecentCalibration('mic',input.label,[]);
    [pathName_thev1,folderName_thev1,fileName_thev1] = mostRecentCalibration('thev',output.label,1);
    [pathName_thev2,folderName_thev2,fileName_thev2] = mostRecentCalibration('thev',output.label,2);
    
    if isempty(fileName_mic) % if there is no microphone correction, then this step is unessary
        return
    end
    dummy = load([pathName_mic,folderName_mic,fileName_mic]); % load the microphone recording
    micCorrection = dummy.micCorrection;

    dummy = load([pathName_thev1,folderName_thev1,fileName_thev1]); % load the thev cal
    t = dummy.t;
    if length(t.micCorrection) ~= 1 % if there is a mic correciton written,
        t.applyMicCorrection = 1;
        t.micCorrection = 1; % write in no mic correction value
        t.calculate; % recalculate thevenin source
        save([pathName_thev1,folderName_thev1,fileName_thev1],'t') % save file (overwrites old value)
    end

    dummy = load([pathName_thev2,folderName_thev2,fileName_thev2]); % load the thev cal
    t = dummy.t;
    if length(t.micCorrection) ~= 1 % if there is a mic correciton written,
        t.applyMicCorrection = 1;
        t.micCorrection = 1;
        t.calculate;
        save([pathName_thev2,folderName_thev2,fileName_thev2],'t')
    end

end
function [] = showInstructions()
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
        set(h1,'CData',imread('refCoupler.jpg'))
    txt = ({['1. Place ER10X probe and Ref mic into brass coupler.'];'';...
        '2. Press OK when ready.'});
    uiwait(msgbox(txt,'Calibration Routine'));
    try
        delete(h)
    catch
    end % end instructions
end

% OLD CODE
