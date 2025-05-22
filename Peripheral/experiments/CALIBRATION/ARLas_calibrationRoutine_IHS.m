function [] = ARLas_calibrationRoutine_IHS(varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ARLas_calibrationRoutine_IHS(varargin);
%
% Calibration routine for Etymotic ER10X. This code to be loaded into ARLas
% as an experiment and run. It does the following:
%   0. Backs up old calibration files (if requested to do so)
%   1. Calculates microphone correction (either using long tube or data saved in an Excel sheet)
%   2. Performs Thevenin source calibration on the 10X with new mic correction.
%
% This version (ER10x) is optimized for calibrations in humans from 100 - 20 kHz.
% Calibrator tube in ER10X is used for Thevenin source calibration.
%
% All of these calibrations are saved in 'C:\myWork\ARLas\Peripheral\calibrations\'
% under the appropriate sub-directory (micCals, or thevCals).
% Each calibration is further saved in a sub-directory named by date in the form 'mm_dd_yyyy'.
%
% NOTE: calibration for long tube (i.e. incident pressure calibration is
% performed separatly using the experiment file ARLas_longTubeCal.m, and
% these files are saved in a separate sub-directory (LTCals).
%
% Author: Shawn S Goodman, PhD
% Auditory Research Lab, Dept. of Comm. Sciences & Disorders, The University of Iowa
% Original Date: October 12, 2019
% Updated: October 12-22, 2019 -- ssg
% Updated: February 4, 2020 -- ssg
% Updated: May 28, 2021 -- ssg
% Updated: June 1, 2021 -- ssg
% Updated: July 5, 2021 -- ssg -- created setting for output limiter
% Updated: September 17, 2021 -- ssg -- updated to run new version of ARLas
%                               micCal (version _v2).
% Updated: September 22, 2021 -- ssg -- added mic correction to thev cal,
%                               assuming mic correction already exits]
% Updated: October 19, 2021 -- ssg -- added third version of micCal, which
%                               uses data provided by Etymotic research rather 
%                               than our own recordings.
% Updated: March 15, 2021 -- ssg -- updated mic calibration routine to
%                               bring instructions in line with previous update on Oct 19.
% Updated: January 29, 2025 -- ssg -- updated for IHS stuff
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    obj = varargin{1};
    V = '28APR2025'; % version 
    

%------ USER MODIFIABLE PARAMETERS ----------------------------------------
%--------------------------------------------------------------------------
    testProbe = 'A';        % set to 'A' or 'B';
    outputLimiterEnabled = 0; % 1 = limiter enabled; 0 = disabled (15 dB higher output)
    
    cavityTemperature = 32; % ~34(C); read this off of the 10X and update as needed
    fmin = 100;     % minimum frequency (Hz) over which to calculate fit
    fmax = 20000;   % maximum frequency (Hz) over which to calculate fit
                            % NOTE: These are different from the stimulus chirp range.
                            % The above values allow you to perform a calibration over a sub-range of the
                            % frequencies that are in the calibration chirp, which always goes from 200-20000 Hz.
                            % In the 10x internal cavity, the calibration falls apart above 22000 Hz.
    doBackup = 0;   % turn on (1) and off (0) backup of calibration files.
    doMicCal = 0;   % turn on (1) and off (0) microphone calibration.
    doThevCal = 1;  % turn on (1) and off (0) Thevenin calibration.
                    
    which10X = '[]'; % Specify which ER10X probe Microphone to "calibrate".
                     % Choose from 'FX','JL', or 'SG'.
                     % Note that this requires the spreadsheet 'MicData.xlsx', 
                     % which should be saved in the folder
                     % C:\myWork\ARLas\Peripheral\experiments\CALIBRATION\
                    
%------ END USER MODIFIABLE PARAMETERS ------------------------------------
%--------------------------------------------------------------------------

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
    params.outputLimiterEnabled = outputLimiterEnabled;
    params.which10X = which10X;
    
    disp(' '); disp(' '); disp(' '); disp(' ')
    alertTxt = {'Starting ER10X Thevenin Calibration Routine.'
         ['  Calibrating ER10X probe ',testProbe,'.']
         '..................................'
        };
    nn = size(alertTxt,1);
    for ii=1:nn
        cprintf([0,0,.4],[alertTxt{ii},'\n']);
    end
    disp(' ')    


% 0. Back up existing calibration files. ----------------------------------
    % Backups will be kept in myWork\ARLas_oldCalibrations\
    if doBackup == 1
        ARLas_calibrationBackup % backup the old calibration files
    end

% 1. Calculate microphone correction.--------------------------------------
    if doMicCal == 1
        alertTxt = {'Performing microphone correction routine.'
            ['  This version uses the values provided by Etymotic.']
            ['  Note that this corrects for magnitude only--no phase correction is provided!']
            '  No actual measurement is neede for this routine; sit back and relax :-)'
            };
        nn = size(alertTxt,1);
        for ii=1:nn
            cprintf([0,0,.4],[alertTxt{ii},'\n']);
        end
        disp(' ')        

        ARLas_micCal_v3(obj,params) % perform calibration
        % Create microphone calibration using Etymotic's provided mic calibrations.
        % The mic calibration excel file was provided by Steve Viranyi on Sept 29,
        % 2021 in an email sent to Shawn, Jeff, Vicky, and Dan.
        % This version (_v3) reads in the data from the excel file. The Excel file
        % contains microphone corrections (amplitude only--no phase corrections!)        

        alertTxt = {'Applying microphone correction to Thevenin Source calculations.'
            '  This process may take a few minutes.'};
        nn = size(alertTxt,1);
        for ii=1:nn
            cprintf([0,0,.4],[alertTxt{ii},'\n']);
        end
        disp(' ')    
        insertMicCorrection(testProbe,obj.map) % put in the mic correction, recalculate, and save
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
        params.applyMicCorrection = 0; % note: thevenin calibration is performed without the mic calibration on!
        ARLas_cavityRecordings_IHS(obj,params) % perform Thevenin calibration

        insertMicCorrection(testProbe,obj.map) % put in the mic correction, recalculate, and save
    end
    
end

% internal functions ------------------------------------------------------
function [] = insertMicCorrection(testProbe,map)
    % put the microphone correction into the Thevenin calibration, re-calculate, and re-save.
    [inputs,outputs] = hardwareSetup; %hardwareSetup; % read in the saved hardware setup
    if strcmp(testProbe,'A')
        input = inputs{1};        % ER10xA microphone
        output = outputs{1};      % ER10xA loudspeakers
    elseif strcmp(testProbe,'B')
        input = inputs{2};        % ER10xB microphone
        output = outputs{2};      % ER10xB loudspeakers
    end
    %inputRef = inputs{4};         % GRAS 1/8" reference microphone
   
    [pathName_mic,folderName_mic,fileName_mic] = mostRecentCalibration('mic',input.label,[],map);
    [pathName_thev1,folderName_thev1,fileName_thev1] = mostRecentCalibration('thev',output.label,1,map);
    [pathName_thev2,folderName_thev2,fileName_thev2] = mostRecentCalibration('thev',output.label,2,map);
    
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

% OLD CODE ----------------------------------------------------------------

%     
% keyboard
% 
% if doMicCal == 1
%     if strcmp(testProbe,'A')
%         playProbe = 'B';
%     elseif strcmp(testProbe,'B')
%         playProbe = 'A';
%     else
%         error('Unrecognized value for testProbe: must be A or B.')
%     end
%     alertTxt = {'Performing calibration to find microphone correction.'
%         ['  Place ER10X probe ',testProbe,' in the end of the long brass tube with a clear connector.']
%         ['  Place ER10X probe ',playProbe,' in the end of the long brass tube with a colored connector.']
%         '  Ensure that the Reference microphone is turned on.'
%         };
%     nn = size(alertTxt,1);
%     for ii=1:nn
%         cprintf([0,0,.4],[alertTxt{ii},'\n']);
%     end
%     disp(' ')    
%     
%     showInstructions
%     params.applyMicCorrection = 0;
%     params.createMicCorrection = 1;
%     %removeMicCorrection(testProbe)
%     ARLas_couplerRecordings_DW10x(obj,params); % make recordings in coupler
%     %ARLas_micCorrection_DW10x(obj,params)
%     
%     alertTxt = {'Applying microphone correction to Thevenin Source calculations.'
%         '  This process may take a few minutes.'
%         '  Please wait patiently (or wait impatiently, if you prefer!)...'
%         };
%     nn = size(alertTxt,1);
%     for ii=1:nn
%         cprintf([0,0,.4],[alertTxt{ii},'\n']);
%     end
%     disp(' ')    
%     insertMicCorrection(testProbe) % put in the mic correction, recalculate, and save
% end
% % 3. Verify Thevenin and microphone calibrations. -------------------------   
% if doVerify == 1
%     alertTxt = {'Performing Coupler verification of Thevenin Source and microphone calibrations.'
%         '  Place ER10X probe in the brass coupler.'
%         '  Ensure that the Reference microphone is also attached and turned on.'
%         };
%     nn = size(alertTxt,1);
%     for ii=1:nn
%         cprintf([0,0,.4],[alertTxt{ii},'\n']);
%     end
%     disp(' ')    
%     showInstructions
%     params.applyMicCorrection = 1;
%     params.createMicCorrection = 0;
%     ARLas_couplerRecordings_DW10x(obj,params);
% end
% % 4. Double check system distortion and noise floor. ----------------------
%     % NOT ACTUALLY IMPLEMENTED HERE
%     alertTxt = {'..................................'
%          'Finished Calibration Routine'
%         };
%     nn = size(alertTxt,1);
%     disp(' ')
%     for ii=1:nn
%         cprintf([0,0,.4],[alertTxt{ii},'\n']);
%     end
%     disp(' '); disp(' ');
%     

