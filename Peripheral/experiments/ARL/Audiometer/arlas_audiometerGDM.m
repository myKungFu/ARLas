classdef arlas_audiometerGDM < handle
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Class definition arlas_audiometerGDM
%
% This program is an audiometer for use with ARLas 
% (Auditory Research Laboratory auditory software)
%
% Required functions:
%    myLogisticReg.m
% NOTE: Objects of this class only work with ARLas version '2025.01.07' or later.
%
% The Audiometer may be called in three different modes:
% RUN MODE 1) To call this program in standard mode, use an experiment file 
%     to call it. Run the experiment file in ARLas. The Experiment file should 
%     contain the following lines of code: 
%       obj = varargin{1}; % get the arlas object
%       probeL = 'A'; % probe that is in LEFT ear; 'A', 'B', or []
%       probeR = 'B'; % probe that is in RIGHT ear; 'A', 'B', or []
%       subjAge = 32; % change this to be the subjects age in years
%       subjID = 'Bill'; % chang this to be the subject's ID (as a string)
%       objA = arlas_audiometerGDM(obj,probeL,probeR,subjAge,subjID);
% RUN MODE 2) To use this program to examine and plot previously collected data, 
%     find the previously saved data and load the saved copy of the object.
%     This will be saved as a variable in a file called XXXXXX_analyzedAUDIO.mat. 
%     The variable has the name objCopy. Call the audiometer using objCopy
%     as a single input argument:
%       arlas_audiometerGDM(audiogramData);
% RUN MODE 3) To call this program in simulation mode (no actual stimuli presented
%     so you don't have to have a sound card, no calibration needed),
%     simply call the program by typing its name at the command prompt with
%     no input arguments: 
%       arlas_audiometerGDM;
%     Data will be saved Here:
%     'C:\myWork\ARLas\Data\HumptyD\todaysDate\AudiophileAudiophile_analysis'
%     Data will be saved in the files named
%       HumptyD_analyzedAUDIO.mat
%       HumptyD_analyzedAUDIO.xlsx
%
% Auditory Research Lab, The University of Iowa
% Deptartment of Communication Sciences & Disorders
% The University of Iowa
% Author: Shawn S. Goodman, PhD
% Date: July 26, 2021
% Last Updated: July 26, 2021
% Last Updated: July 30, 2021 -- ssg
% Last Updated: January 6-20, 2025 -- ssg
% Last Updated: May 8, 2025 -- mek -- added keyboard bindings
%
% To do
% reset buttons when crashes obj.buttonManager(10) to reset buttons
% make run with only left or right channels
% save output voltages
% measure voltages in dB SPL in a coupler
% get real calibration
% figure out response button (auditory and/or visual?)
% do semi-automatic mode?
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
properties (SetAccess = private)
    version = '08MAY2025';
    sep % path delimiter appriate for the current operating system 
    map % struct containing file paths
    arlasObj % object that is the arlas object (not the audiometer object)
    fs % sampling rate
    runMode % running to collect data (1), examining saved data (2), test mode (no actual sound presented) (3)
    
    %----- handles fror graphical user interface ----- -----    
    H % Figure handle ----- 
        guiSize
        
        CONTROL % CONTROL Panel -----
            h % control panel gui
            gui % dimensions of the control panel
            nButtons = 5; % number of buttons
            h_splash
            h_freqDn
            h_freqUp
            h_present
            h_lvlUp
            h_lvlDn
            h_cal
            h_save
            h_respYES
            h_respNO
            
            TEST_EARtxt
            TEST_EAR
            FREQ
            LVL
            %AUTO
            %AUTOtxt
            RESPtxt
            STEP
            STEPtxt
            REF
            REFtxt
            
            h_autoRun
            h_earNow
            h_freqNow 
            h_lvlNow
            h_thdNow
            
        VIEW % VIEW pannel (performance intensity function and tracking -----
            % The following belong in this pane, but they are public properties
            %   h_pif % performance intensity function
            %   h_track % response tracker
            %   h_ag % audiogram
            h_switchPlots
            marker1 % marker (triangle) for pi function
            marker2 % marker (triangle) for progress tracker
            L1 % lines for audiogram
            L2 % lines for audiogram
            markerX % marker (triangle) for audiogram
            markerY % marker (triangle) for audiogram
            LX % thicker line for audiogram
            LY % thicker line for audiogram
            LfHL % line showing dB fHL
            LfHLa % line showing dB fHL for patient's age
            thdMarks % audiometric threshold markers (x or o)
            thdDots % dotted lines connecting thresholds
            AGjunk % handle for audiogram title
            AGjunk2 % handle for audiogram confidence intervals
            PIjunk % handles to plotted aspects of the pi function
            PIjunk2 % handles to text plotted on the pi function
            PRjunk % handles to plotted aspects of the progress tracker
            PRjunk2
            doingSwitch % switching screen views or not
end
properties (SetAccess = public) 
    % visualization/plotting ----------------------------------------------
    h_pif % performance intensity function
    h_track % response tracker
    h_ag % audiogram

    % ear info ------------------------------------------------------------
    ear  % current test ear (L = 1, R = 2)
    
    % frequency info ------------------------------------------------------
    freqs     % row vector of test frequencies
    nFreqs    % total number of test frequencies
    fIndx     % index of current frequency (2 x nFreqs vector; L and R ears)

    % level info ----------------------------------------------------------
    lvlRefs    % level reference: SPL, FPL, IPL
    stepSizes % level step size (5 or 2 dB)
    minOutput % row vector of minimum output levels (1 x nFreqs vector)
    maxOutput % row vector of maximum output levels (1 x nFreqs vector)
    levels    % vector of current possible levels (for active frequency)
    nLevels   % total number of possible levels
    lvlIndx   % index of current level in the levels vector
    
    % threshold values ----------------------------------------------------
    %thd % row vector of minimum output levels (2 x nFreqs vector)
    
    autorunNow % automatically run audiometer
    autorespNow % automaticly get response
    earNow
    freqNow 
    lvlNow
    thdNow
    stepNow
    refNow
    iscNow
    viewNow

    % stimuli ------------------------------------------
    Stimuli
    stimNow
    foldN
    NFolds
    nReps
    zpad
    
    DONE % for total auto sequence
    doneL 
    doneR
    
    nReversals
    descend1
    
    % calibration info
    iscS1L
    iscS2L
    iscS1R
    iscS2R
    probeOutputL
    probeOutputR
    probeInputL
    probeInputR
    clickerInput
    micCorrectionL
    micCorrectionR
    C1L
    C2L
    C1R
    C2R
    probeR
    probeL
    doMicCorrection
    fmin
    fmax
    whichLoudspeaker % most probes have 2 loudspeakers--which one are you using (1 or 2)
    dBConversionL % conversion factor to make FPL into IPL (in dB)
    dBConversionR
    calFinished

    AAA % audiometer results structure. By ear and frequency
    minPresentations % number of stim presentations needed to get three "reversals" in standard Hewson-Westlake paradigm
    maxPresentations
    minReversals % minumim number of reversals to say threshold achieved

    subjAge % age of person being tested (years)
    subjID % id of person being tested
    PRIORS % structure containing priors information for current subject
    priorsPathName % file containing priors data (for all ages)
    priorsFileName % location of priors
    pdPrior % probability distribution of priors
    Q % structure used when finding pi function
    epsilon
    fHL % expected mean thresholds for young, normal hearing 20 year old
    fHLa % expected mean thresholds for the decade of life of the subject
    startingLvls_L % starting test levels for each frequency, given age
    startingLvls_R
    audiogramData
end
methods
    
    function keyPressHandler(obj, event)
        switch event.Key
            case 'uparrow'
                obj.lvlUpManager(obj.h_lvlUp, []);
            case 'downarrow'
                obj.lvlDnManager(obj.h_lvlDn, []);
            case 'rightarrow'
                obj.freqUpManager(obj.h_freqUp, []);
            case 'leftarrow'
                obj.freqDnManager(obj.h_freqDn, []);
            case 'space'
                obj.presentManager(obj.h_present, []);
            case 'y'
                obj.doYes(obj.h_respYES, []);
            case 'n'
                obj.doNo(obj.h_respNO, []);
            % Add more keys if needed here Shawn
        end
    end

    function obj = arlas_audiometerGDM(subjID,subjAge,arlasObj,probeL,probeR,freqs) % initialize object of class initARLas
        if nargin >= 5 % if being called from ARLas (to collect real data)
            obj.runMode = 1;
            obj.sep = filesep; % get the path delimiter appropriate to the operating system being used            
            obj.arlasObj = arlasObj;
            obj.fs = arlasObj.fs; % get the system sampling rate
            obj.probeL = probeL;
            obj.probeR = probeR;
            if nargin == 6
                freqs = sort(freqs(:));
            else % frequencies not specified
                freqs = [0.25 0.5 1 2 4 8 10 11.2 12.5 14 16]' * 1000; % use default frequency set
            end
            obj.whichLoudspeaker = 1;
            obj.calFinished = 0;
        elseif nargin == 1 % is being called to examine a saved session
            audiogramData = subjID;
            obj.runMode = 2;
            obj.sep = filesep; % get the path delimiter appropriate to the operating system being used
            obj.probeL = audiogramData.probeL; %objCopy.probeL;
            obj.probeR = audiogramData.probeR; %objCopy.probeR;
            freqs = audiogramData.freqs;
            obj.AAA = audiogramData.AAA;
            subjAge = audiogramData.subjAge;
            subjID = audiogramData.subjID;
            obj.dBConversionL = audiogramData.dBConversionL;
            obj.dBConversionR = audiogramData.dBConversionR;
            obj.calFinished = 1;
            obj.whichLoudspeaker = audiogramData.whichLoudspeaker;
        elseif nargin == 0 % being called to run in development or practice mode
            obj.runMode = 3;
            obj.sep = filesep; % get the path delimiter appropriate to the operating system being used            
            obj.fs = 96000; % set the system sampling rate
            obj.probeL = 'A';
            obj.probeR = 'B';
            freqs = [0.25 0.5 1 2 4 8 10 11.2 12.5 14 16]' * 1000; % default frequency set
            subjID = 'TestCase';
            subjAge = 58;
            obj.whichLoudspeaker = 1;
            obj.calFinished = 1;
            obj.dBConversionL = ones(size(freqs))*6;
            obj.dBConversionR = ones(size(freqs))*6;
       end
        obj.freqs = freqs(:)';
        obj.nFreqs = length(obj.freqs);
        [~,obj.fIndx] = min(abs(freqs - 1000)); % initialize to frequency closest to 1 kHz
        obj.ear = 1; % initialize to left ear

        obj.lvlRefs = [{'FPL'},{'IPL'}]; % initialize to fpl (forward pressure level)
        obj.stepSizes = [5,2]; % level step size options
        obj.minOutput = zeros(1,obj.nFreqs) -10;  % initialize to minus 10 dB
        obj.maxOutput = ones(1,obj.nFreqs) * 110;  % initialize to 100 dB

        % ----- default client values:
        obj.priorsPathName = 'C:\myWork\ARLas\Peripheral\experiments\ARL\Audiometer\';
        obj.priorsFileName = 'HarpAudioFPL.mat';
        % get dB fHL (HL, but in fpl), defined as mean thresholds for ages 15-25
        obj.setAge(20); % set age for purposes of getting HL (for fpl)
        obj.fHL = obj.PRIORS.mu; % dB HL (but for fpl)
        obj.setAge(subjAge); % now set to subjects age
        obj.fHLa = obj.PRIORS.mu; % dB HL (but for fpl--age adjusted)
        obj.setID(subjID); % set subject ID

        obj.stepNow = obj.stepSizes(1);
        obj.levels = (obj.minOutput(obj.fIndx):obj.stepNow:obj.maxOutput(obj.fIndx))';
        obj.nLevels = length(obj.levels);
        [~,obj.lvlIndx] = min(abs(obj.levels - 50)); % initialize to value to 50. Will change when enter age
        
        obj.autorunNow = 0; % initualize to manual control
        obj.autorespNow = 0; % initialize to manual response
        obj.earNow = 1;
        obj.freqNow = obj.freqs(obj.fIndx); % initialize to closest to 1 kHz

        obj.lvlNow = obj.startingLvls_L(obj.fIndx); %obj.levels(obj.lvlIndx);
        [~,obj.lvlIndx] = min(abs(obj.levels - obj.lvlNow)); % initialize to value to 50. Will change when enter age

        obj.refNow = obj.lvlRefs(1);
        obj.thdNow = [];
        
        obj.DONE = 0; % for total auto sequence
        obj.doneL = zeros(obj.nFreqs,1); % for total auto sequence
        obj.doneR = zeros(obj.nFreqs,1); % for total auto sequence
        if obj.runMode == 1
            obj.getStim;
        end

        % AAA the audiometer results structure ---------------------------------------
        if obj.runMode ~= 2
            obj.AAA.freqs = obj.freqs;
            obj.AAA.nFreqs = obj.nFreqs;
            obj.AAA.levels = obj.levels;
            obj.AAA.DONE = 0;
            obj.doneL = zeros(obj.nFreqs,1);
            obj.doneR = zeros(obj.nFreqs,1);
            ear = 'L';
            for ii=1:obj.nFreqs
                obj.AAA.(ear).(['f_',num2str(obj.freqs(ii))]).X = []; % corresponding stimulus level vector
                obj.AAA.(ear).(['f_',num2str(obj.freqs(ii))]).Y = []; % subject response vector
                obj.AAA.(ear).(['f_',num2str(obj.freqs(ii))]).trialCounter = 0; % total number of trials presented in the track
                obj.AAA.(ear).(['f_',num2str(obj.freqs(ii))]).thd = []; % threshold (50% point on logistic function)
                obj.AAA.(ear).(['f_',num2str(obj.freqs(ii))]).sigma = []; % spread of normal CDF approximation of logistic function
                obj.AAA.(ear).(['f_',num2str(obj.freqs(ii))]).CDF = []; % normal CDF (mu=thd, sigma) at X = B = stimulus levels
                obj.AAA.(ear).(['f_',num2str(obj.freqs(ii))]).B = []; % stimulus levels over which cdf is evaluated
                obj.AAA.(ear).(['f_',num2str(obj.freqs(ii))]).badIndx = []; % indices of identified false positives/negatives
                obj.AAA.(ear).(['f_',num2str(obj.freqs(ii))]).nRejects = []; % number of identified false positives/negatives
                obj.AAA.(ear).(['f_',num2str(obj.freqs(ii))]).prior = []; % prior probability of threshold
            end
            ear = 'R';
            for ii=1:obj.nFreqs
                obj.AAA.(ear).(['f_',num2str(obj.freqs(ii))]).X = []; % corresponding stimulus level vector
                obj.AAA.(ear).(['f_',num2str(obj.freqs(ii))]).Y = []; % subject response vector
                obj.AAA.(ear).(['f_',num2str(obj.freqs(ii))]).trialCounter = 0; % total number of trials presented in the track
                obj.AAA.(ear).(['f_',num2str(obj.freqs(ii))]).thd = []; % threshold (50% point on logistic function)
                obj.AAA.(ear).(['f_',num2str(obj.freqs(ii))]).sigma = []; % spread of normal CDF approximation of logistic function
                obj.AAA.(ear).(['f_',num2str(obj.freqs(ii))]).CDF = []; % normal CDF (mu=thd, sigma) at X = B = stimulus levels
                obj.AAA.(ear).(['f_',num2str(obj.freqs(ii))]).B = []; % stimulus levels over which cdf is evaluated
                obj.AAA.(ear).(['f_',num2str(obj.freqs(ii))]).badIndx = []; % indices of identified false positives/negatives
                obj.AAA.(ear).(['f_',num2str(obj.freqs(ii))]).nRejects = []; % number of identified false positives/negatives
                obj.AAA.(ear).(['f_',num2str(obj.freqs(ii))]).prior = []; % prior probability of threshold
            end
        end
        % -----------------------------------------------------------------
        
        obj.minPresentations = 6; % minimum number of presentations to get threshold
        obj.maxPresentations = 100; % maximum number of presentations to get
        obj.minReversals = 3; % minimum reversals to get threshold
        obj.doingSwitch = 1;
        obj.viewNow = 0; % initialize to pi / progress view (0)
        obj.initGui % initialize the gui
        obj.audiometer_plotPI
        if nargin == 0 | nargin == 1 % simulation or previously collected data
            obj.buttonManager(20) % do not need to run cal first
        end
        if obj.runMode == 2 % examine a saved session, so put in audiogram view
            obj.doSwitch
        end
        if obj.runMode == 1 % if in the data collection mode, set up data to save
            obj.audiogramData.freqs = obj.freqs;     % row vector of test frequencies
            obj.audiogramData.AAA = obj.AAA; % audiometer results structure. By ear and frequency
            obj.audiogramData.subjAge = obj.subjAge; % age of person being tested (years)
            obj.audiogramData.subjID = obj.subjID; % id of person being tested
            obj.audiogramData.PRIORS = obj.PRIORS; % structure containing priors information for current subject
            obj.audiogramData.pdPrior = obj.pdPrior; % probability distribution of priors
        end
    end  
    function abort(varargin) % instructions for aborting when gui closed
        % try
        %     obj = varargin{1};
        % catch
        % end
        try
            obj = varargin{1};
            delete(obj.H)
            delete(obj);
        catch
            while gcf~=1 
                delete(gcf);
            end
            delete(gcf);
        end
    end
    function initGui(varargin) % initialize the arlas graphical user interface
        obj = varargin{1};
        try % delete figure if it already exists
            delete(obj.H)
        catch
        end
        try % create new main figure
            %[left, bottom,width, height]
            obj.guiSize.width = 545; % figure width in pixels
            obj.guiSize.height = 288; %620
            scrsz = get(0,'ScreenSize'); % get the current screen size
            obj.guiSize.left = round(scrsz(4) * .1); % location of left edge
            obj.guiSize.bottom = scrsz(4) * .1; % location of bottom
            overhang = (obj.guiSize.bottom + obj.guiSize.height)-scrsz(4); % adjust so gui doesn't hang off the top of the screen
            if overhang > 0 % if positive value (meaning figure is off the screen)
                overhang = overhang + 0; % give a little extra to account for top of the figure
                %if overhang <= obj.guiSize.bottom % try to fix this problem
                    obj.guiSize.bottom = obj.guiSize.bottom - overhang; % correct for the overhang
                %else % otherwise, do the best you can
                %    obj.guiSize.bottom = 1; % put it as low as possible
                %end
            end
            rect = [obj.guiSize.left, obj.guiSize.bottom, obj.guiSize.width, obj.guiSize.height];
            obj.H = figure('Position',rect,'Color',[1 1 1],'Units','Pixels',...
                'CloseRequestFcn',@obj.abort,'Name',['ARLas audiometer version ',obj.version],...
                'NumberTitle','off','MenuBar','none','Resize','on','Color',[1 1 1],'Tag','ARLas');
            set(obj.H, 'KeyPressFcn', @(src, event) obj.keyPressHandler(event)); %% Key press function
       obj.H.Position = [887.0000   67.4000  544.8000  744.4000];
            % create panels within main gui figure -----
            obj.CONTROL = uipanel('Parent',obj.H,'Title','CONTROL PANEL','FontSize',12,...
                'BackgroundColor','white','Units','Pixels','Position',[10 10 105*5+4 250]); % 140
            obj.VIEW = uipanel('Parent',obj.H,'Title','VIEW','FontSize',12,...
                'BackgroundColor','white','Units','Pixels','Position',[10 290 105*5+4 455]);     
            % Populate control panel -----
            obj.gui.height = 105;
            obj.gui.width = 105;
            obj.h_splash = uicontrol('Parent',obj.CONTROL,'Style','togglebutton',...
                'BackgroundColor',[1 1 1],'Position',[(obj.gui.width*0) 1 obj.gui.width*5 obj.gui.height],...
                'Visible','on','CData',imread('splash.jpg'),'BusyAction','queue','Interruptible','off');
            pause(1.5)
            delete(obj.h_splash)
            obj.h_freqDn = uicontrol('Parent',obj.CONTROL,'Style','togglebutton',...
                'BackgroundColor',[1 1 1],'Position',[(obj.gui.width*0) 1 obj.gui.width obj.gui.height],...
                'Callback',@obj.freqDnManager,'Visible','on','TooltipString','Increase Frequency',...
                'CData',imread('freqDnGray.jpg'),'BusyAction','queue','Interruptible','on',...
                'KeyPressFcn',@obj.keyPress);
            
            obj.h_freqUp = uicontrol('Parent',obj.CONTROL,'Style','togglebutton',...
                'BackgroundColor',[1 1 1],'Position',[(obj.gui.width*1) 1 obj.gui.width obj.gui.height],...
                'Callback',@obj.freqUpManager,'Visible','on','TooltipString','Increase Frequency',...
                'CData',imread('freqUpGray.jpg'),'BusyAction','queue','Interruptible','on');
            obj.h_present = uicontrol('Parent',obj.CONTROL,'Style','togglebutton',...
                'BackgroundColor',[1 1 1],'Position',[(obj.gui.width*2) 1 obj.gui.width obj.gui.height],...
                'Callback',@obj.presentManager,'Visible','on','TooltipString','Present Stimulus',...
                'CData',imread('presentGray.jpg'),'BusyAction','queue','Interruptible','on');
            obj.h_lvlDn = uicontrol('Parent',obj.CONTROL,'Style','togglebutton',...
                'BackgroundColor',[1 1 1],'Position',[(obj.gui.width*3) 1 obj.gui.width obj.gui.height],...
                'Callback',@obj.lvlDnManager,'Visible','on','TooltipString','Decrease Stimulus Level',...
                'CData',imread('lvlDnGray.jpg'),'BusyAction','queue','Interruptible','on');
            obj.h_lvlUp = uicontrol('Parent',obj.CONTROL,'Style','togglebutton',...
                'BackgroundColor',[1 1 1],'Position',[(obj.gui.width*4) 1 obj.gui.width obj.gui.height],...
                'Callback',@obj.lvlUpManager,'Visible','on','TooltipString','Increase Stimulus Level',...
                'CData',imread('lvlUpGray.jpg'),'BusyAction','queue','Interruptible','on');
            % -----            
             % [left,bottom,width,height]       
            % create controls -----
            % view input channel
            obj.FREQ = uicontrol('Parent',obj.CONTROL,'Style','text',...
                'String',[num2str(obj.freqNow),' Hz'],'FontSize',16,'BackgroundColor',...
                'white','HorizontalAlignment','Left','Units','Normalized',...
                'Position',[.13 .55 .3 .1]);            
            obj.LVL = uicontrol('Parent',obj.CONTROL,'Style','text',...
                'String',[num2str(obj.lvlNow),' dB ',char(obj.refNow)],'FontSize',16,'BackgroundColor',...
                'white','HorizontalAlignment','Left','Units','Normalized',...
                'Position',[.72 .55 .3 .1]);

            obj.REFtxt = uicontrol('Parent',obj.CONTROL,'Style','text',...
                'BackgroundColor','white','String','Lvl Reference:',...
                'HorizontalAlignment','Left','Units','Normalized',...
                'Position',[.12 .85 .3 .1],'FontSize',12);
            obj.REF = uicontrol('Parent',obj.CONTROL,'Style','popup',...
                'String',{'FPL','IPL'},'FontSize',12,'BackgroundColor',...
                'white','HorizontalAlignment','Left','Units','Normalized',...
                'Position',[.12 .75 .16 .1],'Value',1,...
                'Callback',@obj.refManager);
            obj.TEST_EARtxt = uicontrol('Parent',obj.CONTROL,'Style','text',...
                'BackgroundColor','white','String','Test Ear:',...
                'HorizontalAlignment','Left','Units','Normalized',...
                'Position',[.45 .85 .3 .1],'FontSize',12);
            obj.TEST_EAR = uicontrol('Parent',obj.CONTROL,'Style','popup',...
                'String',{'Left','Right'},'FontSize',12,'BackgroundColor',...
                'white','HorizontalAlignment','Left','Units','Normalized',...
                'Position',[.42 .75 .16 .1],'Value',obj.earNow,...
                'Callback',@obj.earManager);
            obj.STEPtxt = uicontrol('Parent',obj.CONTROL,'Style','text',...
                'BackgroundColor','white','String','Step Size:',...
                'HorizontalAlignment','Left','Units','Normalized',...
                'Position',[.73 .85 .3 .1],'FontSize',12);
            obj.STEP = uicontrol('Parent',obj.CONTROL,'Style','popup',...
                'String',{'5','2'},'FontSize',12,'BackgroundColor',...
                'white','HorizontalAlignment','Left','Units','Normalized',...
                'Position',[.72 .75 .16 .1],'Value',1,...
                'Callback',@obj.stepManager);

            % obj.AUTOtxt = uicontrol('Parent',obj.CONTROL,'Style','text',...
            %     'BackgroundColor','white','String','Control:',...
            %     'HorizontalAlignment','Left','Units','Normalized',...
            %     'Position',[.45 .59 .16 .1],'FontSize',12);
            obj.RESPtxt = uicontrol('Parent',obj.CONTROL,'Style','text',...
                'BackgroundColor','white','String','Subj Response:',...
                'HorizontalAlignment','Left','Units','Normalized',...
                'Position',[.395 .59 .22 .1],'FontSize',12);
            % obj.AUTO = uicontrol('Parent',obj.CONTROL,'Style','popup',...
            %     'String',[{'Manual','Auto'}],'FontSize',12,'BackgroundColor',...
            %     'white','HorizontalAlignment','Left','Units','Normalized',...
            %     'Position',[.42 .5 .16 .1],'Value',1,...
            %     'Callback',@obj.autoManager);
            obj.h_respYES = uicontrol('Parent',obj.CONTROL,'Style','pushbutton',...
            'BackgroundColor',[.7 .7 .7],'Units','Normalized',...
                  'Position',[.4 .47 .09 .13],'Visible','on','Callback',@obj.doYes,...
                  'TooltipString','Did the subject hear the tone?',...
                  'String','Yes','BusyAction','queue','Interruptible','on');    
            obj.h_respNO = uicontrol('Parent',obj.CONTROL,'Style','pushbutton',...
            'BackgroundColor',[.7 .7 .7],'Units','Normalized',...
                  'Position',[.505 .47 .09 .13],'Visible','on','Callback',@obj.doNo,...
                  'TooltipString','Did the subject hear the tone?',...
                  'String','No','BusyAction','queue','Interruptible','on');    

            obj.h_cal = uicontrol('Parent',obj.CONTROL,'Style','togglebutton',...
                'BackgroundColor',[.7 .7 .7],'Units','Normalized',...
                'Position',[.01     .9   .1   .1],'Callback',@obj.doCal,'Visible','on',...
                'TooltipString','Run Calibration Routine',...
                'String','Cal','BusyAction','queue','Interruptible','on');
            obj.h_save = uicontrol('Parent',obj.CONTROL,'Style','togglebutton',...
                'BackgroundColor',[.7 .7 .7],'Units','Normalized',...
                'Position',[.895     .9   .1   .1],'Callback',@obj.doSave,'Visible','on',...
                'TooltipString','Save Data',...
                'String','Save','BusyAction','queue','Interruptible','on');
            obj.h_switchPlots = uicontrol('Parent',obj.VIEW,'Style','togglebutton',...
                'BackgroundColor',[.7 .7 .7],'Units','Normalized',...
                'Position',[.855     .01   .14   .05],'Callback',@obj.doSwitch,'Visible','on',...
                'TooltipString','Switch from PI/Progress view to Audiogram',...
                'String','Switch View','BusyAction','queue','Interruptible','on','Value',1);
            % -----            
            obj.buttonManager(10)

            obj.h_ag = axes('Parent',obj.VIEW,'Visible','on','Color',[1 1 1],...
                'Units','Normalized','Position',[.125 .15 .85 .8],'XTick',[]);
            axes(obj.h_ag)
            xlabel('Frequency (kHz)','FontSize',12);
            ylabel('Stim Level (dB)','FontSize',12);
            obj.h_ag.Visible = 'off';
            
            obj.h_pif = axes('Parent',obj.VIEW,'Visible','on','Color',[1 1 1],...
                'Units','Normalized','Position',[.125 .6 .85 .35],'XTick',[]);
            axes(obj.h_pif)
            xlabel('Level (dB FPL)','FontSize',12);
            ylabel('Response Probability','FontSize',12);
            obj.h_pif.Visible = 'on';
            
            obj.h_track = axes('Parent',obj.VIEW,'Visible','on','Color',[1 1 1],...
                'Units','Normalized','Position',[.125 .145 .85 .35]); % [.125 .125 .85 .35]
            axes(obj.h_track)
            xlabel('Trial Number','FontSize',12);
            ylabel('Stim Lvl (dB FPL)','FontSize',12);  
            obj.h_track.Visible = 'on';
            
        catch ME
            errorTxt = {'  Issue: Error creating audiometer GUI.'
                 '  Action: None.'
                 '  Location: in arlas_audiometer.initGui.'
                };
            errorMsgARLas(errorTxt);
            obj.printError(ME)
        end
    end    
    function buttonManager(varargin) % change colors and states of buttons on bottom bar
        obj = varargin{1};
        try
            code = varargin{2};
            flickerLen = 0.1; % length of flicker for error notification
            flickerReps = 5; % number of flickers
            % 10s are for arlas_audiometer gui setup
            % 20s are for initializing the system
            % 30s are for pausing 
            % 40s are for loading experiment files
            % 50s are for running
            % 60s are for aborting
            %
            % the 'Tag' field is used to indicate whether a button can be used
            % or not. Will be set to zero when its callback is currently
            % running, as well.
            switch code
                case 10 % initialize arlas_audiometer gui (FUNCTION initGui)
                    set(obj.h_cal,'Value',0,'Tag','on','BackgroundColor',[1 1 0])
                    set(obj.h_save,'Value',0,'Tag','off','BackgroundColor',[.7 .7 .7])
                    set(obj.h_switchPlots,'Value',0,'Tag','off','BackgroundColor',[.7 .7 .7])
                    set(obj.h_respYES,'Value',0,'Tag','off','BackgroundColor',[.7 .7 .7])
                    set(obj.h_respNO,'Value',0,'Tag','off','BackgroundColor',[.7 .7 .7])
                    set(obj.h_freqUp,'Value',0,'Tag','off','CData',imread('freqUpGray.jpg'))
                    set(obj.h_freqDn,'Value',0,'Tag','off','CData',imread('freqDnGray.jpg'))
                    set(obj.h_present,'Value',0,'Tag','off','CData',imread('presentGray.jpg'))
                    set(obj.h_lvlUp,'Value',0,'Tag','off','CData',imread('lvlUpGray.jpg'))
                    set(obj.h_lvlDn,'Value',0,'Tag','off','CData',imread('lvlDnGray.jpg'))
                case 20 % after calibrating
                    set(obj.h_cal,'Value',0,'Tag','on','BackgroundColor',[1 1 0])
                    set(obj.h_save,'Value',0,'Tag','on','BackgroundColor',[1 1 0])
                    set(obj.h_switchPlots,'Value',0,'Tag','on','BackgroundColor',[1 1 0])
                    set(obj.h_respYES,'Value',0,'Tag','off','BackgroundColor',[.7 .7 .7])
                    set(obj.h_respNO,'Value',0,'Tag','off','BackgroundColor',[.7 .7 .7])
                    set(obj.h_freqUp,'Value',0,'Tag','on','CData',imread('freqUpYellow.jpg'))
                    set(obj.h_freqDn,'Value',0,'Tag','on','CData',imread('freqDnYellow.jpg'))
                    set(obj.h_present,'Value',0,'Tag','on','CData',imread('presentYellow.jpg'))
                    set(obj.h_lvlUp,'Value',0,'Tag','on','CData',imread('lvlUpYellow.jpg'))
                    set(obj.h_lvlDn,'Value',0,'Tag','on','CData',imread('lvlDnYellow.jpg'))
                case 21 % -- unsuccessful initialization
                    for ii=1:flickerReps
                        set(obj.h_init,'Value',1,'Tag','off','CData',imread('initRed.jpg'))
                        pause(flickerLen)
                        set(obj.h_init,'Value',1,'Tag','off','CData',imread('initGreen.jpg'))
                        pause(flickerLen)
                    end
                    set(obj.h_init,'Value',0,'Tag','on','CData',imread('initGreen.jpg'))
                    set(obj.h_load,'Value',0,'Tag','off','CData',imread('loadGray.jpg'))
                    set(obj.h_pause,'Value',0,'Tag','off','CData',imread('pauseGray.jpg'))
                    set(obj.h_run,'Value',0,'Tag','off','CData',imread('runGray.jpg'))
                    set(obj.h_abort,'Value',0,'Tag','off','CData',imread('abortGray.jpg'))
                    set(obj.h_respYES,'Value',0,'Tag','off','BackgroundColor',[.7 .7 .7])
                    set(obj.h_respNO,'Value',0,'Tag','off','BackgroundColor',[.7 .7 .7])
                case 22 % -- successful initialization
                    set(obj.h_init,'Value',0,'Tag','on','CData',imread('initYellow.jpg'))
                    set(obj.h_load,'Value',0,'Tag','on','CData',imread('loadYellow.jpg'))
                    set(obj.h_pause,'Value',0,'Tag','off','CData',imread('pauseGray.jpg'))
                    set(obj.h_run,'Value',0,'Tag','off','CData',imread('runGray.jpg'))
                    set(obj.h_abort,'Value',0,'Tag','off','CData',imread('abortGray.jpg'))
                    set(obj.h_respYES,'Value',0,'Tag','off','BackgroundColor',[.7 .7 .7])
                    set(obj.h_respNO,'Value',0,'Tag','off','BackgroundColor',[.7 .7 .7])
                case 30 % Begin pausing the system (FUNCTION pauseExperiment)
                    set(obj.h_init,'Value',0,'Tag','off','CData',imread('initGray.jpg'))
                    set(obj.h_load,'Value',0,'Tag','off','CData',imread('loadGray.jpg'))
                    set(obj.h_pause,'Value',1,'Tag','off','CData',imread('pauseGreen.jpg'))
                    set(obj.h_run,'Value',0,'Tag','off','CData',imread('runGray.jpg'))
                    set(obj.h_abort,'Value',0,'Tag','off','CData',imread('abortGray.jpg'))
                    set(obj.h_respYES,'Value',0,'Tag','off','BackgroundColor',[.7 .7 .7])
                    set(obj.h_respNO,'Value',0,'Tag','off','BackgroundColor',[.7 .7 .7])
                case 31 % -- unsuccessful pause
                    for ii=1:flickerReps
                        set(obj.h_pause,'Value',1,'Tag','off','CData',imread('pauseRed.jpg'))
                        pause(flickerLen)
                        set(obj.h_pause,'Value',1,'Tag','off','CData',imread('pauseGreen.jpg'))
                        pause(flickerLen)
                    end
                    set(obj.h_init,'Value',0,'Tag','on','CData',imread('initGray.jpg'))
                    set(obj.h_load,'Value',0,'Tag','off','CData',imread('loadGray.jpg'))
                    set(obj.h_pause,'Value',0,'Tag','on','CData',imread('pauseYellow.jpg'))
                    set(obj.h_run,'Value',0,'Tag','on','CData',imread('runGreen.jpg'))
                    set(obj.h_abort,'Value',0,'Tag','on','CData',imread('abortRed.jpg'))
                    set(obj.h_respYES,'Value',0,'Tag','off','BackgroundColor',[.7 .7 .7])
                    set(obj.h_respNO,'Value',0,'Tag','off','BackgroundColor',[.7 .7 .7])
                case 32 % -- successful pause ON
                    set(obj.h_init,'Value',0,'Tag','off','CData',imread('initGray.jpg'))
                    set(obj.h_load,'Value',0,'Tag','off','CData',imread('loadGray.jpg'))
                    set(obj.h_pause,'Value',0,'Tag','on','CData',imread('pauseGreen.jpg'))
                    set(obj.h_run,'Value',0,'Tag','off','CData',imread('runGray.jpg'))
                    set(obj.h_abort,'Value',0,'Tag','on','CData',imread('abortRed.jpg'))
                    set(obj.h_respYES,'Value',0,'Tag','off','BackgroundColor',[.7 .7 .7])
                    set(obj.h_respNO,'Value',0,'Tag','off','BackgroundColor',[.7 .7 .7])
                case 33 % -- successful pause OFF
                    set(obj.h_init,'Value',0,'Tag','off','CData',imread('initGray.jpg'))
                    set(obj.h_load,'Value',0,'Tag','off','CData',imread('loadGray.jpg'))
                    set(obj.h_pause,'Value',0,'Tag','on','CData',imread('pauseYellow.jpg'))
                    set(obj.h_run,'Value',0,'Tag','on','CData',imread('runGreen.jpg'))
                    set(obj.h_abort,'Value',0,'Tag','on','CData',imread('abortRed.jpg'))
                    set(obj.h_respYES,'Value',0,'Tag','off','BackgroundColor',[.7 .7 .7])
                    set(obj.h_respNO,'Value',0,'Tag','off','BackgroundColor',[.7 .7 .7])
                case 40 % Begin loading experiment file (FUNCTION loadExperiment)
                    set(obj.h_init,'Value',0,'Tag','off','CData',imread('initGray.jpg'))
                    set(obj.h_load,'Value',1,'Tag','off','CData',imread('loadGreen.jpg'))
                    set(obj.h_pause,'Value',0,'Tag','off','CData',imread('pauseGray.jpg'))
                    set(obj.h_run,'Value',0,'Tag','off','CData',imread('runGray.jpg'))
                    set(obj.h_abort,'Value',0,'Tag','off','CData',imread('abortGray.jpg'))
                    set(obj.h_respYES,'Value',0,'Tag','off','BackgroundColor',[.7 .7 .7])
                    set(obj.h_respNO,'Value',0,'Tag','off','BackgroundColor',[.7 .7 .7])
                case 41 % -- failed to load experiment file
                    for ii=1:flickerReps
                        set(obj.h_load,'Value',1,'Tag','off','CData',imread('loadRed.jpg'))
                        pause(flickerLen)
                        set(obj.h_load,'Value',1,'Tag','off','CData',imread('loadGreen.jpg'))
                        pause(flickerLen)
                    end
                    set(obj.h_init,'Value',0,'Tag','on','CData',imread('initYellow.jpg'))
                    set(obj.h_load,'Value',0,'Tag','on','CData',imread('loadYellow.jpg'))
                    set(obj.h_pause,'Value',0,'Tag','off','CData',imread('pauseGray.jpg'))
                    set(obj.h_run,'Value',0,'Tag','off','CData',imread('runGray.jpg'))
                    set(obj.h_abort,'Value',0,'Tag','off','CData',imread('abortGray.jpg'))
                    set(obj.h_respYES,'Value',0,'Tag','off','BackgroundColor',[.7 .7 .7])
                    set(obj.h_respNO,'Value',0,'Tag','off','BackgroundColor',[.7 .7 .7])
                case 42 % -- successfully loaded experiment file
                    set(obj.h_init,'Value',0,'Tag','on','CData',imread('initYellow.jpg'))
                    set(obj.h_load,'Value',0,'Tag','on','CData',imread('loadedYellow.jpg'))
                    set(obj.h_pause,'Value',0,'Tag','off','CData',imread('pauseGray.jpg'))
                    set(obj.h_run,'Value',0,'Tag','on','CData',imread('runYellow.jpg'))
                    set(obj.h_abort,'Value',0,'Tag','off','CData',imread('abortGray.jpg'))
                    set(obj.h_respYES,'Value',0,'Tag','off','BackgroundColor',[.7 .7 .7])
                    set(obj.h_respNO,'Value',0,'Tag','off','BackgroundColor',[.7 .7 .7])
                case 50 % Begin running an experiment file (FUNCTION runExperiment)
                    set(obj.h_init,'Value',0,'Tag','off','CData',imread('initGray.jpg'))
                    set(obj.h_load,'Value',0,'Tag','off','CData',imread('loadedGray.jpg'))
                    set(obj.h_pause,'Value',0,'Tag','on','CData',imread('pauseYellow.jpg'))
                    set(obj.h_run,'Value',1,'Tag','off','CData',imread('runGreen.jpg'))
                    set(obj.h_abort,'Value',0,'Tag','on','CData',imread('abortRed.jpg'))
                    set(obj.h_respYES,'Value',0,'Tag','off','BackgroundColor',[.7 .7 .7])
                    set(obj.h_respNO,'Value',0,'Tag','off','BackgroundColor',[.7 .7 .7])
                case 51 % -- unsuccessfully ran experiment
                    for ii=1:flickerReps
                        set(obj.h_run,'Value',1,'Tag','off','CData',imread('runRed.jpg'))
                        pause(flickerLen)
                        set(obj.h_run,'Value',1,'Tag','off','CData',imread('runGreen.jpg'))
                        pause(flickerLen)
                    end
                    set(obj.h_init,'Value',0,'Tag','on','CData',imread('initYellow.jpg'))
                    set(obj.h_load,'Value',0,'Tag','on','CData',imread('loadedYellow.jpg'))
                    set(obj.h_pause,'Value',0,'Tag','off','CData',imread('pauseGray.jpg'))
                    set(obj.h_run,'Value',0,'Tag','on','CData',imread('runYellow.jpg'))
                    set(obj.h_abort,'Value',0,'Tag','off','CData',imread('abortGray.jpg'))
                case 52 % -- successfully ran experiment
                    set(obj.h_init,'Value',0,'Tag','on','CData',imread('initYellow.jpg'))
                    set(obj.h_load,'Value',0,'Tag','on','CData',imread('loadedYellow.jpg'))
                    set(obj.h_pause,'Value',0,'Tag','off','CData',imread('pauseGray.jpg'))
                    set(obj.h_run,'Value',0,'Tag','on','CData',imread('runYellow.jpg'))
                    set(obj.h_abort,'Value',0,'Tag','off','CData',imread('abortGray.jpg'))

                case 60 % Begin aborting run (FUNCTION abortExperiment)
                    set(obj.h_init,'Value',0,'Tag','off','CData',imread('initGray.jpg'))
                    set(obj.h_load,'Value',0,'Tag','off','CData',imread('loadedGray.jpg'))
                    set(obj.h_pause,'Value',0,'Tag','off','CData',imread('pauseGray.jpg'))
                    set(obj.h_run,'Value',0,'Tag','off','CData',imread('runGray.jpg'))
                    set(obj.h_abort,'Value',1,'Tag','off','CData',imread('abortRed.jpg'))
                case 61 % -- failed to abort run
                    for ii=1:flickerReps
                        set(obj.h_abort,'Value',1,'Tag','off','CData',imread('abortGray.jpg'))
                        pause(flickerLen)
                        set(obj.h_abort,'Value',1,'Tag','off','CData',imread('abortRed.jpg'))
                        pause(flickerLen)
                    end
                    set(obj.h_init,'Value',0,'Tag','on','CData',imread('initYellow.jpg'))
                    set(obj.h_load,'Value',0,'Tag','on','CData',imread('loadedYellow.jpg'))
                    set(obj.h_pause,'Value',0,'Tag','off','CData',imread('pauseGray.jpg'))
                    set(obj.h_run,'Value',0,'Tag','on','CData',imread('runYellow.jpg'))
                    set(obj.h_abort,'Value',0,'Tag','off','CData',imread('abortGray.jpg'))
                case 62 % -- successfully aborted run
                    set(obj.h_init,'Value',0,'Tag','on','CData',imread('initYellow.jpg'))
                    set(obj.h_load,'Value',0,'Tag','on','CData',imread('loadedYellow.jpg'))
                    set(obj.h_pause,'Value',0,'Tag','off','CData',imread('pauseGray.jpg'))
                    set(obj.h_run,'Value',0,'Tag','on','CData',imread('runYellow.jpg'))
                    set(obj.h_abort,'Value',0,'Tag','off','CData',imread('abortGray.jpg')) 
                otherwise
            end
            pause(.05)
        catch
            errorTxt = {'  Issue: Button management error.'
                 '  Action: None.'
                 '  Location: in arlas.buttonManager.'
                };
            errorMsgARLas(errorTxt);
            obj.printError(ME)
        end
    end

    % define managing functions
    function presentManager(varargin) % present the stimulus
        obj = varargin{1};
        if strcmp(obj.h_present.Tag,'off')
            obj.h_present.Value = 0;
            return
        end
        if obj.runMode == 2 % when examining old sessions, do NOT make new save
            obj.h_present.Value = 0;
            return
        end
        
        set(obj.h_present,'Value',1,'Tag','off','CData',imread('presentGreen.jpg'))
        pause(0.1)
        
        if obj.autorunNow == 1 % AUTO RUN ---------------------------------
            obj.AUTO.BackgroundColor = [0 1 0];
            while obj.DONE == 0 % while entire test is not finished ------
                if obj.earNow == 1
                    while any(obj.doneL == 0) % while left ear is not finished ----
                        for ii=1:obj.nFreqs
                            %obj.trialCounter = 1;
                            while obj.doneL(ii) == 0
                                obj.freqUpManager
                                obj.prepStim
                                obj.applyISC
                                obj.playStim
                                obj.getResponse
                                obj.analyze
                                obj.audiometer_plotPI
                                obj.checkIsFinished
                                if obj.DONE == 1
                                    break
                                else
                                    %obj.trialCounter = obj.trialCounter + 1;
                                end
                            end
                        end
                    end
                end
                if obj.earNow == 2
                    while any(obj.doneL == 0) % while left ear is not finished ----
                        for ii=1:obj.nFreqs
                            %obj.trialCounter = 1;
                            while obj.doneL(ii) == 0
                                obj.freqUpManager
                                obj.prepStim
                                obj.applyISC
                                obj.playStim
                                obj.getResponse
                                obj.analyze
                                obj.audiometer_plotPI
                                obj.checkIsFinished
                                if obj.DONE == 1
                                    break
                                else
                                    %obj.trialCounter = obj.trialCounter + 1;
                                end
                            end
                        end
                    end
                end
            end
            obj.AUTO.BackgroundColor = [1 1 1];
        else % MANUAL RUN -------------------------------------------------
            if obj.runMode ~= 3 % if actually running data (mode 1)
                obj.prepStim
                obj.applyISC
                obj.playStim
                obj.getResponse
                %obj.checkIsFinished 
            else % for demo mode (3) just ask for a response
                obj.getResponse
            end
        end
        
        set(obj.h_present,'Value',0,'Tag','on','CData',imread('presentYellow.jpg'))
        pause(0.1)
    end
    function freqDnManager(varargin) % control which input channel is currenly being viewed
        obj = varargin{1};
        if strcmp(obj.h_freqDn.Tag,'off')
            obj.h_freqDn.Value = 0;
            return
        end
        
        set(obj.h_freqDn,'Value',1,'Tag','off','CData',imread('freqDnGreen.jpg'))
        pause(0.1)
        
        obj.doSave; % always save first

        obj.fIndx = obj.fIndx - 1; % decrement frequency index by 1
        if obj.fIndx < 1
            obj.fIndx = obj.nFreqs;
        end
        obj.freqNow = obj.freqs(obj.fIndx);
        obj.FREQ.String = [num2str(obj.freqNow),' Hz'];
        
        if obj.earNow == 1
            obj.lvlNow = obj.startingLvls_L(obj.fIndx); %obj.levels(obj.lvlIndx);
        else
            obj.lvlNow = obj.startingLvls_R(obj.fIndx);
        end
        [~,obj.lvlIndx] = min(abs(obj.lvlNow - obj.levels));
        obj.LVL.String = [num2str(obj.lvlNow),' dB ',char(obj.refNow)];

        % update the plots with location of new presentation level --------
        if obj.viewNow == 0 % pi and tracking view --
            obj.audiometer_plotPI
        elseif obj.viewNow == 1 % audiogram view --
            obj.audiometer_plotAG
        end

        set(obj.h_freqDn,'Value',0,'Tag','on','CData',imread('freqDnYellow.jpg'))
        pause(0.1)
    end
    function freqUpManager(varargin) % control which input channel is currenly being viewed
        obj = varargin{1};
        if strcmp(obj.h_freqUp.Tag,'off')
            obj.h_freqUp.Value = 0;
            return
        end
        
        set(obj.h_freqUp,'Value',1,'Tag','off','CData',imread('freqUpGreen.jpg'))
        pause(0.1)

        obj.doSave; % always save first

        obj.fIndx = obj.fIndx + 1;
        if obj.fIndx > obj.nFreqs
            obj.fIndx = 1;
        end
        obj.freqNow = obj.freqs(obj.fIndx);
        obj.FREQ.String = [num2str(obj.freqNow),' Hz'];

        if obj.earNow == 1
            obj.lvlNow = obj.startingLvls_L(obj.fIndx); %obj.levels(obj.lvlIndx);
        else
            obj.lvlNow = obj.startingLvls_R(obj.fIndx);
        end
        [~,obj.lvlIndx] = min(abs(obj.lvlNow - obj.levels));
        obj.LVL.String = [num2str(obj.lvlNow),' dB ',char(obj.refNow)];

        % update the plots with location of new presentation level --------
         if obj.viewNow == 0 % pi and tracking view --
            obj.audiometer_plotPI
        elseif obj.viewNow == 1 % audiogram view --
            obj.audiometer_plotAG
        end
        
        set(obj.h_freqUp,'Value',0,'Tag','on','CData',imread('freqUpYellow.jpg'))
        pause(0.1)
    end    
    function lvlDnManager(varargin) % control which input channel is currenly being viewed
        obj = varargin{1};
        if strcmp(obj.h_lvlDn.Tag,'off')
            obj.h_lvlDn.Value = 0;
            return
        end
        if obj.runMode == 2 % when examining old sessions, do NOT make new save
            obj.h_lvlDn.Value = 0;
            return
        end
        
        set(obj.h_lvlDn,'Value',1,'Tag','off','CData',imread('lvlDnGreen.jpg'))
        pause(0.1)

        obj.lvlIndx = obj.lvlIndx - 1;
        if obj.lvlIndx < 1 % can't go lower than lowest level
            obj.lvlIndx = 1;
            
            flickerLen = 0.1; % length of flicker for error notification
            flickerReps = 5; % number of flickers
            for ii=1:flickerReps
                set(obj.h_lvlDn,'Value',1,'Tag','on','CData',imread('lvlDnRed.jpg'))
                pause(flickerLen)
                set(obj.h_lvlDn,'Value',1,'Tag','on','CData',imread('lvlDnGreen.jpg'))
                pause(flickerLen)
            end            
        end
        obj.lvlNow = obj.levels(obj.lvlIndx);
        obj.LVL.String = [num2str(obj.lvlNow),' dB ',char(obj.refNow)];
        
        % update the plots with location of new presentation level --------
        if obj.viewNow == 0 % pi and tracking view --
            obj.audiometer_plotPI
        elseif obj.viewNow == 1 % audiogram view --
            obj.audiometer_plotAG
        end


        % % update the plots
        % if obj.h_switchPlots.Value == 0 % pi and tracking view --
        %     %obj.audiometer_plotPI
        %     h1 = obj.h_pif; % performance intensity function
        %     axes(h1); % plot the performance-intensity function -----
        %     try delete(obj.marker1); catch; end
        %     ymin = -0.1;
        %     obj.marker1 = plot(obj.lvlNow,ymin,'k^','MarkerFaceColor',[0 0 0]);
        %     h2 = obj.h_track; % response tracker        
        %     axes(h2) % plot the response track -----------------
        %     try delete(obj.marker2); catch; end
        %     if obj.earNow == 1
        %         ear = 'L';
        %     elseif obj.earNow == 2
        %         ear = 'R';
        %     end
        %     A = obj.AAA.(ear).(['f_',num2str(obj.freqNow)]);            
        %     RESP = A.RESP;
        %     N = length(RESP);
        %     obj.marker2 = plot(N+1,obj.lvlNow,'k^','MarkerFaceColor',[0 0 0]);  
        % elseif obj.h_switchPlots.Value == 1 % audiogram view --
        %     obj.audiometer_plotAG
        % end
        
        set(obj.h_lvlDn,'Value',0,'Tag','on','CData',imread('lvlDnYellow.jpg'))
        pause(0.1)
    end    
    function lvlUpManager(varargin) % control which input channel is currenly being viewed
        obj = varargin{1};
        if strcmp(obj.h_lvlUp.Tag,'off')
            obj.h_lvlUp.Value = 0;
            return
        end
        if obj.runMode == 2 % when examining old sessions, do NOT make new save
            obj.h_lvlUp.Value = 0;
            return
        end
        
        set(obj.h_lvlUp,'Value',1,'Tag','off','CData',imread('lvlUpGreen.jpg'))
        pause(0.1)

        obj.lvlIndx = obj.lvlIndx + 1;
        if obj.lvlIndx > obj.nLevels % can't go higher than highest level
            obj.lvlIndx = obj.nLevels;
            
            flickerLen = 0.1; % length of flicker for error notification
            flickerReps = 5; % number of flickers
            for ii=1:flickerReps
                set(obj.h_lvlUp,'Value',1,'Tag','on','CData',imread('lvlUpRed.jpg'))
                pause(flickerLen)
                set(obj.h_lvlUp,'Value',1,'Tag','on','CData',imread('lvlUpGreen.jpg'))
                pause(flickerLen)
            end            
        end
        obj.lvlNow = obj.levels(obj.lvlIndx);
        obj.LVL.String = [num2str(obj.lvlNow),' dB ',char(obj.refNow)];
        
        % update the plots with location of new presentation level --------
        if obj.viewNow == 0 % pi and tracking view --
            obj.audiometer_plotPI
        elseif obj.viewNow == 1 % audiogram view --
            obj.audiometer_plotAG
        end
        
        set(obj.h_lvlUp,'Value',0,'Tag','on','CData',imread('lvlUpYellow.jpg'))
        pause(0.1)
    end
    function earManager(varargin) % control which input channel is currenly being viewed
        obj = varargin{1};
        obj.doSave; % always save first
        if obj.TEST_EAR.Value == 1 % left ear
            obj.earNow = 1;
        elseif obj.TEST_EAR.Value == 2 % right ear
            obj.earNow = 2;
        end
        % update the plots with location of new presentation level --------
        if obj.viewNow == 0 % pi and tracking view --
            obj.audiometer_plotPI
        elseif obj.viewNow == 1 % audiogram view --
            obj.audiometer_plotAG
        end
    end    
    function stepManager(varargin) % control level step size
        obj = varargin{1};
        if obj.STEP.Value == 1
            obj.stepNow = obj.stepSizes(1);
        elseif obj.STEP.Value == 2
            obj.stepNow = obj.stepSizes(2);
        elseif obj.STEP.Value == 3
            obj.stepNow = obj.stepSizes(3);            
        end
        obj.levels = (obj.minOutput(obj.fIndx):obj.stepNow:obj.maxOutput(obj.fIndx))';
        obj.nLevels = length(obj.levels);
        [~,obj.lvlIndx] = min(abs(obj.levels - obj.lvlNow));
        obj.lvlNow = obj.levels(obj.lvlIndx);
        obj.LVL.String = [num2str(obj.lvlNow),' dB ',obj.refNow];
        % update the plots
        if obj.viewNow == 0 % pi and tracking view --
            obj.audiometer_plotPI
        elseif obj.viewNow == 1 % audiogram view --
            obj.audiometer_plotAG
        end
    end    
    function refManager(varargin) % control level reference
        obj = varargin{1};
        if obj.REF.Value == 1
            obj.refNow = obj.lvlRefs(1);
        elseif obj.REF.Value == 2
            obj.refNow = obj.lvlRefs(2);
        end

        % update the plots with location of new presentation level --------
        if obj.viewNow == 0 % pi and tracking view --
            obj.audiometer_plotPI
        elseif obj.viewNow == 1 % audiogram view --
            obj.audiometer_plotAG
        end
        % update the level string
        if obj.earNow == 1
            ear = 'L';
            if obj.calFinished == 0
                conv = 0;
            else
                if strcmp(char(obj.refNow),'FPL')
                    conv = 0;
                elseif strcmp(char(obj.refNow),'IPL')
                    conv = obj.dBConversionL(obj.fIndx); % converting from fpl to ipl
                end
            end
        elseif obj.earNow == 2
            ear = 'R';
            if obj.calFinished == 0
                conv = 0;
            else
                if strcmp(char(obj.refNow),'FPL')
                    conv = 0;
                elseif strcmp(char(obj.refNow),'IPL')
                    conv = obj.dBConversionR(obj.fIndx); % converting from fpl to ipl
                end
            end
        end
        obj.LVL.String = [num2str(obj.lvlNow+conv,3),' dB ',char(obj.refNow)];
    end
    function autoManager(varargin) % control whether in auto or manual mode
        obj = varargin{1};
        if obj.AUTO.Value == 1
            obj.autorunNow = 0;
            obj.AUTO.BackgroundColor = [1 1 1];
        elseif obj.AUTO.Value == 2 % turn on the auto sequence
            obj.autorunNow = 1;
        end
    end
    function doCal(varargin) % do in-situ calibration
        obj = varargin{1};
        if strcmp(obj.h_cal.Tag,'off')
            obj.h_cal.Value = 0;
            return
        end
        if obj.runMode == 2 | obj.runMode == 3 % when examining old sessions or simulating, do NOT make new save
            obj.h_cal.Value = 0;
            return
        end
        obj.buttonManager(10)
        set(obj.h_cal,'Value',0,'Tag','off','BackgroundColor',[0 1 0])
        
        ao = obj.arlasObj; % ao is the arlas object (as opposed to the audiometer object)
        probeL = obj.probeL;
        probeR = obj.probeR;
        doMicCorrection = obj.doMicCorrection;

        % specify calibration to use to set stimulus levels
        calType = 'thev'; % use 'thev' for Thevenin source
        targetCalType = char(obj.refNow); % 'fpl', 'spl', 'ipl'.
        ISCL = struct;
        ISCR = struct;
        if ~(isempty(probeL) | isempty(probeR)) % if both ears are being tested
            testProbe = 'Both';
        elseif ~(isempty(probeL))
            testProbe = 'Left';
        elseif ~(isempty(probeR))
            testProbe = 'Right';
        end
        
        if ~(isempty(probeL) | isempty(probeR)) % if both ears are being tested
            testEar = 'Both';
        else
            % Get stimulus ear:
            prompt = {'Enter the test ear (Left or Right)'};
            title = 'Test Ear'; 
            defaultAns = {'Right'};
            numlines = 1;
            answer = inputdlg(prompt,title,numlines,defaultAns);        
            if isempty(answer) % user hit cancel
                return
            elseif strcmp(answer{1},'') % if user left field blank
                return
            else % user put something in the field
                % check for legal input
                testEar = answer{1}; % voltage rms must be positive
                if strcmp(testEar,'left')
                    testEar = 'Left';
                end
                if strcmp(testEar,'right')
                    testEar = 'Right';
                end
                if ~(strcmp(testEar,'Left') | strcmp(testEar,'Right'))
                    disp('Error: Invalid input. Must be Left or Right.')
                    return
                end
            end
        end
        
        % get probe settings
        validateStimLevels = [];
        [probeInputL,probeOutputL,probeInputR,probeOutputR,~] = getProbeSettings(probeL,probeR,validateStimLevels);
        % get most recent calibration files -----
        if isempty(obj)
            tempobj = struct;
            tempobj.map = [];
        else
            tempobj = obj;
        end    
        micCorrectionL = getMicCorrection(probeInputL,doMicCorrection,tempobj.map);
        micCorrectionR = getMicCorrection(probeInputR,doMicCorrection,tempobj.map);
        if doMicCorrection==1
            if ~isempty(probeInputL)
                if micCorrectionL == 1
                    warning('MIC CORRECTION IS NOT BEING APPLIED!')
                end
            end
            if  ~isempty(probeInputR)
                if micCorrectionR == 1
                    warning('MIC CORRECTION IS NOT BEING APPLIED!')
                end
            end
        end
        [C1L,C2L,calPath1L,calPath2L] = getOutputCal(calType,probeOutputL,tempobj.map);
        [C1R,C2R,calPath1R,calPath2R] = getOutputCal(calType,probeOutputR,tempobj.map);
        
        % PERFORM IN-SITU CALIBRATION -----
        fmin = 100;
        fmax = 18000;
        disp('----- Running in-situ calibration -----')
        inSituReps = 6;
        doIndividual = 1;
        doSimultaneous = 0;
        [iscS1L,iscS2L,~] = ARLas_runISC(obj.arlasObj,probeInputL,probeOutputL,calPath1L,calPath2L,calType,fmin,fmax,micCorrectionL,inSituReps,doIndividual,doSimultaneous);
        [iscS1R,iscS2R,~] = ARLas_runISC(obj.arlasObj,probeInputR,probeOutputR,calPath1R,calPath2R,calType,fmin,fmax,micCorrectionR,inSituReps,doIndividual,doSimultaneous);
        
        % these strucures will be used to apply the in-situ calibration
        obj.iscS1L = iscS1L;
        obj.iscS2L = iscS2L;
        obj.iscS1R = iscS1R;
        obj.iscS2R = iscS2R;
        obj.probeInputL = probeInputL;
        obj.probeInputR = probeInputR;
        obj.probeOutputL = probeOutputL;
        obj.probeOutputR = probeOutputR;

        % Get the clicker information (for auto mode only)
        % [inputs,outputs] = hardwareSetup; % read in the hardware setup
        % N = length(inputs);
        % done = 0;
        % counter = 1;
        % while done == 0
        %     label = inputs{counter}.label;
        %     if strcmp(label,'clicker')
        %         done = 1;
        %         obj.clickerInput = inputs{counter};
        %     else
        %         counter = counter + 1;
        %     end
        %     if counter > N
        %         done = 1;
        %         error('Unable to locate clicker. Make sure there is a clicker label in hardware setup.')
        %     end
        % end

        disp('----- Finished in-situ calibration -----')
        % get IPL conversions for left ear
        if obj.whichLoudspeaker == 1
            iscNow = obj.iscS1L;
        else
            iscNow = obj.iscS2L;
        end
        lvlNow = 70;
        for ii=1:obj.nFreqs
            stim = obj.Stimuli(:,ii);
            [~,fplScaling(ii,1),~] = ARLas_applyISC(iscNow,lvlNow,'fpl',obj.freqs(ii),stim(:));
            [~,iplScaling(ii,1),~] = ARLas_applyISC(iscNow,lvlNow,'ipl',obj.freqs(ii),stim(:));
        end
        obj.dBConversionL = 20*log10(abs(fplScaling) ./ abs(iplScaling));
        % get IPL conversions for right ear
        if obj.whichLoudspeaker == 1
            iscNow = obj.iscS1R;
        else
            iscNow = obj.iscS2R;
        end
        lvlNow = 70;
        for ii=1:obj.nFreqs
            stim = obj.Stimuli(:,ii);
            [~,fplScaling(ii,1),~] = ARLas_applyISC(iscNow,lvlNow,'fpl',obj.freqs(ii),stim(:));
            [~,iplScaling(ii,1),~] = ARLas_applyISC(iscNow,lvlNow,'ipl',obj.freqs(ii),stim(:));
        end
        obj.dBConversionR = 20*log10(abs(fplScaling) ./ abs(iplScaling));
        
        if isempty(obj.probeInputL) & isempty(obj.probeInputR)
            obj.calFinished = 0;
            warning('Error: both Left and Right probe inputs are empty.')
        else
            obj.calFinished = 1;
            % calibration info
            obj.audiogramData.iscS1L = obj.iscS1L;
            obj.audiogramData.iscS2L = obj.iscS2L;
            obj.audiogramData.iscS1R = obj.iscS1R;
            obj.audiogramData.iscS2R = obj.iscS2R;
            %obj.audiogramData.probeOutput = obj.probeOutput;
            %obj.audiogramData.probeOutputR = obj.probeOutputR;
            %obj.audiogramData.probeInputL = obj.probeInputL;
            %obj.audiogramData.probeInputR = obj.probeInputR;
            obj.audiogramData.micCorrectionL = obj.micCorrectionL;
            obj.audiogramData.micCorrectionR = obj.micCorrectionR;
            obj.audiogramData.C1L = obj.C1L;
            obj.audiogramData.C2L = obj.C2L;
            obj.audiogramData.C1R = obj.C1R;
            obj.audiogramData.C2R = obj.C2R;
            obj.audiogramData.probeR = obj.probeR;
            obj.audiogramData.probeL = obj.probeL;
            obj.audiogramData.doMicCorrection = obj.doMicCorrection;
            obj.audiogramData.whichLoudspeaker = obj.whichLoudspeaker; % most probes have 2 loudspeakers--which one are you using (1 or 2)
            obj.audiogramData.dBConversionL = obj.dBConversionL; % conversion factor to make FPL into IPL (in dB)
            obj.audiogramData.dBConversionR = obj.dBConversionR;
        end
        obj.buttonManager(20)
    end
    function doSave(varargin) % save the results to excel spreadsheet
        obj = varargin{1};
        if strcmp(obj.h_save.Tag,'off')
            obj.h_save.Value = 0;
            return
        end
        if obj.runMode == 2 % when examining old sessions, do NOT make new save
            obj.h_save.Value = 0;
            return
        end
        obj.buttonManager(10)
        set(obj.h_cal,'Value',0,'Tag','off','BackgroundColor',[.7 .7 .7])
        set(obj.h_save,'Value',0,'Tag','off','BackgroundColor',[0 1 0])
        
        % first extract thresholds
        freqs = obj.AAA.freqs';
        nFreqs = obj.AAA.nFreqs;
        for ii=1:nFreqs
            dummy = obj.AAA.L.(['f_',num2str(freqs(ii))]).thd;
            if ~isempty(dummy)
                ThdL(ii,1) = dummy;
            else
                ThdL(ii,1) = NaN;
            end
            dummy = obj.AAA.R.(['f_',num2str(freqs(ii))]).thd;
            if ~isempty(dummy)
                ThdR(ii,1) = dummy;
            else
                ThdR(ii,1) = NaN;
            end
        end
        
        if obj.runMode == 3 % if running in simulation mode
            dummy = which('arlas'); % get the location of the currently-running verison of arlas.m
            sep = filesep;
            indx = length(dummy); % strip off the file name at the end
            stop = 0;
            while ~stop
                if strcmp(dummy(indx),sep)
                    stop = 1;
                    dummy = dummy(1:indx); % location of arlas
                    base = dummy(1:end-13); % base location is dummy less the 13 characters: Core\classes\
                end
                indx = indx - 1;
            end
            map = pathSetup(base);
            subjectID = 'HumptyD';
            experimentID = 'Audiophile';
            d = datetime('today');
            formatOut = 'dd';
            str1 = datestr(d,formatOut,'local');
            formatOut = 'mmm';
            str2 = datestr(d,formatOut,'local');
            str2 = upper(str2);
            formatOut = 'yyyy';
            str3 = datestr(d,formatOut,'local');
            runDate = [str1,str2,str3];
            savingGrace = [map.data,subjectID,sep,runDate,sep,experimentID];
            ao.savingGrace = savingGrace;
            ao.experimentID = experimentID;
            ao.sep = obj.sep;
            ao.subjectID = subjectID;
        else
            ao = obj.arlasObj;
        end
        try 
            savePath = [ao.savingGrace,ao.experimentID,'_analysis',ao.sep];
            if exist(savePath,'dir') == 0 % if path does not exist
                success = mkdir(savePath);
                if ~success
                    warning('Unable to create new Experiment Directory: data save path')
                end
            end 
            saveFileName = [ao.subjectID,'_analyzedAUDIO.mat'];
            %saveFileName = ARLas_saveName(savePath,saveFileName);
            AUDIO = obj.AAA;
            AUDIO.timeStamp = cellstr(datetime('now'));
            AUDIO.subjectID = ao.subjectID;
            AUDIO.ThdL = ThdL;
            AUDIO.ThdR = ThdR;
            
            %objCopy = obj; 
            %objCopy.arlasObj = [];
            %objCopy.H = [];
            obj.audiogramData.AAA = obj.AAA;
            obj.audiogramData.pdPrior = obj.pdPrior;
            audiogramData = obj.audiogramData;
            %warning off
            save([savePath,saveFileName],'AUDIO','audiogramData')
            %warning on
        catch
            disp('Warning: Audiometric Analysis not saved!')
        end
        
        % write thresholds to excel file here ----------
        try
        saveFileName = [ao.subjectID,'_analyzedAUDIO.xlsx'];
        catch ME
            keyboard
        end
        obj.writeData2Excel(savePath,saveFileName,freqs,ThdL,ThdR)

        obj.buttonManager(20)
    end
    function doSwitch(varargin) % switch between audiogram and pi/tracking views
        obj = varargin{1};
        if strcmp(obj.h_switchPlots.Tag,'off')
            return
        end
        obj.doingSwitch = 1;
        %obj.buttonManager(10)

        set(obj.h_switchPlots,'Tag','on','BackgroundColor',[0 1 0])
        pause(0.01)

        % change the view
        if obj.viewNow == 0 % current view is PI
            obj.viewNow = 1; % new view is audiogram
            obj.h_switchPlots.Value = 1;
            try delete(obj.marker1); catch; end
            try delete(obj.marker2); catch; end
            try delete(obj.PIjunk{7}); catch; end
            try delete(obj.PIjunk{6}); catch; end
            try delete(obj.PIjunk{5}); catch; end
            try delete(obj.PIjunk{4}); catch; end
            try delete(obj.PIjunk{3}); catch; end
            try delete(obj.PIjunk{2}); catch; end
            try delete(obj.PIjunk{1}); catch; end        
            try delete(obj.PIjunk2); catch; end 
            try delete(obj.PRjunk{5}); catch; end
            try delete(obj.PRjunk{4}); catch; end
            try delete(obj.PRjunk{3}); catch; end
            try delete(obj.PRjunk{2}); catch; end
            try delete(obj.PRjunk{1}); catch; end
            try delete(obj.PRjunk2); catch; end

            axes(obj.h_pif) % performance intensity function -----
            axis off            
            axes(obj.h_track) % tracking progress
            axis off

            axes(obj.h_ag) % audiogram -----
            axis on
            %obj.h_ag.Visible = 'on';
            set(obj.h_switchPlots,'Tag','on','BackgroundColor',[1 1 0])
            obj.audiometer_plotAG
        elseif obj.viewNow == 1 % current view is audiogram
            obj.viewNow = 0; % new view is PI
            obj.h_switchPlots.Value = 0;
            % turn off the audiogram stuff
            axes(obj.h_ag) % audiogram ; NOTE: obj.h_ag.Visible = 'off' doesnt work well
            axis off
            try delete(obj.L1); catch; end
            try delete(obj.L2); catch; end
            try delete(obj.markerX); catch; end
            try delete(obj.markerY); catch; end
            try delete(obj.LX); catch; end
            try delete(obj.LY); catch; end
            try delete(obj.LfHL); catch; end
            try delete(obj.LfHLa); catch; end
            try delete(obj.AGjunk{1}); catch; end
            try delete(obj.thdDots); catch; end
            for ii=1:length(obj.freqs)
                try delete(obj.thdMarks{ii}); catch; end
            end
            for ii=1:length(obj.freqs)
                try delete(obj.AGjunk2{ii}); catch; end
            end

            axes(obj.h_pif) % performance intensity function
            axis on            
            %Aobj.h_pif.Visible = 'on';
            axes(obj.h_track) % progress tracking
            axis on
            %obj.h_track.Visible = 'on';
            set(obj.h_switchPlots,'Tag','on','BackgroundColor',[1 1 0])
            obj.audiometer_plotPI
        end
        obj.doingSwitch = 0;
        %obj.buttonManager(20)
    end
    function doYes(varargin)
        obj = varargin{1};
        if strcmp(obj.h_respYES.Tag,'off')
            return
        end
        if obj.earNow == 1
            ear = 'L';
        elseif obj.earNow == 2
            ear = 'R';
        end
        trialCounter = obj.AAA.(ear).(['f_',num2str(obj.freqNow)]).trialCounter+1; % increment counter here
        Y = obj.AAA.(ear).(['f_',num2str(obj.freqNow)]).Y;
        X = obj.AAA.(ear).(['f_',num2str(obj.freqNow)]).X;
        X(trialCounter,1) = obj.lvlNow;
        Y(trialCounter,1) = 1;
        obj.AAA.(ear).(['f_',num2str(obj.freqNow)]).Y = Y;
        obj.AAA.(ear).(['f_',num2str(obj.freqNow)]).X = X;
        obj.AAA.(ear).(['f_',num2str(obj.freqNow)]).trialCounter = trialCounter;

        obj.analyze
        if obj.viewNow == 0 % pi and tracking view --
            obj.audiometer_plotPI
        elseif obj.viewNow == 1 % audiogram view --
            obj.audiometer_plotAG
        end
        set(obj.h_respYES,'Value',0,'Tag','off','BackgroundColor',[.7 .7 .7])
        set(obj.h_respNO,'Value',0,'Tag','off','BackgroundColor',[.7 .7 .7])
    end
    function doNo(varargin)
        obj = varargin{1};
        if strcmp(obj.h_respNO.Tag,'off')
            return
        end
        if obj.earNow == 1
            ear = 'L';
        elseif obj.earNow == 2
            ear = 'R';
        end
        trialCounter = obj.AAA.(ear).(['f_',num2str(obj.freqNow)]).trialCounter+1; % increment counter here
        Y = obj.AAA.(ear).(['f_',num2str(obj.freqNow)]).Y;
        X = obj.AAA.(ear).(['f_',num2str(obj.freqNow)]).X;
        X(trialCounter,1) = obj.lvlNow;
        Y(trialCounter,1) = 0;
        obj.AAA.(ear).(['f_',num2str(obj.freqNow)]).Y = Y;
        obj.AAA.(ear).(['f_',num2str(obj.freqNow)]).X = X;
        obj.AAA.(ear).(['f_',num2str(obj.freqNow)]).trialCounter = trialCounter;

        obj.analyze
        if obj.viewNow == 0 % pi and tracking view --
            obj.audiometer_plotPI
        elseif obj.viewNow == 1 % audiogram view --
            obj.audiometer_plotAG
        end
        set(obj.h_respYES,'Value',0,'Tag','off','BackgroundColor',[.7 .7 .7])
        set(obj.h_respNO,'Value',0,'Tag','off','BackgroundColor',[.7 .7 .7])
    end

    % define supporting functions
    function getStim(varargin) % create the stimuli
        obj = varargin{1};
        % Pulsed tones, per ANSI S3.6-2004
        % plateau duration of and and off 225 ms +/- 35 ms
        % rise fall within 20-50 ms
        %
        % Here, we will use 3 on pulses and 4 off pulses
        % Rise fall will be 20 ms
        % Steady state will be 150 ms
        % [150 20 150 20 150 20 150 20 150 20 150 20 150] =
        %  off    on     off    on     off     on    off   
        %
        % By starting and ending with off cycles, if the system delay cuts
        % a few of those off, won't matter and won't splatter the signal
        % The total signal duration is 1.170 seconds
        
        %len = 1.17; % stimulus length (s)
        if obj.autorespNow == 1
            respLen = 2; % response window length(s)
        else
            respLen = 0; % in manual mode, there is no response window
        end

        wLen = 0.02;
        wN = round(wLen * obj.fs);
        if mod(wN,2)~=0
            wN = wN + 1;
        end
        h = hann(wN*2);
        h1 = h(1:wN);
        h2 = h(wN+1:end);
        ssLen = 0.150;
        ssN = round(ssLen * obj.fs);
        zpad = zeros(ssN,1);
        opad = ones(ssN,1);

        mask = [zpad;h1;opad;h2;zpad;h1;opad;h2;zpad;h1;opad;h2;zpad];

        N = length(mask);
        respN = round(respLen * obj.fs);
        t = (0:1:N-1)'/obj.fs; % time vector
        nSamples = N;

        Stim = zeros(nSamples,obj.nFreqs);
        for ii=1:obj.nFreqs
           Stim(:,ii) = sin(2*pi*obj.freqs(ii)*t) .* mask;  
        end

        zpad = zeros(respN,1);
        Zpad = repmat(zpad,1,obj.nFreqs);
        obj.Stimuli = [Stim;Zpad];        
        
        NN = size(obj.Stimuli,1);
        foldLen = 0.1;
        foldN = round(obj.fs * foldLen);
        nFolds = ceil(NN / foldN);

        if obj.autorespNow == 1
            extraCount = 4; % add in some extra padding to make it work
            disp('This doesnt work yet!')
            keyboard
            zpad = zeros(foldN*extraCount,1);
            qq = 1;
            stim = obj.Stimuli(:,qq);
            
            stim = reshape(stim,foldN,nFolds);
            nReps = size(stim,2);
            NFolds = nFolds + extraCount;
        else
            stim = obj.Stimuli(:,1);
            k = floor(length(stim)/foldN);
            tag = (foldN * (k+1))- length(stim);
            zpad = zeros(tag,1);
            stim = [stim;zpad];
            nFolds = length(stim) / foldN;
            Zpad = repmat(zpad,1,obj.nFreqs);
            Stimuli = obj.Stimuli;
            Stimuli = [Stimuli;Zpad];
            obj.Stimuli = Stimuli;
            stim = reshape(Stimuli(:,1),foldN,nFolds);
            nReps = size(stim,2);
        end
        obj.foldN = foldN;
        obj.NFolds = nFolds;
        obj.nReps = nReps;
        obj.zpad = [];
    end
    function applyISC(varargin) % Apply the in-situ calibration to the stimuli
        obj = varargin{1};
        if obj.earNow == 1 % left ear
            if obj.whichLoudspeaker == 1
                obj.iscNow = obj.iscS1L;
            else
                obj.iscNow = obj.iscS2L;
            end
        elseif obj.earNow == 2 % right ear
            if obj.whichLoudspeaker == 1
                obj.iscNow = obj.iscS1R;
            else
                obj.iscNow = obj.iscS2R;
            end
        end
        stim = obj.stimNow;
        [rows,cols] = size(stim);
        if strcmp(char(obj.refNow),'FPL')
            targetCalType = 'fpl';
        elseif strcmp(char(obj.refNow),'SPL')
            targetCalType = 'spl';
        elseif strcmp(char(obj.refNow),'IPL')
            targetCalType = 'ipl';
        end
        
        [s1,~,~] = ARLas_applyISC(obj.iscNow,obj.lvlNow+3,targetCalType,obj.freqNow,stim(:));
        obj.stimNow = reshape(s1,rows,cols);
    end
    function prepStim(varargin) % get the current stimulus and fold it for presentation efficiency
        obj = varargin{1};
        Stim = obj.Stimuli(:,obj.fIndx);
        obj.stimNow = [Stim;obj.zpad];
        obj.stimNow = reshape(obj.stimNow,obj.foldN,obj.NFolds);
    end
    function playStim(varargin) % present the stimulus to subject
        obj = varargin{1};
        ao = obj.arlasObj;

        % set the play list -----------------------------------------------
        try
        ao.clearPlayList % clear out whatever was used previously
        catch ME
            keyboard
        end
        if obj.earNow == 1 % left ear
            try
                ao.setPlayList(obj.stimNow,obj.probeOutputL.ch(1));
            catch
                ao.setPlayList(obj.stimNow,obj.probeOutputL.ch(2));
            end
        end
        if obj.earNow == 2 % right ear
            try
                ao.setPlayList(obj.stimNow,obj.probeOutputR.ch(1));
            catch
                ao.setPlayList(obj.stimNow,obj.probeOutputR.ch(2));
            end
        end
        % set the rec list ------------------------------------------------
        ao.clearRecList % clear out whatever was used previously 
        %ao.setRecList(obj.clickerInput.ch,obj.clickerInput.label,obj.clickerInput.micSens,obj.clickerInput.gain);
        if obj.earNow == 1 % left ear
            if ~contains(obj.probeInputL.label,'_audio')
                obj.probeInputL.label = [obj.probeInputL.label,'_audio'];
            end
            ao.setRecList(obj.probeInputL.ch,obj.probeInputL.label,obj.probeInputL.micSens,obj.probeInputL.gain);
        end
        if obj.earNow == 2 % right ear
            if ~contains(obj.probeInputR.label,'_audio')
                obj.probeInputR.label = [obj.probeInputR.label,'_audio'];
            end
            ao.setRecList(obj.probeInputR.ch,obj.probeInputR.label,obj.probeInputR.micSens,obj.probeInputR.gain);
        end
        ao.setNReps(obj.nReps); % number of times to play stimulus
        ao.setFilter(0); % note: this is a highpass filter with a 75 Hz cutoff frequency.
        % set the user info -----------------------------------------------
        ao.objPlayrec.userInfo = []; % clear out previous, non-structure array info
        % ao.objPlayrec.userInfo.f = obj.freqNow;
        % ao.objPlayrec.userInfo.fs = obj.fs;
        % ao.objPlayrec.userInfo.audiometer_version = obj.version;
        % ao.objPlayrec.userInfo.calType = obj.refNow;
        % ao.objPlayrec.userInfo.C1L = obj.C1L;
        % ao.objPlayrec.userInfo.C2L = obj.C2L;
        % ao.objPlayrec.userInfo.C1R = obj.C1R;
        % ao.objPlayrec.userInfo.C2R = obj.C2R;
        % ao.objPlayrec.userInfo.iscS1 = obj.iscS1L;
        % ao.objPlayrec.userInfo.iscS2 = obj.iscS2L;
        % ao.objPlayrec.userInfo.iscS1 = obj.iscS1R;
        % ao.objPlayrec.userInfo.iscS2 = obj.iscS2R;
        % ao.objPlayrec.userInfo.probeR = obj.probeR;
        % ao.objPlayrec.userInfo.probeL = obj.probeL;
        if obj.autorunNow % if running on autopilot, 
            pauseLen = rand(1,1)*1;
            pause(pauseLen) % add some randomness to the presentation delay
        end
        ao.objPlayrec.servedNeat = 1; % put arlas into minimalist mode (no recordign, no alighment of input/output offsets--no frills!)
        ao.objPlayrec.run % playback and record -----
        if ao.killRun
            obj.DONE = 1;
            obj.AAA.DONE = 1;
            return
        end
    end
    function getResponse(varargin) % extract the recordings to see if button was pressed (tone was heard)
        obj = varargin{1};
        %ao = obj.arlasObj;
        set(obj.h_respYES,'Value',0,'Tag','on','BackgroundColor',[0 1 0])
        set(obj.h_respNO,'Value',0,'Tag','on','BackgroundColor',[1 0 0])
    end
    function analyze(varargin) % analyze and get next level (if in auto mode)
        obj = varargin{1};

        if obj.earNow == 1
            ear = 'L';
        elseif obj.earNow == 2
            ear = 'R';
        end
        X = obj.AAA.(ear).(['f_',num2str(obj.freqNow)]).X;
        Y = obj.AAA.(ear).(['f_',num2str(obj.freqNow)]).Y;
        trialCounter = obj.AAA.(ear).(['f_',num2str(obj.freqNow)]).trialCounter;
        
        % audiometer performance-intensity function from data
        mu = obj.pdPrior.mu(obj.fIndx);
        sigma = obj.pdPrior.sigma(obj.fIndx);
        pd = makedist('normal','mu',mu,'sigma',sigma);
        OUT = myLogisticReg(X,Y,pd);

        if trialCounter >= 5
            thd = OUT.thd + 5-mod(OUT.thd,5); % set starting levels for testing near threshold
            obj.startingLvls_L(obj.fIndx) = thd;
            obj.startingLvls_R(obj.fIndx) = thd;
            try % also try to set adjacent starting levels to current threshold
                obj.startingLvls_L(obj.fIndx-1) = thd;
            end
            try
                obj.startingLvls_L(obj.fIndx+1) = thd;
            end
            try % also try to set adjacent starting levels to current threshold
                obj.startingLvls_R(obj.fIndx-1) = thd;
            end
            try
                obj.startingLvls_R(obj.fIndx+1) = thd;
            end
        end

        obj.AAA.(ear).(['f_',num2str(obj.freqNow)]).X = X;
        obj.AAA.(ear).(['f_',num2str(obj.freqNow)]).Y = Y;
        obj.AAA.(ear).(['f_',num2str(obj.freqNow)]).trialCounter = trialCounter;
        obj.AAA.(ear).(['f_',num2str(obj.freqNow)]).thd = OUT.thd; % threshold (50% point on logistic function)
        obj.AAA.(ear).(['f_',num2str(obj.freqNow)]).sigma = OUT.sigma; % spread of normal CDF approximation of logistic function
        obj.AAA.(ear).(['f_',num2str(obj.freqNow)]).CDF = OUT.CDF; % normal CDF (mu=thd, sigma) at X = B = stimulus levels
        obj.AAA.(ear).(['f_',num2str(obj.freqNow)]).B = OUT.B; % stimulus levels over which cdf is evaluated
        obj.AAA.(ear).(['f_',num2str(obj.freqNow)]).badIndx = OUT.badIndx; % indices of identified false positives/negatives
        obj.AAA.(ear).(['f_',num2str(obj.freqNow)]).nRejects = OUT.nRejects; % number of identified false positives/negatives
        obj.AAA.(ear).(['f_',num2str(obj.freqNow)]).prior = OUT.prior;
    end
    function countReversals(varargin) % figure out if threshold criterion reached
        obj = varargin{1};
        
        if obj.earNow == 1
            ear = 'L';
        elseif obj.earNow == 2
            ear = 'R';
        end        
        LVL = obj.AAA.(ear).(['f_',num2str(obj.freqs(obj.fIndx))]).LVL; 
        
        counter = 2;
        N = length(LVL);
        peakIndx = [];
        while counter < N-1
            if LVL(counter-1)<LVL(counter) && LVL(counter+1)<LVL(counter)
                peakIndx = [peakIndx;counter];
            end
            counter = counter + 1;
        end
        if isempty(peakIndx)
            nReversals = 0;
        else
            nReversals = length(peakIndx);
        end
        obj.AAA.(ear).(['f_',num2str(obj.freqs(obj.fIndx))]).nReversals = nReversals; 
    end
    
    function audiometer_plotPI(varargin) % plot audiometer PI function and history
        obj = varargin{1};

        if obj.earNow == 1
            ear = 'L';
            if obj.calFinished == 0
                conv = 0;
            else
                if strcmp(char(obj.refNow),'FPL')
                    conv = 0;
                elseif strcmp(char(obj.refNow),'IPL')
                    conv = obj.dBConversionL(obj.fIndx); % converting from fpl to ipl
                end
            end
        elseif obj.earNow == 2
            ear = 'R';
            if obj.calFinished == 0
                conv = 0;
            else
                if strcmp(char(obj.refNow),'FPL')
                    conv = 0;
                elseif strcmp(char(obj.refNow),'IPL')
                    conv = obj.dBConversionR(obj.fIndx); % converting from fpl to ipl
                end
            end
        end
        A = obj.AAA.(ear).(['f_',num2str(obj.freqNow)]); % results structure
        xmin = obj.minOutput(obj.fIndx) + conv;
        xmax = obj.maxOutput(obj.fIndx) + conv;

        X = A.X + conv; % x values = stimulus levels
        Y = A.Y; % y values = yes / no response
        thd = A.thd + conv; 
        trialCounter = A.trialCounter;
        sigma = A.sigma;
        CDF = A.CDF;
        B = A.B + conv;
        badIndx = A.badIndx;
        nRejects = A.nRejects;
        prior = A.prior;

        % plot performance-intensity function ------------------
        h1 = obj.h_pif; % performance intensity function
        h2 = obj.h_track; % response tracker        
        
        ymin = -0.1;
        ymax = 1.05; 
        N = length(Y);
        
        axes(h1); % plot the performance-intensity function -----
        try delete(obj.PIjunk{8}); catch; end
        try delete(obj.PIjunk{7}); catch; end
        try delete(obj.PIjunk{6}); catch; end
        try delete(obj.PIjunk{5}); catch; end
        try delete(obj.PIjunk{4}); catch; end
        try delete(obj.PIjunk{3}); catch; end
        try delete(obj.PIjunk{2}); catch; end
        try delete(obj.PIjunk{1}); catch; end        
        try delete(obj.PIjunk2); catch; end 
        if isempty(X)
            hold off
            obj.PIjunk{2} = plot(0,0,'.w'); % plot individual responses
            hold on
            b_min = obj.minOutput(obj.fIndx) + conv;
            b_max = obj.maxOutput(obj.fIndx) + conv;
            N = 200;
            B = linspace(b_min,b_max,N)' + conv;
            mu = obj.pdPrior.mu(obj.fIndx) + conv;
            sigma = obj.pdPrior.sigma(obj.fIndx);
            pd = makedist('normal','mu',mu,'sigma',sigma);
            prior = pd.cdf(B);
            obj.PIjunk{4} = plot(B,prior,'--','Color',[.6 .6 .6]); % plot prior propability of threshold values
            xlim([b_min,b_max])
            ylim([ymin,ymax])
            thd = pd.median;
            obj.PIjunk{6} = line([thd,thd],[ymin,ymax],'Color',[0 0 0],'LineStyle','--','LineWidth',1); % plot current threshold
            xlabel(['Stim Level (dB ',char(obj.refNow),')'])
            ylabel('CDF')
            grid on
        else
            hold off
            obj.PIjunk{1} = fill([thd-sigma*2,thd+sigma*2,thd+sigma*2,thd-sigma*2],[ymin,ymin,ymax,ymax],[.9 .9 .9],'EdgeColor','none'); % 95% confidence interval
            hold on
            indx = find(Y == 1);
            obj.PIjunk{2} = plot(X(indx),Y(indx)+0.01*randn(size(X(indx))),'og','MarkerSize',4,'LineWidth',1.5); % individual yes (heard) responses
            indx = find(Y == 0);
            obj.PIjunk{3} = plot(X(indx),Y(indx)+0.01*randn(size(X(indx))),'xk','LineWidth',1.5); % individual no(not heard) responses
            obj.PIjunk{4} = plot(B,prior,'--','Color',[.6 .6 .6]); % plot prior propability of threshold values
            if obj.earNow == 1 % left ear
                color = [0 0 1];
            elseif obj.earNow == 2 % right ear
                color = [1 0 0];
            end
            obj.PIjunk{5} = plot(B,CDF,'Color',color); % plot current probability of threshold values
            obj.PIjunk{6} = line([thd,thd],[ymin,ymax],'Color',[0 0 0],'LineStyle','--','LineWidth',1); % plot current threshold
            if nRejects > 0
                for jj=1:length(nRejects)
                    obj.PIjunk2{jj} = plot(X(badIndx),Y(badIndx),'*','Color',[.7 .7 .7],'LineWidth',1.5);
                end
            end
            xlim([xmin xmax])
            ylim([ymin,ymax])
            grid on

            try
                obj.PIjunk{7} = title([obj.subjID,':  ',num2str(obj.freqNow),' Hz Threshold: ',num2str(round(thd)),' dB ',obj.refNow{1},' +/- ',num2str(sigma*2,2),' dB']);
            catch
                obj.PIjunk{7} = title(obj.subjID);
            end
        end
        xlabel(['Stim Level (dB ',char(obj.refNow),')'])
        ylabel('CDF')
        grid on
        pause(0.01)

        try delete(obj.marker1); catch; end
        if obj.earNow == 1
            color = [0 0 1];
        elseif obj.earNow == 2
            color = [1 0 0];
        end
        obj.marker1 = plot(obj.lvlNow+conv,ymin,'k^','MarkerFaceColor',color);

        ylim([ymin,ymax])

        if trialCounter < 5
            xlim([xmin xmax])
        else
            xlim([min(X)-5,max(X)+5])
        end
    
        axes(h2) % plot the response track --------------------------------
        try delete(obj.PRjunk{5}); catch; end
        try delete(obj.PRjunk{4}); catch; end
        try delete(obj.PRjunk{3}); catch; end
        try delete(obj.PRjunk{2}); catch; end
        try delete(obj.PRjunk{1}); catch; end
        try delete(obj.PRjunk2); catch; end

        if isempty(X)
            hold off
            obj.PRjunk{1} = plot(1,0,'.w');
            hold on
            obj.PRjunk{2} = line([1 6],[thd thd],'Color',[0 0 0],'LineWidth',0.5,'LineStyle','--');
            xlabel('Trial Number')
            ylabel(['Stim Level (dB ',char(obj.refNow),')'])
            ymin = thd-pd.sigma*3;
            ymax = thd+pd.sigma*3;
            if ymin < obj.minOutput(obj.fIndx)
                ymin = obj.minOutput(obj.fIndx);
            end
            if ymax > obj.maxOutput(obj.fIndx)
                ymax = obj.maxOutput(obj.fIndx);
            end
            ylim([ymin,ymax])
            xlim([1 6])
            xticks((1:1:6))
            grid on
            try delete(obj.marker2); catch; end
            if obj.earNow == 1
                color = [0 0 1];
            elseif obj.earNow == 2
                color = [1 0 0];
            end
            obj.marker2 = plot(1,obj.lvlNow,'k^','MarkerFaceColor',color);
        else
            hold off
            if N < 6
                xmin = 1;
                xmax = 6;
            else
                xmin = 1;
                xmax = N+1;
            end
            obj.PRjunk{1} = fill([xmin,xmax,xmax,xmin],[thd-sigma*2,thd-sigma*2,thd+sigma*2,thd+sigma*2],[.9 .9 .9],'EdgeColor','none');
            hold on
            if N  < 6
                obj.PRjunk{2} = line([1 6],[thd thd],'Color',[0 0 0],'LineWidth',0.5,'LineStyle','--');
            else
                obj.PRjunk{2} = line([1 N+1],[thd thd],'Color',[0 0 0],'LineWidth',0.5,'LineStyle','--');
            end
            x = (1:N)';
            indx1 = find(Y == 1);
            indx0 = find(Y == 0);
            obj.PRjunk{3} = plot(x(indx1),X(indx1),'go','LineWidth',2);
            obj.PRjunk{4} = plot(x(indx0),X(indx0),'kx','LineWidth',2);
            obj.PRjunk{5} = plot(x,X,'.-','Color',color);

            if nRejects > 0
                for jj=1:length(nRejects)
                    obj.PRjunk2{jj} = plot(badIndx,X(badIndx),'*','Color',[.7 .7 .7],'LineWidth',1.5);
                end
            end

            xlabel('Trial Number')
            ylabel(['Stim Level (dB ',char(obj.refNow),')'])
            if N < 5
                ymin = min([thd-sigma*3,min(X),obj.lvlNow]);
                ymax = max([thd+sigma*3,max(X),obj.lvlNow]);
                if ymin < obj.minOutput(obj.fIndx)
                    ymin = obj.minOutput(obj.fIndx);
                end
                if ymax > obj.maxOutput(obj.fIndx)
                    ymax = obj.maxOutput(obj.fIndx);
                end
            else
                ymin = min([thd-sigma*2,min(X),obj.lvlNow]);
                ymax = max([thd+sigma*2,max(X),obj.lvlNow]);
            end
            if ymin == ymax
                ymin = obj.minOutput(1);
                ymax = obj.maxOutput(1);
            end
            ylim([ymin,ymax])
            if N < 6
                xlim([1 6])
                xticks((1:1:6))
            else
                xlim([1 N+1])
                xticks((1:1:N))
            end

            try delete(obj.marker2); catch; end
            if obj.earNow == 1
                color = [0 0 1];
            elseif obj.earNow == 2
                color = [1 0 0];
            end
            obj.marker2 = plot(N+1,obj.lvlNow+conv,'k^','MarkerFaceColor',color);
            grid on
        end
        
    end
    function audiometer_plotAG(varargin) % plot audiometer Audiogram
        obj = varargin{1};

       if obj.earNow == 1
            ear = 'L';
            if obj.calFinished == 0
                conv = 0;
            else
                if strcmp(char(obj.refNow),'FPL')
                    conv = 0;
                elseif strcmp(char(obj.refNow),'IPL')
                    conv = obj.dBConversionL(obj.fIndx); % converting from fpl to ipl
                end
            end
        elseif obj.earNow == 2
            ear = 'R';
            if obj.calFinished == 0
                conv = 0;
            else
                if strcmp(char(obj.refNow),'FPL')
                    conv = 0;
                elseif strcmp(char(obj.refNow),'IPL')
                    conv = obj.dBConversionR(obj.fIndx); % converting from fpl to ipl
                end
            end
        end        
        
        FF = [250,500,1000,2000,4000,8000,16000,20000]';
        LL = (-10:10:110)'+conv;
        ff = log2(FF);
        nf = length(ff);
        ymin = min(LL);
        ymax = max(LL);
        xmin = min(ff);
        xmax = max(ff);
        
        h1 = obj.h_ag; % audiogram
        axes(h1)
        
        % get the current thresholds
        for ii=1:obj.nFreqs
            dummy = obj.AAA.L.(['f_',num2str(obj.freqs(ii))]).thd + conv;
            dummy2 = obj.AAA.L.(['f_',num2str(obj.freqs(ii))]).sigma;
            if isempty(dummy)
                thdL(ii,1) = NaN;
                sigmaL(ii,1) = NaN;
            else
                thdL(ii,1) = dummy;
                sigmaL(ii,1) = dummy2;
            end
        end
        for ii=1:obj.nFreqs
            dummy = obj.AAA.R.(['f_',num2str(obj.freqs(ii))]).thd + conv;
            dummy2 = obj.AAA.R.(['f_',num2str(obj.freqs(ii))]).sigma;
            if isempty(dummy)
                thdR(ii,1) = NaN;
                sigmaR(ii,1) = NaN;
            else
                thdR(ii,1) = dummy;
                sigmaR(ii,1) = dummy2;
            end
        end 
        
        if obj.doingSwitch == 1 % if switching views, redraw the level and frequency lines   
            hold off
            plot(0,0,'w.')
            ylim([ymin,ymax])
            xlim([xmin,xmax])
            for ii=1:nf
                axes(h1)
                obj.L1(ii) = line([ff(ii),ff(ii)],[ymin,ymax],'Color',[0 0 0],'LineWidth',0.5,'LineStyle','-');
                hold on
            end
            for ii=1:length(LL)
                axes(h1)
                obj.L2(ii) = line([xmin,xmax],[LL(ii),LL(ii)],'Color',[0 0 0],'LineWidth',0.5,'LineStyle','-');
            end
            obj.LfHL = plot(log2(obj.freqs),obj.fHL+conv,'k--'); % line showing fHL (fpl Hearing Level re: 20 year olds)
            obj.LfHLa = plot(log2(obj.freqs),obj.fHLa+conv,'k:'); % line showing fHLa (fpl Hearing Level re: subject's decade of life)
        end
        ylim([ymin,ymax])
        xlim([xmin,xmax])
        try delete(obj.markerX); catch; end
        try delete(obj.markerY); catch; end
        if obj.earNow == 1
            color = [0 0 1];
        elseif obj.earNow == 2
            color = [1 0 0];
        end
        obj.markerX = plot(log2(obj.freqNow),ymin,'k^','MarkerFaceColor',color);
        obj.markerY = plot(xmin,obj.lvlNow+conv,'k^','MarkerFaceColor',color);
        try delete(obj.LX); catch; end
        try delete(obj.LY); catch; end
        if obj.earNow == 1 % ear is left
            obj.LX = line([log2(obj.freqNow),log2(obj.freqNow)],[ymin,ymax],'Color',[0 0 .5],'LineWidth',1,'LineStyle','-');
            obj.LY = line([xmin,xmax],[obj.lvlNow+conv,obj.lvlNow+conv],'Color',[0 0 .5],'LineWidth',1,'LineStyle','-');
        elseif obj.earNow == 2 % ear is right
            obj.LX = line([log2(obj.freqNow),log2(obj.freqNow)],[ymin,ymax],'Color',[.5 0 0],'LineWidth',1,'LineStyle','-');
            obj.LY = line([xmin,xmax],[obj.lvlNow+conv,obj.lvlNow+conv],'Color',[.5 0 0],'LineWidth',1,'LineStyle','-');
        end
        if all(isnan(thdL)) && all(isnan(thdR)) % there is nothing else to plot
        else

            for jj=obj.nFreqs:-1:1
                try delete(obj.thdMarks{jj}); catch; end
            end
            try delete(obj.thdMarks); catch; end
            try delete(obj.thdDots); catch; end
            try delete(obj.AGjunk{1}); catch; end
            for ii=1:length(obj.freqs)
                try delete(obj.AGjunk2{ii}); catch; end
            end
            hold on
            freqs = log2(obj.freqs);
            % plot audiometric thresholds -------------------------------

            if obj.earNow == 1 % ear is left
                for ii=1:obj.nFreqs
                    obj.AGjunk2{ii} = plot([freqs(ii),freqs(ii)],[thdL(ii)-sigmaL(ii)*2,thdL(ii)+sigmaL(ii)*2],'Color',[.7 .7 .7],'LineWidth',3); % 95% confidence interval
                    obj.thdMarks{ii} = plot(freqs(ii),thdL(ii),'bx','MarkerSize',12,'LineWidth',1);
                end
                obj.thdDots = plot(freqs(:)',thdL,'b-','LineWidth',0.5);
            elseif obj.earNow == 2 % ear is right
                for ii=1:obj.nFreqs
                    obj.AGjunk2{ii} = plot([freqs(ii),freqs(ii)],[thdR(ii)-sigmaR(ii)*2,thdR(ii)+sigmaR(ii)*2],'Color',[.7 .7 .7],'LineWidth',3); % 95% confidence interval
                    obj.thdMarks{ii} = plot(freqs(ii),thdR(ii),'ro','MarkerSize',11,'LineWidth',1);
                end
                obj.thdDots = plot(freqs(:)',thdR,'r-','LineWidth',0.5);                
            end
        end    
        xticks(ff)
        xticklabels({'0.25','0.5','1.0','2.0','4.0','8.0','16.0','20'})
        xlabel('Frequency (kHz)','FontSize',12);
        ylabel(['Stim Level (dB ',char(obj.refNow),')'],'FontSize',12);
        if obj.earNow == 1
            thd = thdL(obj.fIndx);
            sigma = sigmaL(obj.fIndx);
        elseif obj.earNow == 2
            thd = thdR(obj.fIndx);
            sigma = sigmaR(obj.fIndx);
        end
        try
            obj.AGjunk{1} = title([obj.subjID,':  ',num2str(obj.freqNow),' Hz Threshold: ',num2str(round(thd)),' dB ',obj.refNow{1},' +/- ',num2str(sigma*2,2),' dB']);
        catch
            obj.AGjunk{1} = title(obj.subjID);
        end
    end
    
    function checkIsFinished(varargin) % see if finished running (auto mode)
        obj = varargin{1};
        
        done = 0;
        minN = obj.minPresentations; % minimum number of presentations required for 3 ascending runs
        maxN = obj.maxPresentations; % maximum number of allowed presentations; if >, assume something wrong
        if obj.earNow == 1
            ear = 'L';
        elseif obj.earNow == 2
            ear = 'R';
        end
        A = obj.AAA.(ear).(['f_',num2str(obj.freqNow)]);
        A.done = 'not finished';
        lvl = A.LVL;   % vector of presentation levels as a function of trial number
        %resp = A.RESP; % vector of responses as a function of trial number
        N = length(lvl);
        obj.countReversals
        nReversals = A.nReversals;
        %nReversals = obj.nReversals;
        if isempty(nReversals)
            nReversals = 0;
        end

        if nReversals >= obj.minReversals
            A.done = 'min number of reversals achieved';
            done = 1;
        elseif N < minN
            A.done = 'min number of trials not achieved';
            done = 0;
        elseif N > maxN
            warning('Number of trials exceeds maximum number')
            done = 1;
            A.done = 'max number of trials exceeded';
        end
        if done == 1
            if obj.earNow == 1
                obj.doneL(obj.fIndx) = 1;
                if ~any(obj.doneL == 0)
                    obj.DONE = 1;
                end
            elseif obj.earNow == 2
                obj.doneR(obj.fIndx) = 1;
                if ~any(obj.doneR == 0)
                    obj.DONE = 1;
                end
            end
        end
        
        
    end
    function writeData2Excel(varargin)
        obj = varargin{1};
        pathName = varargin{2};
        fileName = varargin{3};
        freqs = varargin{4};
        ThdL = varargin{5};
        ThdR = varargin{6};
        
        warning off
        writematrix(' ',[pathName,fileName],'Sheet',1,'Range','A1');
        writematrix('Freq (Hz)',[pathName,fileName],'Sheet',1,'Range','A3');
        writematrix(freqs(:,1),[pathName,fileName],'Sheet',1,'Range','A4');
        writematrix('Thd L',[pathName,fileName],'Sheet',1,'Range','B3');
        writematrix(ThdL,[pathName,fileName],'Sheet',1,'Range','B4');
        writematrix('Thd R',[pathName,fileName],'Sheet',1,'Range','C3');
        writematrix(ThdR,[pathName,fileName],'Sheet',1,'Range','C4');
        warning on
    end
    function keyPress(varargin) % does not currently work
        obj = varargin{1};
        e = varargin{2};
        switch e.Key
            case 'c'
                obj.freqDnManager
            case 30 % up
            case 31 % down
            case 28 % left
            case 29 % right
            case 'rightarrow'
            case 'leftarrow'
            case 'uparrow'
            case 'downarrow'
            case 'space'
            otherwise
                keyboard
        end
        
    end
    function setAge(varargin)
        obj = varargin{1};
        obj.subjAge = varargin{2};
        obj.getPriors
    end
    function setID(varargin)
        obj = varargin{1};
        subjID = varargin{2};
        if ~isa(subjID,'char')
            error('subjID must be a string.')
        end
        obj.subjID = subjID;
    end
    function getPriors(varargin)
        obj = varargin{1};
        % Find the priors for audiogram, based on Sumit's data
        % FPL thresholds for 357 subjects
        load([obj.priorsPathName,obj.priorsFileName])

        % look only over the subset needed for the current audiometer
        nFreqs = length(obj.freqs);
        for ii=1:nFreqs
            [~,freqIndx(ii,1)] = min(abs(FREQS-obj.freqs(ii)));
        end
        % look only at the client age, +/- 5 years in either direction
        age = obj.subjAge;
        ageMin = 15;
        ageMax = 60;
        if age < ageMin
            warning('Subject age < available prior data. Using priors from ',num2str(ageMin),' years.');
        end
        if age > ageMax
            warning('Subject age > available prior data. Using priors from ',num2str(ageMax),' years.');
        end
        ageIndx = find(AGE>age-5 & AGE<age+5);

        nPoints = 5;
        cix = linspace(median(obj.minOutput),median(obj.maxOutput),nPoints)';
        PRIORS.cix = cix;
        PRIORS.nPoints = nPoints;
        for jj=1:nFreqs
            q = AUDIO(ageIndx,jj);
            nanindx = find(isnan(q));
            q(nanindx) = [];
            MU = median(q);
            IQR = iqr(q);
            PRIORS.mu(1,jj) = MU;
            PRIORS.sigma(1,jj) = IQR;
            
            pd = makedist('normal','mu',MU,'sigma',IQR);
            cd = pd.cdf(cix);
            PRIORS.points(:,jj) = cd(:); % prior probabilities at each point
        end
        obj.PRIORS = PRIORS; % this is the fixed starting point
        obj.pdPrior = PRIORS; % this will change as data are collected
        % starting level is round up to nearest 5 dB
        for ii=1:length(obj.pdPrior.mu)
            obj.startingLvls_L(1,ii) = obj.pdPrior.mu(ii) + 5-mod(obj.pdPrior.mu(ii),5); % set starting levels for testing
        end
        obj.startingLvls_R = obj.startingLvls_L;
    end
end
end

% OLD CODE ----------------------------------------------------------------

    % function [C1,C2,calPath1,calPath2] = getOutputCal(obj,calType,probeOutput)
    %     if isempty(probeOutput)
    %         C1 = [];
    %         C2 = [];
    %         calPath1.pathName = [];
    %         calPath1.folderName = [];
    %         calPath1.fileName = [];
    %         calPath2.pathName = [];
    %         calPath2.folderName = [];
    %         calPath2.fileName = [];
    %         return
    %     end
    %     if strcmp(calType,'thev')
    %         [pathName,folderName,fileName] = mostRecentCalibration('thev',probeOutput.label,1);
    %         calPath1.pathName = pathName;
    %         calPath1.folderName = folderName;
    %         calPath1.fileName = fileName;
    %         [pathName,folderName,fileName] = mostRecentCalibration('thev',probeOutput.label,2);
    %         calPath2.pathName = pathName;
    %         calPath2.folderName = folderName;
    %         calPath2.fileName = fileName;
    %     elseif strcmp(calType,'long')
    %         [pathName,folderName,fileName] = mostRecentCalibrationLT('longTube',probeOutput.label,1);
    %         calPath1.pathName = pathName;
    %         calPath1.folderName = folderName;
    %         calPath1.fileName = fileName;
    %         [pathName,folderName,fileName] = mostRecentCalibrationLT('longTube',probeOutput.label,2);
    %         calPath2.pathName = pathName;
    %         calPath2.folderName = folderName;
    %         calPath2.fileName = fileName;
    %     else
    %         error('Unrecognized calibration type. Must be thev or long.')
    %     end
    %     % load the calibration files
    %     dummy = load([calPath1.pathName,calPath1.folderName,calPath1.fileName]);
    %     if strcmp(calType,'thev')
    %         C1 = dummy.t;
    %     elseif strcmp(calType,'long')
    %         C1 = dummy.LTC;
    %     end
    %     dummy = load([calPath2.pathName,calPath2.folderName,calPath1.fileName]);
    %     if strcmp(calType,'thev')
    %         C2 = dummy.t;
    %     elseif strcmp(calType,'long')
    %         C2 = dummy.LTC;
    %     end
    % end

        % function [micCorrection] = getMicCorrection(obj,probeInput,doMicCorrection)
    %     if isempty(probeInput)
    %         micCorrection = [];
    %         return
    %     end
    %     [pathName_mic,folderName_mic,fileName_mic] = mostRecentCalibration('mic',probeInput.label,[]);
    %     if doMicCorrection == 1
    %         if ~isempty(fileName_mic)
    %             dummy = load([pathName_mic,folderName_mic,fileName_mic]);
    %             micCorrection = dummy.micCorrection; % microphone correction
    %         else
    %             micCorrection = 1; % convolving by this changes nothing
    %         end
    %     else
    %         micCorrection = 1; % convolving by this changes nothing
    %     end
    % end
