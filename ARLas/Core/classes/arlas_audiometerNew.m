classdef arlas_audiometerNew < handle
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Class definition arlas_audiometer
% For use with ARLas (Auditory Research Laboratory auditory software)
%
% Auditory Research Lab, The University of Iowa
% Deptartment of Communication Sciences & Disorders
% The University of Iowa
% Author: Shawn S. Goodman, PhD
% Date: July 26, 2021
% Last Updated: July 26, 2021
% Last Updated: July 30, 2021 -- ssg
% Last Updated: January 6-7, 2025 -- ssg
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
properties (SetAccess = private)
    version = '07JAN2025';
    sep % path delimiter appriate for the current operating system 
    map % struct containing file paths
    arlasObj % object that is the arlas object (not the audiometer object)
    fs % sampling rate
    b % lowpass filter coefficients
    
    
    %----- handles for graphical user interface ----- -----    
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
            thdMarks % audiometric threshold markers (x or o)
            thdDots % dotted lines connecting thresholds
            PIjunk % handles to plotted aspects of the pi function
            PIjunk2 % handles to text plotted on the pi function
            PRjunk % handles to plotted aspects of the progress tracker
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
    thd % row vector of minimum output levels (2 x nFreqs vector)
    
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
    %x0
    %trialCounter
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
       
    AAA % audiometer results structure. By ear and frequency
    minPresentations % number of stim presentations needed to get three "reversals" in standard Hewson-Westlake paradigm
    maxPresentations
    minReversals % minumim number of reversals to say threshold achieved
    
end
methods
    function obj = arlas_audiometerNew(arlasObj,probeL,probeR,freqs) % initialize object of class initARLas
        obj.sep = filesep; % get the path delimiter appropriate to the operating system being used            
        obj.arlasObj = arlasObj;
        obj.fs = arlasObj.fs; % get the system sampling rate
        obj.getLpf;
        obj.probeL = probeL;
        obj.probeR = probeR;
        if nargin >3
            freqs = sort(freqs(:));
        else
            freqs = [1 2 4 8 10 11.2 12.5 14 16]' * 1000; % default frequency set
        end
        obj.freqs = freqs(:)';
        obj.nFreqs = length(obj.freqs);
        obj.fIndx = 1; % initialize to lowest frequency
        
        obj.ear = 1; % initialize to left ear

        obj.lvlRefs = [{'FPL'},{'SPL'},{'IPL'}]; % initialize to fpl (forward pressure level)
        obj.stepSizes = [5,2]; % level step size options
        obj.minOutput = zeros(1,obj.nFreqs) -10;  % initialize to minus 10 dB
        obj.maxOutput = ones(1,obj.nFreqs) * 100;  % initialize to 100 dB

        obj.stepNow = obj.stepSizes(1);
        obj.levels = (obj.minOutput(obj.fIndx):obj.stepNow:obj.maxOutput(obj.fIndx))';
        obj.nLevels = length(obj.levels);
        [~,obj.lvlIndx] = min(abs(obj.levels - 80)); % initialize to 80 dB
        
        obj.autorunNow = 0; % initualize to manual control
        obj.autorespNow = 0; % initialize to manual response
        obj.earNow = 1;
        obj.freqNow = obj.freqs(obj.fIndx);
        obj.lvlNow = obj.levels(obj.lvlIndx);
        obj.refNow = obj.lvlRefs(1);
        obj.thdNow = [];
        
        obj.DONE = 0; % for total auto sequence
        obj.doneL = zeros(obj.nFreqs,1); % for total auto sequence
        obj.doneR = zeros(obj.nFreqs,1); % for total auto sequence
        obj.getStim;
        
        % AAA the audiometer results structure ---------------------------------------
        obj.AAA.freqs = obj.freqs;
        obj.AAA.nFreqs = obj.nFreqs;
        obj.AAA.levels = obj.levels;
        obj.AAA.DONE = 0;
        obj.doneL = zeros(obj.nFreqs,1);
        obj.doneR = zeros(obj.nFreqs,1);
        ear = 'L';
        for ii=1:obj.nFreqs
            obj.AAA.(ear).(['f_',num2str(obj.freqs(ii))]).thd = []; % current threshold from pi function
            obj.AAA.(ear).(['f_',num2str(obj.freqs(ii))]).thd2 = []; % current threshold from hewson westlake
            obj.AAA.(ear).(['f_',num2str(obj.freqs(ii))]).ThdEst = []; % vector of threshold estimates (from pi function) as a function of trial number

            obj.AAA.(ear).(['f_',num2str(obj.freqs(ii))]).RESP = [];  % subject response vector
            obj.AAA.(ear).(['f_',num2str(obj.freqs(ii))]).LVL = [];   % corresponding stimulus level vector
            
            obj.AAA.(ear).(['f_',num2str(obj.freqs(ii))]).trialCounter = 0; % total number of trials presented in the track
            obj.AAA.(ear).(['f_',num2str(obj.freqs(ii))]).descend1 = 1; % is in initial descent phase (for auto run)
            obj.AAA.(ear).(['f_',num2str(obj.freqs(ii))]).nTrials = 0; % total number of trials presented in the track
            obj.AAA.(ear).(['f_',num2str(obj.freqs(ii))]).nReversals = 0; % total number of reversals in the track
            obj.AAA.(ear).(['f_',num2str(obj.freqs(ii))]).u = []; % vector of unique levels presented
            obj.AAA.(ear).(['f_',num2str(obj.freqs(ii))]).v = []; % vector of percent percent "yes" responses at each unique level
            obj.AAA.(ear).(['f_',num2str(obj.freqs(ii))]).nn = []; % vector of number of presentations at each unique level
            
            obj.AAA.(ear).(['f_',num2str(obj.freqs(ii))]).x0 = [];    % pi fit x-offset
            obj.AAA.(ear).(['f_',num2str(obj.freqs(ii))]).k = [];     % pi fit slope
            obj.AAA.(ear).(['f_',num2str(obj.freqs(ii))]).rmse = [];  % pi fit root mean square error
            obj.AAA.(ear).(['f_',num2str(obj.freqs(ii))]).piw = []; % prediction interval width around threshold (dB)
            obj.AAA.(ear).(['f_',num2str(obj.freqs(ii))]).cix = []; % x-values (levels) for pi function and ci.
            obj.AAA.(ear).(['f_',num2str(obj.freqs(ii))]).pif = []; % most recent performance-intensity function
            obj.AAA.(ear).(['f_',num2str(obj.freqs(ii))]).ci = []; % most recent 95% confidence intervals around pif
        end
        ear = 'R';
        for ii=1:obj.nFreqs
            obj.AAA.(ear).(['f_',num2str(obj.freqs(ii))]).thd = []; % current threshold from pi function
            obj.AAA.(ear).(['f_',num2str(obj.freqs(ii))]).thd2 = []; % current threshold from hewson westlake
            obj.AAA.(ear).(['f_',num2str(obj.freqs(ii))]).ThdEst = []; % vector of threshold estimates (from pi function) as a function of trial number

            obj.AAA.(ear).(['f_',num2str(obj.freqs(ii))]).RESP = [];  % subject response vector
            obj.AAA.(ear).(['f_',num2str(obj.freqs(ii))]).LVL = [];   % corresponding stimulus level vector
            
            obj.AAA.(ear).(['f_',num2str(obj.freqs(ii))]).trialCounter = 1; % total number of trials presented in the track  
            obj.AAA.(ear).(['f_',num2str(obj.freqs(ii))]).descend1 = 1; % is in initial descent phase (for auto run)
            obj.AAA.(ear).(['f_',num2str(obj.freqs(ii))]).nTrials = 0; % total number of trials presented in the track
            obj.AAA.(ear).(['f_',num2str(obj.freqs(ii))]).nReversals = 0; % total number of reversals in the track
            obj.AAA.(ear).(['f_',num2str(obj.freqs(ii))]).u = []; % vector of unique levels presented
            obj.AAA.(ear).(['f_',num2str(obj.freqs(ii))]).v = []; % vector of percent percent "yes" responses at each unique level
            obj.AAA.(ear).(['f_',num2str(obj.freqs(ii))]).nn = []; % vector of number of presentations at each unique level
            
            obj.AAA.(ear).(['f_',num2str(obj.freqs(ii))]).x0 = [];    % pi fit x-offset
            obj.AAA.(ear).(['f_',num2str(obj.freqs(ii))]).k = [];     % pi fit slope
            obj.AAA.(ear).(['f_',num2str(obj.freqs(ii))]).rmse = [];  % pi fit root mean square error
            obj.AAA.(ear).(['f_',num2str(obj.freqs(ii))]).piw = []; % prediction interval width around threshold (dB)
            obj.AAA.(ear).(['f_',num2str(obj.freqs(ii))]).cix = []; % x-values (levels) for pi function and ci.
            obj.AAA.(ear).(['f_',num2str(obj.freqs(ii))]).pif = []; % most recent performance-intensity function
            obj.AAA.(ear).(['f_',num2str(obj.freqs(ii))]).ci = []; % most recent 95% confidence intervals around pif
        end
        % -----------------------------------------------------------------
        
        obj.minPresentations = 8;
        obj.maxPresentations = 30;
        obj.minReversals = 3;
        obj.doingSwitch = 1;
        obj.viewNow = 1; % initialize to pi / progress view
        obj.initGui % initialize the gui
    end  
    function abort(varargin) % instructions for aborting when gui closed
        try
            obj = varargin{1};
        catch
        end
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
                'String',{'FPL','SPL','IPL'},'FontSize',12,'BackgroundColor',...
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
                                    keyboard
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
            obj.prepStim
            obj.applyISC
            obj.playStim
            obj.getResponse
            %obj.analyze
            % if obj.h_switchPlots.Value == 0 % pi and tracking view --
            %     obj.audiometer_plotPI
            % elseif obj.h_switchPlots.Value == 1 % audiogram view --
            %     obj.audiometer_plotAG
            % end            
            obj.checkIsFinished            
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

        obj.fIndx = obj.fIndx - 1;
        if obj.fIndx < 1
            obj.fIndx = obj.nFreqs;
        end
        obj.freqNow = obj.freqs(obj.fIndx);
        obj.FREQ.String = [num2str(obj.freqNow),' Hz'];
        
        if obj.h_switchPlots.Value == 1 % pi and tracking view --
            obj.audiometer_plotPI
            % h1 = obj.h_pif; % performance intensity function
            % axes(h1); % plot the performance-intensity function -----
            % try delete(obj.marker1); catch; end
            % ymin = -0.1;
            % obj.marker1 = plot(obj.lvlNow,ymin,'k^','MarkerFaceColor',[0 0 0]);
            % h2 = obj.h_track; % response tracker        
            % axes(h2) % plot the response track -----------------
            % try delete(obj.marker2); catch; end
            % if obj.earNow == 1
            %     ear = 'L';
            % elseif obj.earNow == 2
            %     ear = 'R';
            % end
            % A = obj.AAA.(ear).(['f_',num2str(obj.freqNow)]);            
            % RESP = A.RESP;
            % N = length(RESP);
            % obj.marker2 = plot(N+1,obj.lvlNow,'k^','MarkerFaceColor',[0 0 0]);  
        elseif obj.h_switchPlots.Value == 0 % audiogram view --
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
        
        if obj.h_switchPlots.Value == 1 % pi and tracking view --
            obj.audiometer_plotPI
            % h1 = obj.h_pif; % performance intensity function
            % axes(h1); % plot the performance-intensity function -----
            % try delete(obj.marker1); catch; end
            % ymin = -0.1;
            % obj.marker1 = plot(obj.lvlNow,ymin,'k^','MarkerFaceColor',[0 0 0]);
            % h2 = obj.h_track; % response tracker        
            % axes(h2) % plot the response track -----------------
            % try delete(obj.marker2); catch; end
            % if obj.earNow == 1
            %     ear = 'L';
            % elseif obj.earNow == 2
            %     ear = 'R';
            % end
            % A = obj.AAA.(ear).(['f_',num2str(obj.freqNow)]);            
            % RESP = A.RESP;
            % N = length(RESP);
            % obj.marker2 = plot(N+1,obj.lvlNow,'k^','MarkerFaceColor',[0 0 0]);  
        elseif obj.h_switchPlots.Value == 0 % audiogram view --
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
        
        % update the plots
        if obj.h_switchPlots.Value == 0 % pi and tracking view --
            %obj.audiometer_plotPI
            h1 = obj.h_pif; % performance intensity function
            axes(h1); % plot the performance-intensity function -----
            try delete(obj.marker1); catch; end
            ymin = -0.1;
            obj.marker1 = plot(obj.lvlNow,ymin,'k^','MarkerFaceColor',[0 0 0]);
            h2 = obj.h_track; % response tracker        
            axes(h2) % plot the response track -----------------
            try delete(obj.marker2); catch; end
            if obj.earNow == 1
                ear = 'L';
            elseif obj.earNow == 2
                ear = 'R';
            end
            A = obj.AAA.(ear).(['f_',num2str(obj.freqNow)]);            
            RESP = A.RESP;
            N = length(RESP);
            obj.marker2 = plot(N+1,obj.lvlNow,'k^','MarkerFaceColor',[0 0 0]);  
        elseif obj.h_switchPlots.Value == 1 % audiogram view --
            obj.audiometer_plotAG
        end
        
        set(obj.h_lvlDn,'Value',0,'Tag','on','CData',imread('lvlDnYellow.jpg'))
        pause(0.1)
    end    
    function lvlUpManager(varargin) % control which input channel is currenly being viewed
        obj = varargin{1};
        if strcmp(obj.h_lvlUp.Tag,'off')
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
        if obj.h_switchPlots.Value == 1 % pi and tracking view --
            obj.audiometer_plotPI
        elseif obj.h_switchPlots.Value == 0 % audiogram view --
            obj.audiometer_plotAG
        end
        
        % xmin = obj.minOutput(obj.fIndx);
        % xmax = obj.maxOutput(obj.fIndx);
        % h1 = obj.h_pif; % performance intensity function
        % axes(h1)
        % try delete(obj.marker1)
        % catch
        % end
        % obj.marker1 = plot(obj.lvlNow,0,'k^','MarkerFaceColor',[0 0 0]);
        % ylim([-0.1,1.05]);
        % xlim([xmin,xmax]);
        % xlabel('Stim Level (dB FPL)')
        % ylabel('Response Probability')
        % grid on
        % 
        % h2 = obj.h_track; % response tracker   
        % axes(h2)
        % if obj.earNow == 1
        %     ear = 'L';
        % elseif obj.earNow == 2
        %     ear = 'R';
        % end
        % A = obj.AAA.(ear).(['f_',num2str(obj.freqNow)]);
        % RESP = A.RESP;
        % N = length(RESP);
        % try delete(obj.marker2)
        % catch
        % end
        % obj.marker2 = plot(N+1,obj.lvlNow,'k^','MarkerFaceColor',[0 0 0]);
        % xlabel('Trial Number')
        % ylabel('Stim Level')
        % if N < 6
        %     xlim([1 6])
        %     xticks((1:1:6))
        % else
        %     xlim([1 N+1])
        %     xticks((1:1:N))
        % end
        % -----------------------------------------------------------------
        
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
        if obj.h_switchPlots.Value == 1 % pi and tracking view --
            obj.audiometer_plotPI
        elseif obj.h_switchPlots.Value == 0 % audiogram view --
            obj.doingSwitch = 1;
            obj.audiometer_plotAG
            obj.doingSwitch = 0;
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
        if obj.h_switchPlots.Value == 0 % pi and tracking view --
            obj.audiometer_plotPI
        elseif obj.h_switchPlots.Value == 1 % audiogram view --
            obj.audiometer_plotAG
        end
        
        % with location of new presentation level --------
        % xmin = obj.minOutput(obj.fIndx);
        % xmax = obj.maxOutput(obj.fIndx);
        % h1 = obj.h_pif; % performance intensity function
        % axes(h1)
        % try delete(obj.marker1)
        % catch
        % end
        % obj.marker1 = plot(obj.lvlNow,0,'k^','MarkerFaceColor',[0 0 0]);
        % ylim([-0.1,1.05]);
        % xlim([xmin,xmax]);
        % xlabel('Stim Level (dB FPL)')
        % ylabel('Response Probability')
        % grid on
        % 
        % h2 = obj.h_track; % response tracker   
        % axes(h2)
        % if obj.earNow == 1
        %     ear = 'L';
        % elseif obj.earNow == 2
        %     ear = 'R';
        % end
        % A = obj.AAA.(ear).(['f_',num2str(obj.freqNow)]);
        % RESP = A.RESP;
        % N = length(RESP);
        % try delete(obj.marker2)
        % catch
        % end
        % obj.marker2 = plot(N+1,obj.lvlNow,'k^','MarkerFaceColor',[0 0 0]);
        % xlabel('Trial Number')
        % ylabel('Stim Level')
        % if N < 6
        %     xlim([1 6])
        %     xticks((1:1:6))
        % else
        %     xlim([1 N+1])
        %     xticks((1:1:N))
        % end
        % -----------------------------------------------------------------
    end    
    function refManager(varargin) % control level reference
        obj = varargin{1};
        if obj.REF.Value == 1
            obj.refNow = obj.lvlRefs(1);
        elseif obj.REF.Value == 2
            obj.refNow = obj.lvlRefs(2);
        elseif obj.REF.Value == 3
            obj.refNow = obj.lvlRefs(3);
        end
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
        obj.buttonManager(10)
        set(obj.h_cal,'Value',0,'Tag','off','BackgroundColor',[0 1 0])
        
        ao = obj.arlasObj;

        % obj.doMicCorrection = 1;
        % % get most recent calibration files -----
        % micCorrectionL = obj.getMicCorrection(obj.probeInputL,obj.doMicCorrection);
        % micCorrectionR = obj.getMicCorrection(obj.probeInputR,obj.doMicCorrection);
        % [C1L,C2L,calPath1L,calPath2L] = obj.getOutputCal('thev',obj.probeOutputL);
        % [C1R,C2R,calPath1R,calPath2R] = obj.getOutputCal('thev',obj.probeOutputR);
        % 
        % % PERFORM IN-SITU CALIBRATION -----
        % disp('----- Running in-situ calibration -----')
        % inSituReps = 6;
        % doIndividual = 1;
        % doSimultaneous = 0;
        % calType = 'thev';
        % fmin = obj.fmin;
        % fmax = obj.fmax;
        % 
        % if ~isempty(obj.probeInputL)
        %     [obj.iscS1L,obj.iscS2L,iscS12L] = ARLas_runISC(ao,obj.probeInputL,obj.probeOutputL,calPath1L,calPath2L,calType,fmin,fmax,micCorrectionL,inSituReps,doIndividual,doSimultaneous);
        % else
        %     obj.iscS1L = [];
        %     obj.iscS2L = [];
        %     %iscS12L = [];
        % end
        % if ~isempty(obj.probeInputR)
        %     [obj.iscS1R,obj.iscS2R,iscS12R] = ARLas_runISC(ao,obj.probeInputR,obj.probeOutputR,calPath1R,calPath2R,calType,fmin,fmax,micCorrectionR,inSituReps,doIndividual,doSimultaneous);
        % else
        %     obj.iscS1R = [];
        %     obj.iscS2R = [];
        %     %iscS12R = [];
        % end
    

        probeL = obj.probeL;
        probeR = obj.probeR;
        doMicCorrection = 1;

        % specify calibration to use to set stimulus levels
        calType = 'thev'; % use 'thev' for Thevenin source
        targetCalType = 'fpl'; % 'fpl', 'spl', 'ipl'.
        ISCL = struct;
        ISCR = struct;
        % If presenting X dB FPL, this is 3-6 dB higher in dB SPL.
        % In order to make more comparable with other studies that calibrate in
        % dB SPL, can adjust the level accordingly:
        % if strcmp(calType,'thev') 
        %     if strcmp(targetCalType,'fpl')
        %         adjust = -3; % if 55 dB fpl, this is more like 58 dB SPL. Therefore, subtract 3 to make the level 52 dB FPL ~= 55 dB SPL.
        %     else
        %         adjust = 0;
        %     end
        % end
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
        fmax = 20000;
        disp('----- Running in-situ calibration -----')
        inSituReps = 6;
        doIndividual = 1;
        doSimultaneous = 0;
        [iscS1L,iscS2L,iscS12L] = ARLas_runISC(obj.arlasObj,probeInputL,probeOutputL,calPath1L,calPath2L,calType,fmin,fmax,micCorrectionL,inSituReps,doIndividual,doSimultaneous);
        [iscS1R,iscS2R,iscS12R] = ARLas_runISC(obj.arlasObj,probeInputR,probeOutputR,calPath1R,calPath2R,calType,fmin,fmax,micCorrectionR,inSituReps,doIndividual,doSimultaneous);
        
        obj.iscS1L = iscS1L;
        obj.iscS2L = iscS2L;
        obj.iscS1R = iscS1R;
        obj.iscS2R = iscS2R;
        obj.probeInputL = probeInputL;
        obj.probeInputR = probeInputR;
        obj.probeOutputL = probeOutputL;
        obj.probeOutputR = probeOutputR;
        [inputs,outputs] = hardwareSetup; % read in the hardware setup
        N = length(inputs);
        done = 0;
        counter = 1;
        while done == 0
            label = inputs{counter}.label;
            if strcmp(label,'clicker')
                done = 1;
                obj.clickerInput = inputs{counter};
            else
                counter = counter + 1;
            end
            if counter > N
                done = 1;
                error('Unable to locate clicker. Make sure there is a clicker label in hardware setup.')
            end
        end

        disp('----- Finished in-situ calibration -----')

        if isempty(obj.probeInputL) & isempty(obj.probeInputR)
            warning('Error: both Left and Right probe inputs are empty.')
        end
        obj.buttonManager(20)

    end
    function doSave(varargin) % save the results to excel spreadsheet
        obj = varargin{1};
        if strcmp(obj.h_save.Tag,'off')
            obj.h_save.Value = 0;
            return
        end
        obj.buttonManager(10)
        set(obj.h_cal,'Value',0,'Tag','off','BackgroundColor',[.7 .7 .7])
        set(obj.h_save,'Value',0,'Tag','off','BackgroundColor',[0 1 0])
        
        
        %disp('----- Saving audiogram files -----')
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
        
        try 
            ao = obj.arlasObj;
            savePath = [ao.savingGrace,ao.experimentID,'_analysis',ao.sep];
            if exist(savePath,'dir') == 0 % if path does not exist
                success = mkdir(savePath);
                if ~success
                    warning('Unable to create new Experiment Directory: data save path')
                end
            end 
            saveFileName = [ao.subjectID,'_analyzedAUDIO.mat'];
            saveFileName = ARLas_saveName(savePath,saveFileName);
            AUDIO = obj.AAA;
            AUDIO.timeStamp = cellstr(datetime('now'));
            AUDIO.subjectID = ao.subjectID;
            AUDIO.ThdL = ThdL;
            AUDIO.ThdR = ThdR;
            
            save([savePath,saveFileName],'AUDIO')
        catch
            disp('Warning: Audiometric Analysis not saved!')
        end
        
        % write thresholds to excel file here ----------
        %disp('----- Writing Excel files -----')
        saveFileName = [ao.subjectID,'_analyzedAUDIO.xlsx'];
        %saveFileName = ARLas_saveName(savePath,saveFileName);
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
        % when you push the button to the down position, Value is 1
        % we use 1 for the audiogram view
        if obj.h_switchPlots.Value == 1 % switch to pi and tracking view -------------
            set(obj.h_switchPlots,'Tag','on','BackgroundColor',[0 1 0])
            pause(0.01)
            % turn off the audiogram stuff
            axes(obj.h_ag) % audiogram ; NOTE: obj.h_ag.Visible = 'off' doesnt work well
            axis off
            try delete(obj.L1); catch; end
            try delete(obj.L2); catch; end
            try delete(obj.markerX); catch; end
            try delete(obj.markerY); catch; end
            try delete(obj.LX); catch; end
            try delete(obj.LY); catch; end

            axes(obj.h_pif) % performance intensity function
            axis on            
            %Aobj.h_pif.Visible = 'on';
            axes(obj.h_track) % progress tracking
            axis on
            %obj.h_track.Visible = 'on';
            set(obj.h_switchPlots,'Tag','on','BackgroundColor',[1 1 0])
            obj.audiometer_plotPI
        elseif obj.h_switchPlots.Value == 0 % switch to audiogram view ------------------------------------------------
            set(obj.h_switchPlots,'Tag','on','BackgroundColor',[0 1 0])
            pause(0.01)
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
            try delete(obj.PRjunk{4}); catch; end
            try delete(obj.PRjunk{3}); catch; end
            try delete(obj.PRjunk{2}); catch; end
            try delete(obj.PRjunk{1}); catch; end

            axes(obj.h_pif) % performance intensity function -----
            axis off            
            axes(obj.h_track) % tracking progress
            axis off

            axes(obj.h_ag) % audiogram -----
            axis on
            %obj.h_ag.Visible = 'on';
            set(obj.h_switchPlots,'Tag','on','BackgroundColor',[1 1 0])
            obj.audiometer_plotAG
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
        RESP = obj.AAA.(ear).(['f_',num2str(obj.freqNow)]).RESP;
        LVL = obj.AAA.(ear).(['f_',num2str(obj.freqNow)]).LVL;
        LVL(trialCounter,1) = obj.lvlNow;
        RESP(trialCounter,1) = 1;
        obj.AAA.(ear).(['f_',num2str(obj.freqNow)]).RESP = RESP;
        obj.AAA.(ear).(['f_',num2str(obj.freqNow)]).LVL = LVL;
        obj.AAA.(ear).(['f_',num2str(obj.freqNow)]).trialCounter = trialCounter;

        obj.analyze
        if obj.h_switchPlots.Value == 0 % pi and tracking view --
            obj.audiometer_plotPI
        elseif obj.h_switchPlots.Value == 1 % audiogram view --
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
        RESP = obj.AAA.(ear).(['f_',num2str(obj.freqNow)]).RESP;
        LVL = obj.AAA.(ear).(['f_',num2str(obj.freqNow)]).LVL;
        LVL(trialCounter,1) = obj.lvlNow;
        RESP(trialCounter,1) = 0;
        obj.AAA.(ear).(['f_',num2str(obj.freqNow)]).RESP = RESP;
        obj.AAA.(ear).(['f_',num2str(obj.freqNow)]).LVL = LVL;
        obj.AAA.(ear).(['f_',num2str(obj.freqNow)]).trialCounter = trialCounter;

        obj.analyze
        if obj.h_switchPlots.Value == 0 % pi and tracking view --
            obj.audiometer_plotPI
        elseif obj.h_switchPlots.Value == 1 % audiogram view --
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
    function getLpf(varargin) % lowpass filter coefficients
        obj = varargin{1};
        Fpass = 10;
        Fstop = 1.7 * Fpass;    % Stopband Frequency
        Dpass = 0.057501127785;  % Passband Ripple
        Dstop = 0.0001;          % Stopband Attenuation
        flag  = 'scale';         % Sampling Flag
        [N,Wn,BETA,TYPE] = kaiserord([Fpass Fstop]/(obj.fs/2), [1 0], [Dstop Dpass]);
        if mod(N,2)~= 0
            N = N + 1;
        end
        obj.b  = fir1(N, Wn, TYPE, kaiser(N+1, BETA), flag)';
    end
    function applyISC(varargin) % Apply the in-situ calibration to the stimuli
        obj = varargin{1};
        if obj.earNow == 1 % left ear
            try
                obj.iscNow = obj.iscS1L;
            catch
                obj.iscNow = obj.iscS2L;
            end
        elseif obj.earNow == 2 % right ear
            try
                obj.iscNow = obj.iscS1R;
            catch
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
        
        [s1,clicksImpulse,errors1] = ARLas_applyISC(obj.iscNow,obj.lvlNow+3,targetCalType,obj.freqNow,stim(:));
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
        ao.clearPlayList % clear out whatever was used previously
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
        ao.setRecList(obj.clickerInput.ch,obj.clickerInput.label,obj.clickerInput.micSens,obj.clickerInput.gain);
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
        ao.objPlayrec.userInfo.f = obj.freqNow;
        ao.objPlayrec.userInfo.fs = obj.fs;
        ao.objPlayrec.userInfo.audiometer_version = obj.version;
        ao.objPlayrec.userInfo.calType = obj.refNow;
        ao.objPlayrec.userInfo.C1L = obj.C1L;
        ao.objPlayrec.userInfo.C2L = obj.C2L;
        ao.objPlayrec.userInfo.C1R = obj.C1R;
        ao.objPlayrec.userInfo.C2R = obj.C2R;
        ao.objPlayrec.userInfo.iscS1 = obj.iscS1L;
        ao.objPlayrec.userInfo.iscS2 = obj.iscS2L;
        ao.objPlayrec.userInfo.iscS1 = obj.iscS1R;
        ao.objPlayrec.userInfo.iscS2 = obj.iscS2R;
        ao.objPlayrec.userInfo.probeR = obj.probeR;
        ao.objPlayrec.userInfo.probeL = obj.probeL;
        if obj.autorunNow % if running on autopilot, 
            pauseLen = rand(1,1)*1;
            pause(pauseLen) % add some randomness to the presentation delay
        end
        ao.objPlayrec.servedNeat = 1;
        ao.objPlayrec.run % playback and record -----
        if ao.killRun
            obj.DONE = 1;
            obj.AAA.DONE = 1;
            return
        end
    end
    function getResponse(varargin) % extract the recordings to see if button was pressed (tone was heard)
        obj = varargin{1};
        ao = obj.arlasObj;
        if obj.autorespNow == 1
            [header,Data] = ao.retrieveData(['Ch',num2str(obj.clickerInput.ch)]); % retrive recorded data
            q = fastFilter(obj.b,Data(:));
            [indx,nRejects] = AR(q,'moderate',0);
            r = nRejects / length(q);
    
            if obj.earNow == 1
                ear = 'L';
            elseif obj.earNow == 2
                ear = 'R';
            end
            trialCounter = obj.AAA.(ear).(['f_',num2str(obj.freqNow)]).trialCounter+1; % increment counter here
            RESP = obj.AAA.(ear).(['f_',num2str(obj.freqNow)]).RESP;
            LVL = obj.AAA.(ear).(['f_',num2str(obj.freqNow)]).LVL;
            LVL(trialCounter,1) = obj.lvlNow;
    
            if r > 0.0005
                RESP(trialCounter,1) = 1;
            else
                RESP(trialCounter,1) = 0;
            end
            
            if obj.autorunNow == 1 % check to see if still descending (auto run only)
                if obj.AAA.(ear).(['f_',num2str(obj.freqNow)]).descend1 == 1
                    if RESP(trialCounter,1) == 0 % initial decent is ended
                        obj.AAA.(ear).(['f_',num2str(obj.freqNow)]).descend1 = 0;
                    end
                end
            end
            obj.AAA.(ear).(['f_',num2str(obj.freqNow)]).RESP = RESP;
            obj.AAA.(ear).(['f_',num2str(obj.freqNow)]).LVL = LVL;
            obj.AAA.(ear).(['f_',num2str(obj.freqNow)]).trialCounter = trialCounter;
        else
            set(obj.h_respYES,'Value',0,'Tag','on','BackgroundColor',[0 1 0])
            set(obj.h_respNO,'Value',0,'Tag','on','BackgroundColor',[1 0 0])
        end
    end
    function analyze(varargin) % analyze and get next level (if in auto mode)
        obj = varargin{1};

        if obj.earNow == 1
            ear = 'L';
        elseif obj.earNow == 2
            ear = 'R';
        end
        RESP = obj.AAA.(ear).(['f_',num2str(obj.freqNow)]).RESP;
        LVL = obj.AAA.(ear).(['f_',num2str(obj.freqNow)]).LVL;
        x0 = obj.AAA.(ear).(['f_',num2str(obj.freqNow)]).x0;
        k = obj.AAA.(ear).(['f_',num2str(obj.freqNow)]).k;
        rmse = obj.AAA.(ear).(['f_',num2str(obj.freqNow)]).rmse;
        minout = obj.minOutput(obj.fIndx);
        maxout = obj.maxOutput(obj.fIndx);
        trialCounter = obj.AAA.(ear).(['f_',num2str(obj.freqNow)]).trialCounter;
        descend1 = obj.AAA.(ear).(['f_',num2str(obj.freqNow)]).descend1;
        
        % audiometer performance-intensity function
        [x0(trialCounter,1),k(trialCounter,1),ci,cix,rmse(trialCounter,1)] = audiometer_piFit(LVL,RESP,minout,maxout);
        
        % calculate new test stimulus lvl
        if trialCounter >= 5
            obj.countReversals
        end
        
        if obj.autorunNow == 1 % check to see if still descending (auto run only)
            if obj.descend1 ==1 % initial descent has larger steps
                stepSize = 10;
            else
                stepSize = 5;
                obj.STEP.Value = 1;
                obj.stepManager                
            end
        else
            stepSize = 5;
            obj.STEP.Value = 1;
            obj.stepManager
        end
        
        if obj.autorunNow == 1 % automatically calculate new level  (auto run only)        
            xo = LVL(end);
            remainder = mod(xo,stepSize); % force to step size multiples
            if remainder<(stepSize/2)
                xo = xo - remainder;
            else
                xo = xo + (stepSize - remainder);
            end            
            if RESP(trialCounter,1) == 1
                %obj.lvlNow = xo -(2*stepSize);
                obj.lvlDnManager
                pause(0.01)
                obj.lvlDnManager
            else
                %obj.lvlNow = xo + (1*stepSize);
                obj.lvlUpManager
            end
            %if obj.lvlNow < minout
            %    obj.lvlNow = minout;
            %end
            %if obj.lvlNow > maxout
            %    obj.lvlNow = maxout;
            %end
            %keyboard
            
        end
        
        u = unique(LVL); % average all responses from same levels
        for jj=1:length(u)
           indxx = find(LVL==u(jj));
           resp = RESP(indxx);
           v(jj,1) = mean(resp); % averaged response rate
           nn(jj,1) = length(indxx);
        end
        
        % add logic to end if no response        
        %         indxMax = find(u==maxout);
        %         if~isempty(indxMax)
        %             keyboard
        %         end
        %         indxMin = find(u==minout);
        %         if~isempty(indxMin)
        %             keyboard
        %         end
        
        pif = 1./(1+exp(-(k(trialCounter)*(cix-x0(trialCounter))))); % performance intensity function
        [~,indx1] = min(abs(0.5 - ci(:,1)));
        [~,indx2] = min(abs(0.5 - ci(:,2)));
        piw(trialCounter,1) = cix(indx1) - cix(indx2); % 95% prediction interval width at threshold (dB)        
        
        obj.AAA.(ear).(['f_',num2str(obj.freqNow)]).RESP = RESP;
        obj.AAA.(ear).(['f_',num2str(obj.freqNow)]).LVL = LVL;
        obj.AAA.(ear).(['f_',num2str(obj.freqNow)]).trialCounter = trialCounter;
        obj.AAA.(ear).(['f_',num2str(obj.freqNow)]).descend1 = descend1;
        
        obj.AAA.(ear).(['f_',num2str(obj.freqNow)]).thd = x0(end);   % calculated threshold (dB FPL)
        obj.AAA.(ear).(['f_',num2str(obj.freqNow)]).piw = piw; % prediction interval width around threshold (dB)
        obj.AAA.(ear).(['f_',num2str(obj.freqNow)]).cix = cix; % x-values (levels) for pi function and ci.
        obj.AAA.(ear).(['f_',num2str(obj.freqNow)]).pif = pif; % most recent performance-intensity function
        obj.AAA.(ear).(['f_',num2str(obj.freqNow)]).ci = ci; % most recent 95% confidence intervals around pif
        obj.AAA.(ear).(['f_',num2str(obj.freqNow)]).nTrials = trialCounter; % total number of trials presented
        obj.AAA.(ear).(['f_',num2str(obj.freqNow)]).u = u; % vector of unique levels presented
        obj.AAA.(ear).(['f_',num2str(obj.freqNow)]).v = v; % vector of percent percent "yes" responses at each unique level
        obj.AAA.(ear).(['f_',num2str(obj.freqNow)]).nn = nn; % vector of number of presentations at each unique level
        obj.AAA.(ear).(['f_',num2str(obj.freqNow)]).x0 = x0; % vector of threshold estimates as a function of trial number
        obj.AAA.(ear).(['f_',num2str(obj.freqNow)]).k = k;   % vector of slope estimates as a function of trial number
        obj.AAA.(ear).(['f_',num2str(obj.freqNow)]).rmse = rmse; % root mean square error; estimate of goodness of logistic fit
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
    
    function audiometer_plotPI(varargin) % plot audiometer PI function
        obj = varargin{1};

        xmin = obj.minOutput(obj.fIndx);
        xmax = obj.maxOutput(obj.fIndx);
        if obj.earNow == 1
            ear = 'L';
        elseif obj.earNow == 2
            ear = 'R';
        end
        A = obj.AAA.(ear).(['f_',num2str(obj.freqNow)]);
        
        thd = A.thd;   % calculated threshold (dB FPL)
        piw = A.piw; % prediction interval width around threshold (dB)
        cix = A.cix; % x-values (levels) for pi function and ci.
        pif = A.pif; % most recent performance-intensity function
        ci = A.ci; % most recent 95% confidence intervals around pif    
        LVL = A.LVL;   % vector of presentation levels as a function of trial number
        RESP = A.RESP; % vector of responses as a function of trial number
        u = A.u; % vector of unique levels presented
        v = A.v; % vector of percent percent "yes" responses at each unique level
        nn = A.nn; % vector of number of presentations at each unique level

        % plot performance-intensity function ------------------
        h1 = obj.h_pif; % performance intensity function
        h2 = obj.h_track; % response tracker        
        
        ymin = -0.1;
        ymax = 1.05; 
        N = length(RESP);
        
        axes(h1); % plot the performance-intensity function -----
        
        try delete(obj.PIjunk{7}); catch; end
        try delete(obj.PIjunk{6}); catch; end
        try delete(obj.PIjunk{5}); catch; end
        try delete(obj.PIjunk{4}); catch; end
        try delete(obj.PIjunk{3}); catch; end
        try delete(obj.PIjunk{2}); catch; end
        try delete(obj.PIjunk{1}); catch; end        
        try delete(obj.PIjunk2); catch; end 
        if isempty(LVL)
            hold off
            obj.PIjunk{1} = plot(0,0,'.w'); % plot individual responses
            hold on
        else
            hold off
            obj.PIjunk{1} = plot(LVL,RESP+0.01*randn(size(LVL)),'.b'); % plot individual responses
            hold on
            obj.PIjunk{2} = plot(u,v,'ob'); % plot averaged response rates
            if N > 6
                obj.PIjunk{3} = plot(cix,ci(:,1),'g--');
                obj.PIjunk{4} = plot(cix,ci(:,2),'g--');
                obj.PIjunk{5} = line([thd-(piw(end)/2) thd+(piw(end)/2)],[0.5 0.5],'Color',[0 1 0],'LineWidth',2,'LineStyle','-');
            end
            obj.PIjunk{6} = plot(cix,pif,'r'); % plot logistic fit
            try
                obj.PIjunk{7} = plot(thd,0.5,'or'); % plot current threshold
            catch ME
                return
            end
            for kk=1:length(u) % show how many presentations at each level
                obj.PIjunk2{kk} = text(u(kk),-0.05,num2str(nn(kk)), 'Rotation', 0, 'VerticalAlignment', 'middle','FontSize',8); 
            end
            title(['Threshold: ',num2str(round(thd)),' dB FPL +/- ',num2str(piw(end)),' dB'])
        end
        xlabel('Stim Level (dB FPL)')
        ylabel('Response Probability')
        grid on
        pause(0.01)

        try delete(obj.marker1); catch; end
        obj.marker1 = plot(obj.lvlNow,ymin,'k^','MarkerFaceColor',[0 0 0]);
        if N < 6
            xlim([xmin xmax])
            ylim([ymin,ymax])
        else
            xmin2 = thd - 20;
            xmax2 = thd + 20;
            if xmin < xmin
                xmin2 = xmin;
            end
            if xmax2 > xmax
                xmax2 = xmax;
            end
            xlim([xmin2 xmax2])
            ylim([ymin,ymax])
        end

        axes(h2) % plot the response track -----------------
        try delete(obj.PRjunk{4}); catch; end
        try delete(obj.PRjunk{3}); catch; end
        try delete(obj.PRjunk{2}); catch; end
        try delete(obj.PRjunk{1}); catch; end

        if isempty(LVL)
            hold off
            obj.PRjunk{1} = plot(1,0,'.w');
            hold on
        else
            x = (1:N)';
            indx1 = find(RESP == 1);
            indx0 = find(RESP == 0);
            hold off
            obj.PRjunk{1} = plot(x(indx1),LVL(indx1),'go','LineWidth',2);
            hold on
            obj.PRjunk{2} = plot(x(indx0),LVL(indx0),'rx','LineWidth',2);
            obj.PRjunk{3} = plot(x,LVL,'.-b');
        end
        xlabel('Trial Number')
        ylabel('Stim Level')
        if N < 6
            xlim([1 6])
            xticks((1:1:6))
        else
            xlim([1 N+1])
            xticks((1:1:N))
        end
        try delete(obj.marker2); catch; end
        obj.marker2 = plot(N+1,obj.lvlNow,'k^','MarkerFaceColor',[0 0 0]);
        if N > 6
            obj.PRjunk{4} = line([1 N+1],[thd thd],'Color',[0 0 0],'LineWidth',0.5,'LineStyle','--');
        end
        grid on
    end
    function audiometer_plotAG(varargin) % plot audiometer Audiogram
        obj = varargin{1};
        
        FF = [250,500,1000,2000,4000,8000,16000,20000]';
        LL = (-10:10:110)';
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
            dummy = obj.AAA.L.(['f_',num2str(obj.freqs(ii))]).thd;
            if isempty(dummy)
                thdL(ii,1) = NaN;
            else
                thdL(ii,1) = dummy;
            end
        end
        for ii=1:obj.nFreqs
            dummy = obj.AAA.R.(['f_',num2str(obj.freqs(ii))]).thd;
            if isempty(dummy)
                thdR(ii,1) = NaN;
            else
                thdR(ii,1) = dummy;
            end
        end 
        
        if obj.doingSwitch == 1 % if switching views, redraw the level and frequency lines   
            hold off
            plot(0,0,'w.')
            for ii=1:nf
                axes(h1)
                obj.L1(ii) = line([ff(ii),ff(ii)],[ymin,ymax],'Color',[0 0 0],'LineWidth',0.5,'LineStyle','-');
                hold on
            end
            for ii=1:length(LL)
                axes(h1)
                obj.L2(ii) = line([xmin,xmax],[LL(ii),LL(ii)],'Color',[0 0 0],'LineWidth',0.5,'LineStyle','-');
            end
        end
        ylim([ymin,ymax])
        xlim([xmin,xmax])
        try delete(obj.markerX); 
        catch
        end
        try delete(obj.markerY); 
        catch 
        end
        obj.markerX = plot(log2(obj.freqNow),ymin,'k^','MarkerFaceColor',[0 0 0]);
        obj.markerY = plot(xmin,obj.lvlNow,'k^','MarkerFaceColor',[0 0 0]);
        try delete(obj.LX); 
        catch
        end
        try delete(obj.LY); 
        catch 
        end
        if obj.earNow == 1 % ear is left
            obj.LX = line([log2(obj.freqNow),log2(obj.freqNow)],[ymin,ymax],'Color',[0 0 .5],'LineWidth',1,'LineStyle','-');
            obj.LY = line([xmin,xmax],[obj.lvlNow,obj.lvlNow],'Color',[0 0 .5],'LineWidth',1,'LineStyle','-');
        elseif obj.earNow == 2 % ear is right
            obj.LX = line([log2(obj.freqNow),log2(obj.freqNow)],[ymin,ymax],'Color',[.5 0 0],'LineWidth',1,'LineStyle','-');
            obj.LY = line([xmin,xmax],[obj.lvlNow,obj.lvlNow],'Color',[.5 0 0],'LineWidth',1,'LineStyle','-');
        end
        if all(isnan(thdL)) && all(isnan(thdR)) % there is nothing else to plot
        else
            try delete(obj.thdMarks)
            catch
            end
            try delete(obj.thdDots)
            catch
            end
            hold on
            freqs = log2(obj.freqs);
            if obj.earNow == 1
                for ii=1:obj.nFreqs
                    obj.thdMarks = plot(freqs(ii),thdL(ii),'bx','MarkerSize',12,'LineWidth',1);
                end
                obj.thdDots = plot(freqs(:)',thdL,'b:','LineWidth',0.5);
            elseif obj.earNow == 2
                for ii=1:obj.nFreqs
                    obj.thdMarks = plot(freqs(ii),thdR(ii),'ro','MarkerSize',11,'LineWidth',1);
                end
                obj.thdDots = plot(freqs(:)',thdR,'r:','LineWidth',0.5);                
            end
        end    
        xticks(ff)
        xticklabels({'0.25','0.5','1.0','2.0','4.0','8.0','16.0','20'})
        xlabel('Frequency (kHz)','FontSize',12);
        ylabel(['Stim Level (dB ',char(obj.refNow),')'],'FontSize',12);
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
    function [micCorrection] = getMicCorrection(obj,probeInput,doMicCorrection)
        if isempty(probeInput)
            micCorrection = [];
            return
        end
        [pathName_mic,folderName_mic,fileName_mic] = mostRecentCalibration('mic',probeInput.label,[]);
        if doMicCorrection == 1
            if ~isempty(fileName_mic)
                dummy = load([pathName_mic,folderName_mic,fileName_mic]);
                micCorrection = dummy.micCorrection; % microphone correction
            else
                micCorrection = 1; % convolving by this changes nothing
            end
        else
            micCorrection = 1; % convolving by this changes nothing
        end
    end
    function [C1,C2,calPath1,calPath2] = getOutputCal(obj,calType,probeOutput)
        if isempty(probeOutput)
            C1 = [];
            C2 = [];
            calPath1.pathName = [];
            calPath1.folderName = [];
            calPath1.fileName = [];
            calPath2.pathName = [];
            calPath2.folderName = [];
            calPath2.fileName = [];
            return
        end
        if strcmp(calType,'thev')
            [pathName,folderName,fileName] = mostRecentCalibration('thev',probeOutput.label,1);
            calPath1.pathName = pathName;
            calPath1.folderName = folderName;
            calPath1.fileName = fileName;
            [pathName,folderName,fileName] = mostRecentCalibration('thev',probeOutput.label,2);
            calPath2.pathName = pathName;
            calPath2.folderName = folderName;
            calPath2.fileName = fileName;
        elseif strcmp(calType,'long')
            [pathName,folderName,fileName] = mostRecentCalibrationLT('longTube',probeOutput.label,1);
            calPath1.pathName = pathName;
            calPath1.folderName = folderName;
            calPath1.fileName = fileName;
            [pathName,folderName,fileName] = mostRecentCalibrationLT('longTube',probeOutput.label,2);
            calPath2.pathName = pathName;
            calPath2.folderName = folderName;
            calPath2.fileName = fileName;
        else
            error('Unrecognized calibration type. Must be thev or long.')
        end
        % load the calibration files
        dummy = load([calPath1.pathName,calPath1.folderName,calPath1.fileName]);
        if strcmp(calType,'thev')
            C1 = dummy.t;
        elseif strcmp(calType,'long')
            C1 = dummy.LTC;
        end
        dummy = load([calPath2.pathName,calPath2.folderName,calPath1.fileName]);
        if strcmp(calType,'thev')
            C2 = dummy.t;
        elseif strcmp(calType,'long')
            C2 = dummy.LTC;
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
    
end
end
