classdef playrecARLas < handle
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Class definition: playrecARLas
% For use with ARLas (Auditory Research Laboratory auditory software)
%
% Auditory Research Lab, The University of Iowa
% Deptartment of Communication Sciences & Disorders
% The University of Iowa
% Author: Shawn S. Goodman, PhD
% Date: September 13, 2016
% Last Updated: July 24, 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
properties (SetAccess = private)
    arlasVersion = '2017.07.24';
    sep                 % path delimiter appriate for the current operating system 
    map                 % struct containing file paths
    binFileName         % binary file path (full) and name (partial)
    fid                 % file identifier for binary files
    nSamples            % number of samples in the stimulus train
    nColumns            % number of columns in each channel of the stimulus train
    recDuration         % recording duration; -1 sets to length of the stimulus
    pageFiles           % pageFiles files; handles to playrec recordings
    completedPageFiles  % number of curently completed page files
    completedBuffers    % number of buffers written to disk
    writePauseLen       % how long to wait before re-checkng for completed page files
    currentWaveform     % the most currently recorded waveform
    mu                  % the average of the recorded waveforms
    objInit             % initialization object, passed from initARLas.

    H                   % handle for the main gui window, passed from initARLas, passed from arlas.
      VIEW              % handle for the main view window, passed from initARLAS, passed from arlas.
        WAVEFORM        % panel for plotting the waveforms
          CURRENTtxt    % panel for the instantaneous waveforms
          h_current     % plot for instantaneous waveforms
        AVERAGEtxt      % panel for the averaged waveforms
          h_average     % plot for the averaged waveforms.
      SLIDER1           % control the lower bound of time axis
      SLIDER2           % control toe upper bound of time axis
      CHANtxt_in        % text associated with current viewing channel
      CHAN_in           % popup for choosing current viewing channel
        indx_in_now     % current viewing this channel number
        LABEL           % input channel label for currently viewed channel
        AVGtxt          % completed buffers
        AVG             % how many buffers are completed in the current recording
        
    colorScheme         % what colors get plotted. Changes for each channel
    xAxisTxt            % x-axis label: time(s or ms)
    xAxisUnits          % x-axis label units: s, ms, us, etc.
    yAxisTxt            % y-axis label: amplitude (volts or pascals)
    yAxisUnits          % y-axis label units: V, mV, uV, Pa, mPa, uPa, etc.
    time                % time vector for plotting
    
    SOS                 % IIR filter saved as second-order-sections
    G                   % IIR filter gain
    ORDER               % IIR filter order
    extraBuffers        % number of extra buffers to collect to account for system delay
    plotDelay           % number of samples to shift plots to account for system delay
    skippedPages        % skipped page files
    failedRun           % determines whether page files have run out
    deadInTheWater      % detemines whether error already occured and whether to print error msg
    nReps               % number of stimulus repetitions to play/record    
end
properties (SetAccess = public) 
    fs                  % sample rate
    systemDelay         % delay of the system (in samples)
    card2volts          % multiplier to convert sound card units to voltage. Passed from objInit
    id_out              % currently chosen soundcard device ID
    name_out            % currently chosen device name
    id_in               % currently chosen value
    name_in             % currently chosen name
    maxChans_out        % max number of ouptut channels
    maxChans_in         % max number of input channels
    
    playChanList        % list of channels for playback (includes Ch0 for initialization purposes)
    recChanList         % list of channes for recording (includes Ch0 for initialization purposes)
    micSens             % microphone sensitivity for each inputchannel (includes initialization value)
    ampGain             % amplifier gain for each input channel (includes initialization value)
    label               % name label for each input channel (includes initialization value)
    loadingDock         % the stimuli to be used for playback (includes initialzation value)

    chans_out_now       % list of currently used channels (without Ch0 initialization)
    nChans_out_now      % number of currently used output channels
    chans_in_now        % list of currently used channels (without Ch0 initialization)
    nChans_in_now       % number of currently used input channels
    
    %nReps               % number of stimulus repetitions to play/record
    xStart              % sets the start of xlim for viewing the time plots
    xFinish             % sets the end of xlim for viewing the time plots
    
    Data                % the recorded data
    killRun             % stop running when the abort button was pressed
    pauseRun            % pause running when the pause button is pressed
    savedFiles          % list of the file names most recently saved
    savedFilesPath      % the path where these files are saved
    doFilter            % specifies whether or not to apply IIR high-pass filtering (125 Hz cutoff)
    
    userInfo            % Empty variable that can be filled by user with whatever they want.
                        %   will be written to the header file.
end
methods
    function objPlayrec = playrecARLas(objInit) % instantiate object of class initARLas
        try 
            if objInit.obj.ping == 1
                return
            end            
            objPlayrec.sep = filesep; % get the path delimiter appropriate to the operating system being used
            objPlayrec.map = objInit.obj.map; %pathRegistry(objPlayrec);
            objPlayrec.binFileName = [objPlayrec.map.data,'arlasBinary']; % location to write raw binary data
            objPlayrec.recDuration = -1; % recording duration; -1 sets to length of the stimulus
            objPlayrec.id_out = objInit.id_out_now;     % currently chosen ID
            objPlayrec.name_out = objInit.name_out_now;   % currently chosen device name
            objPlayrec.maxChans_out = objInit.chans_out;  % max number of output channels
            objPlayrec.id_in = objInit.id_in_now;     % currently chosen value
            objPlayrec.name_in = objInit.name_in_now;   % currently chosen name
            objPlayrec.maxChans_in = objInit.chans_in;  % max number of input channels 
            objPlayrec.playChanList = 0; % initialize to no channels
            objPlayrec.recChanList = 0;
            objPlayrec.fs = objInit.fs_now;
            objPlayrec.systemDelay = objInit.delay_now;
            objPlayrec.card2volts = objInit.card2volts_now;
            objPlayrec.ampGain = 0;
            objPlayrec.loadingDock.Ch0 = [];
            objPlayrec.label{1,1} = 'Empty';
            objPlayrec.micSens = 1;
            objPlayrec.objInit = objInit;
            objPlayrec.H = objInit.H; % pass handle to main gui
            objPlayrec.VIEW = objInit.VIEW; % pass handle to view panel of main gui
            objPlayrec.indx_in_now = 1; % set default view to first item in recChanList
            objPlayrec.completedBuffers = 0;
            objPlayrec.nReps = 0;
            objPlayrec.setFilter % load the filters in case they need to be used
            objPlayrec.doFilter = 1; % default is to turn on filtering
            objPlayrec.failedRun = 0;
            objPlayrec.deadInTheWater = 0;
        catch ME
            errorTxt = {'  Issue: Unexpected error creating object of class playrecARLas.'
                 '  Action: None.'
                 '  Location: playrecARLas.'
                };
            errorMsgARLas(errorTxt);
            objPlayrec.objInit.obj.buttonManager(51)
            objPlayrec.printError(ME)
        end
    end
    function abort(varargin) % instructions for aborting when gui closed
        try
            objPlayrec = varargin{1};
            objPlayrec.killPlots
            delete(objPlayrec);
        catch ME
            if objPlayrec.deadInTheWater == 1
                return
            else objPlayrec.deadInTheWater = 1;
            end
            errorTxt = {'  Issue: Unexpected error destructing object of class playrecARLas.'
                 '  Action: None.'
                 '  Location: playrecARLas.abort.'
                };
            errorMsgARLas(errorTxt);
            objPlayrec.objInit.obj.buttonManager(51)
            objPlayrec.printError(ME)            
        end
    end
    function initDataPlot(varargin) % create plots associated with playrecARLas
        objPlayrec = varargin{1};
        try
            set(objPlayrec.VIEW,'Title','VIEW: Playrec')
            % create panel ----- % create panels within figure ([left, bottom,width, height])
            objPlayrec.WAVEFORM = uipanel('Parent',objPlayrec.objInit.VIEW,...
                'FontSize',12,'BackgroundColor','white','Units','Normalized','Position',[.02 .2 .96 .79]); % .8
            objPlayrec.CURRENTtxt = uicontrol('Parent',objPlayrec.WAVEFORM,'Style','text',...
                'BackgroundColor','white','String','Current Waveform',...
                'HorizontalAlignment','Left','Units','Normalized','Position',...
                [.4 .92 .27 .06],'FontSize',12);
            objPlayrec.AVERAGEtxt = uicontrol('Parent',objPlayrec.WAVEFORM,'Style','text',...
                'BackgroundColor','white','String','Averaged Waveform',...
                'HorizontalAlignment','Left','Units','Normalized','Position',...
                [.38 .47 .29 .06],'FontSize',12); % [.38 .45 .29 .06]
            % create plots -----
            objPlayrec.h_current = axes('Parent',objPlayrec.WAVEFORM,'Visible','on','Color',[1 1 1],...
                'Units','Normalized','Position',[.125 .6 .85 .35],'XTick',[]);
            %axes(objPlayrec.h_current)
            ylabel('Amplitude (mPa)','FontSize',12);
            objPlayrec.h_average = axes('Parent',objPlayrec.WAVEFORM,'Visible','on','Color',[1 1 1],...
                'Units','Normalized','Position',[.125 .145 .85 .35]); % [.125 .125 .85 .35]
            %axes(objPlayrec.h_average)
            xlabel('Time (ms)','FontSize',12);
            ylabel('Amplitude (mPa)','FontSize',12);
            % [left,bottom,width,height]       
            % create controls -----
            % view input channel
            adjust = 0.075;
            objPlayrec.CHANtxt_in = uicontrol('Parent',objPlayrec.objInit.VIEW,'Style','text',...
                'BackgroundColor','white','String','View Channel:',...
                'HorizontalAlignment','Left','Units','Normalized',...
                'Position',[.06 .063-adjust .3 .1],'FontSize',10);
            objPlayrec.CHAN_in = uicontrol('Parent',objPlayrec.objInit.VIEW,'Style','popup',...
                'String',objPlayrec.recChanList,'FontSize',10,'BackgroundColor',...
                'white','HorizontalAlignment','Left','Units','Normalized',...
                'Position',[.25 .09-adjust .1 .09],'Value',objPlayrec.indx_in_now,...
                'Callback',@objPlayrec.viewManager);
            objPlayrec.LABEL = uicontrol('Parent',objPlayrec.objInit.VIEW,'Style','text',...
                'BackgroundColor','white','String',objPlayrec.label{objPlayrec.indx_in_now},...
                'HorizontalAlignment','Left','Units','Normalized',...
                'Position',[.25 .015-adjust .5 .1],'FontSize',10);
            % number of buffers in average
            objPlayrec.AVGtxt = uicontrol('Parent',objPlayrec.objInit.VIEW,'Style','text',...
                'BackgroundColor','white','String','Completed Buffers:',...
                'HorizontalAlignment','Left','Units','Normalized',...
                'Position',[.45 .063-adjust .5 .1],'FontSize',11);
            objPlayrec.AVG = uicontrol('Parent',objPlayrec.objInit.VIEW,'Style','text',...
                'String',num2str(objPlayrec.completedBuffers),'FontSize',11,'BackgroundColor',...
                'white','HorizontalAlignment','Left','Units','Normalized',...
                'Position',[.75 .063-adjust .2 .1]);
            % time axis slider
            objPlayrec.SLIDER1 = uicontrol('Parent',objPlayrec.objInit.VIEW,'Style','slider',...
                'BackgroundColor','white','Value',0,'Units','Normalized',...
                'Position',[.115 .15 .86 .03],'Callback',@objPlayrec.sliderManager1);
            objPlayrec.SLIDER2 = uicontrol('Parent',objPlayrec.objInit.VIEW,'Style','slider',...
                'BackgroundColor','white','Value',1,'Units','Normalized',...
                'Position',[.115 .12 .86 .03],'Callback',@objPlayrec.sliderManager2);
        catch ME
            if objPlayrec.deadInTheWater == 1
                return
            else objPlayrec.deadInTheWater = 1;
            end
            errorTxt = {'  Issue: Error initializing playrecARLas plots.'
                 '  Action: None.'
                 '  Location: playrecARLas.initDataPlot.'
                };
            errorMsgARLas(errorTxt);
            objPlayrec.objInit.obj.buttonManager(51)
            objPlayrec.printError(ME)            
        end
    end
    function killPlots(varargin) % close all plots associated with playrecARLas
        objPlayrec = varargin{1};
        try
            try set(objPlayrec.VIEW,'Title','VIEW: '); catch; end
            try delete(objPlayrec.WAVEFORM); catch; end
            try delete(objPlayrec.CURRENTtxt); catch; end
            try delete(objPlayrec.AVERAGEtxt); catch; end
            try delete(objPlayrec.h_current); catch; end
            try delete(objPlayrec.h_average); catch; end
            try delete(objPlayrec.CHANtxt_in); catch; end
            try delete(objPlayrec.CHAN_in); catch; end
            try delete(objPlayrec.LABEL); catch; end
            try delete(objPlayrec.AVGtxt); catch; end
            try delete(objPlayrec.AVG); catch; end
            try delete(objPlayrec.SLIDER1); catch; end
            try delete(objPlayrec.SLIDER2); catch; end
        catch ME
            if objPlayrec.deadInTheWater == 1
                return
            else objPlayrec.deadInTheWater = 1;
            end            
            errorTxt = {'  Issue: Error deleting playrecARLas plots.'
                 '  Action: None.'
                 '  Location: playrecARLas.killPlots.'
                };
            errorMsgARLas(errorTxt);
            objPlayrec.objInit.obj.buttonManager(51)
            objPlayrec.printError(ME)            
        end
    end
    
    function viewManager(varargin) % control which input channel is currenly being viewed
        objPlayrec = varargin{1};
        try
            objPlayrec.LABEL.String = objPlayrec.label{objPlayrec.CHAN_in.Value+1};
            objPlayrec.SLIDER1.Value = objPlayrec.xStart(objPlayrec.CHAN_in.Value,1);
            objPlayrec.SLIDER2.Value = objPlayrec.xFinish(objPlayrec.CHAN_in.Value,1);
            objPlayrec.doPlot
        catch ME
            if objPlayrec.deadInTheWater == 1
                return
            else objPlayrec.deadInTheWater = 1;
            end            
            errorTxt = {'  Issue: Error setting which input channel to view.'
                 '  Action: None.'
                 '  Location: in playrecARLas.viewManager.'
                };
            errorMsgARLas(errorTxt);
            objPlayrec.printError(ME)
        end
    end
    function sliderManager1(varargin) % manage the top slider, which controls the lower bound of the x-axes
        objPlayrec = varargin{1};
        try
            % do not let slider1 go to a higher value than slider2
            if objPlayrec.SLIDER1.Value > objPlayrec.SLIDER2.Value
                objPlayrec.SLIDER1.Value = objPlayrec.SLIDER2.Value * .95;
            end
            indx = objPlayrec.CHAN_in.Value;
            objPlayrec.xStart(indx,1) = objPlayrec.SLIDER1.Value;
            objPlayrec.doPlot
        catch ME
            if objPlayrec.deadInTheWater == 1
                return
            else objPlayrec.deadInTheWater = 1;
            end            
            errorTxt = {'  Issue: Error setting slider 1.'
                 '  Action: None.'
                 '  Location: in playrecARLas.sliderManager1.'
                };
            errorMsgARLas(errorTxt);
            objPlayrec.printError(ME)
        end        
    end
    function sliderManager2(varargin) % manage the bottom slider, which controls the upper bound of the x-axes
        objPlayrec = varargin{1};
        try
            % do not let slider2 go to a lower value than slider1
            if objPlayrec.SLIDER2.Value < objPlayrec.SLIDER1.Value
                objPlayrec.SLIDER2.Value = objPlayrec.SLIDER1.Value * 1.05;
            end
            indx = objPlayrec.CHAN_in.Value;
            objPlayrec.xFinish(indx,1) = objPlayrec.SLIDER2.Value;
            objPlayrec.doPlot
        catch ME
            if objPlayrec.deadInTheWater == 1
                return
            else objPlayrec.deadInTheWater = 1;
            end            
            errorTxt = {'  Issue: Error setting slider 2.'
                 '  Action: None.'
                 '  Location: in playrecARLas.sliderManager2.'
                };
            errorMsgARLas(errorTxt);
            objPlayrec.printError(ME)
        end                
    end
        
    function run(varargin) % run playback and record. This is the main subfunction
        objPlayrec = varargin{1};
        try
            objPlayrec.initData
            objPlayrec.initStream
            if objPlayrec.killRun == 1
                objPlayrec.killPlots
                return
            end
            objPlayrec.makePlotVariables
            objPlayrec.queue
            while objPlayrec.completedBuffers < objPlayrec.nReps+objPlayrec.extraBuffers
                objPlayrec.engine
                pause(0.001)
                if objPlayrec.killRun == 1 || objPlayrec.failedRun == 1
                    return                
                end
            end
            objPlayrec.killPlots
        catch ME
            objPlayrec.killRun = 1;
            objPlayrec.killPlots
            if objPlayrec.deadInTheWater == 1
                return
            else objPlayrec.deadInTheWater = 1;
            end            
            errorTxt = {'  Issue: Unexpected error running playrec.'
                 '  Action: None.'
                 '  Location: in playrecARLas.run.'
                };
            errorMsgARLas(errorTxt);
            objPlayrec.objInit.obj.buttonManager(51)
            objPlayrec.printError(ME)
        end
    end
    function engine(varargin) % sub-routine of run. Does a while loop to write data until finished, pause, or abort
        objPlayrec = varargin{1};
        try
            if objPlayrec.pauseRun == 0
                playrec('pause',0); % turn pause off
            end
            while objPlayrec.completedBuffers < objPlayrec.nReps+objPlayrec.extraBuffers
                if objPlayrec.killRun == 1
                    playrec('delPage');
                    objPlayrec.cleanUp % close down any open fids
                    objPlayrec.retrieveData % get the written raw data
                    objPlayrec.saveData % save the data in mat files, one matrix per channel
                    return
                elseif objPlayrec.pauseRun == 1
                    playrec('pause',1); % turn pause on
                    return
                elseif objPlayrec.failedRun == 1
                    playrec('delPage');
                    objPlayrec.cleanUp % close down any open fids
                    objPlayrec.retrieveData % get the written raw data
                    objPlayrec.saveData % save the data in mat files, one matrix per channel
                    return
                else
                    objPlayrec.writeData % stream raw recorded data to disk
                    pause(objPlayrec.writePauseLen)
                    if isempty(objPlayrec.pageFiles) || length(objPlayrec.pageFiles)<1
                        objPlayrec.failedRun = 1;
                    end
                end
            end
            playrec('delPage'); % kill any remaining sound
            objPlayrec.cleanUp % close down any open fids
            objPlayrec.retrieveData % get the written raw data
            objPlayrec.saveData % save the data in mat files, one matrix per channel
        catch ME
            objPlayrec.killRun = 1;
            objPlayrec.killPlots
            if objPlayrec.deadInTheWater == 1
                return
            else objPlayrec.deadInTheWater = 1;
            end            
            errorTxt = {'  Issue: Internal combustion failure ;-).'
                 '  Action: None.'
                 '  Location: in playrecARLas.engine.'
                };
            errorMsgARLas(errorTxt);
            objPlayrec.objInit.obj.buttonManager(51)
            objPlayrec.printError(ME)
        end
    end
    function initData(varargin) % initialize all variables
        try
            objPlayrec = varargin{1};
            objPlayrec.initDataPlot % initialize data plot
            objPlayrec.killRun = 0; % reset: run is not killed
            objPlayrec.pauseRun = 0; % reset: run is not paused
            objPlayrec.completedBuffers = 0; % reset: buffer counter starts at zero
            objPlayrec.Data = []; % there are no data saved
            objPlayrec.yAxisUnits = [];
            objPlayrec.mu = 0;
            objPlayrec.failedRun = 0;
            objPlayrec.pageFiles = [];
            objPlayrec.savedFiles = [];
            objPlayrec.skippedPages = [];
            if isempty(objPlayrec.nReps)
                set(objPlayrec.AVG,'String',[num2str(objPlayrec.completedBuffers),...
                    ' of N'])
            else
                set(objPlayrec.AVG,'String',[num2str(objPlayrec.completedBuffers),...
                    ' of ',num2str(objPlayrec.nReps)])
            end
        catch ME
            objPlayrec.killRun = 1;
            objPlayrec.killPlots
            if objPlayrec.deadInTheWater == 1
                return
            else objPlayrec.deadInTheWater = 1;
            end            
            errorTxt = {'  Issue: Error initializing playrecARLas variables.'
                 '  Action: None.'
                 '  Location: in playrecARLas.initData.'
                };
            errorMsgARLas(errorTxt);
            objPlayrec.objInit.obj.buttonManager(51)
            objPlayrec.printError(ME)
        end
    end
    function initStream(varargin) % initialize streaming of recorded data to disk
        objPlayrec = varargin{1};
        try
            objPlayrec.cleanUp % close down any open fids
            if objPlayrec.recChanList(1) ~= 0
                disp('Error: unexpected leading value for recChanList (0).')
                return
            end
            if length(objPlayrec.recChanList) < 2 || length(objPlayrec.playChanList) < 2
                errorTxt = {'  Issue: playChanList and/or recChanList of length < 1.'
                     '  Fix: Specify at least one value  in the experiment file using obj.setPlayList and obj.setRecList.'
                     '  Location: in playrecARLas.initStream.'
                    };
                errorMsgARLas(errorTxt);
                objPlayrec.objInit.obj.buttonManager(51)
                objPlayrec.killRun = 1;
                return
            end
            objPlayrec.chans_in_now = objPlayrec.recChanList(2:end);
            objPlayrec.chans_out_now = objPlayrec.playChanList(2:end);
            objPlayrec.nChans_in_now = length(objPlayrec.chans_in_now);
            objPlayrec.nChans_out_now = length(objPlayrec.chans_out_now);
            objPlayrec.CHAN_in.String =  objPlayrec.chans_in_now;
            objPlayrec.LABEL.String = objPlayrec.label{objPlayrec.CHAN_in.Value+1}; % new!!       
            for ii=1:objPlayrec.nChans_in_now
                try
                    if exist([objPlayrec.binFileName,'_',num2str(objPlayrec.chans_in_now(ii)),'.bin'],'file') ~= 0
                        delete([objPlayrec.binFileName,'_',num2str(objPlayrec.chans_in_now(ii)),'.bin']); % file id, write append, read, native machine format (little-endian on windows machines)
                    end
                catch ME
                    if objPlayrec.deadInTheWater == 1
                        return
                    else objPlayrec.deadInTheWater = 1;
                    end                    
                    errorTxt = {'  Issue: Error deleting old files.'
                         '  Action: None.'
                         '  Location: in playrecARLas.initStream.'
                        };
                    errorMsgARLas(errorTxt);
                    objPlayrec.objInit.obj.buttonManager(51)
                    objPlayrec.printError(ME)
                end
                if exist([objPlayrec.binFileName,'_',num2str(objPlayrec.chans_in_now(ii)),'.bin'],'file') ~= 0
                    delete([objPlayrec.binFileName,'_',num2str(objPlayrec.chans_in_now(ii)),'.bin']); % file id, write append, read, native machine format (little-endian on windows machines)
                    errorTxt = {'  Issue: Error deleting old files.'
                         '  Action: None.'
                         '  Location: in playrecARLas.initStream.'
                        };
                    errorMsgARLas(errorTxt);
                    objPlayrec.objInit.obj.buttonManager(51)
                end
            end
            try
                for ii=1:objPlayrec.nChans_in_now
                    objPlayrec.fid(1,ii) = fopen([objPlayrec.binFileName,'_',num2str(objPlayrec.chans_in_now(ii)),'.bin'],'a+'); % file id, write append, read, native machine format (little-endian on windows machines)
                end
            catch ME
                if objPlayrec.deadInTheWater == 1
                    return
                else objPlayrec.deadInTheWater = 1;
                end                
                errorTxt = {'  Issue: Error setting fid (file id) using fopen.'
                     '  Action: None.'
                     '  Location: in playrecARLas.initStream.'
                    };
                errorMsgARLas(errorTxt);
                objPlayrec.objInit.obj.buttonManager(51)
                objPlayrec.printError(ME)
            end
            if any(objPlayrec.fid==-1) % If fopen cannot open the file, it returns -1.
                errorTxt = {'  Issue: fopen unable to open one or more files.'
                     '  Action: None.'
                     '  Location: in playrecARLas.initStream.'
                    };
                errorMsgARLas(errorTxt);
                objPlayrec.objInit.obj.buttonManager(51)
                return
            end
        catch ME
            if objPlayrec.deadInTheWater == 1
                return
            else objPlayrec.deadInTheWater = 1;
            end            
            errorTxt = {'  Issue: Non-specific error initializing streaming files.'
                 '  Action: None.'
                 '  Location: in playrecARLas.initStream.'
                };
            errorMsgARLas(errorTxt);
            objPlayrec.objInit.obj.buttonManager(51)
            objPlayrec.printError(ME)
        end
    end
    function makePlotVariables(varargin) % create plotting labels
        objPlayrec = varargin{1};
        try
            for ii=1:objPlayrec.nChans_out_now
                [nSamples(ii,1),nColumns(ii,1)] = size(objPlayrec.loadingDock.(matlab.lang.makeValidName(['Ch',num2str(objPlayrec.chans_out_now(ii))]))); 
            end        
            if sum(diff(nSamples)) ~= 0
                errorTxt = {'  Issue: the number of rows in all stimTrain channels are not the same.'
                     '  Action: Aborting playrec.'
                     '  Location: in playrecARLas.makePlotVariables.'
                    };
                errorMsgARLas(errorTxt);
                objPlayrec.killRun = 1;
                objPlayrec.objInit.obj.killRun = 1;
                objPlayrec.killPlots
                return  
            else
                objPlayrec.nSamples = nSamples(1);
                objPlayrec.nColumns = nColumns;
            end            
            %objPlayrec.nSamples = size(objPlayrec.stimTrain,1); % number of times samples
            objPlayrec.time = (0:1:objPlayrec.nSamples-1)'/objPlayrec.fs; % time vector
            if objPlayrec.time(end) <= 0.5 % if less than half a second...
                objPlayrec.time = objPlayrec.time * 1000; %...express time in ms
                objPlayrec.xAxisUnits = '(ms)';
            else
                 objPlayrec.xAxisUnits = '(s)';
            end
            nChans = length(objPlayrec.chans_in_now); %nPlaybackChans
            objPlayrec.xStart = zeros(nChans,1);
            objPlayrec.xFinish = ones(nChans,1);
            for kk=1:nChans
               if objPlayrec.micSens(kk+1) == 1 % if mic sens = 1,
                   objPlayrec.yAxisUnits{kk} = char('V)'); % assume no microphone being used and plot in volts
               else
                   objPlayrec.yAxisUnits{kk} = char('Pa)'); % otherwise, microphone is being used, and plot in Pa
               end
            end
            objPlayrec.xAxisTxt = 'Time ';
            objPlayrec.yAxisTxt = 'Amplitude ';
            q = colormap(hsv(64)); % get the full hsv color map as a 64 by 3 matrix;
            q(5:15,:) = []; % throw out the yellow colors
            stepSize = floor(length(q) / objPlayrec.maxChans_in);
            indx = (1:1:objPlayrec.maxChans_in) * stepSize;
            objPlayrec.colorScheme = q(indx,:);
        catch ME
            if objPlayrec.deadInTheWater == 1
                return
            else objPlayrec.deadInTheWater = 1;
            end            
            errorTxt = {'  Issue: Error creating variables for plotting.'
                 '  Action: None.'
                 '  Location: in playrecARLas.makePlotVariables.'
                };
            errorMsgARLas(errorTxt);
            objPlayrec.objInit.obj.buttonManager(51)
            objPlayrec.printError(ME)
        end            
    end
    function queue(varargin) % load the queue and start it running
        objPlayrec = varargin{1};
        try % perform initial basic stimulus error checking
            if objPlayrec.nReps < 1 % number of reps must be 1 or more
                errorTxt = {'  Issue: nReps must be >= 1.'
                     '  Fix: In experiment file, make nReps >=2 and re-run.'
                     '  Location: in playrecARLas.queue.'
                    };
                errorMsgARLas(errorTxt);
                objPlayrec.objInit.obj.buttonManager(51)
                objPlayrec.killRun = 1;
                objPlayrec.killPlots
                return
            end
            for jj=1:objPlayrec.nChans_out_now % make sure output stream is longer than 3x filter order
                buffer = objPlayrec.loadingDock.(matlab.lang.makeValidName(['Ch',num2str(objPlayrec.chans_out_now(jj))]));
                if length(buffer) < (objPlayrec.ORDER * 3)
                    errorTxt = {'  Issue: Stimulus length must exceed IIR filter order by a factor of 3.'
                         ['  Fix: Make stimulus longer than ',num2str(objPlayrec.ORDER * 3),' samples.']
                         '  Location: in playrecARLas.queue.'
                        };
                    errorMsgARLas(errorTxt);
                    objPlayrec.objInit.obj.buttonManager(51)
                    objPlayrec.killRun = 1;
                    objPlayrec.killPlots
                    return
                end                
            end
            % make sure that no stimulus amplitudes exceed the max (+/- 1)
            for jj=1:objPlayrec.nChans_out_now % loop across output channels
                buffer = objPlayrec.loadingDock.(matlab.lang.makeValidName(['Ch',num2str(objPlayrec.chans_out_now(jj))]));
                if any(max(abs(buffer)) >= 1)
                    errorTxt = {'  Issue: Stimulus amplitudes cannot be or exceed +/- 1.'
                         '  Fix: Ensure that for all stimuli, max(abs(stimulus))<1 in experiment file.'
                         '  Location: in playrecARLas.queue.'
                        };
                    errorMsgARLas(errorTxt);
                    objPlayrec.objInit.obj.buttonManager(51)
                    objPlayrec.killRun = 1;
                    objPlayrec.killPlots
                    return
                end                        
            end
        catch ME
            if objPlayrec.deadInTheWater == 1
                return
            else objPlayrec.deadInTheWater = 1;
            end            
            errorTxt = {'  Issue: Error in basic stimulus check prior ot starting queue.'
                 '  Action: None.'
                 '  Location: in playrecARLas.queue.'
                };
            errorMsgARLas(errorTxt);
            objPlayrec.objInit.obj.buttonManager(51)
            objPlayrec.printError(ME)            
        end
        try % re-initialize playrec, get current system latencies
            if playrec('isInitialised')
                playrec('reset');
            end
            playrec('init',objPlayrec.fs,objPlayrec.id_out,objPlayrec.id_in);
            [playLatency,playSuggestedLatency] = playrec('getPlayLatency');
            [recLatency,recSuggestedLatency] = playrec('getRecLatency');
            playrec('reset');
            fs = objPlayrec.fs; % sample rate
            pd = objPlayrec.id_out; % play device
            rd = objPlayrec.id_in; % rec device
            pmc = objPlayrec.maxChans_out; % play max channel
            rmc = objPlayrec.maxChans_in; % rec max channel
            fpb = 0; % frames per buffer
            psl = playSuggestedLatency; % play suggested latency
            rsl = recSuggestedLatency; % rec suggested latency
                %objPlayrec.systemDelay = round(fs*playLatency); % note: this
                %was coming up too short; use the following command instead
            preferredDelay = round(fs*(recSuggestedLatency + playSuggestedLatency)); % system delay according to playrec
            if preferredDelay > objPlayrec.systemDelay
                ratio = objPlayrec.systemDelay / preferredDelay;
            else
                ratio = preferredDelay / objPlayrec.systemDelay;
            end
            if ratio > 0.95 % if preferred delay (playrec's estimate) is close to the user's calculated delay
                objPlayrec.systemDelay = preferredDelay; % use playrec's estimate
            else
                if objPlayrec.systemDelay ~= 0
                    disp(' ')
                    disp('ALERT: playrec estimated delay is more than 95% of measured delay!')
                    disp('       Using measured delay, not playrec estimated delay.')
                    disp(' ')
                end
            end
            playrec('init',fs,pd,rd,pmc,rmc,fpb,psl,rsl);
            pause(.001)
            objPlayrec.writePauseLen = (objPlayrec.nSamples / objPlayrec.fs) * 0.5; % estimate delay between each check for completed buffers
            if objPlayrec.nReps == 1
                NN = 1;
            else
                NN = objPlayrec.nReps * 2; % number of buffers to load into queue. Will load twice what was asked for, and then abort when done.
            end
            if objPlayrec.nReps > 1
                objPlayrec.extraBuffers = floor(objPlayrec.systemDelay / objPlayrec.nSamples) + 1; % number of extra buffers needed to account for system delay
                objPlayrec.plotDelay = mod(objPlayrec.systemDelay,objPlayrec.nSamples); % how much delay (based on system delay) to account for when plotting
            else
                objPlayrec.extraBuffers = 0;
                objPlayrec.plotDelay = 0;
            end
        catch ME
            if objPlayrec.deadInTheWater == 1
                return
            else objPlayrec.deadInTheWater = 1;
            end            
            errorTxt = {'  Issue: Error re-initializing playrec prior to building the queue.'
                 '  Action: None.'
                 '  Location: in playrecARLas.queue.'
                };
            errorMsgARLas(errorTxt);
            objPlayrec.objInit.obj.buttonManager(51)
            objPlayrec.printError(ME)
        end
        try % load the queue
            if ~any(objPlayrec.nColumns>1) % there is only one column in each channel--use simplified queue method
                stimTrain = [];
                for jj=1:objPlayrec.nChans_out_now
                    buffer = objPlayrec.loadingDock.(matlab.lang.makeValidName(['Ch',num2str(objPlayrec.chans_out_now(jj))]));
                    stimTrain = [stimTrain,buffer];
                end                            
                for ii=1:NN
                    pageNumList = playrec('playrec',stimTrain,objPlayrec.chans_out_now,objPlayrec.recDuration,objPlayrec.chans_in_now);
                    objPlayrec.pageFiles = [objPlayrec.pageFiles,pageNumList];
                    pause(0.0001)
                end
            else % there is more than one column in each channel--use full method
                Counters = zeros(NN,objPlayrec.nChans_out_now); % make a matrix of indices
                for jj=1:length(objPlayrec.chans_out_now)
                    set = (1:1:objPlayrec.nColumns(jj));
                    nSets = ceil(NN/length(set));
                    counter = repmat(set,1,nSets);
                    counter = counter(:);
                    counter = counter(1:NN,1);
                    Counters(:,jj) = counter;
                end
                for ii=1:NN
                    for jj=1:objPlayrec.nChans_out_now % loop across output channels
                        buffer = objPlayrec.loadingDock.(matlab.lang.makeValidName(['Ch',num2str(objPlayrec.chans_out_now(jj))]));
                        buffer = buffer(:,Counters(ii,jj));
                        stimTrain(:,jj) = buffer;
                    end
                    pageNumList = playrec('playrec',stimTrain,objPlayrec.chans_out_now,objPlayrec.recDuration,objPlayrec.chans_in_now);
                    objPlayrec.pageFiles = [objPlayrec.pageFiles,pageNumList];
                    pause(0.0001)
                end
            end
        catch ME
            objPlayrec.killRun = 1;
            objPlayrec.killPlots            
            if objPlayrec.deadInTheWater == 1
                return
            else objPlayrec.deadInTheWater = 1;
            end            
            errorTxt = {'  Issue: Error queuing playrec.'
                 '  Action: None.'
                 '  Location: in playrecARLas.queue.'
                };
            errorMsgARLas(errorTxt);
            objPlayrec.objInit.obj.buttonManager(51)
            objPlayrec.printError(ME)
        end
    end
    function writeData(varargin) % get the finished page files and write to disk
        objPlayrec = varargin{1};
        try % get completed page files from playrec and determine whether to continue with write process
            n = length(objPlayrec.pageFiles); % total number of page files
            isdone = zeros(n,1); % initialize vector showing finished recordings
            for ii=1:n % find and index which recordings are finished
                isdone(ii,1) = playrec('isFinished',objPlayrec.pageFiles(ii));
            end
            if objPlayrec.mu == 0  % if this is the first recorded buffer set
                if sum(isdone) < 2 && objPlayrec.nReps > 1 % if there are less than 2 buffers
                    return         % wait until there are at least 2 buffers
                end
            end
            objPlayrec.completedPageFiles = sum(isdone); % total number of completed recordings
            if objPlayrec.completedPageFiles == 0
                return
            elseif objPlayrec.completedPageFiles < 0
                errorTxt = {'  Issue: Unexpected page file error: Completed page files < 0.'
                     '  Action: None.'
                     '  Location: in playrecARLas.writeData.'
                    };
                errorMsgARLas(errorTxt);
                objPlayrec.objInit.obj.buttonManager(51)
            end
            skipped = playrec('getSkippedSampleCount'); % check to see if any "glitches" have occurred (according to playrec)
            if objPlayrec.nReps > 1
                if skipped ~= 0 % if glitches have occurred
                    objPlayrec.skippedPages = [objPlayrec.skippedPages;skipped]; % keep track of skipped samples
                    playrec('resetSkippedSampleCount'); % reset counter to zero
                    for ii=1:objPlayrec.completedPageFiles % the current set of files have a glitch, so delete them
                        playrec('delPage',objPlayrec.pageFiles(1)); % clear the current pageFiles
                        objPlayrec.pageFiles(1) = []; % delete finished pageFiles
                    end
                    return
                end
            end
        catch ME
            if objPlayrec.deadInTheWater == 1
                return
            else objPlayrec.deadInTheWater = 1;
            end            
            errorTxt = {'  Issue: Error checking pageFiles.'
                 '  Action: None.'
                 '  Location: in playrecARLas.writeData.'
                };
            errorMsgARLas(errorTxt);
            objPlayrec.objInit.obj.buttonManager(51)
            objPlayrec.printError(ME)
        end
        try % read in the finished recordings
            Recording = zeros(objPlayrec.nSamples*objPlayrec.completedPageFiles,objPlayrec.nChans_in_now); % initialize variable for finished buffers
            indx1 = 1; % initialize indices for placing data
            indx2 = objPlayrec.nSamples;
            for ii=1:objPlayrec.completedPageFiles % retrieve the finished recordings
                dummy = playrec('getRec',objPlayrec.pageFiles(1)); % get recorded data
                Recording(indx1:indx2,:) = dummy;
                playrec('delPage',objPlayrec.pageFiles(1)); % clear the current pageFile
                objPlayrec.pageFiles(1) = []; % delete finished pageFile
                indx1 = indx1 + objPlayrec.nSamples; % increment the indices
                indx2 = indx2 + objPlayrec.nSamples;
            end
            Recording = double(Recording); % convert from single precision (returned by playrec) to double precision 
            Recording = Recording * objPlayrec.card2volts; % convert from card units to a voltage. Done for all channels
            for kk=1:size(Recording,2)
               Recording(:,kk) = Recording(:,kk) / objPlayrec.micSens(kk+1);  % apply mic sensitivity
               Recording(:,kk) = Recording(:,kk) / 10^(objPlayrec.ampGain(kk+1)/20);  % apply input amplifier gain
            end
        catch ME
            if objPlayrec.deadInTheWater == 1
                return
            else objPlayrec.deadInTheWater = 1;
            end            
            errorTxt = {'  Issue: Error retrieving pageFiles.'
                 '  Action: None.'
                 '  Location: in playrecARLas.writeData.'
                };
            errorMsgARLas(errorTxt);
            objPlayrec.objInit.obj.buttonManager(51)
            objPlayrec.printError(ME)
        end
        try % write the finished recordings to disk
            if ~isempty(Recording)
                for ii=1:objPlayrec.nChans_in_now % loop over input channels; write each channel separately
                    fwrite(objPlayrec.fid(ii),Recording(:,ii),'double'); % stream the data to disk
                end
                objPlayrec.completedBuffers = objPlayrec.completedBuffers + objPlayrec.completedPageFiles; % increment the files written to disk counter
            end
        catch ME
            if objPlayrec.deadInTheWater == 1
                return
            else objPlayrec.deadInTheWater = 1;
            end            
            errorTxt = {'  Issue: Error writing recordings to disk. Probable error using fwrite.'
                 '  Action: None.'
                 '  Location: in playrecARLas.writeData.'
                };
            errorMsgARLas(errorTxt);
            objPlayrec.objInit.obj.buttonManager(51)
            objPlayrec.printError(ME)
        end
        try % % Update the waveform plots
            if ~isempty(Recording)
                % Note: Written data are not yet corrected for system delay. 
                %       Here we apply an "artificial" correction that will be
                %       re-done properly later when the data are extracted using objPlayrec.retrieveData()
                startingBuffer = (objPlayrec.completedBuffers - size(Recording,1)/objPlayrec.nSamples) +1;
                if startingBuffer < objPlayrec.extraBuffers-1 %size(Recording,1) <= (objPlayrec.nSamples*objPlayrec.extraBuffers-1)
                    return
                end
                if objPlayrec.doFilter == 1
                    Recording = filtfilt(objPlayrec.SOS,objPlayrec.G,Recording); % zero-phase IIR filtering
                end
                objPlayrec.currentWaveform = Recording(1:objPlayrec.nSamples,:); % plot the first waveform of the set
                objPlayrec.currentWaveform = [objPlayrec.currentWaveform(objPlayrec.plotDelay+1:end,:);objPlayrec.currentWaveform(1:objPlayrec.plotDelay,:)]; % account for system delay
                % Create mean of current Recording set. Note: must do this
                % in a loop to account for multiple input channels
                newBuffer = zeros(objPlayrec.nSamples,objPlayrec.nChans_in_now); 
                muIndx1 = 1;
                muIndx2 = objPlayrec.nSamples;
                start = 1;
                if objPlayrec.mu == 0 % the first set; nothing is currently in the mean
                    if startingBuffer < objPlayrec.extraBuffers
                        muIndx1 = (objPlayrec.extraBuffers-1)*objPlayrec.nSamples;
                        muIndx2 = muIndx1 + objPlayrec.nSamples -1;
                        start = 2;
                    end
                    objPlayrec.mu = Recording(muIndx1:muIndx2,:);
                    objPlayrec.mu = [objPlayrec.mu(objPlayrec.plotDelay+1:end,:);objPlayrec.mu(1:objPlayrec.plotDelay,:)]; % account for system delay
                end
                for ii=start:objPlayrec.completedPageFiles
                    newBuffer = newBuffer + Recording(muIndx1:muIndx2,:);
                    muIndx1 = muIndx1 + objPlayrec.nSamples;
                    muIndx2 = muIndx2 + objPlayrec.nSamples;
                end
                newBuffer = newBuffer / objPlayrec.completedPageFiles;
                newBuffer = [newBuffer(objPlayrec.plotDelay+1:end,:);newBuffer(1:objPlayrec.plotDelay,:)]; % account for system delay
                delta = newBuffer - objPlayrec.mu;
                objPlayrec.mu = objPlayrec.mu + (delta ./ (objPlayrec.completedPageFiles - (start-1)));
                objPlayrec.doPlot % update the plots
                objPlayrec.completedPageFiles = 0; % reset to zero
            end
        catch ME
            if objPlayrec.deadInTheWater == 1
                return
            else objPlayrec.deadInTheWater = 1;
            end            
            errorTxt = {'  Issue: Error updating waveform plots.'
                 '  Action: None.'
                 '  Location: in playrecARLas.writeData.'
                };
            errorMsgARLas(errorTxt);
            objPlayrec.objInit.obj.buttonManager(51)
            objPlayrec.printError(ME)
        end
    end
    function doPlot(varargin) % update the waveform plots
        objPlayrec = varargin{1};
        try % plot current data waveforms, using best scaling
            orig = gcf; % get the current figure; return that figure to top after plotting
            indx = objPlayrec.CHAN_in.Value;
            color = objPlayrec.colorScheme(indx,:); % select the plotting color for the waveform
            % scale y-axis appropriately
            N = length(objPlayrec.time);
            indxLow = ceil(N*objPlayrec.xStart(indx,1))+1;
            indxHigh = floor(N*objPlayrec.xFinish(indx,1));
            if indxLow > indxHigh
                return
            end
            ymax = max(abs(objPlayrec.currentWaveform(indxLow:indxHigh,indx)));
            decPlaces = [-9,-6,-3,0,3,6,9]'; % possible sets of decimal places for engineering notation (pico to giga)
            residual = round(ymax .* 10.^decPlaces); % try several possible engineering exponents
            indx2 = find(residual > 0); % find the one that gives the smallest mantissa
            if isempty(indx2)
                decPlaces = 9; % the best answer is smaller than this, but this is as small as we want to go
            else
                residual = residual(indx2);
                decPlaces = decPlaces(indx2);
                [~,indx2] = min(residual);
                decPlaces = decPlaces(indx2);
            end
            base = objPlayrec.yAxisUnits{indx};
            switch decPlaces
                case 9
                    units = ['(p',base]; % nano
                case 6
                    units = ['(\mu',base]; % micro
                case 3
                    units = ['(m',base]; % milli
                case 0
                    units = ['(',base]; % base
                case -3
                    units = ['(k',base]; % kilo
                case -6
                    units = ['(M',base]; % Mega
                case -9
                    units = ['(G',base]; % giga
                otherwise
            end
            axes(objPlayrec.h_current) % plot current waveform -----
            plot(objPlayrec.time(indxLow:indxHigh),objPlayrec.currentWaveform(indxLow:indxHigh,indx)*10^decPlaces,'Color',color) % plot the currently selected input channel
            xlim([objPlayrec.time(indxLow),objPlayrec.time(indxHigh)])
            ylabel([objPlayrec.yAxisTxt,units],'FontSize',12)
            axes(objPlayrec.h_average) % plot averaged waveform -----
            plot(objPlayrec.time(indxLow:indxHigh),objPlayrec.mu(indxLow:indxHigh,indx)*10^decPlaces,'Color',color) % plot the currently selected input channel
            xlim([objPlayrec.time(indxLow),objPlayrec.time(indxHigh)])
            xlabel([objPlayrec.xAxisTxt,objPlayrec.xAxisUnits],'FontSize',12)
            ylabel([objPlayrec.yAxisTxt,units],'FontSize',12)
            if objPlayrec.completedBuffers > objPlayrec.nReps
                set(objPlayrec.AVG,'String',[num2str(objPlayrec.nReps),' of ',num2str(objPlayrec.nReps)])    
            else
                set(objPlayrec.AVG,'String',[num2str(objPlayrec.completedBuffers),' of ',num2str(objPlayrec.nReps)])
            end
            %pause(0.0001)
            %figure(orig)
        catch ME
            objPlayrec.killRun = 1;
            objPlayrec.killPlots
            if objPlayrec.deadInTheWater == 1
                return
            else objPlayrec.deadInTheWater = 1;
            end            
            if objPlayrec.completedBuffers == 0
                return
            end
            errorTxt = {'  Issue: Error plotting data.'
                 '  Action: None.'
                 '  Location: in playrecARLas.doPlot.'
                };
            errorMsgARLas(errorTxt);
            objPlayrec.objInit.obj.buttonManager(51)
            objPlayrec.printError(ME)
        end
    end
    function retrieveData(varargin) % get the raw data that were written to the disk
        objPlayrec = varargin{1};
        try % open the files written to disk by writeData
            for ii=1:objPlayrec.nChans_in_now % get the file identifiers (fid) for the written binary files
                objPlayrec.fid(1,ii) = fopen([objPlayrec.binFileName,'_',num2str(objPlayrec.chans_in_now(ii)),'.bin'],'r'); % file id, read-only
            end
            if any(objPlayrec.fid == -1)
                errorTxt = {'  Issue: Error retrieving written data. Probable error obtaining fid using fopen.'
                     '  Action: None.'
                     '  Location: in playrecARLas.retrieveData.'
                    };
                errorMsgARLas(errorTxt);
                objPlayrec.objInit.obj.buttonManager(51)
                return
            end            
        catch ME
            if objPlayrec.deadInTheWater == 1
                return
            else objPlayrec.deadInTheWater = 1;
            end            
            errorTxt = {'  Issue: Error retrieving written data.'
                 '  Action: None.'
                 '  Location: in playrecARLas.retrieveData.'
                };
            errorMsgARLas(errorTxt);
            objPlayrec.objInit.obj.buttonManager(51)
            objPlayrec.printError(ME)
        end
        try % read the written data, filter, correct for system delay, and write to arlas 
            for ii=1:objPlayrec.nChans_in_now % loop across input channels
                status = fseek(objPlayrec.fid(1,ii),0,'bof'); % go to beginning of file
                if status == -1 % status is 0 on success and -1 on failure.
                    errorTxt = {'  Issue: Error using fseek to retrieve written data.'
                         '  Action: None.'
                         '  Location: in playrecARLas.retrieveData.'
                        };
                    errorMsgARLas(errorTxt);
                    objPlayrec.objInit.obj.buttonManager(51) 
                    return
                end
                X = fread(objPlayrec.fid(1,ii),inf,'double'); % read in the data
                if objPlayrec.doFilter == 1
                    X = filtfilt(objPlayrec.SOS,objPlayrec.G,X); % zero-phase IIR filtering
                end
                expectedLength = (objPlayrec.completedBuffers) * objPlayrec.nSamples; % number of samples to read in
                X = X(1:expectedLength); % include only the desired number of samples (discard any extras, if they exist)
                if objPlayrec.nReps > 1
                    X = X(objPlayrec.systemDelay+1:end); % cut off the leading system delay
                    X = X(1:end-mod(length(X),objPlayrec.nSamples)); % cut off the resulting trailing edge
                    expectedLength = objPlayrec.nReps * objPlayrec.nSamples; % number of samples after shifting to fix system delay
                    actualLength = floor(length(X)/objPlayrec.nSamples) * objPlayrec.nSamples; % actual length may be different if aborted
                    if actualLength < expectedLength % the case when run aborted early
                        X = X(1:actualLength); % include only the desired number of samples (discard any extras, if they exist)
                        warnTxt = {['  Issue: Run aborted. ',num2str(floor(length(X)/objPlayrec.nSamples)),' of ',num2str(objPlayrec.nReps),' reps saved.']
                             '  Action: None.'
                             '  Location: in playrecARLas.retrieveData.'
                            };
                        warnMsgARLas(warnTxt);
                        objPlayrec.objInit.obj.buttonManager(51)
                        objPlayrec.nReps = floor(length(X)/objPlayrec.nSamples); % update nReps to reflect actual number recorded
                    else
                        X = X(1:expectedLength); % include only the desired number of samples (discard any extras, if they exist)
                    end
                end
                objPlayrec.Data(:,ii) = X; % pass the recorded data to arlas
                if ~isempty(objPlayrec.skippedPages)
                    if length(objPlayrec.skippedPages) > 1
                        warnTxt = {'  Issue: Page files were skipped!.'
                             ['  ',num2str(length(objPlayrec.skippedPages)),' pages were identified and discarded.']
                             '  Action: None.'
                             '  Location: in playrecARLas.retrieveData.'
                            };
                        warnMsgARLas(warnTxt);
                    end
                end
            end
            objPlayrec.cleanUp
        catch ME
            if objPlayrec.deadInTheWater == 1
                return
            else objPlayrec.deadInTheWater = 1;
            end            
            errorTxt = {'  Issue: Error reading saved data.'
                 '  Action: None.'
                 '  Location: in playrecARLas.retrieveData.'
                };
            errorMsgARLas(errorTxt);
            objPlayrec.objInit.obj.buttonManager(51)
            objPlayrec.printError(ME)
        end        
    end
    function saveData(varargin) % save data (header and data) in .mat files to the folders specified by arlas
        objPlayrec = varargin{1};
        try % determine location and filename to use when saving data
            basePath = objPlayrec.map.data; % base location to write data
            expName = objPlayrec.objInit.obj.experimentID;
            subjName = objPlayrec.objInit.obj.subjectID;
            pathName = [basePath,expName,objPlayrec.sep,subjName,objPlayrec.sep];
            if exist(pathName,'dir') ~= 7 % if directory does not exist
                errorTxt = {'  Issue: Expected directory does not exist.'
                         '  Action: aborting data save operation.'
                         '  Location: in playrecARLas.saveData.'
                        };
                errorMsgARLas(errorTxt);
                objPlayrec.objInit.obj.buttonManager(51)                
                return
            end
            % never overwrite previously-written data. Always add an incrementing tag (4 digit number) to the end.
            for ii=1:objPlayrec.nChans_in_now % look for the appropriate counter
                label = objPlayrec.label{1,ii+1}; % part of the saved file name
                fileName = ['Ch',num2str(objPlayrec.chans_in_now(ii)),'_',label];
                d = dir([pathName,'*.mat']); % get list of existing files in the directory where data will be saved
                nFiles = length(d); % number of files to check
                hitCounter = 1;
                if nFiles > 0 % look for all the files with the same base file name (not including incrementing tag)
                    for jj=1:nFiles
                        if verLessThan('matlab', '9.1') % the function 'contains.m' was introduced in 2016b (9.1) version of matlab
                            dummy = strfind(d(jj).name,fileName);
                            if ~isempty(dummy)
                                dummy = 1;
                            end
                        else
                            dummy = contains(d(jj).name,fileName);
                        end
                        if dummy == 1
                            hitIndx(hitCounter,1) = jj;
                            hitCounter = hitCounter + 1;
                        end
                    end
                end
                hitCounter = hitCounter - 1;
                nHits = hitCounter; % number of previous files with same base filename
                if nHits > 0 % if same filenames exist, determine highest value incrementing tag that exists
                    q = [];
                    for jj=1:nHits
                        dummy = d(hitIndx(jj)).name;
                        dummy = dummy(end-7:end-4);
                        try
                            q(jj,1) = str2num(dummy);
                        catch
                        end
                    end
                    counter = max(q) + 1; % take the max incrementing tag value and increment it by 1
                else
                    counter = 1;
                end
                Counter(ii,1) = counter; % search across all recording channels
            end
            counter = max(Counter); % new incrementing tag is max value across channels
            % create a new incrementing tag string
            if counter < 10
                tag = ['000',num2str(counter)];
            elseif counter < 100
                tag = ['00',num2str(counter)];
            elseif counter < 1000
                tag = ['0',num2str(counter)];
            else
                tag = num2str(counter);
            end
        catch ME
            if objPlayrec.deadInTheWater == 1
                return
            else objPlayrec.deadInTheWater = 1;
            end            
            errorTxt = {'  Issue: Error determining save information (location, filename, or incrementing tag).'
                 '  Action: Data not saved.'
                 '  Location: in playrecARLas.saveData.'
                };
            errorMsgARLas(errorTxt);
            objPlayrec.objInit.obj.buttonManager(51)
            objPlayrec.printError(ME)
        end
        try % reshape the data into matrices, prepare the headers, and save
            for ii=1:objPlayrec.nChans_in_now
                data = objPlayrec.Data(:,ii); % each column is one channel
                % reshape into a matrix, with buffers in a column
                nSamples = length(data) / objPlayrec.nReps;
                if isnan(nSamples) % if no data were actually read in
                    return
                end
                if mod(nSamples,1) ~= 0 % if nSamples is a non-integer
                    nSamples = floor(nSamples); % round down
                    data = data(1:nSamples * objPlayrec.nReps); % cut to proper length
                    % and alert user, but continue
                    warnTxt = {'  Issue: Unexpected number of samples in recording.'
                         '  Action: Truncating samples to expected value.'
                         '  Location: in playrecARLas.saveData.'
                        };
                    warnMsgARLas(warnTxt);
                    objPlayrec.objInit.obj.buttonManager(51)
                end
                data = reshape(data,nSamples,objPlayrec.nReps);
                label = objPlayrec.label{1,ii+1};
                fileName = ['Ch',num2str(objPlayrec.chans_in_now(ii)),'_',label];
                fileName = [fileName,'_',tag,'.mat'];
                if exist([pathName,fileName],'file') == 2 % if file already exists
                    warnTxt = {'  Issue: Proposed new file name already exists: about to overwrite saved data!'
                         '  Action: Going into debug mode. Tell me what to do....'
                         '  Location: in playrecARLas.saveData.'
                        };
                    warnMsgARLas(warnTxt);
                    objPlayrec.objInit.obj.buttonManager(51)
                    keyboard
                end
                % header is a structure with important info about the data
                header.timeStamp = datestr(clock, 0); % date and time when this file was created
                header.subjID = objPlayrec.objInit.obj.subjectID;
                header.expID = objPlayrec.objInit.obj.experimentID;
                header.operator = objPlayrec.objInit.obj.operatorID; % person who collected the data
                header.fileName = fileName;
                header.pathName = pathName;
                header.arlasVersion = objPlayrec.objInit.obj.arlasVersion;
                header.fs = objPlayrec.fs; % sampling rate in Hz
                header.systemDelay = objPlayrec.systemDelay;
                header.card2volts = objPlayrec.card2volts;
                header.label = label;
                header.ampGain = objPlayrec.ampGain(ii+1);
                header.micSens = objPlayrec.micSens(ii+1);
                header.nSamples = objPlayrec.nSamples;
                header.nReps = objPlayrec.nReps;
                header.completedBuffers = objPlayrec.completedBuffers;
                header.mu = objPlayrec.mu;
                header.time = objPlayrec.time;
                header.userInfo = objPlayrec.userInfo;
                header.skippedPages = objPlayrec.skippedPages; % lists skipped files
                header.aborted = objPlayrec.killRun; % whether or not the run was aborted using the abort button
                
                save([pathName,fileName],'data','header')
                objPlayrec.savedFiles{ii} = fileName; % tell arlas what the most recent save was (file name)
                objPlayrec.savedFilesPath = pathName; % tell arlas what the most recent save was (location)
            end
        catch ME
            if objPlayrec.deadInTheWater == 1
                return
            else objPlayrec.deadInTheWater = 1;
            end            
            errorTxt = {'  Issue: Error preping and saving data.'
                 '  Action: Data not saved.'
                 '  Location: in playrecARLas.saveData.'
                };
            errorMsgARLas(errorTxt);
            objPlayrec.objInit.obj.buttonManager(51)
            objPlayrec.printError(ME)            
        end
    end
    function cleanUp(~) % close all open file identifiers
        try 
            fids = fopen('all'); % close all files open for read or write access
            if ~isempty(fids)
                for ii=1:length(fids)
                    fclose(fids(ii));
                end
            end
        catch ME
            if objPlayrec.deadInTheWater == 1
                return
            else objPlayrec.deadInTheWater = 1;
            end            
            warnTxt = {'  Issue: Error using fopen or fclose.'
                 '  Action: Check to see if open fid still exist.'
                 '  Location: in playrecARLas.cleanUp.'
                };
            warnMsgARLas(warnTxt);
            objPlayrec.objInit.obj.buttonManager(51)
            objPlayrec.printError(ME)
        end
    end  
    function printError(varargin)
        try
            ME = varargin{2};
            disp(ME.identifier)
            disp(ME.message)
            disp(ME.stack(1).line)
        catch
        end
    end
    function setFilter(varargin) % get IIR lowpass filter. All filters here
        objPlayrec = varargin{1};
        try % get the IIR filer needed for the sampling rate being used
            % are stored as second-order sections, direct form 1, with 
            % fstop = 75 Hz and fpass = 125 Hz.
            if objPlayrec.fs == 44100
                objPlayrec.SOS = [
                       1.000000000000000  -2.000000000000000   1.000000000000000   1.000000000000000  -1.997059754185550   0.997346175646593
                       1.000000000000000  -2.000000000000000   1.000000000000000   1.000000000000000  -1.991838843069538   0.992124515739271
                       1.000000000000000  -2.000000000000000   1.000000000000000   1.000000000000000  -1.986837300802882   0.987122256143531
                       1.000000000000000  -2.000000000000000   1.000000000000000   1.000000000000000  -1.982174292981907   0.982458579546617
                       1.000000000000000  -2.000000000000000   1.000000000000000   1.000000000000000  -1.977959289398182   0.978242971440427
                       1.000000000000000  -2.000000000000000   1.000000000000000   1.000000000000000  -1.974289895038583   0.974573050810499
                       1.000000000000000  -2.000000000000000   1.000000000000000   1.000000000000000  -1.971250043041826   0.971532762833373
                       1.000000000000000  -2.000000000000000   1.000000000000000   1.000000000000000  -1.968908537778030   0.969190921747197
                       1.000000000000000  -2.000000000000000   1.000000000000000   1.000000000000000  -1.967317928279519   0.967600084120962
                       1.000000000000000  -2.000000000000000   1.000000000000000   1.000000000000000  -1.966513689370135   0.966795729866364];
                objPlayrec.G = [ 
                       0.998601482458036
                       0.995990839702202
                       0.993489889236603
                       0.991158218132131
                       0.989050565209652
                       0.987215736462271
                       0.985695701468800
                       0.984524864881307
                       0.983729503100120
                       0.983327354809125
                       1.000000000000000];
                objPlayrec.ORDER = 20;
            elseif objPlayrec.fs == 48000
                objPlayrec.SOS = [
                       1.000000000000000  -2.000000000000000   1.000000000000000   1.000000000000000  -1.997319715515155   0.997561512368863
                       1.000000000000000  -2.000000000000000   1.000000000000000   1.000000000000000  -1.992520797315129   0.992762013208607
                       1.000000000000000  -2.000000000000000   1.000000000000000   1.000000000000000  -1.987921606685138   0.988162265797538
                       1.000000000000000  -2.000000000000000   1.000000000000000   1.000000000000000  -1.983632032144771   0.983872171958423
                       1.000000000000000  -2.000000000000000   1.000000000000000   1.000000000000000  -1.979753184017290   0.979992854255007
                       1.000000000000000  -2.000000000000000   1.000000000000000   1.000000000000000  -1.976375348882584   0.976614610197333
                       1.000000000000000  -2.000000000000000   1.000000000000000   1.000000000000000  -1.973576270910103   0.973815193366610
                       1.000000000000000  -2.000000000000000   1.000000000000000   1.000000000000000  -1.971419754576607   0.971658415963811
                       1.000000000000000  -2.000000000000000   1.000000000000000   1.000000000000000  -1.969954575701692   0.970193059713367
                       1.000000000000000  -2.000000000000000   1.000000000000000   1.000000000000000  -1.969213684450179   0.969452078769066];
                objPlayrec.G = [
                       0.998720306971005
                       0.996320702630934
                       0.994020968120669
                       0.991876051025798
                       0.989936509568074
                       0.988247489769979
                       0.986847866069178
                       0.985769542635105
                       0.985036908853765
                       0.984666440804811
                       1.000000000000000];
                objPlayrec.ORDER = 20;
            elseif objPlayrec.fs == 64000
                objPlayrec.SOS = [
                       1.000000000000000  -2.000000000000000   1.000000000000000   1.000000000000000  -1.998034479828863   0.998170534672852
                       1.000000000000000  -2.000000000000000   1.000000000000000   1.000000000000000  -1.994430660154271   0.994566469598529
                       1.000000000000000  -2.000000000000000   1.000000000000000   1.000000000000000  -1.990972765237639   0.991108339218817
                       1.000000000000000  -2.000000000000000   1.000000000000000   1.000000000000000  -1.987744063904802   0.987879418029686
                       1.000000000000000  -2.000000000000000   1.000000000000000   1.000000000000000  -1.984821521814066   0.984956676930368
                       1.000000000000000  -2.000000000000000   1.000000000000000   1.000000000000000  -1.982274155591861   0.982409137246935
                       1.000000000000000  -2.000000000000000   1.000000000000000   1.000000000000000  -1.980161618967374   0.980296456770656
                       1.000000000000000  -2.000000000000000   1.000000000000000   1.000000000000000  -1.978533028410711   0.978667755316190
                       1.000000000000000  -2.000000000000000   1.000000000000000   1.000000000000000  -1.977426029688228   0.977560681213357
                       1.000000000000000  -2.000000000000000   1.000000000000000   1.000000000000000  -1.976866103027017   0.977000716424308];
                objPlayrec.G = [
                       0.999051253625429
                       0.997249282438200
                       0.995520276114114
                       0.993905870483622
                       0.992444549686108
                       0.991170823209699
                       0.990114518934507
                       0.989300195931725
                       0.988746677725396
                       0.988466704862831
                       1.000000000000000];
                objPlayrec.ORDER = 20;
            elseif objPlayrec.fs == 88200
                objPlayrec.SOS = [ 
                       1.000000000000000  -2.000000000000000   1.000000000000000   1.000000000000000  -1.998600489269353   0.998672145001581
                       1.000000000000000  -2.000000000000000   1.000000000000000   1.000000000000000  -1.995982744128268   0.996054306006598
                       1.000000000000000  -2.000000000000000   1.000000000000000   1.000000000000000  -1.993468559529376   0.993540031266760
                       1.000000000000000  -2.000000000000000   1.000000000000000   1.000000000000000  -1.991118861113119   0.991190248606873
                       1.000000000000000  -2.000000000000000   1.000000000000000   1.000000000000000  -1.988990172526968   0.989061483700947
                       1.000000000000000  -2.000000000000000   1.000000000000000   1.000000000000000  -1.987133354619049   0.987204599220620
                       1.000000000000000  -2.000000000000000   1.000000000000000   1.000000000000000  -1.985592504179986   0.985663693537517
                       1.000000000000000  -2.000000000000000   1.000000000000000   1.000000000000000  -1.984404024797579   0.984475171544613
                       1.000000000000000  -2.000000000000000   1.000000000000000   1.000000000000000  -1.983595877874523   0.983666995647102
                       1.000000000000000  -2.000000000000000   1.000000000000000   1.000000000000000  -1.983187018487897   0.983258121601659];
                objPlayrec.G = [
                       0.999318158567733
                       0.998009262533716
                       0.996752147699034
                       0.995577277429998
                       0.994512914056979
                       0.993584488459917
                       0.992814049429376
                       0.992219799085548
                       0.991815718380406
                       0.991611285022389
                       1.000000000000000];
                objPlayrec.ORDER = 20;
            elseif objPlayrec.fs == 96000
                objPlayrec.SOS = [
                       1.000000000000000  -2.000000000000000   1.000000000000000   1.000000000000000  -1.998719476348764   0.998779964451356
                       1.000000000000000  -2.000000000000000   1.000000000000000   1.000000000000000  -1.996313888395105   0.996374303696360
                       1.000000000000000  -2.000000000000000   1.000000000000000   1.000000000000000  -1.994002986021370   0.994063331386798
                       1.000000000000000  -2.000000000000000   1.000000000000000   1.000000000000000  -1.991842843536058   0.991903123528170
                       1.000000000000000  -2.000000000000000   1.000000000000000   1.000000000000000  -1.989885525275958   0.989945746032910
                       1.000000000000000  -2.000000000000000   1.000000000000000   1.000000000000000  -1.988177914130336   0.988238083209121
                       1.000000000000000  -2.000000000000000   1.000000000000000   1.000000000000000  -1.986760684810996   0.986820810999564
                       1.000000000000000  -2.000000000000000   1.000000000000000   1.000000000000000  -1.985667434835688   0.985727527938764
                       1.000000000000000  -2.000000000000000   1.000000000000000   1.000000000000000  -1.984923982043184   0.984984052646830
                       1.000000000000000  -2.000000000000000   1.000000000000000   1.000000000000000  -1.984547834237828   0.984607893457952];
                objPlayrec.G = [
                       0.999374860200030
                       0.998172048022866
                       0.997016579352042
                       0.995936491766057
                       0.994957817827217
                       0.994103999334864
                       0.993395373952640
                       0.992848740693613
                       0.992477008672503
                       0.992288931923945
                       1.000000000000000];
                objPlayrec.ORDER = 20;
            elseif objPlayrec.fs == 128000
                objPlayrec.SOS = [
                       1.000000000000000  -2.000000000000000   1.000000000000000   1.000000000000000  -1.999050798606797   0.999084828519590
                       1.000000000000000  -2.000000000000000   1.000000000000000   1.000000000000000  -1.997245480027280   0.997279479208070
                       1.000000000000000  -2.000000000000000   1.000000000000000   1.000000000000000  -1.995510196791163   0.995544166432165
                       1.000000000000000  -2.000000000000000   1.000000000000000   1.000000000000000  -1.993887213416425   0.993921155429323
                       1.000000000000000  -2.000000000000000   1.000000000000000   1.000000000000000  -1.992415860445422   0.992449777411427
                       1.000000000000000  -2.000000000000000   1.000000000000000   1.000000000000000  -1.991131628080282   0.991165523184753
                       1.000000000000000  -2.000000000000000   1.000000000000000   1.000000000000000  -1.990065364329789   0.990099241283215
                       1.000000000000000  -2.000000000000000   1.000000000000000   1.000000000000000  -1.989242590437450   0.989276453384766
                       1.000000000000000  -2.000000000000000   1.000000000000000   1.000000000000000  -1.988682943270738   0.988716796691161
                       1.000000000000000  -2.000000000000000   1.000000000000000   1.000000000000000  -1.988399751689874   0.988433600289517];
                objPlayrec.G = [
                       0.999533906781597
                       0.998631239808838
                       0.997763590805832
                       0.996952092211437
                       0.996216409464212
                       0.995574287816259
                       0.995041151403251
                       0.994629760955554
                       0.994349934990475
                       0.994208337994848
                       1.000000000000000];
                objPlayrec.ORDER = 20;
            elseif objPlayrec.fs == 176000
                objPlayrec.SOS = [
                       1.000000000000000  -2.000000000000000   1.000000000000000   1.000000000000000  -1.999316333763183   0.999334335355423
                       1.000000000000000  -2.000000000000000   1.000000000000000   1.000000000000000  -1.998002708504451   0.998020698268973
                       1.000000000000000  -2.000000000000000   1.000000000000000   1.000000000000000  -1.996739433887339   0.996757412277497
                       1.000000000000000  -2.000000000000000   1.000000000000000   1.000000000000000  -1.995557371294714   0.995575339041729
                       1.000000000000000  -2.000000000000000   1.000000000000000   1.000000000000000  -1.994485292569836   0.994503250663989
                       1.000000000000000  -2.000000000000000   1.000000000000000   1.000000000000000  -1.993549204357446   0.993567154023179
                       1.000000000000000  -2.000000000000000   1.000000000000000   1.000000000000000  -1.992771745924802   0.992789688590398
                       1.000000000000000  -2.000000000000000   1.000000000000000   1.000000000000000  -1.992171671580954   0.992189608843556
                       1.000000000000000  -2.000000000000000   1.000000000000000   1.000000000000000  -1.991763426599580   0.991781360186395
                       1.000000000000000  -2.000000000000000   1.000000000000000   1.000000000000000  -1.991556823480408   0.991574755206995];
                objPlayrec.G = [
                       0.999662667279651
                       0.999005851693356
                       0.998374211541209
                       0.997783177584111
                       0.997247135808456
                       0.996779089595156
                       0.996390358628800
                       0.996090320106127
                       0.995886196696494
                       0.995782894671851
                       1.000000000000000];
                objPlayrec.ORDER = 20;
            elseif objPlayrec.fs == 192000
                objPlayrec.SOS = [
                       1.000000000000000  -2.000000000000000   1.000000000000000   1.000000000000000  -1.999374663338452   0.999389790103919
                       1.000000000000000  -2.000000000000000   1.000000000000000   1.000000000000000  -1.998170371650984   0.998185489305083
                       1.000000000000000  -2.000000000000000   1.000000000000000   1.000000000000000  -1.997012115533866   0.997027224424890
                       1.000000000000000  -2.000000000000000   1.000000000000000   1.000000000000000  -1.995928209678654   0.995943310369119
                       1.000000000000000  -2.000000000000000   1.000000000000000   1.000000000000000  -1.994945062440357   0.994960155692578
                       1.000000000000000  -2.000000000000000   1.000000000000000   1.000000000000000  -1.994086553129288   0.994101639886243
                       1.000000000000000  -2.000000000000000   1.000000000000000   1.000000000000000  -1.993373476130980   0.993388557492974
                       1.000000000000000  -2.000000000000000   1.000000000000000   1.000000000000000  -1.992823062421166   0.992838139618869
                       1.000000000000000  -2.000000000000000   1.000000000000000   1.000000000000000  -1.992448587023558   0.992463661388074
                       1.000000000000000  -2.000000000000000   1.000000000000000   1.000000000000000  -1.992259069034386   0.992274141965057];
                objPlayrec.G = [
                       0.999691113360593
                       0.999088965239017
                       0.998509834989689
                       0.997967880011943
                       0.997476304533234
                       0.997047048253883
                       0.996690508405989
                       0.996415300510009
                       0.996228062102908
                       0.996133302749861
                       1.000000000000000];
                objPlayrec.ORDER = 20;
            elseif objPlayrec.fs == 200000
                objPlayrec.SOS = [
                       1.000000000000000  -2.000000000000000   1.000000000000000   1.000000000000000  -1.999400250199571   0.999414191199800
                       1.000000000000000  -2.000000000000000   1.000000000000000   1.000000000000000  -1.998244073129201   0.998258006067880
                       1.000000000000000  -2.000000000000000   1.000000000000000   1.000000000000000  -1.997132039841839   0.997145965026764
                       1.000000000000000  -2.000000000000000   1.000000000000000   1.000000000000000  -1.996091343031063   0.996105260959636
                       1.000000000000000  -2.000000000000000   1.000000000000000   1.000000000000000  -1.995147349147031   0.995161260493520
                       1.000000000000000  -2.000000000000000   1.000000000000000   1.000000000000000  -1.994322999180300   0.994336904778934
                       1.000000000000000  -2.000000000000000   1.000000000000000   1.000000000000000  -1.993638273375844   0.993652174200165
                       1.000000000000000  -2.000000000000000   1.000000000000000   1.000000000000000  -1.993109730180287   0.993123627319293
                       1.000000000000000  -2.000000000000000   1.000000000000000   1.000000000000000  -1.992750127791144   0.992764022422789
                       1.000000000000000  -2.000000000000000   1.000000000000000   1.000000000000000  -1.992568134819156   0.992582028181838];
                objPlayrec.G = [
                       0.999703610349843
                       0.999125519799270
                       0.998569501217151
                       0.998049150997675
                       0.997577152410138
                       0.997164975989808
                       0.996822611894002
                       0.996558339374895
                       0.996378537553483
                       0.996287540750249
                       1.000000000000000];
                objPlayrec.ORDER = 20;
            end
        catch ME
            if objPlayrec.deadInTheWater == 1
                return
            else objPlayrec.deadInTheWater = 1;
            end            
            errorTxt = {'  Issue: Error obtaining IIR filter coefficients.'
                 '  Action: Unable to filter. Aborting run.'
                 '  Location: in playrecARLas.setFilter.'
                };
            errorMsgARLas(errorTxt);
            objPlayrec.objInit.obj.buttonManager(51)
            objPlayrec.printError(ME)                
        end
    end
    function setNReps(varargin) % set the number of playback/record repetitions
        objPlayrec = varargin{1};
        if nargin < 2
            error('setNReps must be given an input argument. Example: obj.setNReps(50);')
        end
        N = varargin{2};
        if N < 1
            error('number of repetitions mube be >= 1')
        end
        N = round(N); % force N to be an integer
        objPlayrec.nReps = N;
    end    
    
end
end