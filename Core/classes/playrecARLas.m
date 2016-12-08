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
% Last Updated: December 7, 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
properties (SetAccess = private)
    sep                 % path delimiter appriate for the current operating system 
    map                 % struct containing file paths
    binFileName         % binary file path (full) and name (partial)
    fid                 % file identifier for binary files
    nSamples            % number of samples in the stimulus train
    nColumns            % number of columns in each channel of the stimulus train
    recDuration         % recording duration; -1 sets to length of the stimulus
    playChanList        % list of channels for playback
    recChanList         % list of channes for recording
    pageFiles           % pageFiles files; handles to playrec recordings
    completedPageFiles  % number of curently completed page files
    completedBuffers    % number of buffers written to disk
    holdingPen          % storage space for buffering
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
    xAxisTxt
    xAxisUnits
    yAxisTxt
    yAxisUnits
    time
end
properties (SetAccess = public) 
    card2volts % multiplier to convert sound card units to voltage. Passed from objInit
    micSens    % microphone sensitivity for each inputchannel
    ampGain    % amplifier gain for each input channel
    label      % name label for each input channel
    
    fs         % currently chosen sample rate
    id_out     % currently chosen soundcard device ID
    name_out   % currently chosen device name
    chans_out  % number of channels for currently chosen device
    id_in      % currently chosen value
    name_in    % currently chosen name
    chans_in   % curently chosen value
    systemDelay % delay of the system (in samples)
    
    stimTrain  % the stimulus to present (older, non-explicit way)
    stimTrainEx % the stimulus to present (newer, explicit way)
    nReps      % number of stimulus repetitions to play/record
    xStart     % sets the start of xlim for viewing the time plots
    xFinish    % sets the end of xlim for viewing the time plots
    
    Data       % the recorded data
    killRun    % stop running when the abort button was pressed
    pauseRun   % pause running when the pause button is pressed
    savedFiles % alist of the file names most recently saved
    savedFilesPath % the path where these files are saved
    
    userInfo   % Empty variable that can be filled by user with whatever they want.
               %   will be written to the header file.
end

methods
    function objPlayrec = playrecARLas(objInit) % initialize object of class initARLas
        objPlayrec.sep = filesep; % get the path delimiter appropriate to the operating system being used
        objPlayrec.map = objInit.obj.map; %pathRegistry(objPlayrec);
        objPlayrec.binFileName = [objPlayrec.map.data,'arlasBinary']; % location to write raw binary data
        objPlayrec.recDuration = -1; % recording duration; -1 sets to length of the stimulus
        objPlayrec.id_out = objInit.id_out_now;     % currently chosen ID
        objPlayrec.name_out = objInit.name_out_now;   % currently chosen device name
        objPlayrec.chans_out = objInit.chans_out_now;  % number of channels for currently chosen device
        objPlayrec.id_in = objInit.id_in_now;     % currently chosen value
        objPlayrec.name_in = objInit.name_in_now;   % currently chosen name
        objPlayrec.chans_in = objInit.chans_in_now;  % curently chosen value 
        if objPlayrec.chans_out == num2str(0)
            objPlayrec.playChanList = 0;
        else
            objPlayrec.playChanList = (1:1:objPlayrec.chans_out);
        end
        if objPlayrec.chans_in == num2str(0)
            objPlayrec.recChanList = 0;
        else
            objPlayrec.recChanList = (1:1:objPlayrec.chans_in);
        end
        objPlayrec.fs = objInit.fs_now;
        objPlayrec.systemDelay = objInit.delay_now;
        objPlayrec.card2volts = objInit.card2volts_now;
        names = fieldnames(objInit.micSens);
        for ii=1:size(names,1)-1
            objPlayrec.micSens(ii,1) = getfield(objInit.micSens,names{ii+1});
            objPlayrec.ampGain(ii,1) = getfield(objInit.ampGain,names{ii+1});
            objPlayrec.label{ii,1} = getfield(objInit.label,names{ii+1});
        end
        for ii=1:objInit.chans_out_now % make stimTrainEx into a structure with one field for each initialized output channel
            objPlayrec.stimTrainEx.(matlab.lang.makeValidName(['Ch',num2str(ii)])) = []; 
        end
        objPlayrec.ampGain = 10.^(objPlayrec.ampGain / 20); % convert to linear multiplier
        objPlayrec.objInit = objInit;
        objPlayrec.H = objInit.H; % pass handle to main gui
        objPlayrec.VIEW = objInit.VIEW; % pass handle to view panel of main gui
        objPlayrec.indx_in_now = 1; % set default view to first item in recChanList
        objPlayrec.initDataPlot % initialize data plot
        objPlayrec.completedBuffers = 0;
    end
    function abort(varargin) % instructions for aborting when gui closed
        objPlayrec = varargin{1};
        objPlayrec.killPlots
        delete(objPlayrec);
    end
    function initDataPlot(varargin)
        objPlayrec = varargin{1};
        set(objPlayrec.VIEW,'Title','VIEW: Playrec')
        %objPlayrec.objInit.makeInvisible
        % create panel -----
        objPlayrec.WAVEFORM = uipanel('Parent',objPlayrec.objInit.VIEW,...
            'FontSize',12,'BackgroundColor','white','Units','Normalized','Position',[.02 .2 .96 .8]);
        objPlayrec.CURRENTtxt = uicontrol('Parent',objPlayrec.WAVEFORM,'Style','text',...
            'BackgroundColor','white','String','Current Waveform',...
            'HorizontalAlignment','Left','Units','Normalized','Position',...
            [.4 .92 .27 .06],'FontSize',12);
        objPlayrec.AVERAGEtxt = uicontrol('Parent',objPlayrec.WAVEFORM,'Style','text',...
            'BackgroundColor','white','String','Averaged Waveform',...
            'HorizontalAlignment','Left','Units','Normalized','Position',...
            [.38 .45 .29 .06],'FontSize',12);
        % create plots -----
        objPlayrec.h_current = axes('Parent',objPlayrec.WAVEFORM,'Visible','on','Color',[1 1 1],...
            'Units','Normalized','Position',[.125 .6 .85 .35],'XTick',[]);
        ylabel('Amplitude (mPa)','FontSize',12);
        objPlayrec.h_average = axes('Parent',objPlayrec.WAVEFORM,'Visible','on','Color',[1 1 1],...
            'Units','Normalized','Position',[.125 .125 .85 .35]);
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
            'Position',[.25 .07-adjust .1 .1],'Value',objPlayrec.indx_in_now,...
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
    end
    function makeVisible(varargin)
        objPlayrec = varargin{1};
        try
            set(objPlayrec.VIEW,'Title','VIEW: Playrec')
            set(objPlayrec.WAVEFORM,'Visible','on')
            set(objPlayrec.CURRENTtxt,'Visible','on')
            set(objPlayrec.AVERAGEtxt,'Visible','on')
            set(objPlayrec.CHAN_in,'Visible','on')
            set(objPlayrec.CHANtxt_in,'Visible','on')
            set(objPlayrec.AVGtxt,'Visible','on')
            set(objPlayrec.AVG,'Visible','on')
            set(objPlayrec.LABEL,'Visible','on')
            set(objPlayrec.SLIDER1,'Visible','on')
            set(objPlayrec.SLIDER2,'Visible','on')
        catch ME
            errorTxt = {'  Issue: Error making visible playrec plots.'
                 '  Action: None.'
                 '  Location: in playrecARLas.makeVisible.'
                };
            errorMsgARLas(errorTxt);
            objPlayrec.printError(ME)
        end
    end
    function makeInvisible(varargin)
        objPlayrec = varargin{1};
        try
            set(objPlayrec.VIEW,'Title','VIEW: ')
            set(objPlayrec.WAVEFORM,'Visible','off')
            set(objPlayrec.CURRENTtxt,'Visible','off')
            set(objPlayrec.AVERAGEtxt,'Visible','off')
            set(objPlayrec.CHAN_in,'Visible','off')
            set(objPlayrec.CHANtxt_in,'Visible','off')
            set(objPlayrec.AVGtxt,'Visible','off')
            set(objPlayrec.AVG,'Visible','off')
            set(objPlayrec.LABEL,'Visible','off')
            set(objPlayrec.SLIDER1,'Visible','off')
            set(objPlayrec.SLIDER2,'Visible','off')
        catch ME
            errorTxt = {'  Issue: Error making invisible playrec plots.'
                 '  Action: None.'
                 '  Location: in playrecARLas.makeInvisible.'
                };
            errorMsgARLas(errorTxt);
            objPlayrec.printError(ME)
        end
    end
    function killPlots(varargin)
        objPlayrec = varargin{1};
        try
            a = findall(gcf);
            counter = 1;
            for ii=1:size(a,1)
                try
                if strcmp(a(ii).Parent.Title,'VIEW: Playrec')
                    indx(counter,1) = ii;
                    counter = counter + 1;
                end
                catch
                end
            end
            for ii=1:size(indx,1)
                delete(a(indx(ii)))
            end        
            set(objPlayrec.VIEW,'Title','VIEW: ')
        catch ME
            errorTxt = {'  Issue: Error deleting playrec plots.'
                 '  Action: None.'
                 '  Location: in playrecARLas.killPlots.'
                };
            errorMsgARLas(errorTxt);
            objPlayrec.printError(ME)
        end
    end
    
    function viewManager(varargin) % control which input channel is currenly being viewed
        objPlayrec = varargin{1};
        try
            objPlayrec.LABEL.String = objPlayrec.label{objPlayrec.CHAN_in.Value};
            objPlayrec.SLIDER1.Value = objPlayrec.xStart(objPlayrec.CHAN_in.Value,1);
            objPlayrec.SLIDER2.Value = objPlayrec.xFinish(objPlayrec.CHAN_in.Value,1);
            objPlayrec.doPlot
        catch ME
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
        % do not let slider1 go to a higher value than slider2
        if objPlayrec.SLIDER1.Value > objPlayrec.SLIDER2.Value
            objPlayrec.SLIDER1.Value = objPlayrec.SLIDER2.Value * .95;
        end
        indx = objPlayrec.recChanList(objPlayrec.CHAN_in.Value); % get the index of the currently-selected input channel        
        objPlayrec.xStart(indx,1) = objPlayrec.SLIDER1.Value;
        objPlayrec.doPlot
    end
    function sliderManager2(varargin) % manage the bottom slider, which controls the upper bound of the x-axes
        objPlayrec = varargin{1};
        % do not let slider2 go to a lower value than slider1
        if objPlayrec.SLIDER2.Value < objPlayrec.SLIDER1.Value
            objPlayrec.SLIDER2.Value = objPlayrec.SLIDER1.Value * 1.05;
        end
        indx = objPlayrec.recChanList(objPlayrec.CHAN_in.Value); % get the index of the currently-selected input channel
        objPlayrec.xFinish(indx,1) = objPlayrec.SLIDER2.Value;
        objPlayrec.doPlot
    end
        
    function run(varargin) % run playback and record. This is the main subfunction
        objPlayrec = varargin{1};
        try
            objPlayrec.makeVisible
            objPlayrec.initData
            objPlayrec.initStream
            if objPlayrec.killRun == 1
                objPlayrec.makeInvisible
                return
            end
            objPlayrec.makePlotVariables
            objPlayrec.queue
            while objPlayrec.completedBuffers < objPlayrec.nReps
                objPlayrec.engine
                pause(0.001)
                if objPlayrec.killRun == 1
                    return                
                end
            end
        catch ME
            errorTxt = {'  Issue: Unexpected error running playrec.'
                 '  Action: None.'
                 '  Location: in playrecARLas.run.'
                };
            errorMsgARLas(errorTxt);
            objPlayrec.objInit.obj.buttonManager(51)
            objPlayrec.printError(ME)
        end
    end
    function engine(varargin)
        objPlayrec = varargin{1};
        if objPlayrec.pauseRun == 0
            playrec('pause',0); % turn pause off
        end
        while objPlayrec.completedBuffers < objPlayrec.nReps
            if objPlayrec.killRun == 1
                playrec('delPage');
                objPlayrec.cleanUp % close down any open fids
                return
            elseif objPlayrec.pauseRun == 1
                playrec('pause',1); % turn pause on
                return
            else
                objPlayrec.writeData % stream raw recorded data to disk
                pause(objPlayrec.writePauseLen)
            end
        end
        playrec('delPage'); % kill any remaining sound
        objPlayrec.cleanUp % close down any open fids
        objPlayrec.retrieveData % get the written raw data
        objPlayrec.saveData % save the data in mat files, one matrix per channel
    end
    function initData(varargin) % initialize all data containers to zero
        objPlayrec = varargin{1};
        objPlayrec.killRun = 0; % reset: run is not killed
        objPlayrec.pauseRun = 0; % reset: run is not paused
        objPlayrec.completedBuffers = 0; % reset: buffer counter starts at zero
        objPlayrec.holdingPen = []; % nothing is in the holding pen
        objPlayrec.Data = []; % there are no data saved
        objPlayrec.yAxisUnits = [];
        objPlayrec.mu = 0;
        objPlayrec.pageFiles = [];
        objPlayrec.savedFiles = [];
        if isempty(objPlayrec.nReps)
            set(objPlayrec.AVG,'String',[num2str(objPlayrec.completedBuffers),...
                ' of N'])
        else
            set(objPlayrec.AVG,'String',[num2str(objPlayrec.completedBuffers),...
                ' of ',num2str(objPlayrec.nReps)])
        end            
    end
    function initStream(varargin) % initialize streaming recorded data to disk
        objPlayrec = varargin{1};
        objPlayrec.cleanUp % close down any open fids
        if objPlayrec.chans_out == 0
            errorTxt = {'  Issue: Number of output channels must be > 0.'
                 '  Suggested Action: If no output is desired, make output stimulus a vector of zeros of desired length.'
                 '  Location: in playrecARLas.initStream.'
                };
            errorMsgARLas(errorTxt);
            objPlayrec.objInit.obj.buttonManager(51)
            objPlayrec.killRun = 1;
            return
        end
        if objPlayrec.chans_in == 0
            errorTxt = {'  Issue: Number of input channels must be > 0.'
                 '  Action: None.'
                 '  Location: in playrecARLas.initStream.'
                };
            errorMsgARLas(errorTxt);
            objPlayrec.objInit.obj.buttonManager(51)
            objPlayrec.killRun = 1;
            return
        end
        for ii=1:objPlayrec.chans_in
            try
                if exist([objPlayrec.binFileName,'_',num2str(ii),'.bin'],'file') ~= 0
                    delete([objPlayrec.binFileName,'_',num2str(ii),'.bin']); % file id, write append, read, native machine format (little-endian on windows machines)
                end
            catch ME
                errorTxt = {'  Issue: Error deleting old files.'
                     '  Action: None.'
                     '  Location: in playrecARLas.initStream.'
                    };
                errorMsgARLas(errorTxt);
                objPlayrec.objInit.obj.buttonManager(51)
                objPlayrec.printError(ME)
            end
        end
        try
            for ii=1:objPlayrec.chans_in
                objPlayrec.fid(1,ii) = fopen([objPlayrec.binFileName,'_',num2str(ii),'.bin'],'a+'); % file id, write append, read, native machine format (little-endian on windows machines)
            end
        catch ME
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
    end
    function makePlotVariables(varargin) % create plotting labels
        objPlayrec = varargin{1};
        try
            if ~isempty(objPlayrec.stimTrain) % old, non-explicit method is being used
                [nRows,nCols] = size(objPlayrec.stimTrain);
                for ii=1:objPlayrec.objInit.chans_out_now % loop across output channels
                    if ii <= nCols % where possible, load the old format (stimTrain) into the new format (stimTrainEx)
                        objPlayrec.stimTrainEx.(matlab.lang.makeValidName(['Ch',num2str(ii)])) = objPlayrec.stimTrain(:,ii);
                    else
                        objPlayrec.stimTrainEx.(matlab.lang.makeValidName(['Ch',num2str(ii)])) = zeros(nRows,1);
                    end
                end                       
            end
        catch
        end
        
        try
            for ii=1:objPlayrec.objInit.chans_out_now 
                [nSamples(ii,1),nColumns(ii,1)] = size(objPlayrec.stimTrainEx.(matlab.lang.makeValidName(['Ch',num2str(ii)]))); 
            end        
            if sum(diff(nSamples)) ~= 0
                errorTxt = {'  Issue: the number of rows in all stimTrainEx channels are not the same.'
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
            %objPlayrec.nSamples = size(objPlayrec.stimTrainEx,1); % number of times samples
            objPlayrec.time = (0:1:objPlayrec.nSamples-1)'/objPlayrec.fs; % time vector
            if objPlayrec.time(end) <= 0.5 % if less than half a second...
                objPlayrec.time = objPlayrec.time * 1000; %...express time in ms
                objPlayrec.xAxisUnits = '(ms)';
            else
                 objPlayrec.xAxisUnits = '(s)';
            end
            
            nChans = objPlayrec.chans_in;
            objPlayrec.xStart = zeros(nChans,1);
            objPlayrec.xFinish = ones(nChans,1);
            for kk=1:nChans
               if objPlayrec.micSens(kk) == 1 % if mic sens = 1,
                   objPlayrec.yAxisUnits{kk} = char('V)'); % assume no microphone being used and plot in volts
               else
                   objPlayrec.yAxisUnits{kk} = char('Pa)'); % otherwise, microphone is being used, and plot in Pa
               end
            end
            objPlayrec.xAxisTxt = 'Time ';
            objPlayrec.yAxisTxt = 'Amplitude ';
            q = colormap(hsv(64)); % get the full hsv color map as a 64 by 3 matrix;
            q(5:15,:) = []; % throw out the yellow colors
            stepSize = floor(length(q) / objPlayrec.chans_in);
            indx = (1:1:objPlayrec.chans_in) * stepSize;
            objPlayrec.colorScheme = q(indx,:);
        catch ME
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
        objPlayrec.writePauseLen = (objPlayrec.nSamples / objPlayrec.fs) * 0.5; % estimate delay between each check for completed buffers
        NN = objPlayrec.nReps * 2; % number of buffers to load into queue. Load twice what was asked for, and then abort when done.
        if playrec('isInitialised')
            playrec('reset')
        end
        try % re-initialize playrec
            playrec('init',objPlayrec.fs,objPlayrec.id_out,objPlayrec.id_in)
            pause(.01)
        catch ME
            errorTxt = {'  Issue: Error initializing playrec.'
                 '  Action: None.'
                 '  Location: in playrecARLas.queue.'
                };
            errorMsgARLas(errorTxt);
            objPlayrec.objInit.obj.buttonManager(51)
            objPlayrec.printError(ME)
        end
        % main loop is here -----
        try
            if ~any(objPlayrec.nColumns>1) % there is only one column in each channel--use simplified queue method
                stimTrainEx = [];
                for jj=1:objPlayrec.objInit.chans_out_now
                    buffer = objPlayrec.stimTrainEx.(matlab.lang.makeValidName(['Ch',num2str(jj)])); 
                    stimTrainEx = [stimTrainEx,buffer];
                end                            
                for ii=1:NN
                    pageNumList = playrec('playrec',stimTrainEx,objPlayrec.playChanList,objPlayrec.recDuration,objPlayrec.recChanList);
                    objPlayrec.pageFiles = [objPlayrec.pageFiles,pageNumList];
                end
            else % otherwise use full method ----------------
                Counters = zeros(NN,objPlayrec.objInit.chans_out_now); % make a matrix of indices
                for jj=1:objPlayrec.objInit.chans_out_now
                    set = (1:1:objPlayrec.nColumns(jj));
                    nSets = ceil(NN/length(set));
                    counter = repmat(set,1,nSets);
                    counter = counter(:);
                    counter = counter(1:NN,1);
                    Counters(:,jj) = counter;
                end
                for ii=1:NN
                    for jj=1:objPlayrec.objInit.chans_out_now % loop across output channels
                        buffer = objPlayrec.stimTrainEx.(matlab.lang.makeValidName(['Ch',num2str(jj)]));
                        buffer = buffer(:,Counters(ii,jj));
                        stimTrainEx(:,jj) = buffer;
                    end
                    pageNumList = playrec('playrec',stimTrainEx,objPlayrec.playChanList,objPlayrec.recDuration,objPlayrec.recChanList);
                    objPlayrec.pageFiles = [objPlayrec.pageFiles,pageNumList];
                end
            end
        catch ME
            errorTxt = {'  Issue: Error queuing playrec.'
                 '  Action: None.'
                 '  Location: in playrecARLas.initStream.'
                };
            errorMsgARLas(errorTxt);
            objPlayrec.objInit.obj.buttonManager(51)
            objPlayrec.printError(ME)
        end
    end
    function writeData(varargin) % get the finished page files and write to disk
        objPlayrec = varargin{1};
        try 
            n = length(objPlayrec.pageFiles); % number of page files file
            isdone = zeros(n,1); % initialize vector showing finished recordings
            for ii=1:n % find which recordings are finished
                isdone(ii,1) = playrec('isFinished',objPlayrec.pageFiles(ii));
            end
            objPlayrec.completedPageFiles = sum(isdone); % total number of completed recordings
            if objPlayrec.completedPageFiles < 1
                %disp('Stuck!')
                return
            elseif objPlayrec.completedPageFiles < 0
                errorTxt = {'  Issue: Unexpected page file error.'
                     '  Action: None.'
                     '  Location: in playrecARLas.writeData.'
                    };
                errorMsgARLas(errorTxt);
                objPlayrec.objInit.obj.buttonManager(51)
            end
            Recording = zeros(objPlayrec.nSamples*objPlayrec.completedPageFiles,objPlayrec.chans_in); % initialize variable for finished buffers
            indx1 = 1;
            indx2 = objPlayrec.nSamples;
            for ii=1:objPlayrec.completedPageFiles % retrieve the finished recordings
                dummy = playrec('getRec',objPlayrec.pageFiles(1)); % get recorded data
                if size(dummy,1) ~= objPlayrec.nSamples
                    if objPlayrec.killRun == 1
                        return
                    else
                        %objPlayrec.killRun
                        disp('Warning in playrecARLas.writeData: unexpected recording length!')
                    end
                end
                Recording(indx1:indx2,:) = dummy;
                playrec('delPage',objPlayrec.pageFiles(1)); % clear the current pageFiles
                objPlayrec.pageFiles(1) = []; % delete finished pageFiles
                indx1 = indx1 + objPlayrec.nSamples; % increment the indices
                indx2 = indx2 + objPlayrec.nSamples; 
            end
        catch  ME
            errorTxt = {'  Issue: Error reading pageFiles.'
                 '  Action: None.'
                 '  Location: in playrecARLas.writeData.'
                };
            errorMsgARLas(errorTxt);
            objPlayrec.objInit.obj.buttonManager(51)
            objPlayrec.printError(ME)
        end
        Recording = double(Recording); % convert from single precision (returned by playrec) to double precision 
        Recording = Recording * objPlayrec.card2volts; % convert from card units to a voltage. Done for all channels
        for kk=1:size(Recording,2)
           Recording(:,kk) = Recording(:,kk) / objPlayrec.micSens(kk);  % apply mic sensitivity
           Recording(:,kk) = Recording(:,kk) / objPlayrec.ampGain(kk);  % apply input amplifier gain
        end
        try
            if isempty(objPlayrec.holdingPen) && objPlayrec.systemDelay > 0 % the holding pen is therefore empty; fill it
                Recording(1:objPlayrec.systemDelay,:) = []; % discard the first part of the recording due to the delay
                if objPlayrec.completedPageFiles == 1
                    indx = 0;
                else
                    indx = ((objPlayrec.completedPageFiles - 1) * objPlayrec.nSamples) + objPlayrec.systemDelay;
                end
                objPlayrec.holdingPen = Recording(indx+1:end,:); % set aside this portion for the next write operation
                Recording(indx+1:end,:) = []; % discard the portion in the holding pen
                objPlayrec.completedPageFiles = objPlayrec.completedPageFiles - 1;
            else % the holding pen is occupied (or checking 0 samples of delay)
                Recording = [objPlayrec.holdingPen;Recording]; % append the portion waiting in the holding pen
                indx = (objPlayrec.completedPageFiles * objPlayrec.nSamples);
                objPlayrec.holdingPen = Recording(indx+1:end,:); % set aside this portion for the next write operation
                Recording(indx+1:end,:) = []; % discard the portion now in the holding pen
            end
        catch ME
            errorTxt = {'  Issue: Unexpected error writing data.'
                 '  Action: None.'
                 '  Location: in playrecARLas.writeData.'
                };
            errorMsgARLas(errorTxt);
            objPlayrec.objInit.obj.buttonManager(51)
            objPlayrec.printError(ME)
        end
        try
            if ~isempty(Recording)
                for ii=1:objPlayrec.chans_in % loop over input channels; write each channel separately
                    fwrite(objPlayrec.fid(ii),Recording(:,ii),'double'); % stream the data to disk
                end
                objPlayrec.completedBuffers = objPlayrec.completedBuffers + objPlayrec.completedPageFiles; % increment the files written to disk counter
            end
        catch ME
            errorTxt = {'  Issue: Error using fwrite.'
                 '  Action: None.'
                 '  Location: in playrecARLas.initStream.'
                };
            errorMsgARLas(errorTxt);
            objPlayrec.objInit.obj.buttonManager(51)
            objPlayrec.printError(ME)
        end
        if ~isempty(Recording)
            try
                newBuffer = zeros(objPlayrec.nSamples,objPlayrec.chans_in); % create mean of current Recording set
                indx1 = 1;
                indx2 = objPlayrec.nSamples;
                objPlayrec.currentWaveform = Recording(indx1:indx2,:);
                if isempty(objPlayrec.currentWaveform)
                    disp('gotcha!')
                    keyboard
                end
                for ii=1:objPlayrec.completedPageFiles
                    newBuffer = newBuffer + Recording(indx1:indx2,:);
                    indx1 = indx1 + objPlayrec.nSamples;
                    indx2 = indx2 + objPlayrec.nSamples;
                end
                newBuffer = newBuffer / objPlayrec.completedPageFiles;
                delta = newBuffer - objPlayrec.mu;
                objPlayrec.mu = objPlayrec.mu + (delta ./ objPlayrec.completedBuffers);
                objPlayrec.doPlot % update the plots
                objPlayrec.completedPageFiles = 0; % reset to zero
            catch ME
                errorTxt = {'  Issue: Error updating running mean.'
                     '  Action: None.'
                     '  Location: in playrecARLas.writeData.'
                    };
                errorMsgARLas(errorTxt);
                objPlayrec.objInit.obj.buttonManager(51)
                objPlayrec.printError(ME)
            end
        end
    end
    function retrieveData(varargin) % get the data from disk
        objPlayrec = varargin{1};
        try
            for ii=1:objPlayrec.chans_in
                objPlayrec.fid(1,ii) = fopen([objPlayrec.binFileName,'_',num2str(ii),'.bin'],'r'); % file id, read-only
            end
        catch ME
            errorTxt = {'  Issue: Error obtaining fid.'
                 '  Action: None.'
                 '  Location: in playrecARLas.retrieveData.'
                };
            errorMsgARLas(errorTxt);
            objPlayrec.objInit.obj.buttonManager(51)
            objPlayrec.printError(ME)
        end
        if any(objPlayrec.fid == -1)
            errorTxt = {'  Issue: Error opening fid file.'
                 '  Action: None.'
                 '  Location: in playrecARLas.retrieveData.'
                };
            errorMsgARLas(errorTxt);
            objPlayrec.objInit.obj.buttonManager(51)
            objPlayrec.printError(ME)
        end
        try
            for ii=1:objPlayrec.chans_in
                status = fseek(objPlayrec.fid(1,ii),0,'bof');
                X = fread(objPlayrec.fid(1,ii),inf,'double');
                expectedLength = objPlayrec.nReps * objPlayrec.nSamples;
                %disp(['Extra Samples = ',num2str(length(X)-expectedLength)])
                X = X(1:expectedLength);
                objPlayrec.Data(:,ii) = X;
                %= reshape(X,objPlayrec.nSamples,objPlayrec.nReps);
            end
            objPlayrec.cleanUp
        catch ME
            errorTxt = {'  Issue: Error getting fid.'
                 '  Action: None.'
                 '  Location: in playrecARLas.retrieveData.'
                };
            errorMsgARLas(errorTxt);
            objPlayrec.objInit.obj.buttonManager(51)
            objPlayrec.printError(ME)
        end        
    end
    function doPlot(varargin) % update the waveform plots
        objPlayrec = varargin{1};
        try
            indx = objPlayrec.recChanList(objPlayrec.CHAN_in.Value); % get the index of the currently-selected input channel
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
                decPlaces = 9; % the best answer is smalle than this, but this is as small as we want to go
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
                objPlayrec.completedBuffers = objPlayrec.nReps;
            end
            set(objPlayrec.AVG,'String',[num2str(objPlayrec.completedBuffers),' of ',num2str(objPlayrec.nReps)])
            pause(0.001)
        catch ME
            errorTxt = {'  Issue: Error plotting data.'
                 '  Action: None.'
                 '  Location: in playrecARLas.doPlot.'
                };
            errorMsgARLas(errorTxt);
            objPlayrec.objInit.obj.buttonManager(51)
            objPlayrec.printError(ME)
        end
    end
    function saveData(varargin) % save the incremental data
        objPlayrec = varargin{1};
        try
            basePath = objPlayrec.map.data;
            expName = objPlayrec.objInit.obj.experimentID;
            subjName = objPlayrec.objInit.obj.subjectID;
            pathName = [basePath,expName,objPlayrec.sep,subjName,objPlayrec.sep];
            if exist(pathName,'dir') ~= 7 % if directory does not exist
                errorTxt = {'  Issue: Expected data file does not exist.'
                         '  Action: aborting data save operation.'
                         '  Location: in playrecARLas.saveData.'
                        };
                errorMsgARLas(errorTxt);
                objPlayrec.objInit.obj.buttonManager(51)                
                return
            end
            
            Data = objPlayrec.Data; % read in all the data
            for ii=1:objPlayrec.chans_in
                data = Data(:,ii); % each column is one channel
                % reshape into a matrix, with buffers in a column
                nSamples = length(data) / objPlayrec.nReps;
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
                                
                label = objPlayrec.objInit.label.(matlab.lang.makeValidName(['Ch',num2str(ii)])); 
                fileName = ['Ch',num2str(ii),'_',label];
                % check to make sure not overwriting another file with the same name
                ok = 0;
                counter = 1;
                while ~ok
                    newFileName = [fileName,'_',num2str(counter),'.mat'];
                    if exist([pathName,newFileName],'file') == 2 % if file already exists
                        counter = counter + 1;
                    else
                        ok = 1;
                        fileName = newFileName;
                    end
                end
                % header is a structure with important info about the data
                header.timeStamp = datestr(clock, 0); % date and time when this file was created
                header.subjID = objPlayrec.objInit.obj.subjectID;
                header.expID = objPlayrec.objInit.obj.experimentID;
                header.operator = objPlayrec.objInit.obj.operatorID; % person who collected the data
                header.arlasVersion = objPlayrec.objInit.obj.arlasVersion;
                header.fs = objPlayrec.fs; % sampling rate in Hz
                header.systemDelay = objPlayrec.systemDelay;
                header.card2volts = objPlayrec.card2volts;
                header.label = objPlayrec.label{ii};
                header.ampGain = objPlayrec.ampGain(ii);
                header.micSens = objPlayrec.micSens(ii);
                header.nSamples = objPlayrec.nSamples;
                header.nReps = objPlayrec.nReps;
                header.completedBuffers = objPlayrec.completedBuffers;
                header.mu = objPlayrec.mu;
                header.time = objPlayrec.time;
                header.userInfo = objPlayrec.userInfo;
                
                save([pathName,fileName],'data','header')
                objPlayrec.savedFiles{ii} = fileName;
                objPlayrec.savedFilesPath = pathName;
            end
        catch ME
            errorTxt = {'  Issue: Error saving data.'
                 '  Action: Data not saved.'
                 '  Location: in playrecARLas.saveData.'
                };
            errorMsgARLas(errorTxt);
            objPlayrec.objInit.obj.buttonManager(51)
            objPlayrec.printError(ME)
        end
    end
    
    function cleanUp(~) % close all open file identifiers
        fids = fopen('all'); % close all files open for read or write access
        if ~isempty(fids)
            for ii=1:length(fids)
                fclose(fids(ii));
            end
        end
    end  
    function printError(varargin)
        objPlayrec = varargin{1};
        ME = varargin{2};
        disp(ME.identifier)
        disp(ME.message)
        disp(ME.stack(1).line)            
    end
    
end
end