classdef initARLas < handle
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Class definition: initARLas
% For use with ARLas (Auditory Research Laboratory auditory software)
%
% Auditory Research Lab, The University of Iowa
% Deptartment of Communication Sciences & Disorders
% The University of Iowas
% Author: Shawn S. Goodman, PhD
% Date: August 26, 2016
% Last Updated: July 20, 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
properties (SetAccess = private)
    arlasVersion = '2017.07.20';
    obj         % passed from arlas
    initPath    % passed from arlas
    initFile    % passed from arlas
    stimTrain   % created internally
    H % ----- ----- handles for graphical user interface ----- -----
    VIEW
     SYSTEM % SYSTEM Panel -----
        OStxt % operating system label
        OS    % operating system type
        HOSTtxt % HOST API, text label
        HOST  % HOST API, popup menu
        SAMPLINGtxt % sampling rate label
        SAMPLING    % currently selected sampling rate
        DELAYtxt    % system delay_now, text label
        DELAY       % system delay_now (playback and record delay_now) for currently selected setup
        CARD2Vtxt
        CARD2V
     INPUT % INPUT panel -----
        IDtxt_in % input device idenifying numbers, text label
        ID_in % identifying numbers for possible input devices (numeric), popup menu
        NAMEtxt_in
        NAME_in
        CHANStxt_in
        CHANS_in        
     OUTPUT % OUTPUT Panel -----
        IDtxt_out % output device idenifying numbers, text label
        ID_out % identifying numbers for possible output devices (numeric), popup menu
        NAMEtxt_out % output device name, text label
        NAME_out % name of currently selected output device, text label
        CHANStxt_out % number of possible output channels, text label
        CHANS_out % number of channels availabel for currently selected device, text label
     % User Buttons -----
     GETDELAY
     GETC2V
     SAVECONFIG % save current configuration
end
properties (SetAccess = public)
    % SYSTEM Data Structures
    usingSavedValues % whether or not previously-saved values are being used
    devs            % list of devices and properties provided by initialization call to playrec
    os              % the operating system currently being used
    fs              % list of possible sampling rates for current selections
    fs_now          % currently chosen sample rate
    delay_now       % system delay (playback and record) for current setup
    card2volts_now  % multiplier for conversion from sound card to actual volts
    
    HostAPI         % list of all output HOST APIs (as cell array)
    hostAPI         % short list of unique output HOST APIs
    hostAPI_now     % currently chosen HOST API
    hostAPIindices  % indices of chosen hostAPI_now in the large list HostAPI)

    % INPUT Data Structures
    deviceIndx_in   % full list of indices of output devices
    indx_in_now     % value of curently chosen device
    deviceID_in     % list of all input device ID options (numeric)
    id_in_now       % currently chosen value
    name_in_now     % currently chosen name
    chans_in        % possible number of channels to use
    
    % OUTPUT Data Structures
    deviceIndx_out % list of indices of output devices, given host API
    indx_out_now   % curently chosen output device
    deviceID_out   % list of output device IDs, given host API
    id_out_now     % currently chosen ID
    name_out_now   % currently chosen device name
    chans_out      % possible number of channels to use
end
methods
    function objInit = initARLas(obj) % initialize object of class initARLas
        objInit.obj = obj; % object of class arlas
        if obj.ping == 1
            return
        end
        objInit.initFile = obj.initFile;
        objInit.initPath = obj.initPath;
        objInit.usingSavedValues = 0; % defaults to not using saved values
        
        initText = 'Unknown';
        objInit.os = computer;
        objInit.devs = playrec('getDevices'); % read in devices
        objInit.fs = [44100,48000,64000,88200,96000,128000,176000,192000,200000]; % possible sampling rates
        objInit.fs_now = initText;
        objInit.delay_now = 0; % default system delay_now is 0 samples
        objInit.card2volts_now = 1; % default card to volts is a multiplier of 1  (test value: 9.949;)
        %
        objInit.HostAPI = initText;
        objInit.hostAPI = initText;
        objInit.hostAPI_now = initText;
        objInit.hostAPIindices = initText;
        objInit.hostAPI_now = initText;
        %
        objInit.deviceIndx_in = 1;
        objInit.indx_in_now = 1;
        objInit.deviceID_in = initText;
        objInit.id_in_now = initText;
        objInit.name_in_now = initText;
        objInit.chans_in = initText;
        %objInit.chans_in_now = initText;
        %
        objInit.deviceIndx_out = 1;
        objInit.indx_out_now = 1;
        objInit.deviceID_out = initText;
        objInit.id_out_now = initText;
        objInit.name_out_now = initText;
        objInit.chans_out = initText;
        %objInit.chans_out_now = initText;
        %
        objInit.H = obj.H;
        objInit.VIEW = obj.VIEW;
    end  
    function abort(varargin) % instructions for aborting when gui closed
        delete(objInit);
    end
    function initGui(varargin) % initialize gui
        objInit = varargin{1};
        % create panels within figure ([left, bottom,width, height])
        objInit.SYSTEM = uipanel('Parent',objInit.VIEW,'Title','SYSTEM','FontSize',12,...
            'BackgroundColor','white','Units','Normalized',...
            'Position',[.02 .445 .96 .27]); % .425
        objInit.INPUT = uipanel('Parent',objInit.VIEW,'Title','INPUT','FontSize',12,...
            'BackgroundColor','white','Units','Normalized',...
            'Position',[.02 .15 .46 .27]); % .25 height
        objInit.OUTPUT = uipanel('Parent',objInit.VIEW,'Title','OUTPUT','FontSize',12,...
            'BackgroundColor','white','Units','Normalized',...
            'Position',[.52 .15 .46 .27]); 
        % SYSTEM ----- populate system panel
        spacing = .25;
        top = .8;
        height = 0.17; % .15
        objInit.OStxt = uicontrol('Parent',objInit.SYSTEM,'Style','text','BackgroundColor','white',...
            'String','Computer:         ','HorizontalAlignment','Left','Units','Normalized',...
            'ToolTipString','Type of computer on which MATLAB is executing.',...
            'Position',[.01 top-(.3*spacing) .5 height],'FontSize',10);
        objInit.OS = uicontrol('Parent',objInit.SYSTEM,'Style','text',...
            'String',objInit.os,'FontSize',10,'BackgroundColor','white',...
            'HorizontalAlignment','Left','Units','Normalized',...
            'Position',[.2 top-(.3*spacing) .5 height],'Value',2);
        objInit.HOSTtxt = uicontrol('Parent',objInit.SYSTEM,'Style','text',...
            'BackgroundColor','white','String','Host API:         ','HorizontalAlignment','Left',...
            'ToolTipString','Application Programming Interface.','Units','Normalized',...
            'Position',[.01 top-(1.5*spacing) .5 height],'FontSize',10);
        objInit.HOST = uicontrol('Parent',objInit.SYSTEM,'Style','popup',...
            'String',objInit.hostAPI,'FontSize',10,'BackgroundColor','white',...
            'HorizontalAlignment','Left','Units','Normalized',...
            'Position',[.2 top-(1.4*spacing) .25 height],'Value',objInit.indx_out_now,...
            'Callback',{@objInit.hostManager,objInit.devs});
        objInit.SAMPLINGtxt = uicontrol('Parent',objInit.SYSTEM,'Style','text','BackgroundColor','white',...
             'String','Sampling Rate:    ','HorizontalAlignment','Left',...
             'ToolTipString','Sampling frequency in Hz.','Units','Normalized',...
             'Position',[.01 top-(2.7*spacing) .5 height],'FontSize',10);
        objInit.SAMPLING = uicontrol('Parent',objInit.SYSTEM,'Style','popup',...
             'String',objInit.fs,'FontSize',10,'Units','Normalized',...
             'Position',[.2 top-(2.7*spacing)+.025 .25 height],'Value',1,...
             'Callback',{@objInit.samplingManager,objInit.devs});
         objInit.DELAYtxt = uicontrol('Parent',objInit.SYSTEM,'Style','text','BackgroundColor','white',...
             'String','System Delay:    ','HorizontalAlignment','Left','Units','Normalized',...
             'Position',[.5 top-(1.5*spacing) .5 height],'FontSize',10,'ToolTipString','Sample delay between playback and recording.');
        objInit.DELAY = uicontrol('Parent',objInit.SYSTEM,'Style','edit',...
             'String',objInit.delay_now,'FontSize',10,'BackgroundColor','white','HorizontalAlignment',...
             'Left','Callback',@objInit.delayManager,'Units','Normalized',...
             'Position',[.7 top-(1.6*spacing) .2 .18]);
         objInit.CARD2Vtxt = uicontrol('Parent',objInit.SYSTEM,'Style','text','BackgroundColor','white',...
             'String','Card to Volts:    ','HorizontalAlignment','Left',...
             'ToolTipString','Multiplier to convert sound card units to voltage.','Units','Normalized',....
             'Position',[.5 top-(2.7*spacing) .5 height],'FontSize',10);
        objInit.CARD2V = uicontrol('Parent',objInit.SYSTEM,'Style','edit',...
            'String',objInit.card2volts_now,'FontSize',10,'BackgroundColor','white',...
            'HorizontalAlignment','Left','Callback',@objInit.c2vManager,'Units','Normalized',...
            'Position',[.7 top-(2.8*spacing) .2 .18]);
        % OUTPUT ----- populate output panel
        spacing = .25;
        top = .8;
        height = 0.15;
        objInit.IDtxt_out = uicontrol('Parent',objInit.OUTPUT,'Style','text',...
            'BackgroundColor','white','String','Device ID:        ',...
            'ToolTipString','ID number used to refer to outupt devices.',...
            'HorizontalAlignment','Left','Units','Normalized',...
            'Position',[.01 top-(.3*spacing) .5 height],'FontSize',10); % [.01 .8 .5 .1]
        objInit.ID_out = uicontrol('Parent',objInit.OUTPUT,'Style','popup',...
            'String',objInit.id_out_now,'FontSize',10,'BackgroundColor',...
            'white','HorizontalAlignment','Left','Units','Normalized',...
            'Position',[.5 top-(.10*spacing) .2 height],'Value',objInit.indx_out_now,...
            'Callback',{@objInit.outputIdManager,objInit.devs});
        objInit.NAMEtxt_out = uicontrol('Parent',objInit.OUTPUT,'Style',...
            'text','BackgroundColor','white','String','Device Name:      ',...
            'ToolTipString','Name of currently selected output device.',...
            'HorizontalAlignment','Left','Units','Normalized',...
            'Position',[.01 top-(1.5*spacing) .5 height],'FontSize',10);       
        objInit.NAME_out = uicontrol('Parent',objInit.OUTPUT,'Style','text',...
            'String',objInit.name_out_now,'HorizontalAlignment','Left',...
            'FontSize',10,'BackgroundColor','white','Units','Normalized',...
            'Position',[.5 top-(1.5*spacing) .5 height]);
        objInit.CHANStxt_out = uicontrol('Parent',objInit.OUTPUT,'Style','text',...
            'BackgroundColor','white','String','Num. of Channels:       ',...
            'ToolTipString','Number of desired output channels (up to maximum supported by the device).',... 
            'HorizontalAlignment','Left','Units','Normalized',...
            'Position',[.01 top-(2.7*spacing) .5 height],'FontSize',10);
        objInit.CHANS_out = uicontrol('Parent',objInit.OUTPUT,'Style','text',...
            'String',objInit.chans_out,...
            'FontSize',10,'BackgroundColor','white','HorizontalAlignment','Left',...
            'Units','Normalized','Position',[.5 top-(2.7*spacing) .47 height]);
        % INPUT ----- populate input panel
        objInit.IDtxt_in = uicontrol('Parent',objInit.INPUT,'Style','text','BackgroundColor','white',...
            'ToolTipString','ID number used to refer to input devices.',...
            'String','Device ID:        ','HorizontalAlignment','Left','Units','Normalized',...
            'Position',[.01 top-(.3*spacing) .5 height],'FontSize',10);
        objInit.ID_in = uicontrol('Parent',objInit.INPUT,'Style','popup',...
             'String',objInit.id_in_now,'FontSize',10,'BackgroundColor','white',...
             'HorizontalAlignment','Left','Units','Normalized',...
             'Position',[.5 top-(.10*spacing) .2 height],'Value',objInit.indx_in_now,...
             'Callback',{@objInit.inputIdManager,objInit.devs});
        objInit.NAMEtxt_in = uicontrol('Parent',objInit.INPUT,'Style','text','BackgroundColor','white',...
             'ToolTipString','Name of currently selected input device.',...
             'String','Device Name:      ','HorizontalAlignment','Left',...
             'Units','Normalized','Position',[.01 top-(1.5*spacing) .5 height],'FontSize',10);       
        objInit.NAME_in = uicontrol('Parent',objInit.INPUT,'Style','text',...
             'String',objInit.name_in_now,'FontSize',10,'BackgroundColor','white',...
             'HorizontalAlignment','Left',...
             'Units','Normalized','Position',[.5 top-(1.5*spacing) .5 height],'Value',2);
        objInit.CHANStxt_in = uicontrol('Parent',objInit.INPUT,'Style','text','BackgroundColor','white',...
            'ToolTipString','Number of desired input channels (up to maximum supported by the device).',...  
            'String','Num. of Channels:       ','HorizontalAlignment','Left',...
             'Units','Normalized','Position',[.01 top-(2.7*spacing) .5 height],'FontSize',10);
        objInit.CHANS_in = uicontrol('Parent',objInit.INPUT,'Style','text',...
             'String',objInit.chans_in,'FontSize',10,'BackgroundColor','white',...
             'HorizontalAlignment','Left','Units','Normalized','Value',1,...
             'Position',[.5 top-(2.7*spacing) .47 height]);         
        objInit.GETDELAY = uicontrol('Parent',objInit.VIEW,'String','Get System Delay',...
            'Interruptible','off','BusyAction','cancel','Units','Normalized',...
            'ToolTipString','Measure system delay. Connect channel 1 output to channel 1 input, then press button.',... 
            'Position',[.02 .015 .3 .1],'Callback',@objInit.getDelay);
        objInit.GETC2V = uicontrol('Parent',objInit.VIEW,'String','Get Card to Volts',...
            'Interruptible','off','BusyAction','cancel','Units','Normalized',...
            'ToolTipString','Find conversion from card units to volts. use a multimeter and direct out/in connection.',... 
            'Position',[.35 .015 .3 .1],'Callback',@objInit.getCard2V);
        objInit.SAVECONFIG = uicontrol('Parent',objInit.VIEW,'String','Save Configuration',...
            'Interruptible','off','BusyAction','cancel','Units','Normalized',...
            'ToolTipString','Save current configuration for future use.',...             
            'Position',[.68 .015  .3 .1],'Callback',@objInit.saveConfig);
    end
    function makeVisible(varargin)
        objInit = varargin{1};
        try
            set(objInit.SYSTEM,'Visible','on')
            set(objInit.OUTPUT,'Visible','on')
            set(objInit.INPUT,'Visible','on')        
            set(objInit.GETDELAY,'Visible','on')
            set(objInit.GETC2V,'Visible','on')
            set(objInit.SAVECONFIG,'Visible','on')
        catch
        end
    end
    function makeInvisible(varargin)
        objInit = varargin{1};
        try
            set(objInit.SYSTEM,'Visible','off')
            set(objInit.OUTPUT,'Visible','off')
            set(objInit.INPUT,'Visible','off')        
            set(objInit.GETDELAY,'Visible','off')
            set(objInit.GETC2V,'Visible','off')
            set(objInit.SAVECONFIG,'Visible','off')
        catch
        end
    end
    function killPlots(varargin)
        %objInit = varargin{1};
        try
            a = findall(gcf);
            counter = 1;
            for ii=1:size(a,1)
                try
                if strcmp(a(ii).Parent.Title,'VIEW: Initialize')
                    indx(counter,1) = ii;
                    counter = counter + 1;
                end
                catch
                end
            end
            for ii=1:size(indx,1)
                delete(a(indx(ii)))
            end
        catch ME
%             errorTxt = {'  Issue: Error deleting init plots.'
%                  '  Action: None.'
%                  '  Location: in initARLas.killPlots.'
%                 };
%             errorMsgARLas(errorTxt);
%             objInit.printError(ME)
        end
    end
    
    function trySavedValues(varargin) % try to use saved values
        objInit = varargin{1};
        dummy = varargin{2};
        default = dummy.default;
        currentDevs = playrec('getDevices'); % get the current system setup
        % find the current locations of the desired host API
        N = size(currentDevs,2);
        match = zeros(N,1);
        for ii=1:N
            q = currentDevs(ii).hostAPI;
            match(ii,1) = strcmp(q,strtrim(default.hostAPI_now));
        end
        indxAPI = find(match); % indices of desired host API
        for ii=1:N
            q = currentDevs(ii).name;
            match(ii,1) = strcmp(q,strtrim(default.name_out_now));
        end
        indxOutName = find(match); % indices of desired device name
        for ii=1:N
            q = currentDevs(ii).name;
            match(ii,1) = strcmp(q,strtrim(default.name_in_now));
        end
        indxInName = find(match); % indices of desired device name
        % the current values to use are the intersection of host API and device name
        indx1 = intersect(indxAPI,indxOutName);
        indx2 = intersect(indxAPI,indxInName);
        indx = intersect(indx1,indx2);
        indx = indx(1);
        objInit.hostAPI_now = default.hostAPI_now; % this stays the same
        objInit.HostAPI = currentDevs(indx).name; % get the current values
        objInit.hostAPI = default.hostAPI; % this stays the same
        objInit.hostAPIindices = indxAPI; % the new indices
        objInit.deviceIndx_out = indxAPI;
        objInit.indx_out_now = indxOutName;
        for ii=1:length(indxAPI)
            jar(ii,1) = currentDevs(indxAPI(ii)).deviceID;
        end
        objInit.deviceID_out = jar;
        objInit.id_out_now = currentDevs(indxOutName).deviceID;
        objInit.name_out_now = default.name_out_now;
        objInit.chans_out = currentDevs(indxOutName).outputChans;
        %objInit.chans_out_now = objInit.chans_out;
        objInit.deviceIndx_in = indxAPI;
        objInit.indx_in_now = indxInName;
        for ii=1:length(indxAPI)
            jar(ii,1) = currentDevs(indxAPI(ii)).deviceID;
        end
        objInit.deviceID_in = jar;                    
        objInit.id_in_now = currentDevs(indxInName).deviceID;
        objInit.name_in_now = default.name_in_now;
        objInit.chans_in = currentDevs(indxInName).inputChans;
        objInit.os = default.os;
        objInit.fs = default.fs;
        objInit.fs_now = default.fs_now;
        objInit.delay_now = default.delay_now;
        objInit.card2volts_now = default.card2volts_now;
        objInit.devs = currentDevs;
        %objInit.initChannels
        objInit.usingSavedValues = 1;
        objInit.updateGui % update the gui with these new values
    end
    function discover(varargin) % discover current system configuration
        objInit = varargin{1};
        try
            objInit.devs = playrec('getDevices'); % read in devices
            N = size(objInit.devs,2); % total number of recognized sound devices on the machine
            objInit.HostAPI = [];
            for ii=1:N
               dummy = objInit.devs(ii).hostAPI;
               objInit.HostAPI = char(objInit.HostAPI,dummy);
            end
            [~,indx] = unique(char(objInit.HostAPI),'rows');
            % there is a known issue using playrec with Windows WDM-KS: it
            % causes playrect to freeze when doing a reset. Therefore, this
            % host API is disallowed in the following lines 
            %-----
            indx2 = [];
            for ii=1:length(indx)
                dummy = objInit.HostAPI(indx(ii),:);
                if strcmp(deblank(dummy),'Windows WDM-KS') == 0
                    indx2 = [indx2;indx(ii)];
                end        
            end
            indx = indx2;
            %-----        
            objInit.hostAPI =  objInit.HostAPI(indx,:); % list of unique HOST APIs
            objInit.hostAPI(1,:) = []; % first row is empty; discard
            objInit.usingSavedValues = 0; % now we are using new,unsaved values.
            objInit.hostManager
        catch ME
            errorTxt = {'  Issue: Unexpected error discovering new system devices.'
                 '  Action: None.'
                 '  Location: initARLas.discover.'
                };
            errorMsgARLas(errorTxt);
            objInit.obj.buttonManager(21)
            objInit.printError(ME)
        end
    end
    function hostManager(varargin) % manage selection of Host API
        % depending on which particular HOST API is selected, set other allowable options.
        objInit = varargin{1}; % get the object
        try
            try
                value = varargin{2}.Value; % get which API name was selected
            catch
                value = 1; % if none found, default to 1
            end
            set(objInit.HOST,'Value',value)
            objInit.hostAPI_now = objInit.hostAPI(value,:); % currently chosen value
            N = length(objInit.devs);
            indxHost = zeros(N,1);
            counter = 0;
            for ii=1:N % find indices of devices that match the selected host API
                counter = counter + 1;
                indxHost(ii,1) = strcmp(objInit.devs(1,ii).hostAPI,deblank(objInit.hostAPI_now));
            end
            indxHost = indxHost(1:counter,:);
            I = find(indxHost);
            n = length(I);
            objInit.deviceIndx_out = [];
            objInit.deviceID_out = [];
            objInit.deviceIndx_in = [];
            objInit.deviceID_in = [];        
            counter_out = 1;
            counter_in = 1;
            for ii=1:n 
                if objInit.devs(I(ii)).outputChans > 0; % find allowable outputs
                    objInit.deviceIndx_out(counter_out,1) = I(ii); % the index of output devices
                    objInit.deviceID_out(counter_out,1) = objInit.devs(I(ii)).deviceID; % the ID number of output devices
                    counter_out = counter_out + 1;
                end
                if objInit.devs(I(ii)).inputChans > 0; % find allowable inputs
                    objInit.deviceIndx_in(counter_in,1) = I(ii); % the index of output devices
                    objInit.deviceID_in(counter_in,1) = objInit.devs(I(ii)).deviceID; % the ID number of output devices
                    counter_in = counter_in + 1;
                end
            end
            value = 1; % default the new possible output and input IDs to first on the list
            set(objInit.ID_out,'Value',value)
            set(objInit.ID_in,'Value',value)
            objInit.indx_out_now = objInit.deviceIndx_out(value);
            objInit.indx_in_now = objInit.deviceIndx_in(value);
            % update associated values
            objInit.id_out_now = objInit.devs(objInit.indx_out_now).deviceID;
            objInit.name_out_now = objInit.devs(objInit.indx_out_now).name;
            objInit.chans_out = objInit.devs(objInit.indx_out_now).outputChans;
            %objInit.chans_out_now = objInit.chans_out;
            objInit.fs_now = objInit.devs(objInit.indx_out_now).defaultSampleRate;
            [~,value] = min(abs(objInit.fs_now-objInit.fs));
            set(objInit.SAMPLING,'Value',value)
            objInit.id_in_now = objInit.devs(objInit.indx_in_now).deviceID;
            objInit.name_in_now = objInit.devs(objInit.indx_in_now).name;
            objInit.chans_in = objInit.devs(objInit.indx_in_now).inputChans; 
            %objInit.chans_in_now = objInit.chans_in;
            %objInit.initChannels % initialize channel values
            objInit.samplingManager
        catch ME
            errorTxt = {'  Issue: Unexpected error selecting Host API.'
                 '  Action: None.'
                 '  Location: initARLas.hostManager.'
                };
            errorMsgARLas(errorTxt);
            objInit.obj.buttonManager(21)
            objInit.printError(ME)
        end
    end
    function samplingManager(varargin) % manage selection of Sampling Rate
        % check to see if selected sampling rate works
        objInit = varargin{1}; % get the object
        try
            objInit.fs_now = objInit.fs(get(objInit.SAMPLING,'value')); % get the current sampling rate value
            if playrec('isInitialised')
                playrec('reset')
            end
            try % try initializing with current values
                playrec('init',objInit.fs_now,objInit.id_out_now,objInit.id_in_now)
                pause(.01)
                playrec('reset')
                pause(.01)
                set(objInit.SAMPLING,'BackgroundColor',[1 1 1])
                set(objInit.ID_out,'BackgroundColor',[1 1 1])
                set(objInit.ID_in,'BackgroundColor',[1 1 1])
            catch ME
               errorTxt = {'  Issue: selected sampling rate is not compatible with the other selections.'
                     '  Suggested Action: Unknown.'
                     '  Location: initARLas.samplingManager.'
                    };
                warnMsgARLas(errorTxt);
                set(objInit.SAMPLING,'BackgroundColor',[1 0 0])
                objInit.obj.buttonManager(21)
                objInit.printError(ME)
            end
            objInit.updateGui
        catch ME
            errorTxt = {'  Issue: Unexpected error selecting Sampling Rate.'
                 '  Action: None.'
                 '  Location: initARLas.samplingManager.'
                };
            errorMsgARLas(errorTxt);
            objInit.obj.buttonManager(21)
            objInit.printError(ME)            
        end
    end
    function inputIdManager(varargin) % manage Device ID for input
        objInit = varargin{1};
        try
            try
                value = varargin{2}.Value; % get the index value from the popup menu
            catch
                value = 1; % if none found, default to 1
            end        
            objInit.indx_in_now = objInit.deviceIndx_in(value);
            objInit.id_in_now = objInit.devs(objInit.indx_in_now).deviceID;
            objInit.name_in_now = objInit.devs(objInit.indx_in_now).name;
            objInit.chans_in = objInit.devs(objInit.indx_in_now).inputChans;
            %objInit.chans_in_now = objInit.chans_in;
            objInit.ID_in.Value = value;
            set(objInit.NAME_in,'String',objInit.name_in_now);
            set(objInit.CHANS_in,'String',num2str(objInit.chans_in))
            % also change output channel to match
            objInit.indx_out_now = objInit.deviceIndx_out(value); 
            objInit.id_out_now = objInit.devs(objInit.indx_out_now).deviceID;
            objInit.name_out_now = objInit.devs(objInit.indx_out_now).name;
            objInit.chans_out = objInit.devs(objInit.indx_out_now).outputChans;
            %objInit.chans_out_now = objInit.chans_out;
            objInit.ID_out.Value = value;       
            set(objInit.NAME_out,'String',objInit.name_out_now);
            set(objInit.CHANS_out,'String',num2str(objInit.chans_out))
            objInit.initChannels % initialize channel values
            objInit.samplingManager
        catch ME
            errorTxt = {'  Issue: Unexpected error selecting Input Device ID.'
                 '  Action: None.'
                 '  Location: initARLas.inputIdManager.'
                };
            errorMsgARLas(errorTxt);
            objInit.obj.buttonManager(21)
            objInit.printError(ME)
        end
    end
    function outputIdManager(varargin) % manage Device ID for output
        objInit = varargin{1}; % get the object
        try
            value = varargin{2}.Value; % get the index value from the popup menu
        catch
            value = 1; % if none found, default to 1
        end
        objInit.indx_out_now = objInit.deviceIndx_out(value); 
        objInit.id_out_now = objInit.devs(objInit.indx_out_now).deviceID;
        objInit.name_out_now = objInit.devs(objInit.indx_out_now).name;
        objInit.chans_out = objInit.devs(objInit.indx_out_now).outputChans; 
        %objInit.chans_out_now = objInit.chans_out;
        objInit.ID_out.Value = value;
        set(objInit.NAME_out,'String',objInit.name_out_now);
        set(objInit.CHANS_out,'String',num2str(objInit.chans_out))
        % also change input channel to match
        objInit.indx_in_now = objInit.deviceIndx_in(value); 
        objInit.id_in_now = objInit.devs(objInit.indx_in_now).deviceID;
        objInit.name_in_now = objInit.devs(objInit.indx_in_now).name;
        objInit.chans_in = objInit.devs(objInit.indx_in_now).inputChans;
        %objInit.chans_in_now = objInit.chans_in;
        objInit.ID_in.Value = value;     
        set(objInit.NAME_in,'String',objInit.name_in_now);
        set(objInit.CHANS_in,'String',num2str(objInit.chans_in))
        objInit.fs_now = objInit.devs(objInit.indx_out_now).defaultSampleRate;
        [~,value] = min(abs(objInit.fs_now-objInit.fs));
        set(objInit.SAMPLING,'Value',value)
        objInit.initChannels % initialize channel values
        objInit.samplingManager
    end
    
    function delayManager(varargin) % manage setting of system delay
        objInit = varargin{1}; % get the object
        value = objInit.DELAY.String; % get the current value
        ok = 1; % assume the input is okay, unless an error is found
        try % check the current input value
            try 
                value = str2num(value); 
            catch
            end
            if isempty(value) % do not allow empty values
                ok = 0;
                errorTxt = {'  Issue: System Delay cannot be empty.'
                     '  Action: Defaulting to a value of 0.'
                     '  Location: initARLas.delayManager.'
                    };
                errorMsgARLas(errorTxt);
                objInit.obj.buttonManager(21)
                value = 0;
            elseif isa(value,'numeric') == 0 % do not allow non-numeric values
                ok = 0;
                errorTxt = {'  Issue: System Delay must be a numeric value.'
                     '  Action: Defaulting to a value of 0.'
                     '  Location: initARLas.delayManager.'
                    };
                errorMsgARLas(errorTxt);
                objInit.obj.buttonManager(21)
                value = 0;
            elseif value < 0 % allow only positive numeric values
                ok = 0;
                errorTxt = {'  Issue: System Delay must be non-negative.'
                     '  Action: Defaulting to a value of 0.'
                     '  Location: initARLas.delayManager.'
                    };
                errorMsgARLas(errorTxt);
                objInit.obj.buttonManager(21)
                value = 0;
            elseif mod(value,1) ~= 0 % allow only integer values
                ok = 0;
                errorTxt = {'  Issue: System Delay must be an integer.'
                     '  Action: Forcing to an integer value.'
                     '  Location: initARLas.delayManager.'
                    };
                errorMsgARLas(errorTxt);
                objInit.obj.buttonManager(21)
                value = round(value);
            end
            objInit.delay_now = value;
            objInit.DELAY.String = objInit.delay_now;
            if ok == 0
                objInit.DELAY.BackgroundColor = [1 0 0]; % make background red
            else
                objInit.DELAY.BackgroundColor = [1 1 1]; % make background white
            end
        catch ME
            errorTxt = {'  Issue: Unexpected error setting System Delay value.'
                 '  Action: None.'
                 '  Location: initARLas.delayManager.'
                };
            errorMsgARLas(errorTxt);
            objInit.obj.buttonManager(21)
            objInit.printError(ME)
        end
    end
    function c2vManager(varargin) % manage conversion from card units to voltage
        objInit = varargin{1}; % get the object
        value = objInit.CARD2V.String; % get the current value
        ok = 1; % assume the input is okay, unless an error is found
        try % check the current input value
            try 
                value = str2num(value); 
            catch
            end
            if isempty(value) % do not allow empty values
                ok = 0;
                errorTxt = {'  Issue: Card to Volts cannot be empty.'
                     '  Action: Defaulting to a value of 1.'
                     '  Location: initARLas.c2vManager.'
                    };
                errorMsgARLas(errorTxt);
                objInit.obj.buttonManager(21)
                value = 1;
            elseif isa(value,'numeric') == 0 % do not allow non-numeric values
                ok = 0;
                errorTxt = {'  Issue: Card to Volts must be a numeric value.'
                     '  Action: Defaulting to a value of 1.'
                     '  Location: initARLas.c2vManager.'
                    };
                errorMsgARLas(errorTxt);
                objInit.obj.buttonManager(21)
                value = 1;
            elseif value <= 0 % allow only positive numeric values
                ok = 0;
                errorTxt = {'  Issue: Card to Volts must be > 0.'
                     '  Action: Defaulting to a value of 1.'
                     '  Location: initARLas.c2vManager.'
                    };
                errorMsgARLas(errorTxt);
                objInit.obj.buttonManager(21)
                value = 1;
            end
            objInit.card2volts_now = value;
            objInit.CARD2V.String = objInit.card2volts_now;
            if ok == 0
                objInit.CARD2V.BackgroundColor = [1 0 0]; % make background red
            else
                objInit.CARD2V.BackgroundColor = [1 1 1]; % make background white
            end
        catch ME
            errorTxt = {'  Issue: Unexpected error setting Card to Volts value.'
                 '  Action: None.'
                 '  Location: initARLas.c2vManager.'
                };
            errorMsgARLas(errorTxt);
            objInit.obj.buttonManager(21)
            objInit.printError(ME)
        end
    end
    
    function getCard2V(varargin) % calculate multiplier to convert card units to volts
        objInit = varargin{1};
        txt = ({'Calculate your Card to Voltage multiplier by running the experiment file ARLas_getCard2Volts.m.';'';...
                'Click on the Load Button and select the file.';'';...
                'Before you run, be sure to update the first two lines of code to use the correct input and output channels on your sound card.'});
        choice = msgbox(txt,'Get Card to Volts');            
    end
    function getDelay(varargin) % try running the currently-chosen values to see if they work
        objInit = varargin{1};
        txt = ({'Calculate your system delay by running the experiment file ARLas_getDelay.m.';'';...
                'Click on the Load Button and select the file.';'';...
                'Before you run, be sure to update the first two lines of code to use the correct input and output channels on your sound card.'});
        choice = msgbox(txt,'Get System Delay');     
    end
    function saveConfig(varargin)
        objInit = varargin{1};
        try
            [fileName,pathName] = uiputfile([objInit.initPath,'newInitialization.mat'],'Save initialization settings.');
            if fileName == 0 % user chose not to create a file
                dlgtxt = {'You chose NOT to save the current system values.' ' ' 'Are you sure you want to do this?'};
                choice = questdlg(dlgtxt,'Save Initialization','Yes','No','No');
                if strcmp(choice,'Yes') == 1
                    return
                else % give one more chance
                    [fileName,pathName] = uiputfile([objInit.initPath,'newInitialization.mat'],'Save initialization settings.');
                end
            end
            if isempty(fileName) % user chose not to create a file
                return
            end
            default.devs = objInit.devs;
            default.os = objInit.os;
            default.fs = objInit.fs;
            default.fs_now = objInit.fs_now;
            default.delay_now = objInit.delay_now;
            default.card2volts_now = objInit.card2volts_now;
            default.HostAPI = objInit.HostAPI;
            default.hostAPI = objInit.hostAPI;
            default.hostAPIindices = objInit.hostAPIindices;
            default.hostAPI_now = objInit.hostAPI_now;
            %
            default.deviceIndx_in = objInit.deviceIndx_in;
            default.indx_in_now = objInit.indx_in_now;
            default.deviceID_in = objInit.deviceID_in;
            default.id_in_now = objInit.id_in_now;
            default.name_in_now = objInit.name_in_now;
            default.chans_in = objInit.chans_in;
            %
            default.deviceIndx_out = objInit.deviceIndx_out;
            default.indx_out_now = objInit.indx_out_now;
            default.deviceID_out = objInit.deviceID_out;
            default.id_out_now = objInit.id_out_now;
            default.name_out_now = objInit.name_out_now;
            default.chans_out = objInit.chans_out;
            %default.label = objInit.label;
            %default.ampGain = objInit.ampGain;
            %default.micSens = objInit.micSens;
            if exist(pathName,'dir') == 0 % if correct folder does not exist
                success = mkdir(pathName); % try to create it
                if success ~= 1
                    errorTxt = {'  Issue: Error creating the directory for saving initialization.'
                         '  Action: No default values saved.'
                         '  Location: initARLas.saveInitialization.'
                        };
                    errorMsgARLas(errorTxt);
                    objInit.obj.buttonManager(21)
                    return
                end
            end
            try
                save([pathName,fileName],'default');
                alertTxt = {'  Issue: Current initialization saved.'
                     '  Action: None required.'
                     '  Location: initARLas.saveInitialization.'
                    };
                alertMsgARLas(alertTxt);
            catch
            end
            objInit.obj.buttonManager(22)
        catch ME
            errorTxt = {'  Issue: Error saving the current initialization.'
                 '  Action: None.'
                 '  Location: initARLas.saveInitialization.'
                };
            errorMsgARLas(errorTxt);
            objInit.obj.buttonManager(21)
            objInit.printError(ME)
        end
    end
    
    function updateGui(varargin) % update the graphical user interface
        objInit = varargin{1};
        try
            if isempty(objInit.H)
                objInit.initGui
            end
            if objInit.usingSavedValues == 1
                ground = [1 1 1]; 
            else
                ground = [1 .9 .9];
            end
            set(objInit.HOST,'String',objInit.hostAPI,'Value',get(objInit.HOST,'Value'));
            [~,indx] = min(abs(objInit.fs_now - objInit.fs));
            set(objInit.SAMPLING,'Value',indx)
            set(objInit.SAMPLING,'String',objInit.fs,'Value',get(objInit.SAMPLING,'Value'));
            set(objInit.DELAY,'String',num2str(objInit.delay_now),'BackgroundColor',ground)
            set(objInit.CARD2V,'String',num2str(objInit.card2volts_now),'BackgroundColor',ground)
            % OUTPUT -----
            set(objInit.ID_out,'String',num2cell(objInit.deviceID_out),'Value',get(objInit.ID_out,'Value'))
            set(objInit.NAME_out,'String',objInit.name_out_now)
            set(objInit.CHANS_out,'String',num2cell(objInit.chans_out),'BackgroundColor',[1 1 1])
            % INPUT -----
            set(objInit.ID_in,'String',num2cell(objInit.deviceID_in),'Value',get(objInit.ID_in,'Value'))
            set(objInit.NAME_in,'String',objInit.name_in_now)
           set(objInit.CHANS_in,'String',num2cell(objInit.chans_in),'BackgroundColor',[1 1 1]) 
        catch ME
            errorTxt = {'  Issue: Unexpected error updating the initialization GUI.'
                 '  Action: None.'
                 '  Location: initARLas.updateGui.'
                };
            errorMsgARLas(errorTxt);
            objInit.obj.buttonManager(21)
            objInit.printError(ME)
        end
    end
    function printError(varargin)
        objInit = varargin{1};
        ME = varargin{2};
        disp(ME.identifier)
        disp(ME.message)
        disp(ME.stack(1).line)            
    end
end
end