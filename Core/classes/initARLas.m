classdef initARLas < handle
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Class definition: initARLas
% For use with ARLas (Auditory Research Laboratory auditory software)
%
% Auditory Research Lab, The University of Iowa
% Deptartment of Communication Sciences & Disorders
% The University of Iowa
% Author: Shawn S. Goodman, PhD
% Date: August 26, 2016
% Last Updated: November 17, 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
properties (SetAccess = private)
    obj      % passed from arlas
    initPath % passed from arlas
    initFile % passed from arlas
    stimTrain % created internally

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
     OUTPUT % OUTPUT Panel -----
        IDtxt_out % output device idenifying numbers, text label
        ID_out % identifying numbers for possible output devices (numeric), popup menu
        NAMEtxt_out % output device name, text label
        NAME_out % name of currently selected output device, text label
        CHANStxt_out % number of possible output channels, text label
        CHANS_out % number of channels availabel for currently selected device, text label
     INPUT % INPUT panel -----
        IDtxt_in % input device idenifying numbers, text label
        ID_in % identifying numbers for possible input devices (numeric), popup menu
        NAMEtxt_in
        NAME_in
        CHANStxt_in
        CHANS_in
        
        SPECIFYtxt % specify channel input
        SPECIFY
        LABELtxt
        LABEL
        AMPGAINtxt
        AMPGAIN
        MICSENStxt
        MICSENS
        
     % User Buttons -----
     GETDELAY
     GETC2V
     SAVECONFIG % save current configuration

     callingCard % used to figure out which button was the last one called before an error
end
properties (SetAccess = public)
    usingSavedValues % whether or not previously-saved values are being used

    % SYSTEM Data Structures
    devs       % list of devices and properties provided by initialization call to playrec
    os         % the operating system being used
    fs     % list of possible sampling rates for current selections
    fs_now % currently chosen sample rate
    delay_now  % system delay_now (playback and record) for current setup
    card2volts_now % multiplier for conversion from sound card to actual volts
    
    HostAPI    % list of all output HOST APIs (as cell array)
    hostAPI   % short list of unique output HOST APIs
    hostAPIindices  % indices of chosen hostAPI_now in the large list HostAPI)
    hostAPI_now    % currently chosen HOST API
    
    % OUTPUT Data Structures
    deviceIndx_out % list of indices of output devices, given host API
    indx_out_now   % curently chosen output device
    deviceID_out   % list of output device IDs, given host API
    id_out_now     % currently chosen ID
    name_out_now   % currently chosen device name
    chans_out      % possible number of channels to use
    chans_out_now  % number of channels for currently chosen device
    
    % INPUT Data Structures
    deviceIndx_in % full list of indices of output devices
    indx_in_now   % curently chosen value
    deviceID_in   % list of all output device ID options (numeric)
    id_in_now     % currently chosen value
    name_in_now   % currently chosen name
    chans_in      % possible number of channels to use
    chans_in_now  % curently chosen value
    
    specify       % list of input channels to specify label, gain, and sens for
    label         % list of user-defined labels describing the input of each channel
    label_now
    ampGain       % list of amplifier gain setting for each channel
    ampGain_now
    micSens       % list of microphone sensitivity for each channel
    micSens_now
end

methods
    function objInit = initARLas(obj) % initialize object of class initARLas
        objInit.obj = obj; % object of class arlas
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
        objInit.deviceIndx_out = 1;
        objInit.indx_out_now = 1;
        objInit.deviceID_out = initText;
        objInit.id_out_now = initText;
        objInit.name_out_now = initText;
        objInit.chans_out = initText;
        objInit.chans_out_now = initText;
        %
        objInit.deviceIndx_in = 1;
        objInit.indx_in_now = 1;
        objInit.deviceID_in = initText;
        objInit.id_in_now = initText;
        objInit.name_in_now = initText;
        objInit.chans_in = initText;
        objInit.chans_in_now = initText;
        objInit.specify = 0;

        objInit.label.Ch0 = 'Not Assigned'; % default
        objInit.ampGain.Ch0 = 0; % set defaults to 0 dB    
        objInit.micSens.Ch0 = 1; % set defaults to 1 V/Pa
        objInit.label_now = objInit.label.Ch0;
        objInit.ampGain_now = objInit.ampGain.Ch0;
        objInit.micSens_now = objInit.micSens.Ch0;
        %
        objInit.H = obj.H;
        objInit.VIEW = obj.VIEW;
    end  
    function abort(varargin) % instructions for aborting when gui closed
        %objInit = pickObj(varargin);
        %q = findobj('Type','fig'); % find any open figures
        %delete(q); % and delete them
        delete(objInit);
    end
    function initGui(varargin) % initialize gui
        objInit = varargin{1};
        %set(objInit.VIEW,'Title','VIEW: Initialize')
        % create panels within figure
        %[left, bottom,width, height]
        objInit.SYSTEM = uipanel('Parent',objInit.VIEW,'Title','SYSTEM','FontSize',12,...
            'BackgroundColor','white','Units','Normalized',...
            'Position',[.02 .68 .96 .29]);
        objInit.INPUT = uipanel('Parent',objInit.VIEW,'Title','INPUT','FontSize',12,...
            'BackgroundColor','white','Units','Normalized',...
            'Position',[.02 .13 .46 .53]);
        objInit.OUTPUT = uipanel('Parent',objInit.VIEW,'Title','OUTPUT','FontSize',12,...
            'BackgroundColor','white','Units','Normalized',...
            'Position',[.52 .13 .46 .53]); 
        % SYSTEM ----- populate system panel
        spacing = .3;
        top = .87;
        height = 0.125;
        objInit.OStxt = uicontrol('Parent',objInit.SYSTEM,'Style','text','BackgroundColor','white',...
            'String','Computer:         ','HorizontalAlignment','Left','Units','Normalized',...
            'ToolTipString','Type of computer on which MATLAB is executing.',...
            'Position',[.01 top-(.2*spacing) .5 height],'FontSize',10);
        objInit.OS = uicontrol('Parent',objInit.SYSTEM,'Style','text',...
            'String',objInit.os,'FontSize',10,'BackgroundColor','white',...
            'HorizontalAlignment','Left','Units','Normalized',...
            'Position',[.2 top-(.2*spacing) .5 height],'Value',2);
        objInit.HOSTtxt = uicontrol('Parent',objInit.SYSTEM,'Style','text',...
            'BackgroundColor','white','String','Host API:         ','HorizontalAlignment','Left',...
            'ToolTipString','Application Programming Interface.','Units','Normalized',...
            'Position',[.01 top-(1.15*spacing) .5 height],'FontSize',10);
        objInit.HOST = uicontrol('Parent',objInit.SYSTEM,'Style','popup',...
            'String',objInit.hostAPI,'FontSize',10,'BackgroundColor','white',...
            'HorizontalAlignment','Left','Units','Normalized',...
            'Position',[.2 top-(1*spacing) .25 height],'Value',objInit.indx_out_now,...
            'Callback',{@objInit.hostManager,objInit.devs});
        objInit.SAMPLINGtxt = uicontrol('Parent',objInit.SYSTEM,'Style','text','BackgroundColor','white',...
             'String','Sampling Rate:    ','HorizontalAlignment','Left',...
             'ToolTipString','Sampling frequency in Hz.','Units','Normalized',...
             'Position',[.01 top-(2*spacing) .5 height],'FontSize',10);
        objInit.SAMPLING = uicontrol('Parent',objInit.SYSTEM,'Style','popup',...
             'String',objInit.fs,'FontSize',10,'Units','Normalized',...
             'Position',[.2 top-(2*spacing)+.025 .25 height],'Value',1,...
             'Callback',{@objInit.samplingManager,objInit.devs});
         objInit.DELAYtxt = uicontrol('Parent',objInit.SYSTEM,'Style','text','BackgroundColor','white',...
             'String','System Delay:    ','HorizontalAlignment','Left','Units','Normalized',...
             'Position',[.5 top-(1.15*spacing) .5 height],'FontSize',10,'ToolTipString','Sample delay between playback and recording.');
        objInit.DELAY = uicontrol('Parent',objInit.SYSTEM,'Style','edit',...
             'String',objInit.delay_now,'FontSize',10,'BackgroundColor','white','HorizontalAlignment',...
             'Left','Callback',@objInit.delayManager,'Units','Normalized',...
             'Position',[.7 top-(1.25*spacing) .2 .18]);
         objInit.CARD2Vtxt = uicontrol('Parent',objInit.SYSTEM,'Style','text','BackgroundColor','white',...
             'String','Card to Volts:    ','HorizontalAlignment','Left',...
             'ToolTipString','Multiplier to convert sound card units to voltage.','Units','Normalized',....
             'Position',[.5 top-(2*spacing) .5 height],'FontSize',10);
        objInit.CARD2V = uicontrol('Parent',objInit.SYSTEM,'Style','edit',...
            'String',objInit.card2volts_now,'FontSize',10,'BackgroundColor','white',...
            'HorizontalAlignment','Left','Callback',@objInit.c2vManager,'Units','Normalized',...
            'Position',[.7 top-(2.125*spacing) .2 .18]);
        % OUTPUT ----- populate output panel
        spacing = .1;
        top = .9;
        height = 0.125;
        objInit.IDtxt_out = uicontrol('Parent',objInit.OUTPUT,'Style','text',...
            'BackgroundColor','white','String','Device ID:        ',...
            'ToolTipString','ID number used to refer to outupt devices.',...
            'HorizontalAlignment','Left','Units','Normalized',...
            'Position',[.01 top-(.3*spacing) .5 .1],'FontSize',10); % [.01 .8 .5 .1]
        objInit.ID_out = uicontrol('Parent',objInit.OUTPUT,'Style','popup',...
            'String',objInit.id_out_now,'FontSize',10,'BackgroundColor',...
            'white','HorizontalAlignment','Left','Units','Normalized',...
            'Position',[.5 top-(.10*spacing) .2 .1],'Value',objInit.indx_out_now,...
            'Callback',{@objInit.outputIdManager,objInit.devs});
        objInit.NAMEtxt_out = uicontrol('Parent',objInit.OUTPUT,'Style',...
            'text','BackgroundColor','white','String','Device Name:      ',...
            'ToolTipString','Name of currently selected output device.',...
            'HorizontalAlignment','Left','Units','Normalized',...
            'Position',[.01 top-(1.5*spacing) .5 .1],'FontSize',10);       
        objInit.NAME_out = uicontrol('Parent',objInit.OUTPUT,'Style','text',...
            'String',objInit.name_out_now,'HorizontalAlignment','Left',...
            'FontSize',10,'BackgroundColor','white','Units','Normalized',...
            'Position',[.5 top-(1.5*spacing) .5 .1]);
        objInit.CHANStxt_out = uicontrol('Parent',objInit.OUTPUT,'Style','text',...
            'BackgroundColor','white','String','Num. of Channels:       ',...
            'ToolTipString','Number of desired output channels (up to maximum supported by the device).',... 
            'HorizontalAlignment','Left','Units','Normalized',...
            'Position',[.01 top-(2.7*spacing) .5 .1],'FontSize',10);
        objInit.CHANS_out = uicontrol('Parent',objInit.OUTPUT,'Style','popup',...
            'String',objInit.chans_out,...
            'FontSize',10,'BackgroundColor','white','HorizontalAlignment','Left',...
            'Units','Normalized','Callback',@objInit.chansOutManager,...
            'Position',[.5 top-(2.6*spacing) .2 .1]);
        % INPUT ----- populate input panel
        objInit.IDtxt_in = uicontrol('Parent',objInit.INPUT,'Style','text','BackgroundColor','white',...
            'ToolTipString','ID number used to refer to input devices.',...
            'String','Device ID:        ','HorizontalAlignment','Left','Units','Normalized',...
            'Position',[.01 top-(.3*spacing) .5 .1],'FontSize',10);
        objInit.ID_in = uicontrol('Parent',objInit.INPUT,'Style','popup',...
             'String',objInit.id_in_now,'FontSize',10,'BackgroundColor','white',...
             'HorizontalAlignment','Left','Units','Normalized',...
             'Position',[.5 top-(.10*spacing) .2 .1],'Value',objInit.indx_in_now,...
             'Callback',{@objInit.inputIdManager,objInit.devs});
        objInit.NAMEtxt_in = uicontrol('Parent',objInit.INPUT,'Style','text','BackgroundColor','white',...
             'ToolTipString','Name of currently selected input device.',...
             'String','Device Name:      ','HorizontalAlignment','Left',...
             'Units','Normalized','Position',[.01 top-(1.5*spacing) .5 .1],'FontSize',10);       
        objInit.NAME_in = uicontrol('Parent',objInit.INPUT,'Style','text',...
             'String',objInit.name_in_now,'FontSize',10,'BackgroundColor','white',...
             'HorizontalAlignment','Left',...
             'Units','Normalized','Position',[.5 top-(1.5*spacing) .5 .1],'Value',2);
        objInit.CHANStxt_in = uicontrol('Parent',objInit.INPUT,'Style','text','BackgroundColor','white',...
            'ToolTipString','Number of desired input channels (up to maximum supported by the device).',...  
            'String','Num. of Channels:       ','HorizontalAlignment','Left',...
             'Units','Normalized','Position',[.01 top-(2.7*spacing) .5 .1],'FontSize',10);
        objInit.CHANS_in = uicontrol('Parent',objInit.INPUT,'Style','popup',...
             'String',objInit.chans_in,'FontSize',10,'BackgroundColor','white',...
             'HorizontalAlignment','Left','Units','Normalized','Value',1,...
             'Position',[.5 top-(2.6*spacing) .2 .1],'Callback',@objInit.chansInManager);
        %[left,bottom,width,height]
        objInit.SPECIFYtxt = uicontrol('Parent',objInit.INPUT,'Style','text','BackgroundColor','white',...
             'ToolTipString','Use to help assign label, gain, and sensitivity parameters to each channel.',... 
             'String','Specify Input for Channel:       ','HorizontalAlignment','Left',...
             'Units','Normalized','Position',[.01 top-(4.6*spacing) .65 .1],'FontSize',10);
        objInit.SPECIFY = uicontrol('Parent',objInit.INPUT,'Style','popup',...
             'String',objInit.specify,'FontSize',10,'BackgroundColor','white',...
             'HorizontalAlignment','Left','Units','Normalized','Value',1,...
             'Position',[.7 top-(4.5*spacing) .25 .1],'Callback',@objInit.inputManager);
         objInit.LABELtxt = uicontrol('Parent',objInit.INPUT,'Style','text','BackgroundColor','white',...
             'ToolTipString','User-defined alpha-numeric name. This is for ease of use and is not required.',... 
             'String','Channel Label:       ','HorizontalAlignment','Left',...
             'Units','Normalized','Position',[.1 top-(6.0*spacing) .6 .1],'FontSize',10);
        objInit.LABEL = uicontrol('Parent',objInit.INPUT,'Style','edit',...
             'String',objInit.label_now,'FontSize',10,'BackgroundColor','white',...
             'HorizontalAlignment','Left','Units','Normalized','Value',1,...
             'Position',[.7 top-(5.9*spacing) .25 .1],'Callback',@objInit.labelManager);         
         objInit.AMPGAINtxt = uicontrol('Parent',objInit.INPUT,'Style','text','BackgroundColor','white',...
             'ToolTipString','Total amplifier gain prior to ADC on the sound card.',... 
             'String','Amplifier Gain (dB):       ','HorizontalAlignment','Left',...
             'Units','Normalized','Position',[.1 top-(7.4*spacing) .6 .1],'FontSize',10);
        objInit.AMPGAIN = uicontrol('Parent',objInit.INPUT,'Style','edit',...
             'String',objInit.ampGain_now,'FontSize',10,'BackgroundColor','white',...
             'HorizontalAlignment','Left','Units','Normalized','Value',1,...
             'Position',[.7 top-(7.3*spacing) .25 .1],'Callback',@objInit.ampGainManager);         
       objInit.MICSENStxt = uicontrol('Parent',objInit.INPUT,'Style','text','BackgroundColor','white',...
             'ToolTipString','Microphone sensitivity. If = 1, will assume no mic is being used.',... 
             'String','Mic Sensitivity (V/Pa):       ','HorizontalAlignment','Left',...
             'Units','Normalized','Position',[.1 top-(8.7*spacing) .6 .1],'FontSize',10);
        objInit.MICSENS = uicontrol('Parent',objInit.INPUT,'Style','edit',...
             'String',objInit.micSens_now,'FontSize',10,'BackgroundColor','white',...
             'HorizontalAlignment','Left','Units','Normalized','Value',1,...
             'Position',[.7 top-(8.6*spacing) .25 .1],'Callback',@objInit.micSensManager);         
        % create buttons at bottom of figure
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
        objInit = varargin{1};
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
            objInit.chans_out = (0:1:objInit.devs(objInit.indx_out_now).outputChans)';
            set(objInit.CHANS_out,'String',objInit.chans_out)
            set(objInit.CHANS_out,'Value',value)
            objInit.chans_out_now = objInit.CHANS_out.String(value);
            objInit.fs_now = objInit.devs(objInit.indx_out_now).defaultSampleRate;
            [~,value] = min(abs(objInit.fs_now-objInit.fs));
            set(objInit.SAMPLING,'Value',value)
            objInit.id_in_now = objInit.devs(objInit.indx_in_now).deviceID;
            objInit.name_in_now = objInit.devs(objInit.indx_in_now).name;
            objInit.chans_in = (0:1:objInit.devs(objInit.indx_in_now).inputChans)';
            set(objInit.CHANS_in,'String',objInit.chans_in)
            set(objInit.CHANS_in,'Value',value)
            objInit.chans_in_now = objInit.CHANS_in.String(value);
            objInit.callingCard = 'hostManager';
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
                objInit.callingCard = 'samplingManager';
            catch ME
                if isempty(objInit.callingCard)|| strcmp(objInit.callingCard,'samplingManager')
                    errorTxt = {'  Issue: selected sampling rate is not compatible with the other selections.'
                         '  Suggested Action: Unknown.'
                         '  Location: initARLas.samplingManager.'
                        };
                    warnMsgARLas(errorTxt);
                    set(objInit.SAMPLING,'BackgroundColor',[1 0 0])
                elseif strcmp(objInit.callingCard,'inputIDManager') || strcmp(objInit.callingCard,'outputIDManager')
                    errorTxt = {'  Issue: The selected device IDs are not compatible with each other.'
                         '  Suggested Action: Change one of the device IDs so that they match.'
                         '  Location: initARLas.samplingManager.'
                        };
                    warnMsgARLas(errorTxt);
                    set(objInit.ID_out,'BackgroundColor',[1 0 0])
                    set(objInit.ID_in,'BackgroundColor',[1 0 0])
                else 
                    errorTxt = {'  Issue: Unexpected error selecting Sampling Rate.'
                         '  Action: None.'
                         '  Location: initARLas.samplingManager.'
                        };
                    errorMsgARLas(errorTxt);
                    objInit.obj.buttonManager(21)
                    objInit.callingCard
                end
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
            objInit.chans_in = (0:1:objInit.devs(objInit.indx_in_now).inputChans)';
            set(objInit.CHANS_in,'String',objInit.chans_in)
            set(objInit.CHANS_in,'Value',1)
            objInit.CHANS_in.BackgroundColor = [1 0 0];
            objInit.chans_in_now = objInit.CHANS_in.String(1);
            objInit.callingCard = 'inputIDManager';
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
    function chansInManager(varargin) % manage setting number of input channels
        try
            objInit = varargin{1};
            objInit.chans_in_now = str2double(objInit.CHANS_in.String(objInit.CHANS_in.Value));
            if objInit.chans_in_now < 1 
                objInit.CHANS_in.BackgroundColor = [1 .9 .9];
            else
                objInit.CHANS_in.BackgroundColor = [1 1 1]; 
            end
            % update channel specification options ----------------------------
            names = fieldnames(objInit.micSens);
            sizeDiff = size(names,1) - (objInit.chans_in_now + 1);
            if sizeDiff < 0 % build new fields
                for ii=size(names,1):objInit.chans_in_now
                    objInit.label.(matlab.lang.makeValidName(['Ch',num2str(ii)]))= 'Not Assigned';    
                    objInit.ampGain.(matlab.lang.makeValidName(['Ch',num2str(ii)]))= 0; % set defaults to 0 dB    
                    objInit.micSens.(matlab.lang.makeValidName(['Ch',num2str(ii)]))= 1; % set defaults to 1 V/Pa
                end
            elseif sizeDiff > 0 % delete old fields
                for ii=size(names,1):-1:objInit.chans_in_now+2
                    objInit.label = rmfield(objInit.label,names{ii});
                    objInit.ampGain = rmfield(objInit.ampGain,names{ii});
                    objInit.micSens = rmfield(objInit.micSens,names{ii});
                end
            end
            objInit.specify = (0:1:objInit.chans_in_now)';
            set(objInit.SPECIFY,'String',objInit.specify)
            set(objInit.SPECIFY,'Value',1)
        catch ME
            errorTxt = {'  Issue: Unexpected error selecting Number of Input Channels.'
                 '  Action: None.'
                 '  Location: initARLas.chansInManager.'
                };
            errorMsgARLas(errorTxt);
            objInit.obj.buttonManager(21)
            objInit.printError(ME)
        end
    end
    function inputManager(varargin) % specify input label, gain and sens for channels
        objInit = varargin{1};
        try
            indx = objInit.SPECIFY.Value;
            objInit.label_now = objInit.label.(matlab.lang.makeValidName(['Ch',num2str(indx-1)]));    
            objInit.ampGain_now = objInit.ampGain.(matlab.lang.makeValidName(['Ch',num2str(indx-1)]));
            objInit.micSens_now = objInit.micSens.(matlab.lang.makeValidName(['Ch',num2str(indx-1)]));
            objInit.LABEL.String = objInit.label_now;
            objInit.AMPGAIN.String = objInit.ampGain_now;
            objInit.MICSENS.String = objInit.micSens_now;
        catch ME
            errorTxt = {'  Issue: Unexpected error Specifying Input for Channel.'
                 '  Action: None.'
                 '  Location: initARLas.inputManager.'
                };
            errorMsgARLas(errorTxt);
            objInit.obj.buttonManager(21)
            objInit.printError(ME)
        end
    end
    function labelManager(varargin) % specify a label for each channel
        objInit = varargin{1};
        value = objInit.LABEL.String; % get the current value
        ok = 1; % assume the input is okay, unless an error is found
        try % check the current input value
            indx = objInit.SPECIFY.Value;
            try 
                value = num2str(value); 
            catch
            end
            if indx == 1
                ok = 0;
                errorTxt = {'  Issue: Channel Label cannot be specified for Ch 0 (no input).'
                     '  Action: Defaulting to a value of Not Assigned.'
                     '  Location: initARLas.labelManager.'
                    };
                errorMsgARLas(errorTxt);
                objInit.obj.buttonManager(21)
                value = 'Not Assigned';
            elseif isempty(value) % do not allow empty values
                ok = 0;
                errorTxt = {'  Issue: Channel Label cannot be empty.'
                     '  Action: Defaulting to a value of Not Assigned.'
                     '  Location: initARLas.labelManager.'
                    };
                errorMsgARLas(errorTxt);
                objInit.obj.buttonManager(21)
                value = 'Not Assigned';
            elseif isa(value,'numeric') == 1 % do not allow numeric values
                ok = 0;
                errorTxt = {'  Issue: Channel Label must be a string.'
                     '  Action: Defaulting to a value of Not Assigned.'
                     '  Location: initARLas.labelManager.'
                    };
                errorMsgARLas(errorTxt);
                objInit.obj.buttonManager(21)
                value = 'Not Assigned';
            end
            objInit.label_now = value;
            indx = objInit.SPECIFY.Value;
            objInit.label.(matlab.lang.makeValidName(['Ch',num2str(indx-1)])) = objInit.label_now;
            objInit.LABEL.String = objInit.label_now;
            if ok == 0
                if indx > 1
                    objInit.LABEL.BackgroundColor = [1 0 0]; % make background red
                else % this isn't needing to be changed, so don't alert user to change it
                    objInit.LABEL.BackgroundColor = [1 1 1]; % make background white
                end
            else
                objInit.LABEL.BackgroundColor = [1 1 1]; % make background white
            end
        catch ME
            errorTxt = {'  Issue: Unexpected error setting Channel Label value.'
                 '  Action: None.'
                 '  Location: initARLas.labelManager.'
                };
            errorMsgARLas(errorTxt);
            objInit.obj.buttonManager(21)
            objInit.printError(ME)
        end        
    end
    function ampGainManager(varargin)  % specify gain for each channel
        objInit = varargin{1};
        value = objInit.AMPGAIN.String; % get the current value
        ok = 1; % assume the input is okay, unless an error is found
        try % check the current input value
            indx = objInit.SPECIFY.Value;
            try 
                value = str2num(value); 
            catch
            end
            if indx == 1
                ok = 0;
                errorTxt = {'  Issue: Amplifier Gain cannot be specified for Ch 0 (no input).'
                     '  Action: Defaulting to a value of 0 dB.'
                     '  Location: initARLas.ampGainManager.'
                    };
                errorMsgARLas(errorTxt);
                objInit.obj.buttonManager(21)
                value = 0;
            elseif isempty(value) % do not allow empty values
                ok = 0;
                errorTxt = {'  Issue: Amplifier Gain cannot be empty.'
                     '  Action: Defaulting to a value of 0 dB.'
                     '  Location: initARLas.ampGainManager.'
                    };
                errorMsgARLas(errorTxt);
                objInit.obj.buttonManager(21)
                value = 0;
            elseif isa(value,'numeric') == 0 % do not allow non-numeric values
                ok = 0;
                errorTxt = {'  Issue: Amplifier Gain must be a numeric value.'
                     '  Action: Defaulting to a value of 0 dB.'
                     '  Location: initARLas.ampGainManager.'
                    };
                errorMsgARLas(errorTxt);
                objInit.obj.buttonManager(21)
                value = 0;
            end
            objInit.ampGain_now = value;
            indx = objInit.SPECIFY.Value;
            objInit.ampGain.(matlab.lang.makeValidName(['Ch',num2str(indx-1)])) = objInit.ampGain_now;
            objInit.AMPGAIN.String = objInit.ampGain_now;
            if ok == 0
                if indx > 1
                    objInit.AMPGAIN.BackgroundColor = [1 0 0]; % make background red
                else
                    objInit.AMPGAIN.BackgroundColor = [1 1 1]; % make background white
                end
            else
                objInit.AMPGAIN.BackgroundColor = [1 1 1]; % make background white
            end
        catch ME
            errorTxt = {'  Issue: Unexpected error setting Amplifier Gain value.'
                 '  Action: None.'
                 '  Location: initARLas.ampGainManager.'
                };
            errorMsgARLas(errorTxt);
            objInit.obj.buttonManager(21)
            objInit.printError(ME)
        end        
    end
    function micSensManager(varargin) % specify microphone sensitivity for each channel
        objInit = varargin{1};
        value = objInit.MICSENS.String; % get the current value
        ok = 1; % assume the input is okay, unless an error is found
        try % check the current input value
            indx = objInit.SPECIFY.Value;
            try 
                value = str2num(value); 
            catch
            end
            if indx == 1
                ok = 0;
                errorTxt = {'  Issue: Mic Sensitivity cannot be specified for Ch 0 (no input).'
                     '  Action: Defaulting to a value of 1 (unity).'
                     '  Location: initARLas.micSensManager.'
                    };
                errorMsgARLas(errorTxt);
                objInit.obj.buttonManager(21)
                value = 1;
            elseif isempty(value) % do not allow empty values
                ok = 0;
                errorTxt = {'  Issue: Mic Sensitivity cannot be empty.'
                     '  Action: Defaulting to a value of 1 (unity).'
                     '  Location: initARLas.micSensManager.'
                    };
                errorMsgARLas(errorTxt);
                objInit.obj.buttonManager(21)
                value = 1;
            elseif isa(value,'numeric') == 0 % do not allow non-numeric values
                ok = 0;
                errorTxt = {'  Issue: Mic Sensitivity must be a numeric value.'
                     '  Action: Defaulting to a value of 1 (unity).'
                     '  Location: initARLas.micSensManager.'
                    };
                errorMsgARLas(errorTxt);
                objInit.obj.buttonManager(21)
                value = 1;
            end
            objInit.micSens_now = value;
            indx = objInit.SPECIFY.Value;
            objInit.micSens.(matlab.lang.makeValidName(['Ch',num2str(indx-1)])) = objInit.micSens_now;
            objInit.MICSENS.String = objInit.micSens_now;
            if ok == 0
                if indx > 1
                    objInit.MICSENS.BackgroundColor = [1 0 0]; % make background red
                else
                    objInit.MICSENS.BackgroundColor = [1 1 1]; % make background white
                end
            else
                objInit.MICSENS.BackgroundColor = [1 1 1]; % make background white
            end
        catch ME
            errorTxt = {'  Issue: Unexpected error setting Mic Sensitivity value.'
                 '  Action: None.'
                 '  Location: initARLas.micSensManager.'
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
        
        objInit.chans_out = (0:1:objInit.devs(objInit.indx_out_now).outputChans)';
        set(objInit.CHANS_out,'String',objInit.chans_out)
        set(objInit.CHANS_out,'Value',1)
        objInit.chans_out_now = objInit.CHANS_out.String(1);
        objInit.CHANS_out.BackgroundColor = [1 0 0];
        
        objInit.fs_now = objInit.devs(objInit.indx_out_now).defaultSampleRate;
        [~,value] = min(abs(objInit.fs_now-objInit.fs));
        
        set(objInit.SAMPLING,'Value',value)
        objInit.callingCard = 'outputIDManager';
        objInit.samplingManager
    end
    function chansOutManager(varargin) % manage setting number of output channels
        objInit = varargin{1};
        objInit.chans_out_now = str2double(objInit.CHANS_out.String(objInit.CHANS_out.Value));
        if objInit.chans_out_now < 1 
            objInit.CHANS_out.BackgroundColor = [1 .9 .9]; % make this one pink because not truly and error, just a warning
        else
            objInit.CHANS_out.BackgroundColor = [1 1 1]; 
        end
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
            objInit.callingCard = 'discover';
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
    function getCard2V(varargin) % calculate multiplier to convert card units to volts
        objInit = varargin{1};
        try
            txt = ({'This routine will calculate a multiplier to convert card units to voltage.';'';...
                    'A 100 Hz pure tone will be played through all the output channels designated by "Num. of Channels".';'';...
                    'The tone will be played for 10 seconds.';'';...
                    'Use a multimeter to measure the open circuit AC voltage on one of the output channels.';'';...
                    'Write this value down. You will enter it in the next step.'});
            choice = questdlg(txt,'Get Card to Volts','Continue','Cancel','Continue');
            if strcmp(choice,'Continue')
                % do nothing; continue
            elseif strcmp(choice,'Cancel')
                return
            else % if user shut down box using the x
                return
            end
            oldC2V = objInit.card2volts_now; % hold aside the current value
            objInit.card2volts_now = 1; % set to 1 in order to measure
            set(objInit.CARD2V,'String',num2str(objInit.card2volts_now))
            objPlayrec = playrecARLas(objInit);
            
            % make a sinusoidal stimulus
            len = 0.5; % stimulus length
            nSamples = round(len * objInit.fs_now);
            if mod(nSamples,2) ~= 0
                nSamples = nSamples + 1;
            end
            t = (0:1:nSamples-1)'/objInit.fs_now; % time in seconds
            f = 100; % frequency in Hz
            a = .5; % amplitude (full out is 1)
            stimulus = a * cos(2*pi*f*t);
            Stimulus = repmat(stimulus,1,objInit.chans_out_now);
            objPlayrec.stimTrain = Stimulus; % load the stimulus
            objPlayrec.nReps = 10;
            objInit.makeInvisible
            %objInit.killPlots
            objInit.obj.buttonManager(20)
            try
                objPlayrec.run
            catch ME
                errorTxt = {'  Issue: Unexpected error occurred while calculating the card to volts conversion.'
                     '  Action: None.'
                     '  Location: initARLas.getCard2V.'
                    };
                errorMsgARLas(errorTxt);
                objInit.obj.buttonManager(21)
                objInit.printError(ME)
                return
            end
            objInit.obj.buttonManager(22)
            %objPlayrec.killPlots
            objPlayrec.makeInvisible
            objInit.makeVisible

            prompt = {'Enter the measured voltage (in Volts): '};
            title = 'Open Circuit Voltage'; 
            answer = inputdlg(prompt,title);
            if isempty(answer) % user hit cancel
                return
            elseif strcmp(answer{1},'') % if user left field blank
                return
            else % user put something in the field
                openCircuitVoltage = abs(str2num(answer{1})); % voltage rms must be positive 
            end
            if isempty(openCircuitVoltage) % but that something was not numeric
                errorTxt = {'  Issue: Illegal input for voltage; possibly non-numeric value.'
                     '  Action: getCard2V terminated prematurely.'
                     '  Location: initARLas.getCard2V.'
                    };
                errorMsgARLas(errorTxt);
                objInit.obj.buttonManager(21)
                return
            end
            txt = ({'Reconnect the output to the input channel designated by "Specify Input for Channel"..';'';...
                    'A 10-second, 100 Hz pure tone will again be played.';''});
            choice = questdlg(txt,'Get Card to Volts','Continue','Cancel','Continue');
            if strcmp(choice,'Continue')
                % do nothing; continue
            elseif strcmp(choice,'Cancel')
                return
            else % if user shut down box using the x
                return
            end
            objInit.makeInvisible
            %objInit.killPlots
            set(objInit.VIEW,'Title','VIEW: Initialize')
            objInit.obj.buttonManager(20)
            try
                objPlayrec.run
            catch ME
                errorTxt = {'  Issue: Unexpected error occurred while calculating the card to volts conversion.'
                     '  Action: None.'
                     '  Location: initARLas.getCard2V.'
                    };
                errorMsgARLas(errorTxt);
                objInit.obj.buttonManager(21)
                objInit.printError(ME)
                return
            end
            objInit.obj.buttonManager(22)
            objPlayrec.makeInvisible
            %objPlayrec.killPlots
            objInit.makeVisible
            if objPlayrec.killRun == 0
                Data = objPlayrec.Data;
                indx = objInit.specify(objInit.SPECIFY.Value);
                if indx == 0 % can't calculate delay from no channel
                    indx = 1;
                end
                data = Data(:,indx);
                data = reshape(data,nSamples,objPlayrec.nReps);
                m = (mean(data,2));
                % note that multimeter voltage is assumed to be RMS, not peak!
                closedCircuitVoltage = sqrt(mean(m.^2)); % therefore, this is the correct comparison
                % take into account current input settings
                closedCircuitVoltage = closedCircuitVoltage / objPlayrec.card2volts * objPlayrec.micSens(indx) * objPlayrec.ampGain(indx);
                cardMultiplier = openCircuitVoltage / closedCircuitVoltage; 
                objInit.card2volts_now = cardMultiplier; 
                set(objInit.CARD2V,'String',num2str(objInit.card2volts_now))
                if oldC2V ~= objInit.card2volts_now
                    alertTxt = {'  Issue: New card2volts value does not match the old value.'
                         ['  Action: Replaced old delay value (',num2str(oldC2V),') with new one (',num2str(objInit.card2volts_now),').']
                         '  Recommended Action: Save new configuration using the "Save Configuration" button.'
                         '  Location: initARLas.getCard2V.'
                        };
                    alertMsgARLas(alertTxt);
                    objInit.obj.buttonManager(21)
                end
            else
                warnTxt = {'  Issue: killRun called during objPlayrec.run.'
                     '  Action: Card2V process terminated early.'
                     '  Location: initARLas.getCard2V.'
                    };
                warnMsgARLas(warnTxt);
                objInit.obj.buttonManager(21)
            end
        objInit.obj.buttonManager(22)
        catch ME
            errorTxt = {'  Issue: Unexpected error calculating the card to volts conversion.'
                 '  Action: None.'
                 '  Location: initARLas.getCard2V.'
                };
            errorMsgARLas(errorTxt);
            objInit.obj.buttonManager(21)
            objInit.printError(ME)
        end
    end
    function getDelay(varargin) % try running the currently-chosen values to see if they work
        objInit = varargin{1};
        if (objInit.chans_in_now < 1) || (objInit.chans_out_now < 1)
            errorTxt = {'  Issue: Cannot calculate delay when Num of Channels is zero.'
                 '  Action: Aborting getCard2V.'
                 '  Location: initARLas.getCard2V.'
                };
            errorMsgARLas(errorTxt);
            objInit.obj.buttonManager(21)
            return
        end
        try
            txt = ({'This routine will calculate your system delay.';'';...
                    'Clicks will be played through all the output channels designated by "Num. of Channels".';'';...
                    'Input will be measured through the channel designated by "Specify Input for Channel".';'';...
                    'Direct elecrical connection between output and input is recommended.'});
            choice = questdlg(txt,'Get System Delay','Continue','Cancel','Continue');
            if strcmp(choice,'Continue')
                % do nothing; continue
            elseif strcmp(choice,'Cancel')
                return
            else % if user shut down box using the x
                return
            end
            objInit.makeInvisible
            
            oldDelay = objInit.delay_now; % hold aside the current value
            objInit.delay_now = 0; % set to zero in order to measure
            set(objInit.DELAY,'String',num2str(objInit.delay_now))
            objPlayrec = playrecARLas(objInit);
            
            len = 0.5; % stimulus length
            nSamples = round(len * objInit.fs_now);
            objInit.stimTrain = zeros(nSamples,objInit.chans_out_now);
            objInit.stimTrain(1,:) = .01;
            objPlayrec.stimTrain = objInit.stimTrain;
            objPlayrec.nReps = 10;
            objPlayrec.systemDelay = objInit.delay_now;
            
            objInit.makeInvisible
            set(objInit.VIEW,'Title','VIEW: Initialize')
            objInit.obj.buttonManager(20)
            try
                objPlayrec.run
            catch ME
                errorTxt = {'  Issue: Unexpected error occurred while calculating the system delay.'
                     '  Action: None.'
                     '  Location: initARLas.getDelay.'
                    };
                errorMsgARLas(errorTxt);
                objInit.obj.buttonManager(21)
                objInit.printError(ME)
                return
            end
            objInit.obj.buttonManager(22)
            objPlayrec.makeInvisible
            objInit.makeVisible
            if objPlayrec.killRun == 0
                Data = objPlayrec.Data;
                indx = objInit.specify(objInit.SPECIFY.Value);
                if indx == 0 % can't calculate delay from no channel
                    indx = 1;
                end
                data = Data(:,indx);
                data = reshape(data,nSamples,objPlayrec.nReps);
                m = (mean(data,2));
                [~,maxyIndx] = max(abs(m));
                objInit.delay_now = maxyIndx-1; 
                set(objInit.DELAY,'String',num2str(objInit.delay_now))
                if oldDelay ~= objInit.delay_now
                    alertTxt = {'  Issue: New system delay value does not match the old value.'
                         ['  Action: Replaced old delay value (',num2str(oldDelay),') with new one (',num2str(objInit.delay_now),').']
                         '  Recommended Action: Save new configuration using the "Save Configuration" button.'
                         '  Location: initARLas.getDelay.'
                        };
                    alertMsgARLas(alertTxt);
                    objInit.obj.buttonManager(21)
                end
            else
                warnTxt = {'  Issue: killRun called during objPlayrec.run.'
                     '  Action: getDelay process terminated early.'
                     '  Location: initARLas.getDelay.'
                    };
                warnMsgARLas(warnTxt);
                objInit.obj.buttonManager(21)
            end
        catch ME
            errorTxt = {'  Issue: Unexpected error occurred while calculating the system delay.'
                 '  Action: None.'
                 '  Location: initARLas.getDelay.'
                };
            errorMsgARLas(errorTxt);
            objInit.obj.buttonManager(21)
            objInit.printError(ME)
        end
    end
    function saveConfig(varargin)
        objInit = varargin{1};
        objInit.callingCard = 'saveInitialization';
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
            default.deviceIndx_out = objInit.deviceIndx_out;
            default.indx_out_now = objInit.indx_out_now;
            default.deviceID_out = objInit.deviceID_out;
            default.id_out_now = objInit.id_out_now;
            default.name_out_now = objInit.name_out_now;
            default.chans_out = objInit.chans_out;
            default.chans_out_now = objInit.chans_out_now;
            default.deviceIndx_in = objInit.deviceIndx_in;
            default.indx_in_now = objInit.indx_in_now;
            default.deviceID_in = objInit.deviceID_in;
            default.id_in_now = objInit.id_in_now;
            default.name_in_now = objInit.name_in_now;
            default.chans_in = objInit.chans_in;
            default.chans_in_now = objInit.chans_in_now;
            default.specify = objInit.specify;
            default.label = objInit.label;
            default.ampGain = objInit.ampGain;
            default.micSens = objInit.micSens;
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
                alertTxt = {'  Issue: Current initializaiton saved.'
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
            if isa(objInit.chans_out_now,'numeric') == 0
                try
                objInit.chans_out_now = str2num(objInit.chans_out_now);
                catch
                    objInit.chans_out_now = 0;
                end
            end
            [dummy,indx] = min(abs(objInit.chans_out_now - objInit.chans_out));
            set(objInit.CHANS_out,'String',num2cell(objInit.chans_out),'Value',indx,'BackgroundColor',ground)
            % INPUT -----
            set(objInit.ID_in,'String',num2cell(objInit.deviceID_in),'Value',get(objInit.ID_in,'Value'))
            set(objInit.NAME_in,'String',objInit.name_in_now)
            if isa(objInit.chans_in_now,'numeric') == 0
                try
                objInit.chans_in_now = str2num(objInit.chans_in_now);
                catch
                    objInit.chans_in_now = 0;
                end
            end
            [dummy,indx] = min(abs(objInit.chans_in_now - objInit.chans_in));
            set(objInit.CHANS_in,'String',num2cell(objInit.chans_in),'Value',indx,'BackgroundColor',ground)
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