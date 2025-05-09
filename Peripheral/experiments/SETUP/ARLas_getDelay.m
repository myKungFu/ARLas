function [] = ARLas_getDelay(varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ARLas_getDelay(varargin)
%
% Measure sample delay between card output and input
%
% NOTES: The original method for doing this was to make a direct electrical
% connection between an input and an output channel. This worked fine until
% Thevenin source-based FPL calibration was introduced. If using FPL
% calibration, the system delay MUST include the entire "system" you are
% using, i.e. computer, sound card, power amp (if using one), and probe-microphone
% system. Therefore, the new and preferred method for getting the delay is
% to present the sound through a probe tip that is held up in a "free
% field" away from any nearby highly reflective objects.
%
% Auditory Research Lab, The University of Iowa
% Deptartment of Communication Sciences & Disorders
% The University of Iowa
% Author: Shawn S. Goodman, PhD
% Date: January 9, 2017
% Last Updated: August 15, 2017
% Last Updated: May 28, 2021 -- ssg -- implemented the new method for
%                                      soundfield delay and added notes.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    obj = varargin{1}; % get the arlas object
    V = '05.28.2021'; % this is the current version number of this program

%------ USER MODIFIABLE PARAMETERS ----------------------------------------
    probeName = 'A'; % 'A' or 'B'
    % Choose the probe microphone and loudspeaker channel to use.
    % It shouldn't matter which channels you are running on; they should
    % all yield the same answer. The important thing is to use a set of
    % configured hardware/software that is the same as what you will be
    % using to make actual measurements.
%--------------------------------------------------------------------------

    [inputs,outputs] = hardwareSetup;
    if strcmp(probeName,'A')
        probeInput = inputs{1}; % 1         % ER10xA microphone
        probeOutput = outputs{1}; % 1      % ER10xA loudspeakers
    elseif strcmp(probeName,'B')
        probeInput = inputs{1};         % ER10xB microphone
        probeOutput = outputs{1};       % ER10xB loudspeakers
    else
        error('Unrecognized value assigned to variable probeName.')
    end


    % 1) CREATE THE STIMULUS
    % --------------------------------------------------
    %                                                                                                                                                                                                                                                 
    fs = obj.fs; % get the system sampling rate
    len = 0.5; % desired stimulus length (s)
        nSamples = round(len * fs); % number of samples in stimulus
        if mod(nSamples,2) ~= 0 % force to be even number
            nSamples = nSamples + 1;
        end
    stimulus = zeros(nSamples,1); % create a vector of zeros
    stimulus(1,1) = 0.5; % make a single impulse

    txt = ({'This routine will calculate your system delay.';'';...
            ['Clicks will be played through Output Channel ',num2str(probeOutput.ch(1)),'.'];'';...
            ['Measurements will be made through Input Channel ',num2str(probeInput.ch),'.'];'';...
            'Place the probe assembly attached to these channels in the freefield away from nearby reflective objects.'});
    choice = questdlg(txt,'Get System Delay','Continue','Cancel','Continue');
    if strcmp(choice,'Continue')
        % do nothing; continue
    elseif strcmp(choice,'Cancel')
        return
    else % if user shut down box using the x
        return
    end

    % 2) LOAD THE STIMULUS ----------------------------------------------------

    nReps = 10; % number of click repetitions to play
    obj.clearRecList % clear out whatever was used previously 
    obj.setRecList(probeInput.ch,probeInput.label,probeInput.micSens,probeInput.gain);
    obj.setNReps(nReps); % number of times to play stimulus
    obj.setFilter(1); % note: this is a highpass filter with a 75 Hz cutoff frequency.

    obj.clearPlayList % clear out whatever was used previously
    obj.setPlayList(stimulus,probeOutput.ch(1));

    obj.objPlayrec.userInfo = []; % clear out previous, non-structure array info
    obj.objPlayrec.userInfo.fs = fs;
    obj.objPlayrec.userInfo.getDelay_version = V;

    delay_nowOrig = obj.objInit.delay_now;
    systemDelayOrig = obj.objPlayrec.systemDelay;
    obj.objInit.delay_now = 0; % set to zero in order to measure
    obj.objPlayrec.systemDelay = 0; % set to zero in order to measure
    
    % 3) PLAYBACK & RECORD ----------------------------------------------------
    obj.objPlayrec.run % run the stimulus
        if obj.killRun
           return
        end    

    % return to original values
    obj.objInit.delay_now = delay_nowOrig;
    obj.objPlayrec.systemDelay = systemDelayOrig;
        
    % 4) RETRIEVE DATA ----------------------------------------------------
    [header0,data0] = obj.retrieveData(['Ch',num2str(probeInput.ch)]); % retrive recorded data

    m = (mean(data0,2));
    [~,maxyIndx] = max(abs(m));
    systemDelay = maxyIndx-1;
    %obj.objInit.delay_now = systemDelay; 
    
    figure; plt(m)
    
    disp(['Your calculated system delay is: ',num2str(systemDelay),' samples, at ',num2str(obj.fs),' Hz sampling rate.'])

    txt = ({['Your calculated system delay is ',num2str(systemDelay),' samples'];'';...
            ['at the current sampling rate (',num2str(obj.fs),' Hz).'];'';...
            'Re-initialize ARLas (press the initialize button) and type this number in the System Delay box.';'';...
            'Then be sure to save your initialization configuration before continuing further.'});
    choice = msgbox(txt,'Get System Delay');

end % end of experiment file

% OLD CODE ----------------------------------------------------------------

    % OLD VERSION --- No longer used!
    %outputChannel = 7;
    %inputChannel = 7;
    % Make a direct electrical connnection between these two channels.
    % Then run the experiment file.
