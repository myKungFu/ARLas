function [] = ARLas_getCard2Volts(varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ARLas__getCard2Volts(varargin)
%
% Calculate the multiplier to convert card units to voltage.
% Requires electrical measurement using a voltmeter or oscilloscope.
%
% Auditory Research Lab, The University of Iowa
% Deptartment of Communication Sciences & Disorders
% The University of Iowa
% Author: Shawn S. Goodman, PhD
% Date: January 9, 2017
% Updated: January 9, 2017
% Updated: August 15, 2017  Updated a few lines to make compatible with
%                           current ARLas version
% Updated: August 24, 2017  Fixed a couple more bugs: Changed tone
%                           frequency to 250 Hz (was 100) so below cutoff of filter (125 Hz);
%                           Also, turned filter off for good measure; not
%                           need for a direct in/out electrical connection
%                           anyway.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%--------------------- SET THESE VALUE BEFORE RUNNING! --------------------
outputChannel = 8;
inputChannel = 8;
% Make a direct electrical connnection between these two channels.
% Then run the experiment file.
%--------------------------------------------------------------------------

obj = varargin{1}; % get the arlas object

% 1) CREATE THE STIMULUS --------------------------------------------------
fs = obj.fs; % get the system sampling rate
len = 1; % desired stimulus length (s)
nReps = 10; % number of stimulus repetitions
    nSamples = round(len * fs); % number of samples in stimulus
    if mod(nSamples,2) ~= 0 % force to be even number
        nSamples = nSamples + 1;
    end
t = (0:1:nSamples-1)'/obj.fs; % time in seconds
f = 250; % frequency in Hz
a = .5; % amplitude (full out is 1)
stimulus = a * cos(2*pi*f*t);
            
txt = ({'This routine will calculate a multiplier to convert card units to voltage.';'';...
            ['A 250 Hz pure tone will be played through Output Channel ',num2str(outputChannel),'.'];'';...
            'The tone will be played for 10 seconds.';'';...
            'Use a voltmeter to measure the open circuit AC voltage of the channel.';'';...
            'Write this value down. You will enter it in the next step.'});
    choice = questdlg(txt,'Get Card to Volts','Continue','Cancel','Continue');
if strcmp(choice,'Continue')
    % do nothing; continue
elseif strcmp(choice,'Cancel')
    return
else % if user shut down box using the x
    return
end

% 2) LOAD THE STIMULUS ----------------------------------------------------
obj.clearRecList % clear out whatever was used previously 
obj.clearPlayList % clear out whatever was used previously

% Load Output:
obj.setPlayList(stimulus,outputChannel);
obj.setNReps(nReps); % number of times to play stimulus
obj.setFilter(1);

% Load Input:
label = 'noInput';
micSens = 1;
gain = 0;
obj.setRecList(inputChannel,label,micSens,gain);

objInit.card2volts_now = 1; % set to 1 in order to measure
obj.objPlayrec.card2volts = 1; % set to 1 in order to measure

% 3) PLAYBACK & RECORD ----------------------------------------------------
obj.objPlayrec.run % run the stimulus
if obj.killRun
   return
end    

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

txt = ({['Connect Output Channel ',num2str(outputChannel),' directly to Input Channel ',num2str(inputChannel),'.'];'';...
        'A 10-second, 100 Hz pure tone will again be played.';''});
choice = questdlg(txt,'Get Card to Volts','Continue','Cancel','Continue');
if strcmp(choice,'Continue')
    % do nothing; continue
elseif strcmp(choice,'Cancel')
    return
else % if user shut down box using the x
    return
end

obj.objPlayrec.run % run the stimulus
if obj.killRun
   return
end    

% 4) RETRIEVE DATA ----------------------------------------------------
[header1,data1] = obj.retrieveData(['Ch',num2str(inputChannel)]);

m = (mean(data1,2));
% note that multimeter voltage is assumed to be RMS, not peak!
closedCircuitVoltage = sqrt(mean(m.^2)); % therefore, this is the correct comparison
% take into account current input settings
closedCircuitVoltage = closedCircuitVoltage / obj.objPlayrec.card2volts * micSens * 10^(gain/20);
cardMultiplier = openCircuitVoltage / closedCircuitVoltage; 
disp(['Your calculated card to volts is: ',num2str(cardMultiplier),'.'])

txt = ({['Your calculated Card to Volts is ',num2str(cardMultiplier),'.'];'';...
        'Re-initialize ARLas and type this number in the Card to Volts box.';'';...
        'Then be sure to save your initialization configuration before continuing further.'});
choice = msgbox(txt,'Get Card to Volts');
            
end % end of experiment file

% OLD CODE
% if isempty(openCircuitVoltage) % but that something was not numeric
%     errorTxt = {'  Issue: Illegal input for voltage; possibly non-numeric value.'
%          '  Action: getCard2V terminated prematurely.'
%          '  Location: initARLas.getCard2V.'
%         };
%     errorMsgARLas(errorTxt);
%     objInit.obj.buttonManager(21)
%     return
% end
