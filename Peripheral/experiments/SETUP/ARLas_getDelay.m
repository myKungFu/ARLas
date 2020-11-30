function [] = ARLas_getDelay(varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ARLas_getDelay(varargin)
%
% Measure sample delay between card output and input
%
% Auditory Research Lab, The University of Iowa
% Deptartment of Communication Sciences & Disorders
% The University of Iowa
% Author: Shawn S. Goodman, PhD
% Date: January 9, 2017
% Last Updated: August 15, 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%--------------------- SET THESE VALUES BEFORE RUNNING! --------------------
outputChannel = 3; %1;
inputChannel = 3; %1;
% Make a direct electrical connnection between these two channels.
% Then run the experiment file.
%--------------------------------------------------------------------------

obj = varargin{1}; % get the arlas object

% 1) CREATE THE STIMULUS --------------------------------------------------
fs = obj.fs; % get the system sampling rate
len = 0.5; % desired stimulus length (s)
    nSamples = round(len * fs); % number of samples in stimulus
    if mod(nSamples,2) ~= 0 % force to be even number
        nSamples = nSamples + 1;
    end
stimulus = zeros(nSamples,1); % create a vector of zeros
stimulus(1,1) = 0.25; % make a single impulse

txt = ({'This routine will calculate your system delay.';'';...
        ['Clicks will be played through Output Channel ',num2str(outputChannel),'.'];'';...
        ['Measurements will be made through Input Channel ',num2str(inputChannel),'.'];'';...
        'A direct elecrical connection between output and input is recommended.'});
choice = questdlg(txt,'Get System Delay','Continue','Cancel','Continue');
if strcmp(choice,'Continue')
    % do nothing; continue
elseif strcmp(choice,'Cancel')
    return
else % if user shut down box using the x
    return
end

% 2) LOAD THE STIMULUS ----------------------------------------------------

% Load Output:
obj.setPlayList(stimulus,outputChannel);

% Load Input:
label = 'directInput';
micSens = 1;
gain = 0;
obj.setRecList(inputChannel,label,micSens,gain);


obj.setNReps(10); % number of times to play stimulus
obj.objInit.delay_now = 0; % set to zero in order to measure
obj.objPlayrec.systemDelay = 0; % set to zero in order to measure

% 3) PLAYBACK & RECORD ----------------------------------------------------
obj.objPlayrec.run % run the stimulus
    if obj.killRun
       return
    end    

% 4) RETRIEVE DATA ----------------------------------------------------
[header0,data0] = obj.retrieveData(['Ch',num2str(inputChannel)]);

m = (mean(data0,2));
[~,maxyIndx] = max(abs(m));
systemDelay = maxyIndx-1;
%obj.objInit.delay_now = systemDelay; 
disp(['Your calculated system delay is: ',num2str(systemDelay),' samples, at ',num2str(obj.fs),' Hz sampling.'])

txt = ({['Your calculated system delay is ',num2str(systemDelay),' samples'];'';...
        ['at the current sampling rate (',num2str(obj.fs),' Hz).'];'';...
        'Re-initialize ARLas and type this number in the System Delay box.';'';...
        'Then be sure to save your initialization configuration before continuing further.'});
choice = msgbox(txt,'Get System Delay');

end % end of experiment file