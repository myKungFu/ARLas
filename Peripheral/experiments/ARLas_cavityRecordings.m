function [] = ARLas_cavityRecordings(varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ARLas_cavityRecordings(varargin);
%
% Make recordings from brass tubes for calculation of Thevenin source calibration.
%
% Auditory Research Lab, The University of Iowa
% Deptartment of Communication Sciences & Disorders
% The University of Iowa
% Author: Shawn Goodman
% Date: June 29, 2015;
% Last Updated: January 16, 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
obj = varargin{1};

% IMPORTANT! UPDATE THE FOLLOWING LINES OF CODE BEFORE RUNNING!! -----------------------------------------
cavityLengths = [1.85, 2.56, 4, 5.4, 8.3]; % intended cavity lengths (cm);
cavityDiameter = 0.8; % cavity diameter (cm)
cavityTemperature = 24.5; % (C); 25 = room temp; 37 = body temp
fmin = 200; % in Hz
fmax = 8000; % in Hz
%extra = 0.5; % 0.5 accounts for the standard rubber ER10B+tip, assuming you insert it into the tube up to the cap.
extra = 1.68; %1.86; %2.047; % 1.86 accounts for the standard foam ER10C tip, assuming you insert it into the tube as far as possible

% ----- Specify Input:
input.label = 'ER10C'; % label that will be used for naming saved files 
input.micSens = 0.05; % microphone sensitivity in V/Pa
input.gain = 20; % amplifier gain in dB 
input.ch = 1; % channel on the sound card through which the recording is coming

% ----- Specify Output:
output.label = {'ER10C_Ch1','ER10C_Ch2'}; % label to identify the transducer being calibrated
output.gain = [20,20]; % gain in dB; this value not used, but just documents the physical configuration of the probe, 
output.ch = [1,2]; % channel on the sound card through which the playback is sent
output.receiver = [1,2]; % which receiver on the ER10C is being used
% -------------------------------------------------------------------------------------------------------

stimulus = getStimulus(obj.fs,fmin,fmax); % get stimulus
nCavities = length(cavityLengths); % loop over number of cavity lengths specified
nChannelsOut = length(output.ch); % number of output channnels being calibrated
for ii=1:nChannelsOut
    thevRecordings.(matlab.lang.makeValidName(output.label{ii})) = zeros(size(stimulus,1),nCavities);
end
for jj=1:nCavities % loop over cavities
    showInstructions(cavityLengths,extra,jj,nCavities); % show the user the lengths to use
    for kk=1:nChannelsOut % loop over output channels; test output channels one at a time
        % LOAD THE CHANNELS ----------------------------------------------------
        obj.clearRecList % clear any previously used recording list
        obj.setRecList(input.ch,input.label,input.micSens,input.gain) % load the recording info for ARLas to use
        obj.clearPlayList % clear any previously used playback list
        obj.setPlayList(stimulus,output.ch(kk)) % load the currently tested output channel

        % PLAYBACK & RECORD ----------------------------------------------------
        obj.objPlayrec.nReps = 32; %128; % number of times to play stimulus
        obj.objPlayrec.run % run the stimulus
        ok = obj.checkForErrors;
        if ~ok
           return
        end   
        
        % RETRIEVE DATA ----------------------------------------------------
        [header,data] = obj.retrieveData(input.label); % get raw data
        d = cleanData(data,obj.fs);
        thevRecordings.(matlab.lang.makeValidName(output.label{kk}))(:,jj) = d;
    end
end

% save recording data -------- --------------------------------------------
recordingParams.fs = obj.fs; % sampling rate (Hz)
recordingParams.cavityTemperature = cavityTemperature; % degrees C
recordingParams.cavityDiameter = cavityDiameter; % (cm)
recordingParams.cavityLengths_nominal = cavityLengths;
recordingParams.fmin = fmin;
recordingParams.fmax = fmax;
recordingParams.timeStamp = datestr(clock); % time stamp for when this calibration was done
recordingParams.delay_now = obj.objPlayrec.systemDelay;
recordingParams.card2volts_now = obj.objPlayrec.card2volts; % got to here!!!
recordingParams.output = output;
recordingParams.input = input;
recordingParams.thevCalPathName = obj.map.calibrations;
sep = filesep;
recordingParams.cavityRecordingsPathName = [obj.map.calibrations,'cavityRecordings',sep];

for kk=1:nChannelsOut
    savedFiles.(matlab.lang.makeValidName(output.label{kk})) = saveData(recordingParams,stimulus,thevRecordings.(matlab.lang.makeValidName(output.label{kk})),output.ch(kk));
end

try % Try calculating thevenin source. Note that this can also be done manually, later on.
   for kk=1:nChannelsOut 
        dummy = load(savedFiles.(matlab.lang.makeValidName(output.label{kk})));
        t = thevCal(dummy.recordingParams,dummy.stimulus,dummy.recordings);
        t.calculate
        t.saveThev
   end
catch ME
    
end
end


% internal functions ------------------------------------------------------
function [] = showInstructions(cavityLengths,extra,jj,nCavities)
    gui.height = 143;
    gui.width = 296;
    scrsz = get(0,'ScreenSize'); % get the current screen size
    gui.left = round(scrsz(4) * .05);
    gui.bottom = scrsz(4) * .1;    
    %[left, bottom,width, height]
    rect = [gui.left,gui.bottom+125,gui.width,gui.height]; % give instructions
        h = figure('Position',rect);
        set(h,'Name',' ','NumberTitle','off','MenuBar','none','Resize','off','Color',[1 1 1])
        h1 = uicontrol('Parent',h,'Style','togglebutton','BackgroundColor',[1 1 1],...
            'Position',[1 1 296 194],'Visible','on');
        set(h1,'CData',imread('thevSetup.jpg'))
    if jj < nCavities
        txt = ({['1. Place probe into tube of length = ',num2str(cavityLengths(jj)+extra),' cm.'];'';...
            ['2. (Next length will be ',num2str(cavityLengths(jj+1)+extra),' cm.)'];'';...
            '3. Press OK when ready.'});
    else
        txt = ({['1. Place probe into tube of length = ',num2str(cavityLengths(jj)+extra),' cm.'];'';...
            '2. Press OK when ready.'});
    end
    uiwait(msgbox(txt,'Thevenin Calibration','modal'));
    try
        delete(h)
    catch
    end % end instructions
end

function [X] = cleanData(X,fs)
    [rows,cols] = size(X);
    X = X(:);
    b = bpf(fs); % get bandpass filter
    X = fastFilter(b,X);
    X = reshape(X,rows,cols);
    X = X(:,2:end-1);
    % artifact reject
    tolerance = 'moderate'; doPlot = 0;
    rms = sqrt(mean(X.^2,1));
    [indx,nRejects] = AR(rms,tolerance,doPlot);
    disp([num2str(nRejects),' artifacts'])
    Xtract = AR_engine(X,indx);
    X = mean(X,2); % return mean
    % ramp ends
    rampLen = 0.001; % ramp length in ms
    rampN = round(fs * rampLen); % ramp in samples
    h = hann(rampN*2);
    h = h(1:rampN);
    X(1:rampN,:) = X(1:rampN,:) .* h;
end

function [savedFile] = saveData(recordingParams,stimulus,recordings,channel)
    pathName = recordingParams.cavityRecordingsPathName;
    if ~exist(pathName,'dir')
        mkdir(pathName)
        addpath(genpath(pathName))
    end
    try
        [fileName,pathName] = uiputfile([pathName,'newFPLcalibration_Ch',num2str(channel),'.mat'],'Save FPL Calibration.');
        if fileName == 0 % user chose not to create a file
            dlgtxt = {'You chose NOT to save the current calibration recordings.' ' ' 'Are you sure you want to do this?'};
            choice = questdlg(dlgtxt,'Save Initialization','Yes','No','No');
            if strcmp(choice,'Yes') == 1
                return
            else % give one more chance
                [fileName,pathName] = uiputfile([pathName,'newFPLcalibration_Ch',num2str(channel),'.mat'],'Save FPL Calibration.');
            end
        end
        savedFile = fileName;
        if isempty(fileName) % user chose not to create a file
            return
        end

        if exist(pathName,'dir') == 0 % if correct folder does not exist
            success = mkdir(pathName); % try to create it
            if success ~= 1
                errorTxt = {'  Issue: Error creating the directory for saving cavity recordings.'
                     '  Action: Cavity recordings NOT saved.'
                     '  Location: cavityRecordings.'
                    };
                errorMsgARLas(errorTxt);
                return
            end
        end
        try % save file here ----------------------------------------------
            recordingParams.fileName = fileName;
            save([pathName,fileName],'recordingParams','stimulus','recordings');
        catch
            alertTxt = {'  Issue: Error saving cavity recordings.'
                 '  Action: Cavity recordings NOT saved.'
                 '  Location: cavityRecordings.'
                };
            alertMsgARLas(alertTxt);
        end
    catch ME
        errorTxt = {'  Issue: Error saving cavity recordings.'
             '  Action: Cavity recordings NOT saved..'
             '  Location: cavityRecordings.'
            };
        errorMsgARLas(errorTxt);
        objInit.printError(ME)
    end
end

function [stimulus] = getStimulus(samplingRate,fmin,fmax)
    len = 0.075; % stimulus length in seconds
    N = round(len * samplingRate); % total number of samples
    if mod(N,2)~=0
        N = N + 1;
    end
    rampCycles = 4; % ramps this number of cycles of the stimulus frequency at onset and offset 
    rampN_on = round((rampCycles / fmin) * samplingRate);
    if mod(rampN_on,2)~=0
        rampN_on = rampN_on + 1;
    end
    rampN_off = round((rampCycles / fmax) * samplingRate);
    if mod(rampN_off,2)~=0
        rampN_off = rampN_off + 1;
    end
    edgeExtra = 0.1;
    fmin = fmin-(edgeExtra*fmin);
    fmax = fmax+(edgeExtra*fmax);
    f = [fmin,fmax]; % frequency minimum and maximum (Hz)
    stepSize = (f(2)-f(1))/(N-1); % frequency step size (Hz)
    F = (f(1):stepSize:f(end))';
    fOn = ones(rampN_on,1)*f(1);
    fOff = ones(rampN_off,1)*f(end);
    F = [fOn;F;fOff];
    N = length(F);
    phi = 0;
    for jj = 1:N
        p(jj,1) = cos(phi);
        phi = phi + ((2*pi) / (1/F(jj) * samplingRate));
    end
    rampOn = hann(rampN_on*2);
    rampOn = rampOn(1:rampN_on);
    rampOff = hann(rampN_off*2);
    rampOff = rampOn(1:rampN_off);
    p(1:rampN_on) = p(1:rampN_on) .* rampOn;
    p = flipud(p);
    p(1:rampN_off) = p(1:rampN_off) .* rampOff;
    p = flipud(p);
    if fmax > 10000 % increase high-frequency output above 10 kHz
        [dummy,indx1] = min(abs(F-10000));
        h = hann((N-indx1)*2);
        h = h(1:N-indx1);
        dBincrease = 20;
        multiplier = 10^(dBincrease/20);
        h = h * multiplier;
        amplitude = ones(N,1);
        amplitude(indx1+1:end) = amplitude(indx1+1:end) + h;
        p = p .* amplitude;
        p = p / (max(abs(p)));
        segmt = p * .95;
    else
        desiredGain_dB = -30;
        multiplier = 10.^(desiredGain_dB / 20);
        segmt = p * multiplier;
    end
    padLen = 0.005;
    padN = round(padLen * samplingRate);
    pad = zeros(padN,1);
    Input = [pad;segmt;pad];
    stimulus = Input;
end

function b = bpf(Fs)
    % Fs = Sampling Frequency
    Fstop1 = 50;             % First Stopband Frequency
    Fpass1 = 200;             % First Passband Frequency
    Fpass2 = 18000;           % Second Passband Frequency
    Fstop2 = 19000;           % Second Stopband Frequency
    Dstop1 = 0.0001;          % First Stopband Attenuation
    Dpass  = 0.057501127785;  % Passband Ripple
    Dstop2 = 0.0001;          % Second Stopband Attenuation
    flag   = 'scale';         % Sampling Flag
    [N,Wn,BETA,TYPE] = kaiserord([Fstop1 Fpass1 Fpass2 Fstop2]/(Fs/2), [0 ...
                                 1 0], [Dstop1 Dpass Dstop2]);
    b  = fir1(N, Wn, TYPE, kaiser(N+1, BETA), flag)';
end