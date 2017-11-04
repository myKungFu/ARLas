function [] = ARLas_cavityRecordings10x(varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ARLas_cavityRecordings10x(varargin);
%
% Make recordings from ER10x calibrator for calculation of Thevenin source
% parameters.
%
% Auditory Research Lab, The University of Iowa
% Deptartment of Communication Sciences & Disorders
% The University of Iowa
% Author: Shawn Goodman
% Date: June 29, 2015
% Updated: January 31, 2017
% Updated: May 8, 2017  Expanded bandpass filter highcut frequency
% Updated: June 13, 2017 Updates for ER10x system
% Updated: September 8, 2017
% Updated: October 23 - November 4, 2017 Many updates, including filtering and stimuli.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
obj = varargin{1};

% ---  IMPORTANT! UPDATE THE FOLLOWING LINES OF CODE BEFORE RUNNING!!   ---
% -------------------------------------------------------------------------
[inputs,outputs] = hardwareSetup; % read in the saved hardware setup
input = inputs{1}; % choose the probe to use
output = outputs{1}; % choose the probe to use
% -------------------------------------------------------------------------

cavityLengths = [2.9; 3.6; 4.15; 5.2; 6.9]; % cavity lengths in cm (sammi)
    %cavityLengths = [3.15; 3.79; 4.48; 5.532; 7.021]; % cavity lengths in cm (shawn)
% Account for the extra length on tip when inserted into the tube
% Note: use red rubber ER10x tip, fully insert into calibrator cavity
extra = -.4; % (cm); this value is added to the desired length

cavityDiameter = 0.8; % ER10x cavity diameter (cm)
cavityTemperature = 30.5; % (C); read this off of the 10X and update as needed
nReps = 48; % number of stimulus repetitions; recommend ~48

%pathName_micCorr = 'C:\myWork\ARLas\Peripheral\calibrations\micCals\'; % location of mic calibration files

%--------------------------------------------------------------------------
fmin = 200; % in Hz, minimum frequency tested. Note: this is hardcoded later, so don't change this value!
fmax = 20000; % in Hz, maximum frequency tested. Note: this is hardcoded later, so don't change this value!
fs = obj.fs; % sampling rate (Hz)
[stimulus,F,cut1,cut2] = getStimulus(fs); % get logarithmic chirp
fminStim = min(F); % minimum frequency in the stimulus
fmaxStim = max(F); % maximum frequency in the stimulus
stimulus = adjustOutput(F,stimulus); % adjust for better SNR
stimGain_dB = -5; % reduce so you don't overdrive the system; 
                  % -5 dB works for ER10X with 20 dB input gain and no
                  % output gain. If you use amplifiers, may need to reduce
                  % the output by more
stimulus = stimulus * 10^(stimGain_dB/20);
cavityLengths = cavityLengths(:)'; % force these to be a column vector
nCavities = length(cavityLengths); % loop over number of cavity lengths specified
nChannelsOut = length(output.ch); % number of output channnels being calibrated (typically 2)
for jj=1:nCavities % loop over number of cavities
    showInstructions(cavityLengths,extra,jj);
    for kk=1:nChannelsOut % loop over output channels; test output channels one at a time
        % LOAD THE STIMULUS ----------------------------------------------------
        obj.clearRecList % clear the previously used recording list
        obj.setRecList(input.ch,input.label,input.micSens,input.gain); % load the recording info for ARLas to use
        obj.setNReps(nReps); % number of times to play stimulus
        obj.setFilter(1); % turn on default highpass filter

        % PLAYBACK & RECORD ----------------------------------------------------
        obj.clearPlayList % clear the previously used playback list
        obj.setPlayList(stimulus,output.ch(kk)); % load the currently tested output channel
        obj.objPlayrec.run % run the stimulus
        if obj.killRun
           return
        end    

        % RETRIEVE DATA ----------------------------------------------------
        [header,data] = obj.retrieveData(input.label(kk)); % get raw data
        
        if kk == 1
            [thevRecording1(:,jj),stimf] = cleanData(data,obj.fs,stimulus,header,cut1,cut2);
        elseif kk == 2
            [thevRecording2(:,jj),stimf] = cleanData(data,obj.fs,stimulus,header,cut1,cut2);
        end
    end
end

% save recording data -------- --------------------------------------------
recordingParams.fs = obj.fs; % sampling rate (Hz)
recordingParams.cavityTemperature = cavityTemperature; % degrees C
recordingParams.cavityDiameter = cavityDiameter; % (cm)
recordingParams.cavityLengths_nominal = cavityLengths;
recordingParams.fmin = fmin; % calibrate over this range: fmin-fmax
recordingParams.fmax = fmax;
recordingParams.timeStamp = datestr(clock); % time stamp for when this calibration was done
recordingParams.delay_now = obj.objInit.delay_now;
recordingParams.card2volts_now = obj.objInit.card2volts_now;
recordingParams.hostAPI_now = obj.objInit.hostAPI_now;
recordingParams.input = input;
recordingParams.output = output;
recordingParams.thevCalPathName = obj.objPlayrec.savedFilesPath;
recordingParams.cavityRecordingsPathName = obj.objPlayrec.savedFilesPath;
recordingParams.fminStim = fminStim; % minimum frequency in calibration stimulus
recordingParams.fmaxStim = fmaxStim; % maximum frequency in calibration stimulus
recordingParams.stimOrig = stimulus;
recordingParams.stimFilt = stimf;
stimulus = stimf;

channel = 1;
savedFile1 = saveData(recordingParams,stimulus,thevRecording1,channel,recordingParams.cavityRecordingsPathName);
if nChannelsOut == 2
    channel = 2;
    savedFile2 = saveData(recordingParams,stimulus,thevRecording2,channel,recordingParams.cavityRecordingsPathName);
end

load(savedFile1)
t = thevCal(recordingParams,stimulus,recordings);
t.calculate
t.saveThev
if nChannelsOut == 2
    load(savedFile2)
    t = thevCal(recordingParams,stimulus,recordings);
    t.calculate
    t.saveThev
end
end


% internal functions ------------------------------------------------------
function [stim,F,cut1,cut2] = getStimulus(fs)
disp('tic')
    fmin = 100; % starting frequency (hz)
        % find ending frequency:
        %   2 samples = (half a cycle / f) * fs;
        %   2 = (0.5/x)*fs
        %   0.5*fs/2 = x
        %   2400 = x;
    fmax = 0.5*fs/2;
    dur = 0.125; % desired chirp duration is 125 ms
    N = round(fs * dur); % number of samples in recording
    X = (0:1:N-1)';
    X = X / max(X);
    % formula = fmin*exp(-a*X)
    % X(end) = 1, so the value of a is as follows:
    % fmax = 100*exp(-aX)
    % fmax / 100 = exp(-a)
    % -log(fmax / 100) = a
    a = -log(fmax / fmin);
    F = 100*exp(-a*X);
    lead = round(0.5/F(1)*fs); % lead samples for half a cycle
    lag = round(0.25/F(end)*fs); % lead samples for a quarter cycle
    F = [ones(lead,1)*F(1);F;ones(lag,1)*F(end)]; % instantateous frequencies
    % Rotate the phase:
    % At each time sample, time has advanced by the sampling period (1/fs).
    % At each time sample, the phase has rotated by the frequency (f) times
    % the sampling period: phase advance (in cycles) = f * (1/fs), or f/fs.
    % To make this a radian change, multiply by 2pi: phi = (2*pi*f)/fs;
    % Finally, add the newly accumulated phase to the previous phase value.
    stim = zeros(length(F),1);
    phi = 0; % starting phase is zero
    for jj = 1:length(F)
        stim(jj,1) = sin(phi);
        phi = phi + ((2*pi*F(jj))/fs);
    end

    rampN = round(lead*2); % ramp stimulus on
    h = hann(rampN*2);
    h = h(1:rampN);
    stim(1:rampN) = stim(1:rampN) .* h;
    rampN = round(lead*.5); % ramp stimulus off
    h = hann(rampN*2);
    h = h(1:rampN);
    stim = flipud(stim);
    stim(1:rampN) = stim(1:rampN) .* h;
    stim = flipud(stim);
    % zero pad stimulus with 5 ms on both sides. This is to allow for
    % filtering later. Define cut1 and cut2 so you can remove padding after
    % recording
    extra = round(fs*0.002);
    cut2 = length(stim) + extra;
    padLen = 0.050; % zero-pad the stimulus with 5 ms on both sides
    padN = round(padLen * fs);
    cut1 = padN - extra;
    cut2 = cut2 + padN;
    pad = zeros(padN,1);
    stim = [pad;stim;pad];
    F = [pad;F;pad];
    stim = stim / max(abs(stim));  % rescale so 1 is max out
    stim = stim * .99;
end

function [stim] = adjustOutput(F,stim)
   % adjust the output amplitudes to give approximately flat output across
    % the frequency range. These values based on "average" 10x driver
    % response, re: the supplied operating instructions manual from ER.
    % (plus a few emperical tweaks)
    fff =    [100,200,400,600,1000,2000,4000,6000,10000,12000,16000,20000,30000]'; % frequency (Hz)
    output = [ 72, 55, 50, 50,  40,  40,  40,  40,   50,   60,   65,   65,   65]';
    relativeOut = (output - max(output));
    a = zeros(length(stim),1);
    for jj=1:length(stim)
        a(jj,1) = interp1(fff,relativeOut,F(jj),'pchip'); % relative output amplitude (dB)
    end
    a = 10.^(a/20); % put a into linear units
    stim = a .* stim; % amplitude adjusted stimulus
    stim = stim / max(abs(stim));
    stim = stim * .99;
end

function [] = showInstructions(cavityLengths,extra,jj)
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
    txt = ({['1. Place probe into tube of length = ',(num2str((cavityLengths(jj)+extra)*10)),' mm.'];'';...
        '2. Press OK when ready.'});
    uiwait(msgbox(txt,'Thevenin Calibration','modal'));
    try
        delete(h)
    catch
    end % end instructions
end

function [X,stim] = cleanData(X,fs,stim,header,cut1,cut2)
    [rows,cols] = size(X); % orignal size of the recording
    % it is important to do the same processing to the stimulus as to the
    % recorded waveforms.
    stim = repmat(stim,1,cols+2); % replicate stimulus to matrix the same size as recording
    stim = stim(:); % reshape into a long vector
    X = X(:);
    stim = filtfilt(header.SOS,header.G,stim); % apply the playrec highpass filter to stim (same as was for recording)
    [SOS,G] = getLPiir(fs); % get coefficients for a lowpass filter
    stim = filtfilt(SOS,G,stim); % apply the lowpass filter
    X = filtfilt(SOS,G,X);
    X = reshape(X,rows,cols); % reshape into matrices
    stim = reshape(stim,rows,cols+2);
    stim = stim(:,2:end-1);
    tolerance = 'moderate'; doPlot = 0; % perform artifact rejectioin
    rms = sqrt(mean(X.^2,1));
    [indx,nRejects] = AR(rms,tolerance,doPlot);
    disp([num2str(nRejects),' artifacts'])
    X = AR_engine(X,indx);
    stim = AR_engine(stim,indx);
    X = mean(X,2); % return mean    
    stim = mean(stim,2);
    % cut off extra zero padding on the ends
    X = X(cut1:cut2,:);
    stim = stim(cut1:cut2);
    % apply a hann window to the edges
    % fmin = 100; the period is 1/100 = 10 ms. use this for ramp on
    cut1 = round(0.01 * fs);
    % fmax = 20000; however, this period is too short for a window.
    % Instead, use 3 ms (which is actually about 60 periods of 20kHz)
    cut2 = round(0.003 * fs);
    h = hann(cut1*2);
    h = h(1:cut1);
    H = repmat(h,1,size(X,2));
    X(1:cut1,:) = X(1:cut1,:) .* H;
    stim(1:cut1) = stim(1:cut1) .* h;

    h = hann(cut2*2);
    h = h(1:cut2);
    H = repmat(h,1,size(X,2));
    X = flipud(X);
    stim = flipud(stim);
    X(1:cut2,:) = X(1:cut2,:) .* H;
    stim(1:cut2) = stim(1:cut2) .* h;
    X = flipud(X);
    stim = flipud(stim);
end

function [savedFile] = saveData(recordingParams,stimulus,recordings,channel,pathName)
    if nargin < 5
        pathName = recordingParams.cavityRecordingsPathName;
    end
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

function [SOS,G] = getLPiir(fs)
% get low-pass IIR filter coefficients. filter is saved as
% second-order-sections (SOS) and gain (G). Filters are saved for 3
% sampling rates: 44100, 48000, and 96000.
% filters created using Matlab's Filter Designer tool.
% Fpass = 20,000; Fstop = 22000; Apass = 1; Astop = 60; 
% Using a Butterworth fitler and matching stopband exactly, minimum order.
% filter orders are: 57 (96 kHz), 11 (48 kHz), and 3 (44100).
    G = [];
    SOS = [];
    if fs == 96000
      G = [0.366595645667690
       0.348498030887723
       0.332174418047852
       0.317419737340339
       0.304058580128791
       0.291940298698328
       0.280935025159668
       0.270930419417913
       0.261828999009000
       0.253545936085770
       0.246007231598790
       0.239148195722775
       0.232912178259390
       0.227249504156104
       0.222116578200356
       0.217475129962374
       0.213291575607033
       0.209536477606638
       0.206184086914749
       0.203211954998985
       0.200600605426882
       0.198333256568337
       0.196395588510764
       0.194775548549138
       0.193463190668446
       0.192450545324948
       0.191731516591057
       0.191301804385652
       0.437217165824475
       1.000000000000000];
     s1 = [1.000000000000000
       1.000000000000000
       1.000000000000000
       1.000000000000000
       1.000000000000000
       1.000000000000000
       1.000000000000000
       1.000000000000000
       1.000000000000000
       1.000000000000000
       1.000000000000000
       1.000000000000000
       1.000000000000000
       1.000000000000000
       1.000000000000000
       1.000000000000000
       1.000000000000000
       1.000000000000000
       1.000000000000000
       1.000000000000000
       1.000000000000000
       1.000000000000000
       1.000000000000000
       1.000000000000000
       1.000000000000000
       1.000000000000000
       1.000000000000000
       1.000000000000000
       1.000000000000000];
     s2 = [   2.000000000000000
       2.000000000000000
       2.000000000000000
       2.000000000000000
       2.000000000000000
       2.000000000000000
       2.000000000000000
       2.000000000000000
       2.000000000000000
       2.000000000000000
       2.000000000000000
       2.000000000000000
       2.000000000000000
       2.000000000000000
       2.000000000000000
       2.000000000000000
       2.000000000000000
       2.000000000000000
       2.000000000000000
       2.000000000000000
       2.000000000000000
       2.000000000000000
       2.000000000000000
       2.000000000000000
       2.000000000000000
       2.000000000000000
       2.000000000000000
       2.000000000000000
       1.000000000000000];
     s3 = [ 1.000000000000000
       1.000000000000000
       1.000000000000000
       1.000000000000000
       1.000000000000000
       1.000000000000000
       1.000000000000000
       1.000000000000000
       1.000000000000000
       1.000000000000000
       1.000000000000000
       1.000000000000000
       1.000000000000000
       1.000000000000000
       1.000000000000000
       1.000000000000000
       1.000000000000000
       1.000000000000000
       1.000000000000000
       1.000000000000000
       1.000000000000000
       1.000000000000000
       1.000000000000000
       1.000000000000000
       1.000000000000000
       1.000000000000000
       1.000000000000000
       1.000000000000000
                       0];
     s4 = [ 1.000000000000000
       1.000000000000000
       1.000000000000000
       1.000000000000000
       1.000000000000000
       1.000000000000000
       1.000000000000000
       1.000000000000000
       1.000000000000000
       1.000000000000000
       1.000000000000000
       1.000000000000000
       1.000000000000000
       1.000000000000000
       1.000000000000000
       1.000000000000000
       1.000000000000000
       1.000000000000000
       1.000000000000000
       1.000000000000000
       1.000000000000000
       1.000000000000000
       1.000000000000000
       1.000000000000000
       1.000000000000000
       1.000000000000000
       1.000000000000000
       1.000000000000000
       1.000000000000000];
     s5 = [ -0.481608120584464
      -0.457832720237396
      -0.436387881506046
      -0.417004197795273
      -0.399451229314789
      -0.383531065468300
      -0.369073094763771
      -0.355929732518702
      -0.343972912983102
      -0.333091195171439
      -0.323187364227121
      -0.314176435111331
      -0.305983984693643
      -0.298544753309978
      -0.291801468571436
      -0.285703853422153
      -0.280207787731658
      -0.275274598502727
      -0.270870458410870
      -0.266965876119739
      -0.263535264833237
      -0.260556578001090
      -0.258011003108008
      -0.255882706139917
      -0.254158620707033
      -0.252828276971464
      -0.251883666523229
      -0.251319140211792
      -0.125565668351050];
    s6 = [ 0.947990703255224
       0.851824843788288
       0.765085553697455
       0.686683147156629
       0.615685549829951
       0.551292260261611
       0.492813195402444
       0.439651410190355
       0.391288909019103
       0.347274939514517
       0.307216290622282
       0.270769218002433
       0.237632697731205
       0.207542769934396
       0.180267781372859
       0.155604373271651
       0.133374090159791
       0.113420508929278
       0.095606806069865
       0.079813696115680
       0.065937686540765
       0.053889604274438
       0.043593357151063
       0.034984900336470
       0.028011383380817
       0.022630458271254
       0.018809732887457
       0.016526357754398
                       0];
       SOS = [s1,s2,s3,s4,s5,s6];
    elseif fs == 48000
        G = [0.884115938346670
           0.790003589306179
           0.722566286746239
           0.677543695777765
           0.651781204401333
           0.802122266303265
           1.000000000000000];
        s1 = [1.000000000000000
           1.000000000000000
           1.000000000000000
           1.000000000000000
           1.000000000000000
           1.000000000000000]
        s2 = [2.000000000000000
           2.000000000000000
           2.000000000000000
           2.000000000000000
           2.000000000000000
           1.000000000000000];
        s3 = [1.000000000000000
           1.000000000000000
           1.000000000000000
           1.000000000000000
           1.000000000000000
                           0];
        s4 = [1.000000000000000
           1.000000000000000
           1.000000000000000
           1.000000000000000
           1.000000000000000
           1.000000000000000];
        s5 = [1.660622051331147
           1.483852200974787
           1.357185700738689
           1.272620425838064
           1.224231114638748
           0.604244532606530];
        s6 = [0.875841702055535
           0.676162156249930
           0.533079446246268
           0.437554357272998
           0.382893702966584
                           0];
        SOS = [s1,s2,s3,s4,s5,s6];    
    elseif fs == 44100
        G = [0.964424474001281
           0.965605975484490
           1.000000000000000];
        s1 = [1.000000000000000
           1.000000000000000];
        s2 = [2.000000000000000
           1.000000000000000];
        s3 = [1.000000000000000
                           0];
        s4 = [1.000000000000000
           1.000000000000000];
        s5 = [1.926401776975857
           0.931211950968979];
        s6 = [0.931296119029265
                           0];
        SOS = [s1,s2,s3,s4,s5,s6];
    else
        error('The only supported sampling rates are 44100, 48000, and 96000.')
    end
end
% 
% 
% 
% s1 = [1.000000000000000
%    1.000000000000000
%    1.000000000000000
%    1.000000000000000
%    1.000000000000000
%    1.000000000000000
%    1.000000000000000
%    1.000000000000000
%    1.000000000000000
%    1.000000000000000
%    1.000000000000000
%    1.000000000000000
%    1.000000000000000
%    1.000000000000000
%    1.000000000000000
%    1.000000000000000
%    1.000000000000000
%    1.000000000000000
%    1.000000000000000
%    1.000000000000000
%    1.000000000000000
%    1.000000000000000
%    1.000000000000000
%    1.000000000000000
%    1.000000000000000
%    1.000000000000000
%    1.000000000000000
%    1.000000000000000];
%   
% s2 = [2.000000000000000
%    2.000000000000000
%    2.000000000000000
%    2.000000000000000
%    2.000000000000000
%    2.000000000000000
%    2.000000000000000
%    2.000000000000000
%    2.000000000000000
%    2.000000000000000
%    2.000000000000000
%    2.000000000000000
%    2.000000000000000
%    2.000000000000000
%    2.000000000000000
%    2.000000000000000
%    2.000000000000000
%    2.000000000000000
%    2.000000000000000
%    2.000000000000000
%    2.000000000000000
%    2.000000000000000
%    2.000000000000000
%    2.000000000000000
%    2.000000000000000
%    2.000000000000000
%    2.000000000000000
%    2.000000000000000];
% s3 = [1.000000000000000
%    1.000000000000000
%    1.000000000000000
%    1.000000000000000
%    1.000000000000000
%    1.000000000000000
%    1.000000000000000
%    1.000000000000000
%    1.000000000000000
%    1.000000000000000
%    1.000000000000000
%    1.000000000000000
%    1.000000000000000
%    1.000000000000000
%    1.000000000000000
%    1.000000000000000
%    1.000000000000000
%    1.000000000000000
%    1.000000000000000
%    1.000000000000000
%    1.000000000000000
%    1.000000000000000
%    1.000000000000000
%    1.000000000000000
%    1.000000000000000
%    1.000000000000000
%    1.000000000000000
%    1.000000000000000];
% s4 = [1.000000000000000
%    1.000000000000000
%    1.000000000000000
%    1.000000000000000
%    1.000000000000000
%    1.000000000000000
%    1.000000000000000
%    1.000000000000000
%    1.000000000000000
%    1.000000000000000
%    1.000000000000000
%    1.000000000000000
%    1.000000000000000
%    1.000000000000000
%    1.000000000000000
%    1.000000000000000
%    1.000000000000000
%    1.000000000000000
%    1.000000000000000
%    1.000000000000000
%    1.000000000000000
%    1.000000000000000
%    1.000000000000000
%    1.000000000000000
%    1.000000000000000
%    1.000000000000000
%    1.000000000000000
%    1.000000000000000];
% s5 = [-0.481959636952425
%   -0.457777107877400
%   -0.436003737902187
%   -0.416356720316669
%   -0.398594778970428
%   -0.382511206019014
%   -0.367928224476489
%   -0.354692397743366
%   -0.342670871832577
%   -0.331748283932745
%   -0.321824207331955
%   -0.312811030541027
%   -0.304632189860955
%   -0.297220691216657
%   -0.290517869996617
%   -0.284472347763927
%   -0.279039152688589
%   -0.274178976883205
%   -0.269857548875725
%   -0.266045103506594
%   -0.262715934810984
%   -0.259848020107523
%   -0.257422705693709
%   -0.255424446347227
%   -0.253840592332429
%   -0.252661218876758
%   -0.251878994164639
%   -0.251489082839072];
% s6 = [0.947090773667534
%    0.849394668774468
%    0.761431435881228
%    0.682058551229021
%    0.610301272266044
%    0.545324510520333
%    0.486410056618010
%    0.432938035025320
%    0.384371722845500
%    0.340245060582214
%    0.300152329770248
%    0.263739584748647
%    0.230697512329133
%    0.200755460084652
%    0.173676426168846
%    0.149252844485316
%    0.127303031281827
%    0.107668184826835
%    0.090209850233428
%    0.074807777872585
%    0.061358117041661
%    0.049771897303251
%    0.039973758711731
%    0.031900899412694
%    0.025502215160733
%    0.020737610413401
%    0.017577465033732
%    0.016002244441899];
% SOS = [s1,s2,s3,s4,s5,s6];
% G = [0.366282784178777
%    0.347904390224267
%    0.331356924494760
%    0.316425457728088
%    0.302926623323904
%    0.290703326125330
%    0.279620458035380
%    0.269561409320488
%    0.260425212753231
%    0.252124194162367
%    0.244582030609573
%    0.237732138551905
%    0.231516330617045
%    0.225883692216999
%    0.220789639043057
%    0.216195124180347
%    0.212065969648310
%    0.208372301985908
%    0.205088075339426
%    0.202190668591498
%    0.199660545557669
%    0.197480969298932
%    0.195637763254505
%    0.194119113266367
%    0.192915405707076
%    0.192019097884161
%    0.191424617717273
%    0.191128290400707
%    1.000000000000000];
% end
% 
