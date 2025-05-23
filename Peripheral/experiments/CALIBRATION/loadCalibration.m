function [LC,iscCheckHandleA,iscCheckHandleB,OK] = loadCalibration(varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [LC,iscCheckHandleA,iscCheckHandleB,OK] = loadCalibration(obj,doProbeFitCheck,iscCheckHandleA,iscCheckHandleB);
%
% Perform Thevenin-based load calibration.
% Save the results as well as return them as output arguments.
%
% Use LC = getMostRecentLoadCalibration(obj); to get the most current load
% calibration.
%
% Author: Shawn Goodman, PhD
% Auditory Research Lab, the University of Iowa
% Original Date: May 20, 2025
% Last Updated: May 20, 2025 -- ssg
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    obj = varargin{1}; % get the arlas object
    V = '20MAY2025'; % this is the current version number of this program

%------ USER MODIFIABLE PARAMETERS ----------------------------------------
%--------------------------------------------------------------------------

    % specify which ER10X probes to use    
    probeA = 'A'; % set to 'A' if available, or or [] if not
    probeB = []; % set to 'B' if available, or or [] if not
                  
    % additional parameters that user usually doesn't need to mess with
    doMicCorrection = 1; % turn on and off microphone correction 

    loadCalibrationsPath = [obj.map.data,obj.subjectID,'\','loadCals\'];
    % where to store load calibrations for each subject
    tempResultsPath = [obj.map.data,'TEMP\']; % where to store the results for each subject temporarily, prior to upload


%------ END USER MODIFIABLE PARAMETERS ------------------------------------
%--------------------------------------------------------------------------

    if nargin < 2 % if running program directly
        doProbeFitCheck = 1; 
        LC = [];
        iscCheckHandleA = [];
        iscCheckHandleB = [];
    else % otherwise these should usually be passed as inputs
        doProbeFitCheck = varargin{2}; 
        LC = varargin{3};
        iscCheckHandleA = varargin{4};
        iscCheckHandleB = varargin{5};
    end
    OK = 1; % initialize OK to 1 (probe fit ok)

    if exist(tempResultsPath,'dir') == 0 % if path does not exist
        success = mkdir(tempResultsPath);
        addpath(tempResultsPath)
        if ~success
            warning('Unable to create tempResults folder: data save path')
        end
    end 

    disp(' ')
    disp('----- Starting load calibration -----')

    if isempty(LC) % if a previous load calibration is not being used
        LC.loadCalibrationsPath = loadCalibrationsPath;
        LC.tempResultsPath = tempResultsPath;
        % get new probe settings
        [inputs,outputs] = hardwareSetup; % read in the hardware setup
        % For Probe A
        if isempty(probeA)
            probeInputA = []; 
            probeOutputA = [];
        elseif strcmp(probeA,'A')
            probeInputA = inputs{1};   % IHS_A microphone
            probeOutputA = outputs{1}; % IHS_A loudspeakers
        else
            error('Unexpected value for probe A')
        end
        % For Probe B
        if isempty(probeB)
            probeInputB = []; 
            probeOutputB = [];
        elseif strcmp(probeB,'B')
            probeInputB = inputs{2};   % IHS_A microphone
            probeOutputB = outputs{2}; % IHS_A loudspeakers
        else
            error('Unexpected value for probe B')
        end
    
        % Get which ear is being tested (regardless of which probe is being used)
        if ~(isempty(probeA) | isempty(probeB)) % if both ears are being tested
             testEar = 'Both';
        else
            testEar = getTestEar; % Get stimulus ear:
        end
                
        % Get most recent calibration files:
        % Get microphone calibrations
        micCorrectionA = getMicCorrection(probeInputA,doMicCorrection,obj.map);
        micCorrectionB = getMicCorrection(probeInputB,doMicCorrection,obj.map);
        if doMicCorrection==1
            if ~isempty(probeInputA)
                if micCorrectionA == 1
                    warning('PROBE A MIC CORRECTION IS NOT BEING APPLIED!')
                end
            end
            if  ~isempty(probeInputB)
                if micCorrectionB == 1
                    warning('PROBE B MIC CORRECTION IS NOT BEING APPLIED!')
                end
            end
        end
        % Get Thevenin source calibrations
        [CA1,CA2,calPathA1,calPathA2,fminA,fmaxA] = getThevSource(probeOutputA,obj.map);
        [CB1,CB2,calPathB1,calPathB2,fminB,fmaxB] = getThevSource(probeOutputB,obj.map);
        fmin = max([fminA,fminB]);
        fmax = max([fmaxA,fmaxB]);

        % package the data for saving and returning
        LC.testEar = testEar;
        LC.probeA = probeA;
        LC.probeB = probeB;
        LC.doMicCorrection = doMicCorrection;
        LC.probeInputA = probeInputA;
        LC.probeOutputA = probeOutputA;
        LC.probeInputB = probeInputB;
        LC.probeOutputB = probeOutputB;
        LC.micCorrectionA = micCorrectionA;
        LC.micCorrectionB = micCorrectionB;
        LC.CA1 = CA1;
        LC.CA2 = CA2;
        LC.calPathA1 = calPathA1;
        LC.calPathA2 = calPathA2;
        LC.CB1 = CB1;
        LC.CB2 = CB2;
        LC.calPathB1 = calPathB1;
        LC.calPathB2 = calPathB2;
        LC.fmin = fmin;
        LC.fmax = fmax;
    end

    % PERFORM IN-SITU (LOAD) CALIBRATION ----------------------------------
    inSituReps = 6; % number of times to repeat the calibration chirp
    doIndividual = 1; % do individual loudspeakers (i.e., not played together simultaneously)
    doSimultaneous = 0; % do simultaneous loudspeaker presentation
    [iscSA1,iscSA2,iscSA12] = ARLas_runISC(obj,LC.probeInputA,LC.probeOutputA,LC.calPathA1,LC.calPathA2,'thev',LC.fmin,LC.fmax,LC.micCorrectionA,inSituReps,doIndividual,doSimultaneous);
    [iscSB1,iscSB2,iscSB12] = ARLas_runISC(obj,LC.probeInputB,LC.probeOutputB,LC.calPathB1,LC.calPathB2,'thev',LC.fmin,LC.fmax,LC.micCorrectionB,inSituReps,doIndividual,doSimultaneous);
    % (overwrite these LC fields with new data if they already exist)
    LC.iscSA1 = iscSA1;
    LC.iscSA2 = iscSA2;
    LC.iscSB1 = iscSB1;
    LC.iscSB2 = iscSB2;
    LC.timeStamp = datetime;
    
    if ~exist(loadCalibrationsPath,'dir')
        mkdir(loadCalibrationsPath)
        addpath(loadCalibrationsPath)
    end
    fileName = datestr(now, 'dd-mm-yyyy-HH-MM-ss');
    save([loadCalibrationsPath,fileName,'.mat'],'LC')


% ----- 
% insert files for fPTC here
fff = LC.iscSA1.freq;
fofpl = LC.iscSA1.fpl;
FFF = [125,200,300,400,500,600,700,800,900,1000,1500,2000,3000,4000,5000,6000,8000,10000,12000,14000,16000,18000,20000]';
for ii=1:length(FFF)
    FO(ii,1) = interp1(fff,fofpl,FFF(ii),'pchip');    
end
FO(end) = FO(end-1);
indx = find(FFF==1000);
M = [FFF,round(FO)];
PTCpath = 'C:\Users\Public\PSYCHOACOUSTICS\Headphones\';
PTCname = obj.subjectID;
writematrix(M,[PTCpath,PTCname,'.txt'],'Delimiter','space')

txt = num2str(round(FO(indx)));
cprintf([1,0.5,0],['Value for fPTC','\n']);
cprintf([1,0.5,0],[txt,'\n']);

% ----

    if doProbeFitCheck == 1
        [iscCheckHandleA,iscCheckHandleB,OK] = probeFitCheck(obj,LC,iscCheckHandleA,iscCheckHandleB);
    end

    disp('----- Finished load calibration -----')
    disp(' ')

end

% INTERNAL FUNCTIONS ------------------------------------------------------
function  [testEar] = getTestEar(probeOutput)
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
function [micCorrection] = getMicCorrection(probeInput,doMicCorrection,map)
    if isempty(probeInput)
        micCorrection = [];
        return
    end
    [pathName_mic,folderName_mic,fileName_mic] = mostRecentCalibration('mic',probeInput.label,probeInput.ch,map);
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
function [C1,C2,calPath1,calPath2,fmin,fmax] = getThevSource(probeOutput,map)
    C1 = []; C2 = []; calPath1 = []; calPath2 = []; fmin = []; fmax = [];
    if isempty(probeOutput)
        return
    end
    % Assume there are two recievers for each probe.
    % Get location of most recent Thev source files: Receiver 1
    [pathName,folderName,fileName] = mostRecentCalibration('thev',probeOutput.label,1,map);
    calPath1.pathName = pathName;
    calPath1.folderName = folderName;
    calPath1.fileName = fileName;
    dummy = load([calPath1.pathName,calPath1.folderName,calPath1.fileName]);
    C1 = dummy.t;

    % Get location of most recent Thev source files: Receiver 2
    [pathName,folderName,fileName] = mostRecentCalibration('thev',probeOutput.label,2,map);
    calPath2.pathName = pathName;
    calPath2.folderName = folderName;
    calPath2.fileName = fileName;
    dummy = load([calPath2.pathName,calPath2.folderName,calPath2.fileName]);
    C2 = dummy.t;

    fmin = C1.fmin;
    fmax = C1.fmax;
end

% OLD CODE ----------------------------------------------------------------

    % % specify calibration to use to set stimulus levels
    % calType = 'thev'; % use 'thev' for Thevenin source
    % targetCalType = 'fpl'; % 'fpl', 'spl', 'ipl'.

    % % Define which PROBE is being used (regardless of which ear it is placed in)
    % if ~(isempty(probeA) | isempty(probeB)) % if both ears are being tested
    %     testProbe = 'Both';
    % elseif ~(isempty(probeA))
    %     testProbe = 'Left';
    % elseif ~(isempty(probeB))
    %     testProbe = 'Right';
    % end

    % Define which EAR the probe is placed in (regardless of probe used).
    % Note that if both probes are used, it is assumed that A=Left and
    % B=Right
    % ISCA.(['Rec',num2str(1)]).(['S1']) = iscSA1;
    % ISCA.(['Rec',num2str(1)]).(['S2']) = iscSA2;
    % ISCB.(['Rec',num2str(1)]).(['S1']) = iscSB1;
    % ISCB.(['Rec',num2str(1)]).(['S2']) = iscSB2;            
