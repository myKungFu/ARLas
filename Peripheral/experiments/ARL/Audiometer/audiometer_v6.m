function [] = audiometer_v6(varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% audiometer_v6(varargin)
%
% Measure pure-tone behavioral thresolds (air conduction), version 6
% This version uses Thevenin source calibration.
%
% Author: Shawn Goodman, PhD
% Date: July 9, 2021
% Last Updated: July 9, 2021 -- ssg -- airport leaving USF
% Last Updated: July 22, 2021 -- ssg
% Last Updated: July 29, 2021 -- ssg -- v2
% Last Updated: July 30, 2021 -- ssg -- v3 -- interfacing with arlas_audiometer object
% Last Updated: August 1, 2021 -- ssg -- v4 -- allowing loading of old data
%                                              and "virtual" running option
% Last Updated: October 14, 2024 -- ssg -- v6 making this work!
%
% things to do: make audiometer object save itself periodically
% No resonse button
% plot audiogram button
% new calibration routine
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    obj = varargin{1}; % get the arlas object
    V = '14OCT2024'; % this is the current version number of this program


    %------ USER MODIFIABLE PARAMETERS ------------------------------------

    % specify post-hoc analysis, if using ---------------------------------
        % basePath = 'D:\myWork\ARLas\Data\'; % <-- be sure to end this one with a \ character, but not the others below it!
        % SUBJ = 'FX666';        % subject ID
        % DATE = '23AUG2021';     % date data collected
        % EXPRUN = 'audio_run1'; % experiment name and run number
    %--------------------------------------------------------------------------
    
    
    %--------------------------------------------------------------------------

    disp(' ')
    disp('----- Starting behavioral threshold experiment -----')
    
    % hardware setup -----
    probeL = 'A'; % probe that is in left ear;
    probeR = []; % probe that is in right ear;
    [inputs,outputs] = hardwareSetup; % read in the hardware setup
    if strcmp(probeR,'A')
        probeInputR = inputs{1};         % ER10xA microphone
        probeOutputR = outputs{1};       % ER10xA loudspeakers
    elseif strcmp(probeR,'B')
        probeInputR = inputs{2};         % ER10xB microphone
        probeOutputR = outputs{2};       % ER10xB loudspeakers
    elseif isempty(probeR)
        probeInputR = [];
        probeOutputR = [];
    else
        error('Unrecognized probe value. Must be A or B.')
    end
    if strcmp(probeL,'A')
        probeInputL = inputs{1};         % ER10xA microphone
        probeOutputL = outputs{1};       % ER10xA loudspeakers
    elseif strcmp(probeL,'B')
        probeInputL = inputs{2};         % ER10xB microphone
        probeOutputL = outputs{2};       % ER10xB loudspeakers
    elseif isempty(probeL)
        probeInputL = [];
        probeOutputL = [];
    else
        error('Unrecognized probe value. Must be A or B.')
    end
    clickerInput = inputs{6};

    
     if isempty(obj)
        runPostAnalysis = 1;
        disp('----- Starting audiometer post-hoc analysis -----')
    else
        runPostAnalysis = 0;
    end
    
    
    if runPostAnalysis == 1 % running post-hoc analyses on previously recorded data
        getSavedData(basePath,SUBJ,DATE,EXPRUN);
    else % actually collecting data
        testFreqs =  [1 2 4 8 10 12.5 14 16]' * 1000;
        aHandle = arlas_audiometer(obj,testFreqs);
        aHandle.audiometerProbes(probeInputL,probeOutputL,probeInputR,probeOutputR,clickerInput)
        aHandle.simulatePatient = 0;
    end
    
keyboard    
% can use this to turn back on all buttons to active if they get stuck:
%aHandle.buttonManager(20)
 


% % Save figures ----------
%     keyboard
%     disp('----- Saving data -----')
%     try
%         if ~isempty(probeInputL)
%             savePath = obj.objPlayrec.savedFilesPath;
%             figureFileName = [obj.subjectID,'_dpoaePrimariesL.fig'];
%             figureFileName = ARLas_saveName(savePath,figureFileName);
%             savefig(h1L,[savePath,figureFileName])
%             figureFileName = [obj.subjectID,'_dpoaePrimariesL.bmp'];
%             figureFileName = ARLas_saveName(savePath,figureFileName);    
%             saveas(h1L,[savePath,figureFileName])
% 
%             figureFileName = [obj.subjectID,'_dpoaeL.fig'];
%             figureFileName = ARLas_saveName(savePath,figureFileName);
%             savefig(h2L,[savePath,figureFileName])
%             figureFileName = [obj.subjectID,'_dpoaeL.bmp'];
%             figureFileName = ARLas_saveName(savePath,figureFileName);    
%             saveas(h2L,[savePath,figureFileName])
%         end
%     catch
%     end
% 
%     try
%         if ~isempty(probeInputR)
%             savePath = obj.objPlayrec.savedFilesPath;
%             figureFileName = [obj.subjectID,'_dpoaePrimariesR.fig'];
%             figureFileName = ARLas_saveName(savePath,figureFileName);
%             savefig(h1R,[savePath,figureFileName])
%             figureFileName = [obj.subjectID,'_dpoaePrimariesR.bmp'];
%             figureFileName = ARLas_saveName(savePath,figureFileName);    
%             saveas(h1R,[savePath,figureFileName])
% 
%             figureFileName = [obj.subjectID,'_dpoaeR.fig'];
%             figureFileName = ARLas_saveName(savePath,figureFileName);
%             savefig(h2R,[savePath,figureFileName])
%             figureFileName = [obj.subjectID,'_dpoaeR.bmp'];
%             figureFileName = ARLas_saveName(savePath,figureFileName);    
%             saveas(h2R,[savePath,figureFileName])
%         end
%     catch
%     end
%     
%     % Save analyses ------------------------
%     try 
%         if ~isempty(probeInputL)
%             savePath = obj.objPlayrec.savedFilesPath;
%             saveFileName = [obj.subjectID,'_analyzedDPOAE_L.mat'];
%             saveFileName = ARLas_saveName(savePath,saveFileName);
%             save([savePath,saveFileName],'DPOAE_L')
%         end
%     catch
%         disp('Warning: DPOAE_L Analysis not saved!')
%     end
%     
%     try 
%         if ~isempty(probeInputR)
%             savePath = obj.objPlayrec.savedFilesPath;
%             saveFileName = [obj.subjectID,'_analyzedDPOAE_R.mat'];
%             saveFileName = ARLas_saveName(savePath,saveFileName);
%             save([savePath,saveFileName],'DPOAE_R')
%         end
%     catch
%         disp('Warning: DPOAE_R Analysis not saved!')
%     end
% 
% 
%     disp('----- Finished with DPOAE experiment -----')
%     disp(' ')
%     disp(' ')
    
end

% INTERNAL FUNCTIONS ------------------------------------------------------
function [] = getSavedData(basePath,SUBJ,DATE,EXPRUN)      
    os = '\';
    one = [basePath,SUBJ,os];
    two = [SUBJ,'_',DATE,os];
    three = [SUBJ,'_',DATE,'_',EXPRUN,os];
    four = [SUBJ,'_',DATE,'_',EXPRUN,'_analysis',os];
    postPath = [one,two,three];
    postPathSave = [one,two,three,four];

    % Make sure the analysis directory exists -----
    if exist(postPathSave,'dir') ~= 7 % check to make sure that the backup directory exists
        success = mkdir(postPathSave); % if not, try to create it
        if success ~= 1
            error('Backup directory does not exist and could not be created. Aborting program.')
        end
        addpath(genpath(postPathSave))
    end
    
    d = dir([postPath,'Ch',num2str(electrodeInput.ch),'_',electrodeInput.label,'*.mat']); % elctrical files directory
    dA = dir([postPath,'Ch',num2str(probeInput.ch),'_',probeInput.label,'*.mat']); % acoustical files directory
    nLevels = size(dE,1); % will be used to set iterations of the loop
    fs = 96000;

    obj.subjectID = SUBJ;
    
    dummy = load([postPath,dA(1).name]);
    header = dummy.header;
    if contains('B',header.label)
        ear = 'R';
    else
        ear = 'L';
    end
    
end
