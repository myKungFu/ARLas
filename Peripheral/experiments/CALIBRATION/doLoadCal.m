function [] = doLoadCal(varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% doLoadCal(varargin)
%
% Perform load (typically ear canal) calibration using Thevenin source methods.
% This code plays a broadband chirp and calculates load impedance.
% Checks for low-frequency leaks.
% Plots the ear canal transfer function.
% Plots wideband reflectance.
% Asks user if okay to continue.
%
% If this function is called by another experiment, returns the results.
% Otherwise, if called directly ("stand-alone"), saves the results.

% Author: Shawn Goodman, PhD
% Auditory Research Lab, the University of Iowa
% Original Date: April 23, 2025
% Last Updated: April 23, 2025 -- ssg
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    obj = varargin{1}; % get the arlas object
    V = '23APR2025'; % this is the current version number of this program

%------ USER MODIFIABLE PARAMETERS ----------------------------------------
%--------------------------------------------------------------------------

    inSituReps = 6;


    % specify which ER10X probes to use    
    probeL = 'A'; % probe that is in LEFT ear; 'A', 'B', or []
    probeR = []; % probe that is in RIGHT ear; 'A', 'B', or []
                  
    fmin = 100;
    fmax = 18000;

    % additional parameters that user usually doesn't need to mess with
    doMicCorrection = 0; % turn on and off microphone correction 
    doPreCalCheck = 1; % turn on and off the pre-test cal check plotting

%------ END USER MODIFIABLE PARAMETERS ------------------------------------
%--------------------------------------------------------------------------


    calType = 'thev'; % use 'thev' for Thevenin source
    targetCalType = 'fpl'; % 'fpl', 'spl', 'ipl'.
    ISCL = struct; % initialize structures for left ear
    ISCR = struct; % initialize structures for right ear
    
    if ~(isempty(probeL) | isempty(probeR)) % if both ears are being tested
        testProbe = 'Both';
    elseif ~(isempty(probeL))
        testProbe = 'Left';
    elseif ~(isempty(probeR))
        testProbe = 'Right';
    end
    if ~(isempty(probeL) | isempty(probeR)) % if both ears are being tested
         testEar = 'Both';
    else
        testEar = getStimulusEar;
    end
    % get probe settings
    [probeInputL,probeOutputL,probeInputR,probeOutputR,~] = getProbeSettings(probeL,probeR,[]);
    
    % get most recent calibration files -----------------------------------
    % microphone corrections
    micCorrectionL = getMicCorrection(probeInputL,doMicCorrection,obj.map);
    micCorrectionR = getMicCorrection(probeInputR,doMicCorrection,obj.map);
    if doMicCorrection==1
        if ~isempty(probeInputL)
            if micCorrectionL == 1
                warning('MIC CORRECTION IS NOT BEING APPLIED!')
            end
        end
        if  ~isempty(probeInputR)
            if micCorrectionR == 1
                warning('MIC CORRECTION IS NOT BEING APPLIED!')
            end
        end
    end
    % Thevenin source calibrations
    [C1L,C2L,calPath1L,calPath2L] = getOutputCal(calType,probeOutputL,obj.map);
    [C1R,C2R,calPath1R,calPath2R] = getOutputCal(calType,probeOutputR,obj.map);
    
    % PERFORM IN-SITU (LOAD) CALIBRATION ----------------------------------
    disp('----- Running in-situ calibration -----')
    doIndividual = 1;
    doSimultaneous = 0;
    [iscS1L,iscS2L,iscS12L] = ARLas_runISC(obj,probeInputL,probeOutputL,calPath1L,calPath2L,calType,fmin,fmax,micCorrectionL,inSituReps,doIndividual,doSimultaneous);
    [iscS1R,iscS2R,iscS12R] = ARLas_runISC(obj,probeInputR,probeOutputR,calPath1R,calPath2R,calType,fmin,fmax,micCorrectionR,inSituReps,doIndividual,doSimultaneous);
    ISCL.(['S1']) = iscS1L; % for left probe, loudspeaker 1
    ISCL.(['S2']) = iscS2L; % for left probe, loudspeaker 2
    ISCR.(['S1']) = iscS1R; % for right probe, loudspeaker 1
    ISCR.(['S2']) = iscS2R; % for right probe, loudspeaker 2

    % RUN LOAD CALIBRATION CHECK ------------------------------------------
    if doPreCalCheck == 1
        try
            if ~isempty(probeOutputL)
                labelL = probeOutputL.label;
                dL = dir([obj.savingGrace,'*',labelL,'_inSituCal_*']);
            else
                dL = [];
            end
            if ~isempty(probeOutputR)
                labelR = probeOutputR.label;
                dR = dir([obj.savingGrace,'*',labelR,'_inSituCal_*']);
            else
                dR = [];
            end
            iscCheckHandleL = [];
            if ~isempty(probeOutputL)
                dummy = load([obj.savingGrace,dL(1).name]); % first calibration chirp (loudspeaker 1)
                processData(dummy.data,dummy.header,ISCL.S1);
                %iscCheckHandleL = inSituCheck(probeOutputL.label,probeOutputL.ch(1),dummy.header,dummy.data,obj.fs,1,iscCheckHandleL);
                dummy = load([obj.savingGrace,dL(2).name]); % first calibration chirp (loudspeaker 2)
                processData(dummy.data,dummy.header,ISCL.S2);
                %iscCheckHandleL = inSituCheck(probeOutputL.label,probeOutputL.ch(2),dummy.header,dummy.data,obj.fs,2,iscCheckHandleL);
            end
            iscCheckHandleR = [];
            if ~isempty(probeOutputR)
                dummy = load([obj.savingGrace,dR(1).name]); % first calibration chirp (loudspeaker 1)
                iscCheckHandleR = inSituCheck(probeOutputR.label,probeOutputR.ch(1),dummy.header,dummy.data,obj.fs,1,iscCheckHandleR);
                dummy = load([obj.savingGrace,dR(2).name]); % first calibration chirp (loudspeaker 2)
                iscCheckHandleR = inSituCheck(probeOutputR.label,probeOutputR.ch(2),dummy.header,dummy.data,obj.fs,2,iscCheckHandleR);
            end
        catch ME
            disp('Start of measurement in-situ calibration check failed!')
        end
        answer = questdlg('Peaks at low and high frequencies should be between 1 and 2 Pa. Do NOT close figures!',...
            'In-situ Calibration Check',...
	        'Continue','Abort','Continue');
        switch answer
            case 'Continue'
            case 'Abort'
                return
        end
    end
    % END PRE-TEST IN-SITU CALIBRATION CHECK ------------------------------
        

end

% INTERNAL FUNCTIONS ------------------------------------------------------
function [] = processData(data,header,isc)

    if mod(size(data,1),2)~=0
        data = data(1:end-1,:);
    end
    [pl,Pl,phi,other,wf] = ARLas_convertPL(data,isc);
    other.canalLength
    other.diam

    figure
    plt(Pl.f,Pl.PRR)

keyboard
    [frequency,signal,noiseFloor] = ARLas_fda(data,fs,0.00002);
    plot(frequency/1000,signal,'b')
    hold on
    plot(frequency/1000,noiseFloor,'k')
    hold on
    xlim([.1 20])
    %ylim([-2.5 2.5])
    xlabel('Frequency (kHz)')
    ylabel('Magnitude (dB SPL)')
    title([probeLabel,'   Channel ',num2str(probeChannel)])

    
    % Plot all the waveforms stretched out into a single vector -----------
    data = data(:);
    N = length(data);
    t = (0:1:N-1)'/fs;
    
    if whichProbe == 1
        subplot(3,2,3)
    elseif whichProbe == 2
        subplot(3,2,4)
    end
    plot(t,data)
    hold on
    line([0 t(end)],[1 1],'LineStyle','--','LineWidth',0.5,'Color',[0 0 0])
    line([0 t(end)],[-1 -1],'LineStyle','--','LineWidth',0.5,'Color',[0 0 0])
    line([0 t(end)],[2 2],'LineStyle',':','LineWidth',0.5,'Color',[0 0 0])
    line([0 t(end)],[-2 -2],'LineStyle',':','LineWidth',0.5,'Color',[0 0 0])
    xlim([0 t(end)])
    ylim([-2.5 2.5])
    xlabel('Time (s)')
    ylabel('Amplitude (Pa)')
    title([probeLabel,'   Channel ',num2str(probeChannel)])
end
function [micCorrection] = getMicCorrection(probeInput,doMicCorrection,map)
    if isempty(probeInput)
        micCorrection = [];
        return
    end
    [pathName_mic,folderName_mic,fileName_mic] = mostRecentCalibration('mic',probeInput.label,[],map);
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
function [testEar] = getStimulusEar()
            testEar = [];
            % Get stimulus ear (regardless of which probe channel is being used) from user:
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

        
% OLD CODE ----------------------------------------------------------------
% % Save in situ cal figures ------------------------------------------------
% %   Note: for this program, the analyzed data and figures are saved by the
% %   analysis program, rather than here.
% 
%     disp('----- Saving figures -----')
%     try
%         if ~isempty(probeInputL)
%             savePath = [obj.savingGrace,obj.experimentID,'_analysis',obj.sep];
%             if exist(savePath,'dir') == 0 % if path does not exist
%                 success = mkdir(savePath);
%                 if ~success
%                     warning('Unable to create new Experiment Directory: data save path')
%                 end
%             end 
%             figureFileName = [obj.subjectID,'_inSituCal_',DPOAE_L.ear,'.fig'];
%             figureFileName = ARLas_saveName(savePath,figureFileName);
%             savefig(iscCheckHandleL,[savePath,figureFileName])
%             figureFileName = [obj.subjectID,'_inSituCal_',DPOAE_L.ear,'.bmp'];
%             figureFileName = ARLas_saveName(savePath,figureFileName);    
%             saveas(iscCheckHandleL,[savePath,figureFileName])
%         end
%     catch
%         disp('Warning: One or more figures for DPOAE_L not saved!')        
%     end
% 
%     try
%         if ~isempty(probeInputR)
%             savePath = [obj.savingGrace,obj.experimentID,'_analysis',obj.sep];
%             if exist(savePath,'dir') == 0 % if path does not exist
%                 success = mkdir(savePath);
%                 if ~success
%                     warning('Unable to create new Experiment Directory: data save path')
%                 end
%             end 
%             figureFileName = [obj.subjectID,'_inSituCal_',DPOAE_R.ear,'.fig'];
%             figureFileName = ARLas_saveName(savePath,figureFileName);
%             savefig(iscCheckHandleR,[savePath,figureFileName])
%             figureFileName = [obj.subjectID,'_inSituCal_',DPOAE_R.ear,'.bmp'];
%             figureFileName = ARLas_saveName(savePath,figureFileName);    
%             saveas(iscCheckHandleR,[savePath,figureFileName])
%         end
%     catch
%        disp('Warning: One or more figures for DPOAE_R not saved!')
%     end
% 
% % NOTE: Excel results not saved! Use the .mat files to get saved analyses
%     toc/60
%     disp('----- Finished with Swept DPOAE Fixed L1 experiment -----')
%     disp(' ')
%     disp(' ')
