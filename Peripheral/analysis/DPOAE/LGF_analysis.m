function [LGF,h1,h2,h3] = LGF_analysis(varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [LGF,h1,h2,h3] = LGF_analysis(varargin);
%
% Post-Hoc analysis of data collected using LGF.m
% If no input arguments, the program will use the information at the bottom of 
% USER MODIFIABLE PARAMETERS in order to do the analysis.
%
% LGF = a structure containing the data and analysis
% h1, h2, h3 = handles to figures
% All output arguments are saved by this function.
%
% Author: Shawn Goodman, PhD
% Auditory Research Lab, the University of Iowa
% Date: Original Date: July 19, 2024
% Last Updated: April 21, 2025 -- ssg
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%--------------------------------------------------------------------------
%------ USER MODIFIABLE PARAMETERS ----------------------------------------

    % specify data location for post-hoc analysis, if using:

        basePath = 'C:\myWork2\ARLas\Data\'; % <-- be sure to end this one with a \ character, but not the others below it!
        
        SUBJ = 'dummy';        % subject ID
        DATE = '21APR2025';    % date data collected
        EXPRUN = 'dummy_run2'; % experiment name and run number
              
%------ END USER MODIFIABLE PARAMETERS ------------------------------------
%--------------------------------------------------------------------------

    disp('----- Starting LGF post-hoc analysis -----')
    
    if isempty(varargin) % if getting previously saved data
        disp(' ')
        [d,~,~,~,postPath,postPathSave,~,~,~,~] = getSavedData(basePath,SUBJ,DATE,EXPRUN);
        dummy = load([postPath,d(1).name]);
        header = dummy.header;
        Data = dummy.data;
        clear dummy
    else % analyzing data just collected
        header = varargin{1};
        Data = varargin{2};
        obj = varargin{3};
        os = filesep;
        postPathSave = [header.pathName,header.expID,'_analysis',os];
        % Make sure the analysis directory exists -----
        if exist(postPathSave,'dir') ~= 7 % check to make sure that the backup directory exists
            success = mkdir(postPathSave); % if not, try to create it
            if success ~= 1
                %keyboard
                error('Backup directory does not exist and could not be created. Aborting program.')
            end
            addpath(genpath(postPathSave))
        end
    end

    % set up analysis structure
    LGF.subjID = header.subjID;
    try LGF.ear = header.userInfo.testEar;
    catch
        LGF.ear = 'left';
        disp('Test ear not a part of header.userInfo. Using default of left.')
    end
    LGF.postHocTimeStamp = cellstr(datetime('now'));
    try LGF.timeStamp = header.timeStamp;
    catch
        LGF.timeStamp = LGF.postHocTimeStamp;
    end
    LGF.fs = header.fs;
    LGF.f2min = header.userInfo.fmin;
    LGF.f2max = header.userInfo.fmax;
    LGF.nSweeps = header.userInfo.nSweeps;
    LGF.calType = header.userInfo.calType;
    LGF.targetCalType = header.userInfo.targetCalType;

    % figure out which micCorrection to use
    try
        testProbe = obj.objPlayrec.userInfo.testProbe;
        if strcmp(testProbe,'Left') | strcmp(testProbe,'left')
            LGF.micCorrection = header.userInfo.micCorrectionL;
        elseif strcmp(testProbe,'Right') | strcmp(testProbe,'right')
            LGF.micCorrection = header.userInfo.micCorrectionR;
        else
            warning('Unable to choose test probe.')
        end
    catch % for backwards compatibility before vRT3
        LGF.micCorrection = header.userInfo.micCorrectionL;
        LGF.probe = 'Left';
        if isempty(LGF.micCorrection)
            LGF.micCorrection = header.userInfo.micCorrectionR;
            LGF.probe = 'Right';
        end
    end

    f1 = header.userInfo.f1; % f1 (Hz)
    f2 = header.userInfo.f2; % f2 (Hz)
    fdp = header.userInfo.fdp; % fdp (Hz)
    L1 = header.userInfo.L1; % measured L1 (dB SPL) -- large, 8-second vector
    L2 = header.userInfo.L2; % measured L2 (dB SPL) -- large 8-second vector
    fs = header.userInfo.fs; % sampling rate (samples / sec)
    nSweeps = header.userInfo.nSweeps; % number of sweeps
    LGF.f1 = f1;
    LGF.f2 = f2;
    LGF.fdp = fdp;
    LGF.L1 = L1;
    LGF.L2 = L2;
    LGF.fs = fs;
    LGF.nSweeps = nSweeps;
    
    % bandpass filter from 0.4-20 kHz
    Rows = length(Data(:))/nSweeps;
    Data = reshape(Data,Rows,nSweeps);
    b = getbpf; % fir filter coefficients                 
    [Rows,Cols] = size(Data);
    Data = Data(:);
    Data = fastFilter(b,Data);
    Data = reshape(Data,Rows,Cols);

    % apply mic correction
    Data = applyMicCorrection(Data,LGF.micCorrection);
    % before applying EPL, do a check to make sure no problems:
    isc = header.userInfo.iscS1L; % in-situ calibration info
    LGF.isc = isc;
    ok = 1;
    if isc.diam > 3 | isc.diam < .1
        ok = 0;
    end
    if isc.length > 4 | isc.length < .25
        ok = 0;
    end
    if ok == 1
        % apply EPL correction
        LGF.eplCorrected = 1;
        [pl,Pl,phi,other,wf] = ARLas_convertPL(Data,isc);
    else
        LGF.eplCorrected = 0;
    end
    

    % -------------------------------------------
    try
        Lmin = header.userInfo.L2min;
        Lmax = header.userInfo.L2max;
    catch
        Lmin = 0; % minimum L2 is assumed to be zero dB FPL (20 uPa)
        Lmax = max(L2dB); % maximum L2 in dB FPL
    end

% Analysis is done here ---------------------------------------------------    
    
    x = 20*log10(L2/.00002); % L2 levels in dB FPL
    stepSize = 1; % setp size is 1 dB
    [signal,nf,snr,targets,zbar] = sweptLSF(x,Data,fdp,Lmin,Lmax,stepSize,fs);
    LGF.signal = signal;
    LGF.targets = targets;
    LGF.nf = nf;
    LGF.Enf = mean(nf); % expected noise floor (mean value)

    % doing the analysis is dependent on having good enough SNR
    maxSNR = sort(snr);
    maxSNR = median(maxSNR(end-4:end)); % take mean of the four largest SNR values
    snrCriterion = 10;
    LGF.snrCriterion = snrCriterion;
    LGF.maxSNR = maxSNR;
    if maxSNR > snrCriterion % if you have enough SNR to fit
        [~,indxMax] = max(signal); % location of maximum in the signal
        L2max = targets(indxMax); % L2 at the maximum DP amplitude
        L2cut = L2max - 3; % go 3 dB down from maximum for fitting
        [~,indxCut] = min(abs(L2cut - targets)); % cut down so fitting only done below peak
        [jar,h1,h2,h3] = ricianLGF_fit2(signal(1:indxCut),nf(1:indxCut),targets(1:indxCut));

        LGF.L2max = L2max;
        LGF.L2cut = L2cut;
        LGF.indxCut = indxCut;
        LGF.jar = jar;
    else % don't fit if you don't have enough SNR
        jar = struct;
        LGF.jar = jar;
        h2 = []; % can't make these plots if no fitting done
        h3 = []; 
        h1 = figure;
        plot(targets,signal,'b')
        hold on
        plot(targets,nf,'k')
        xlim([targets(1),targets(end)])
        xlabel('Target L2 (dB FPL)')
        ylabel('Ldp (dB SPL)')
        grid on
    end


% Save Analysis -----------------------------------------------------------
    % Save analysis -----------------------------------
    disp('----- Saving analysis files -----')
    savePath = postPathSave;
    saveFileName = [LGF.subjID,'_analyzedLGF_',LGF.ear,'.mat'];
    save([savePath,saveFileName],'LGF')

    % Save figures ------------------------------------
    disp('----- Saving figures -----')
    savePath = postPathSave;
    figureFileName = [LGF.subjID,'_LGF_',LGF.ear,'.fig'];
    savefig(h1,[savePath,figureFileName])
    figureFileName = [LGF.subjID,'_LGF_',LGF.ear,'.bmp'];
    saveas(h1,[savePath,figureFileName])

    if maxSNR > snrCriterion
        savePath = postPathSave;
        figureFileName = [LGF.subjID,'_LGF_LL_',LGF.ear,'.fig'];
        savefig(h2,[savePath,figureFileName])
        figureFileName = [LGF.subjID,'_LGF_LL_',LGF.ear,'.bmp'];
        saveas(h2,[savePath,figureFileName])
    
        savePath = postPathSave;
        figureFileName = [LGF.subjID,'_LGF_LL_',LGF.ear,'.fig'];
        savefig(h3,[savePath,figureFileName])
        figureFileName = [LGF.subjID,'_LGF_LL_',LGF.ear,'.bmp'];
        saveas(h3,[savePath,figureFileName])
    end

    disp('----- Finished with LGF Swept L2 post-hoc analysis -----')
    disp(' ')
    disp(' ')

  
end
% INTERNAL FUNCTIONS ------------------------------------------------------
function [d,iscS1L,iscS2L,iscS12L,postPath,postPathSave,testEar,fmin,fmax,fs] = getSavedData(basePath,SUBJ,DATE,EXPRUN)      
    os = filesep; %'\';
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
    % Get the saved data
    
    % dL = dir([postPath,'Ch','*.mat']); % get all the files in the folder
    % totalFilesL = size(dL,1);
    
    d = dir([postPath,'*.mat']);
    nFiles = size(d,1);
    counter = 1;
    for ii=1:nFiles
        name = d(ii).name;
        if ~contains(name,'inSituCal')
            d2(counter) = d(ii);
            counter = counter + 1;
        end
    end
    d = d2;
    clear d2
    % now order these by frequency and by l2 (not alphabetically)
    nFiles = size(d,2);
    % for ii=1:nFiles % load each file
    %     dummy = load([postPath,d(ii).name]);
    %     Lvl(ii,1) = dummy.header.userInfo.targetL2; % l2
    %     Frq(ii,1) = dummy.header.userInfo.F2(1,1);  % f2
    % end
    % Freqs = unique(Frq);
    % nFreqs = length(Freqs);
    % Levels = unique(Lvl);
    % nLevels = length(Levels);
    % for ii=1:nLevels
    %     for jj=1:nFreqs
    %         L = Levels(ii);
    %         F = Freqs(jj);
    %         indxii = Lvl==L;
    %         indxjj = Frq==F;
    %         INDX(ii,jj) = find(indxii.*indxjj);
    %         % INDX is a matrix of indices showing which index of d to read
    %         % in such a way that it corresponds to the correct level and
    %         % frequency.
    %     end
    % end
    % % just make a choice here to use L (this is for jeff)
    dL = d;
    dR = [];

    
    % try
    %     %dL = dir([postPath,'Ch',num2str(probeInputL.ch),'_',probeInputL.label,'_dpoae_*.mat']);
    %     %dL = dir([postPath,'Ch',num2str(probeInputL.ch),'_',probeInputL.label,'*.mat']);
    %     dL = dir([postPath,'*.mat']);
    %     nFiles = size(dL,1);
    %     counter = 1;
    %     for ii=1:nFiles
    %         name = dL(ii).name;
    %         if ~contains(name,'inSituCal')
    %             dL2(counter) = dL(ii);
    %             counter = counter + 1;
    %         end
    %     end
    %     dL = dL2;
    %     clear dL2
    % catch
    %     dL = [];
    % end
    % try
    %     %dR = dir([postPath,'Ch',num2str(probeInputR.ch),'_',probeInputR.label,'_dpoae_*.mat']);
    %     %dR = dir([postPath,'Ch',num2str(probeInputR.ch),'_',probeInputR.label,'*.mat']);
    %     dR = dir([postPath,'*.mat']);
    %     nFiles = size(dR,1);
    %     counter = 1;
    %     for ii=1:nFiles
    %         name = dR(ii).name;
    %         if ~contains(name,'inSituCal')
    %             dR2(counter) = dR(ii);
    %             counter = counter + 1;
    %         end
    %     end
    %     dR = dR2;
    %     clear dR2
    % catch
    %     dR = [];
    % end
    % nLevels = max([size(dL,1),size(dR,1)]); %size(dL,1); % will be used to set iterations of the loop
    fs = 96000; % sampling rate (Hz)

    % get the saved in-situ calibrations
    if ~isempty(dL)
        dummy = load([postPath,dL(1).name]);
        q = dummy.header;
        %C1L = q.userInfo.C2L;
        %C2L = q.userInfo.C2L;
        iscS1L = q.userInfo.iscS1L;
        iscS2L = q.userInfo.iscS2L;        
        iscS12L = [];
    else
        iscS1L = [];
        iscS2L = [];        
        iscS12L = [];
    end

    if ~isempty(dR)
        dummy = load([postPath,dR(1).name]);
        q = dummy.header;
        %C1R = q.userInfo.C2R;
        %C2R = q.userInfo.C2R;
        iscS1R = q.userInfo.iscS1R;
        iscS2R = q.userInfo.iscS2R;        
        iscS12R = [];
    else
        iscS1R = [];
        iscS2R = [];        
        iscS12R = [];
    end

    obj.subjectID = SUBJ;
    obj.timeStamp = q.timeStamp;
    obj.fs = q.userInfo.fs;
    fmin = q.userInfo.fmin;
    fmax = q.userInfo.fmax;

    if ~isempty(dR) & ~isempty(dL)
        testEar = 'Both';
    elseif ~isempty(dR)
        testEar = 'Right';
    elseif ~isempty(dL)
        testEar = 'Left';
    end

end
function b = getbpf()
    Fs = 96000;
    Fstop1 = 100;             % First Stopband Frequency 
    Fpass1 = 400;            % First Passband Frequency -->for a low freq of 750, the dp is 450 Hz
    Fpass2 = 20000;           % Second Passband Frequency
    Fstop2 = 24000;           % Second Stopband Frequency
    Dstop1 = 0.001;           % First Stopband Attenuation
    Dpass  = 0.057501127785;  % Passband Ripple
    Dstop2 = 0.0001;          % Second Stopband Attenuation
    flag   = 'scale';         % Sampling Flag
    % Calculate the order from the parameters using KAISERORD.
    [N,Wn,BETA,TYPE] = kaiserord([Fstop1 Fpass1 Fpass2 Fstop2]/(Fs/2), [0 ...
                                 1 0], [Dstop1 Dpass Dstop2]);
    if mod(N,2)~= 0
        N = N + 1;
    end
    b  = fir1(N, Wn, TYPE, kaiser(N+1, BETA), flag)';
end

% OLD CODE ----------------------------------------------------------------
    % % Write analysis to spreadsheet
    % disp('----- Writing Excel files -----')
    % savePath = postPathSave;
    % saveFileName = [LGF.subjID,'_dpSweptLGF_',LGF.ear,'.xls'];
    % %saveFileName = ARLas_saveName(savePath,saveFileName);
    % saveFileName = [saveFileName,'x'];
    % DATA = [LGF.targetL2,LGF.Ldp,LGF.Ndp];
    % writematrix(DATA,[savePath,saveFileName],'Sheet',1,'Range','A2')
    % writematrix('L2',[savePath,saveFileName],'Sheet',1,'Range','A1');
    % writematrix('Ldp',[savePath,saveFileName],'Sheet',1,'Range','B1');
    % writematrix('Ndp',[savePath,saveFileName],'Sheet',1,'Range','C1');
    % writematrix('bHat',[savePath,saveFileName],'Sheet',1,'Range','D1');
    % writematrix('mHat',[savePath,saveFileName],'Sheet',1,'Range','E1');
    % if ~isnan(egg.fit.bHat)
    %     writematrix(LGF.egg.fit.bHat,[savePath,saveFileName],'Sheet',1,'Range','D2');
    %     writematrix(LGF.egg.fit.bHat,[savePath,saveFileName],'Sheet',1,'Range','E2');
    %     writematrix('f2',[savePath,saveFileName],'Sheet',1,'Range','F1');
    %     writematrix(LGF.f2(1),[savePath,saveFileName],'Sheet',1,'Range','F2');
    %     if nFreqs == 2
    %         warning off
    %         DATA = [DPOAE_2.targetL2,DPOAE_2.Ldp,DPOAE_2.Ndp];
    %         writematrix(DATA,[savePath,saveFileName],'Sheet',2,'Range','A2')
    %         writematrix('L2',[savePath,saveFileName],'Sheet',2,'Range','A1');
    %         writematrix('Ldp',[savePath,saveFileName],'Sheet',2,'Range','B1');
    %         writematrix('Ndp',[savePath,saveFileName],'Sheet',2,'Range','C1');
    %         writematrix('bHat',[savePath,saveFileName],'Sheet',2,'Range','D1');
    %         writematrix('mHat',[savePath,saveFileName],'Sheet',2,'Range','E1');
    %         writematrix(DPOAE_2.egg.fit.bHat,[savePath,saveFileName],'Sheet',2,'Range','D2');
    %         writematrix(DPOAE_2.egg.fit.bHat,[savePath,saveFileName],'Sheet',2,'Range','E2');
    %         writematrix('f2',[savePath,saveFileName],'Sheet',2,'Range','F1');
    %         writematrix(DPOAE_2.f2(1),[savePath,saveFileName],'Sheet',2,'Range','F2');
    %         warning on
    %     end
    % else
    % 
    % end
