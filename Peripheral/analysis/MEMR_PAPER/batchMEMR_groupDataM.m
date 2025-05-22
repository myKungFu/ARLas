function [] = batchMEMR_groupDataM()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% batchMEMR_groupData
% 
% Batch analysis of group MEMR data 
%
% Author: Shawn Goodman & Ehsan Khalili
% Date: January 23, 2024
% Last Updated: January 23, 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    parentDrive = 'D'; % use 'C' or 'D'
    dataPathName = [parentDrive,':\MEMR_MepaniAnalysis\']; % location of raw data
    savePath = [parentDrive,':\MEMR_GroupDataM\']; % where to save analyzed data

    d = dir([dataPathName,'*.mat']); % read in file directory information
    nFiles = size(d,1); % number of subject folders
    counter = 1;
    for ii=1:nFiles
        disp(['Analyzing file ',num2str(ii),' of ',num2str(nFiles)])
        fileName = d(ii).name;
        %runNumber = str2num(fileName(10));
        dummy = load([dataPathName,fileName]); 
        MEMR_mem = dummy.MEMR_mem;
        

       % plot (MEMR_mem.d1)
        % plot (MEMR_mem.freq, MEMR_mem.d2);
        % pause
% try 
        % extract data as a row vector
        X(ii,1) = 0; %max(MEMR_mem.trend);
        X(ii,2) = 0; % MEMR_mem.peakTime;
        X(ii,3) = MEMR_mem.peakAmp;
        X(ii,4) = 0;
        try 
            X(ii,5) = MEMR_mem.thdOnsetTime;
        catch
            X(ii,5) = NaN;
        end
        try
            X(ii,6) = MEMR_mem.thdOffsetTime;
        catch
            X(ii,6) = NaN;
        end
        X(ii,7) = MEMR_mem.thdOnsetLvl;
        X(ii,8) = MEMR_mem.thdOffsetLvl;
        X(ii,9) = MEMR_mem.hysteresis;
        X(ii,10) = MEMR_mem.slopeUp;
        X(ii,11) = MEMR_mem.slopeDn;
        try
        X(ii,12) = MEMR_mem.thd;
        catch 
        X(ii,12) = NaN;
        end
        try     
        X(ii,13) = MEMR_mem.thdAmp;
        catch
        X(ii,13) = NaN;
        end
        X(ii,14) = NaN;%median(MEMR_mem.pSPL.'); % peak spl of clicks
        X(ii,15) = NaN; %max(MEMR_mem.RMS).'; % rms elicitor levels (max)
% catch ME
    % keyboard
% end
        % X(1,14) = MEMR_mem.elicitorLevel;
        % X(1,15) = MEMR_mem.elicitor;
        % X(1,16) = MEMR_mem.RMS;
        % X(1,17) = MEMR_mem.RMSC;
        % X(1,18) = MEMR_mem.RMST;
        % X(1,19) = MEMR_mem.pSPL;
        % X(1,20) = MEMR_mem.rmsSPL;
        % 
        % if runNumber == 1
        %     X1(counter,:) = X;
        % elseif runNumber == 2
        %     X2(counter,:) = X;
        % elseif runNumber == 3
        %     X3(counter,:) = X;
        % elseif runNumber == 4
        %     X4(counter,:) = X;
        % end

        % if mod(ii,4)==0
        %     counter = counter + 1;
        % end
    end
 
    trend = X(:,1);
    pTime = X(:,2);
    pAmp = X(:,3);
    % delay = X(:,1);
    thdOnsetTime = X(:,4);
    thdOffsetTime = X(:,6);
    thdOnsetLvl = X(:,7);
    thdOffsetLvl = X(:,8);
    hysteresis = X(:,9);
    slopeUp = X(:,10);
    slopeDn = X(:,11);
    thd = X(:,12);
    thdAmp = X(:,13);

    % pSPL = X(:,14);
    % rms = X(:,15);

    saveName = [savePath,'MEMR_groupData.mat'];
    save(saveName,'X','trend','pTime','pAmp',...
        'thdOnsetTime','thdOffsetTime','thdOnsetLvl','thdOffsetLvl','hysteresis','slopeUp','slopeDn','thd','thdAmp')

end