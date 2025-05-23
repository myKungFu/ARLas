function [] = batchMEMR_groupData515()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% batchMEMR_groupData
% 
% Batch analysis of group MEMR data 
%
% Author: Shawn Goodman & Ehsan Khalili
% Date: January 23, 2024
% Last Updated: January 23, 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    parentDrive = 'C'; % use 'C' or 'D'
    dataPathName = [parentDrive,':\myWork\ARLas\MEMR_Analysis_Paper\']; % location of raw data
    savePath = [parentDrive,':\myWork\ARLas\MEMR_Data_Paper\']; % where to save analyzed data

    d = dir([dataPathName,'*.mat']); % read in file directory information
    nFiles = size(d,1); % number of subject folders
    counter = 1;
    for ii=1:nFiles
        disp(['Analyzing file ',num2str(ii),' of ',num2str(nFiles)])
        fileName = d(ii).name;
        runNumber = str2num(fileName(10));
        dummy = load([dataPathName,fileName]); 
        MEMR_mem = dummy.MEMR_mem;

        % extract data as a row vector
        X(1,1) = max(MEMR_mem.trend);
        X(1,2) = MEMR_mem.peakTime;
        X(1,3) = MEMR_mem.peakAmp;
        X(1,4) = MEMR_mem.delay;
        X(1,5) = MEMR_mem.thdOnsetTime;
        X(1,6) = MEMR_mem.thdOffsetTime;
        X(1,7) = MEMR_mem.thdOnsetLvl;
        X(1,8) = MEMR_mem.thdOffsetLvl;
        X(1,9) = MEMR_mem.hysteresis;
        X(1,10) = MEMR_mem.slopeUp;
        X(1,11) = MEMR_mem.slopeDn;
        X(1,12) = MEMR_mem.thd;
        X(1,13) = MEMR_mem.thdAmp;

        X(1,14) = MEMR_mem.pSPL; % peak spl of clicks
        X(1,15) = max(MEMR_mem.RMS); % rms elicitor levels (max)
        % X(1,14) = MEMR_mem.elicitorLevel;
        % X(1,15) = MEMR_mem.elicitor;
        % X(1,16) = MEMR_mem.RMS;
        % X(1,17) = MEMR_mem.RMSC;
        % X(1,18) = MEMR_mem.RMST;
        % X(1,19) = MEMR_mem.pSPL;
        % X(1,20) = MEMR_mem.rmsSPL;

        if runNumber == 1
            X1(counter,:) = X;
        elseif runNumber == 2
            X2(counter,:) = X;
        elseif runNumber == 3
            X3(counter,:) = X;
        elseif runNumber == 4
            X4(counter,:) = X;
        end

        if mod(ii,4)==0
            counter = counter + 1;
        end
    end
 
    trend = [X1(:,1),X2(:,1),X3(:,1),X4(:,1)];
    pTime = [X1(:,2),X2(:,2),X3(:,2),X4(:,2)];
    pAmp = [X1(:,3),X2(:,3),X3(:,3),X4(:,3)];
    delay = [X1(:,4),X2(:,4),X3(:,4),X4(:,4)];
    thdOnsetTime = [X1(:,5),X2(:,5),X3(:,5),X4(:,5)];
    thdOffsetTime = [X1(:,6),X2(:,6),X3(:,6),X4(:,6)];
    thdOnsetLvl = [X1(:,7),X2(:,7),X3(:,7),X4(:,7)];
    thdOffsetLvl = [X1(:,8),X2(:,8),X3(:,8),X4(:,8)];
    hysteresis = [X1(:,9),X2(:,9),X3(:,9),X4(:,9)];
    slopeUp = [X1(:,10),X2(:,10),X3(:,10),X4(:,10)];
    slopeDn = [X1(:,11),X2(:,11),X3(:,11),X4(:,11)];
    thd = [X1(:,12),X2(:,12),X3(:,12),X4(:,12)];
    thdAmp = [X1(:,13),X2(:,13),X3(:,13),X4(:,13)];

    pSPL = [X1(:,14),X2(:,14),X3(:,14),X4(:,14)];
    rms = [X1(:,15),X2(:,15),X3(:,15),X4(:,15)];

    saveName = [savePath,'MEMR_groupData515.mat'];
    save(saveName,'X1','X2','X3','X4','trend','pTime','pAmp','delay',...
        'thdOnsetTime','thdOffsetTime','thdOnsetLvl','thdOffsetLvl','hysteresis','slopeUp','slopeDn','thd','thdAmp')

end