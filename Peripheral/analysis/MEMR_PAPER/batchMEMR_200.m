function [] = batchMEMR_200()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
% Batch analysis of MEMR data using analyzeMEMR_v14.m.
%analyzeMEMR_part2_v1.m
%
% Author: Shawn Goodman & Ehsan Khalili
% Date: November 13, 2023
% Last Updated: November 13, 2023
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    parentDrive = 'D'; % use 'C' or 'D'
    dataPathName = [parentDrive,':\MEMR_DATA\']; % location of raw data
    fileNameL = 'Ch3_ER10xA_memr_0001.mat';
    fileNameR = 'Ch4_ER10xB_memr_0001.mat';
    savePath = [parentDrive,':\MEMR_Analysis200\']; % where to save analyzed data

    d = dir(dataPathName); % read in file directory information
    nSubjects = size(d,1)-2; % number of subject folders
    
    for ii=1:nSubjects
        disp(['Analyzing subject ',num2str(ii),' of ',num2str(nSubjects)])
        folderName = d(ii+2).name;
        d2 = dir([dataPathName,folderName,'\*run*']); % read in file directory information
        nRuns = size(d2,1);
        for jj=1:nRuns
            disp(['  Analyzing run ',num2str(jj),' of ',num2str(nRuns)])
            runNumber = jj;
            dummy = load([dataPathName,folderName,'\',d2(jj).name,'\',fileNameL]); % the clicks
            header = dummy.header;
            Clicks = dummy.data;
            clear dummy
            dummy = load([dataPathName,folderName,'\',d2(jj).name,'\',fileNameR]); % the noise
            headerN = dummy.header;
            Noise = dummy.data;
            clear dummy 

            [MEMR_inc,MEMR_mem,h1] = analyzeMEMR_v200(header,Clicks,headerN,Noise,folderName,runNumber);
            
            saveName = [folderName,'_Run',num2str(jj),'_Analysis1.mat'];
            save([savePath,saveName],'MEMR_inc','MEMR_mem')
            saveName = [folderName,'_Run',num2str(jj),'_Analysis1.bmp'];
            saveas(h1,[savePath,saveName],'bmp')

            %keyboard
            pause(0.01)
            c
        end

    end

end