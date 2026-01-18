function [] = MEMR_Fig6()
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


    parentDrive = 'C'; % use 'C' or 'D'
    dataPathName = [parentDrive,':\myWork\ARLas\MEMR_DATA\']; % location of raw data
    fileNameL = 'Ch3_ER10xA_memr_0001.mat';
    fileNameR = 'Ch4_ER10xB_memr_0001.mat';
    savePath = [parentDrive,':\myWork\ARLas\MEMR_Analysis_Final_15\']; % where to save analyzed data

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

        
        h1 = figure;
        sp1 = subplot(2,1,1);
        sp2 = subplot(2,1,2);
        sp1.Position = [0.1300    0.7230    0.7750    0.2020];
        sp2.Position = [0.1300    0.1100    0.7750    0.4924];
      axes(sp1)
        ax = gca;
        ax.FontSize = 11; 
        plot(MEMR_mem.timeTrend,20*log10(MEMR_mem.trend),'b')
        hold on
         plot(MEMR_mem.timeTrend,20*log10(MEMR_inc.trend),'b--')
        Q = repmat(20*log10(MEMR_mem.d1),1,15);
        Q = Q(:);
        plot(MEMR_mem.timeTrend,20*log10(MEMR_mem.trend)+Q,'k')
        xlabel('Time (s)','FontSize',11)
        ylabel('Change (dB)','FontSize',11)
        title([subjectName,'    Run# ',num2str(runNumber)])
        grid on
            hold on 
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