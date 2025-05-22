function [] = comap()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% comap
% compare mapani and our test
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    %time = (0:.05:8-0.05)'; % this is for our swept test
    Lmin = 40; % y-axis offset for rising slope
    Lmax = 110; % y-axis offset for falling slope
    m = 17.5; % slope (110 dB - 40 dB / 4 seconds)
    %L1 = m*time + Lmin; % for swept test, level as a function of time (rising slope)
    %L2 = -m*(time-4) + Lmax; % for swept test, level as a function of time (falling slope)
    % time for the discrete (mapani) test as a function of level:
    %l1 = (37.5:5:92.5)';
    %l2 = (90:-5:35)';
    l1 = (40:10:110)';
    l2 = (100:-10:40)';
    T1 = (l1 - Lmin) / m;
    T2 = (l2 - 4*m - Lmax) / -m;
    tt = [T1;T2]; % plotting time

    parentDrive = 'D'; % use 'C' or 'D'
    dataPathName = [parentDrive,':\MEMR_AnalysisM\']; % location of raw data
    fileNameL = 'Ch3_ER10xA_memr_0001.mat';
    fileNameR = 'Ch4_ER10xB_memr_0001.mat';
    savePath = [parentDrive,':\MEMR_AnalysisM\']; % where to save analyzed data

    d = dir(dataPathName); % read in file directory information
    nSubjects = size(d,1)-2; % number of subject folders

    for ii=1:nSubjects
        disp(['subject ',num2str(ii)])
        fileName = d(ii+2).name;
  
        figure (ii)
        load([dataPathName,fileName])
        plot (MEMR_mem.timeMepani, 20*log10(MEMR_mem.d1),'b',LineWidth=2)
        
        hold on
    end


    parentDrive = 'D'; % use 'C' or 'D'
    dataPathName = [parentDrive,':\MEMR_AnalysisComp\']; % location of raw data
    fileNameL = 'Ch3_ER10xA_memr_0001.mat';
    fileNameR = 'Ch4_ER10xB_memr_0001.mat';
    savePath = [parentDrive,':\MEMR_AnalysisM\']; % where to save analyzed data

    d = dir(dataPathName); % read in file directory information
    nSubjects = size(d,1)-2; % number of subject folders
    n =0;
    for ii=1:4:nSubjects
        n = n + 1;
        disp(['subject ',num2str(n)])
        fileName1 = d(ii+2).name;
        fileName2 = d(ii+3).name;
        fileName3 = d(ii+4).name;
        fileName4 = d(ii+5).name;

        figure (n)
        load([dataPathName,fileName1])
        plot (MEMR_mem.t-MEMR_mem.delay, 20*log10(MEMR_mem.d1), 'r')
        hold on
        load([dataPathName,fileName2])
        plot (MEMR_mem.t-MEMR_mem.delay, 20*log10(MEMR_mem.d1), 'r')
        load([dataPathName,fileName3])
        plot (MEMR_mem.t-MEMR_mem.delay, 20*log10(MEMR_mem.d1), 'r')
        load([dataPathName,fileName4])
        plot (MEMR_mem.t-MEMR_mem.delay, 20*log10(MEMR_mem.d1), 'r')
        title (['subject ',num2str(n)])
        xticks(tt)
        xticklabels(num2str([l1;l2]))
        xlabel('Elicitor Level (dB SPL)')
        ylabel('Total Change (dB)')
        xlim([tt(1),tt(end)])
    end

FolderName = 'D:\data\';   % Your destination folder
FigList = findobj(allchild(0), 'flat', 'Type', 'figure');
for iFig = 1:length(FigList)
  FigHandle = FigList(iFig);
  FigName   = get(FigHandle, 'Name');
  savefig(FigHandle, fullfile(FolderName, FigName, '.fig'));
end

end



% OLD CODE ----------------------------------------------------------------
% plot(T1,l1,'r')
% hold on
% plot(T2,l2,'r')
% plot(time,L1,'b:')
% plot(time,L2,'b:')

% convert to db for y axis
% plot all of ours (each a separate red line
% mepani thick blue line (plot last so on top of ours)
% tltle showing subject number