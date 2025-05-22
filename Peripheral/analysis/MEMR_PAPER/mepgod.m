function [] = mepgod()
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
    mepgod1 = zeros(30,4,160);
    mepgod2 = zeros(30,24);
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
         if fileName == 'MEM07_Mepani_Analysis1.mat'
            mepgod2(ii, :) = NaN;
         elseif   fileName == 'MEM10_Mepani_Analysis1.mat'
            mepgod2(ii, :) = NaN;
              elseif   fileName == 'MEM11_Mepani_Analysis1.mat'
            mepgod2(ii, :) = NaN;
              elseif   fileName == 'MEM16_Mepani_Analysis1.mat'
            mepgod2(ii, :) = NaN;
              elseif   fileName == 'MEM22_Mepani_Analysis1.mat'
            mepgod2(ii, :) = NaN;
              elseif   fileName == 'MEM26_Mepani_Analysis1.mat'
            mepgod2(ii, :) = NaN;
              elseif   fileName == 'MEM29_Mepani_Analysis1.mat'
            mepgod2(ii, :) = NaN;
              elseif   fileName == 'MEM30_Mepani_Analysis1.mat'
            mepgod2(ii, :) = NaN;
              elseif   fileName == 'MEM34_Mepani_Analysis1.mat'
            mepgod2(ii, :) = NaN;
              elseif   fileName == 'MEM36_Mepani_Analysis1.mat'
            mepgod2(ii, :) = NaN;
         else
             mepgod2(ii, :) = 20*log10(MEMR_mem.d1);
         end
        plot (MEMR_mem.timeMepani, 20*log10(MEMR_mem.d1),'b',LineWidth=2)
        
        hold on
    end


    meparentDrive = 'D'; % use 'C' or 'D'
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
        
        mepgod1(n,1, :) = 20*log10(MEMR_mem.d1);
        hold on
        
        load([dataPathName,fileName2])
        plot (MEMR_mem.t-MEMR_mem.delay, 20*log10(MEMR_mem.d1), 'r')
        mepgod1(n,2, :) = 20*log10(MEMR_mem.d1);
        load([dataPathName,fileName3])
        plot (MEMR_mem.t-MEMR_mem.delay, 20*log10(MEMR_mem.d1), 'r')
        mepgod1(n,3, :) = 20*log10(MEMR_mem.d1);
        load([dataPathName,fileName4])
        plot (MEMR_mem.t-MEMR_mem.delay, 20*log10(MEMR_mem.d1), 'r')
        mepgod1(n,4, :) = 20*log10(MEMR_mem.d1);
        title (['subject ',num2str(n)])
        xticks(tt)
        xticklabels(num2str([l1;l2]))
        xlabel('Elicitor Level (dB SPL)')
        ylabel('Total Change (dB)')
        xlim([tt(1),tt(end)])
    end
    mepgod1 = squeeze(mepgod1);
    median_mepgod1 = median(mepgod1, 2);
median_mepgod1 = squeeze(median_mepgod1);
n = 0;    
    for ii=1:4:nSubjects
        n = n + 1;
        disp(['subject ',num2str(n)])
        fileName1 = d(ii+2).name;
        load([dataPathName,fileName1])
        figure (n)
        plot (MEMR_mem.t-MEMR_mem.delay, (median_mepgod1(n,:)), 'g')

    end
for participant = 1:size(mepgod2, 1)
    for i = 1:12

        if i <= 10 && all(mepgod2(participant, i:i+4) < mepgod2(participant, i+1:i+5))
            cutoff_index = i;
            mepgod2(participant, 1:cutoff_index-1) = NaN;
            break; 
        end
    end
   
    for i = 13:24
        if i <= 23 && mepgod2(participant, i) < mepgod2(participant, i+1)
            continue;
        end
       
        if i <= 23 && mepgod2(participant, i) > mepgod2(participant, i+1)
            if mepgod2(participant, i+1) < mepgod2(participant, i+2)
                cutoff_index = i+1;
                mepgod2(participant, cutoff_index+1:end) = NaN;
            end
            break; 
        end
    end
end

differences = zeros(30, 24);
for participant = 1:30
    differences(participant, 1:12) = median_mepgod1(participant, 4:4:48) - mepgod2(participant, 1:12);
    differences(participant, 13:24) = median_mepgod1(participant, 106:4:150) - mepgod2(participant, 13:24);
end


figure;
% differences = 20*log10(differences);
boxplot(differences);
xlabel('Elicitor Level (dB SPL)');
 l1 = (44:4:92)';
    l2 = (88:-4:48)';
    T1 = (l1 - Lmin) / m;
    T2 = (l2 - 4*m - Lmax) / -m;
xticks(1:24)
xticklabels(num2str([l1;l2]))
ylabel('Difference in Values');
title('Discrete vs Sweeping');

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