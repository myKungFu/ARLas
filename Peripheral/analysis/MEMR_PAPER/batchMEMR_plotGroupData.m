function [] = batchMEMR_plotGroupData()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% batchMEMR_plotGroupData
% 
% Batch analysis of group MEMR data 
%
% Author: Shawn Goodman & Ehsan Khalili
% Date: January 24, 2024
% Last Updated: January 24, 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    parentDrive = 'D'; % use 'C' or 'D'
    dataPathName = [parentDrive,':\MEMR_Analysis1\']; % location of raw data
    savePath = [parentDrive,':\MEMR_GroupData\']; % where to save analyzed data

    d = dir([dataPathName,'*.mat']); % read in file directory information
    nFiles = size(d,1); % number of subject folders
    counter = 1;
    MEMR_memTotal = 0 ;
    for ii=1:nFiles
        disp(['Analyzing file ',num2str(ii),' of ',num2str(nFiles)])
        fileName = d(ii).name;
        runNumber = str2num(fileName(10));
        dummy = load([dataPathName,fileName]); 
        MEMR_mem = dummy.MEMR_mem;

        % extract data as a row vector

        nd1 = 20*log10(MEMR_mem.d1);
        nd1 = nd1 ./ max(nd1);
        ND1(:,ii) = nd1;

        figure(10)
        hold on
        plot(MEMR_mem.t,nd1,'Color',[.7 .7 .7])
        pause(0.01)


        % figure(11)
        % hold on
        % if mod(ii,5)==0
        %     hold off
        % end
        % plot(MEMR_mem.trend,'Color',[.7 .7 .7])
        % if mod(ii,4)==0
        %    pause
        % else
        %    pause(0.01)
        % end
        % 


        if mod(ii,4)==0
            counter = counter + 1;
        end
    end
nd1 = mean(ND1,2);



figure(10)
xlabel('Time (s)')
ylabel('Change (dB)')
plot(MEMR_mem.t,nd1,'k','LineWidth',2)

load('D:\MEMR_GroupData\MEMR_groupData.mat')


x = 20*log10(pAmp(:));
q5 = quantile(x,.5);
q25 = quantile(x,.25);
q75 = quantile(x,.75);
iqr = q75-q25;
whLow = max([min(20*log10(pAmp(:))),q25-1.5*iqr]);
whHigh = min([max(20*log10(pAmp(:))),q75+1.5*iqr]);

%plot([4,4],[q25,q75],'LineWidth',5,'Color',[0 1 0])
%plot([4,4],[whLow,whHigh],'LineWidth',2,'Color',[0 1 0])
%plot([4,4],[q5,q5],'.r','MarkerSize',14)


second =([zeros(9,1);flipud(nd1(85:end))]);
first = nd1(1:85);
figure
plot(first)
hold on
plot(second)


 keyboard

end