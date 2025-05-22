function [] = mepgod_vSSG()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% mepgod_vSSG
% comap
% compare Mepani and our test
%
% Authors: Shawn Goodman & Ehsan Khalili
% Date: May 21, 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    doDelayCorrection = 1; % apply delay correction (=1) or not (=0)
    doPlotResiduals = 0; % add a subplot of residuals (Delta values)
    %plotGood = 1; % plot good ones (=1) or bad ones (=0)
    parentDrive = 'D'; % use 'C' or 'D'
    dataPathNameM = [parentDrive,':\MEMR_AnalysisM\']; % location of raw data for mepani
    dataPathName = [parentDrive,':\MEMR_Analysis1\']; % location of raw data for mepani
    savePath = [parentDrive,':\Mepgod\']; % where to save analyzed data for mepani

    badNames = {'07','10','11','16','22','26','29','30','34','36'}; % these have no good data, so exlcude them

    % initialize output
    Delta_dB = [];
    TimeComps = [];
    DELTA = [];
    TIME = [];

    dM = dir([dataPathNameM,'*.mat']); % read in file directory information
    nSubjects = size(dM,1); % number of subject folders
    for ii=1:nSubjects
        disp(['Processing subject ',num2str(ii)])
        fileNameM = dM(ii).name;
        nameTag = fileNameM(4:5);
        if ~any(strcmp(nameTag,badNames)) % if the current name is on the "good list", process it
            
            % read in the mepani data; there is one run of this
            dummy = load([dataPathNameM,fileNameM]);
            tM = dummy.MEMR_mem.timeMepani;
            d1M = dummy.MEMR_mem.d1;
            d1M = 20*log10(d1M);

            
            % read in the goodman data; there are 4 runs of this
            nameTagGood = fileNameM(1:5); % look for this file name in the goodman data set
            d = dir([dataPathName,[nameTagGood,'*.mat']]); % read in file directory information
            nRuns = size(d,1); % number of runs
            for jj=1:nRuns
                fileName = d(jj).name;
                dummy = load([dataPathName,fileName]);
                if jj==1
                    tTemplate = dummy.MEMR_mem.t; % template time delay values
                end
                t = dummy.MEMR_mem.t - dummy.MEMR_mem.delay; % time delay corrected for delay
                d1 = dummy.MEMR_mem.d1; % total change
                if doDelayCorrection == 1
                    % use spline interpolation to put corrected onto the template time vector
                    pp = spline(t,d1);
                    d1 = ppval(pp,tTemplate);
                    ThdOn(jj,1) = dummy.MEMR_mem.thdOnsetTime - dummy.MEMR_mem.delay; % onset threshold time
                    ThdOff(jj,1) = dummy.MEMR_mem.thdOffsetTime - dummy.MEMR_mem.delay; % offset threshold time
                else
                    ThdOn(jj,1) = dummy.MEMR_mem.thdOnsetTime; % onset threshold time
                    ThdOff(jj,1) = dummy.MEMR_mem.thdOffsetTime; % offset threshold time
                end
                D1(:,jj) = 20*log10(d1);
            end
            % take the median values for comparison purposes
            d1 = median(D1,2);
            thdOn = median(ThdOn);
            thdOff = median(ThdOff);
            % use a cubic spline to enable finding any time value of the goodman data 
            pp = spline(tTemplate,d1);

            % compare only the values above threshold in both cases
            % find the mepani values that are above threshold
            if mod(size(d1M,2),2)~=0
                keyboard
            else
                halfN = size(d1M,2)/2;
            end
            % process the lower half -----
            tt = tM(1:halfN); % half the time vector
            dd1m = d1M(1:halfN); % half the total change vector
            indxt = find(tt>thdOn); % indices of times to compare
            MepaniComps = dd1m(indxt); % comparison values from mepani
            GoodmanComps = ppval(pp,tt(indxt)); % comparison values from goodman
            delta_dB = GoodmanComps(:) - MepaniComps(:); % dB difference (Mepani is reference)
            timeComps = tt(indxt);
            % save this for group comparison
            Delta_dB = [Delta_dB;delta_dB];
            TimeComps = [TimeComps;timeComps];
            % process the upper half -----
            tt = tM(halfN+1:end); % half the time vector
            dd1m = d1M(halfN+1:end); % half the total change vector
            indxt = find(tt<thdOff); % indices of times to compare
            MepaniComps = dd1m(indxt); % comparison values from mepani
            GoodmanComps = ppval(pp,tt(indxt)); % comparison values from goodman
            delta_dB = GoodmanComps(:) - MepaniComps(:); % dB difference (Mepani is reference)
            timeComps = tt(indxt);
            % save this for group comparison
            Delta_dB = [Delta_dB;delta_dB];
            TimeComps = [TimeComps;timeComps];

            DELTA = [DELTA;Delta_dB];
            TIME = [TIME;TimeComps];


            % plotting ----------------------------------------------------
            h = figure;
            if doPlotResiduals == 1
                sp1 = subplot(2,1,1);
                sp2 = subplot(2,1,2);
                sp1.Position = [0.1300    0.4381    0.7750    0.4869];
                sp2.Position = [0.1300    0.1024    0.7750    0.1833];
                axes(sp1)
            end
            plot (tTemplate,D1,'Color',[0.6 0.6 0.6],'LineWidth',0.5)
            hold on
            plot(tTemplate,d1,'r','LineWidth',2)
            plot (tM(1:halfN),d1M(1:halfN),'b','LineWidth',1)
          
            p = plot (tM(1:halfN),d1M(1:halfN),'o','Color',[1 1 1],'LineWidth',1);
                p.MarkerFaceColor = [1 1 1];
                p.MarkerSize = 6;
                p.MarkerEdgeColor = [0 0 1];
            
            plot (tM(halfN+1:end),d1M(halfN+1:end),'b','LineWidth',1)

            p = plot (tM(halfN+1:end),d1M(halfN+1:end),'o','Color',[1 1 1],'LineWidth',1);
                p.MarkerFaceColor = [1 1 1];
                p.MarkerSize = 6;
                p.MarkerEdgeColor = [0 0 1];

            xmin = 0;
            xmax = 8;
            ymin = 0;
            ymax = 4;
            if doPlotResiduals == 1
                for jj=1:nRuns
                    line([ThdOn(jj),ThdOn(jj)],[ymin,ymax],'LineWidth',0.5,'LineStyle',':','Color',[0 0 0])
                    line([ThdOff(jj),ThdOff(jj)],[ymin,ymax],'LineWidth',0.5,'LineStyle',':','Color',[0 0 0])
                end
            end
            line([thdOn,thdOn],[ymin,ymax],'LineWidth',1,'LineStyle','--','Color',[0 0 0])
            line([thdOff,thdOff],[ymin,ymax],'LineWidth',1,'LineStyle','--','Color',[0 0 0])
            grid off
            xlim([xmin,xmax])
            ylim([ymin,ymax])
            xlabel('Time (s)')
            ylabel('Total Change (dB)')
            title(nameTagGood)
            if doPlotResiduals == 1
                axes(sp2)
                plot(TimeComps,Delta_dB,'*-')
                xlabel('Time (s)','FontSize',12)
                ylabel('Delta (dB)','FontSize',12)
                grid off
                xlim([xmin,xmax])
                ylim([-1,1])
            end
            
            saveName = ['Participant ',num2str(ii),'_Mepgod1.bmp'];
            saveas(gcf,[savePath,saveName],'bmp')
            c
            %keyboard
            
            Delta_dB = [];
            TimeComps = [];
            
        end


% added 5/20/2024
q1 = d1(:,1);
t1 = tTemplate';
[qmax,imax] = max(q1);
criterion = 0.5;
[~,indx1] = min(abs(q1(1:imax)-(qmax-criterion)));
[~,indx2] = min(abs(q1(imax:end)-(qmax-criterion)));
indx2 = indx2 + imax;
qq = q1(indx1:indx2);
tt = t1(indx1:indx2);
pp = polyfit(tt,qq,2);
yy = polyval(pp,tt);
resid = sqrt(mean((yy-qq).^2)); % Quantify residuals to look at, check which one has higher and try to remove it !

dydt = gradient(yy)./gradient(tt);
dydt2 = gradient(dydt)./gradient(tt);
A2(:,ii) = median(dydt2);
A1(:,ii) = polyder(polyder(pp));






    end


    figure
    boxplot(DELTA,TIME)
    hold on
    ymin = -1;
    ymax = 2;
    xmin = 0;
    xmax = 12;
    line([xmin,xmax],[0,0],'LineWidth',0.5,'LineStyle','--','Color',[0 0 0])
    xlim([xmin,xmax])
    ylim([ymin,ymax])
    mean(DELTA)
    std(DELTA)
    if strcmp(doDelayCorrection,'0')
        title('Group Data WITHOUT Delay Correction')
    else
        title('Group Data WITH Delay Correction')
    end
    ylabel('Goodman - Mepani (dB)')
            
    
    saveName = ['GD',num2str(ii),'_Mepgod1.bmp'];
            saveas(gcf,[savePath,saveName],'bmp')
keyboard
end

% OLD CODE ----------------------------------------------------------------

