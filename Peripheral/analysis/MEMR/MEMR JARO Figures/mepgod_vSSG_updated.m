function [] = mepgod_vSSG_updated()

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% mepgod_vSSG
% compare Mepani and our test
%
% Authors: Shawn Goodman & Ehsan Khalili
% Date: Jan 18, 2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

FS = 14; % font size
TS = 12; % tick size

doDelayCorrection = 1; % turn delay correction on/off
doPlotResiduals   = 0; % add the little residual subplot or not

parentDrive  = 'D'; % switch to C if needed
dataPathNameM = [parentDrive,':\myWork\Data\MEMR_DATA\MEMR_AnalysisM_original\'];
dataPathName  = [parentDrive,':\myWork\Data\MEMR_DATA\MEMR_Analysis_Fix1\'];
savePath      = [parentDrive,':\myWork\Data\MEMR_DATA\Mepgod_Paper\'];

badNames = {}; % skip these (no usable data)

Delta_dB   = [];
TimeComps  = [];
DELTA      = [];
TIME       = [];

dM = dir([dataPathNameM,'*.mat']);
nSubjects = size(dM,1);

for ii = 1:nSubjects
    disp(['Processing subject ',num2str(ii)])
    fileNameM = dM(ii).name;
    nameTag = fileNameM(4:5);

    if ~any(strcmp(nameTag,badNames))

        dummy = load([dataPathNameM,fileNameM]);
        tM = dummy.MEMR_mem.timeMepani;
        d1M = dummy.MEMR_mem.d1;
        d1M = 20*log10(d1M);

        nameTagGood = fileNameM(1:5);
        d = dir([dataPathName,[nameTagGood,'*.mat']]);
        nRuns = size(d,1);

        for jj = 1:nRuns
            fileName = d(jj).name;
            dummy = load([dataPathName,fileName]);

            if jj == 1
                tTemplate = dummy.MEMR_mem.t;
            end

            t = dummy.MEMR_mem.t - dummy.MEMR_mem.delay;
            d1 = dummy.MEMR_mem.d1;

            if doDelayCorrection == 1
                pp = spline(t,d1);
                d1 = ppval(pp,tTemplate);
                ThdOn(jj,1)  = dummy.MEMR_mem.thdOnsetTime  - dummy.MEMR_mem.delay;
                ThdOff(jj,1) = dummy.MEMR_mem.thdOffsetTime - dummy.MEMR_mem.delay;
            else
                ThdOn(jj,1)  = dummy.MEMR_mem.thdOnsetTime;
                ThdOff(jj,1) = dummy.MEMR_mem.thdOffsetTime;
            end

            D1(:,jj) = 20*log10(d1);
        end

        d1     = median(D1,2);
        thdOn  = median(ThdOn);
        thdOff = median(ThdOff);

        pp = spline(tTemplate,d1);

        if mod(size(d1M,2),2)~=0
            keyboard
        else
            halfN = size(d1M,2)/2;
        end

        tt   = tM(1:halfN);
        dd1m = d1M(1:halfN);
        indxt = find(tt > thdOn);

        MepaniComps = dd1m(indxt);
        GoodmanComps = ppval(pp,tt(indxt));
        delta_dB = GoodmanComps(:) - MepaniComps(:);
        timeComps = tt(indxt);

        Delta_dB  = [Delta_dB;delta_dB];
        TimeComps = [TimeComps;timeComps];

        tt   = tM(halfN+1:end);
        dd1m = d1M(halfN+1:end);
        indxt = find(tt < thdOff);

        MepaniComps = dd1m(indxt);
        GoodmanComps = ppval(pp,tt(indxt));
        delta_dB = GoodmanComps(:) - MepaniComps(:);
        timeComps = tt(indxt);

        Delta_dB  = [Delta_dB;delta_dB];
        TimeComps = [TimeComps;timeComps];

        DELTA = [DELTA;Delta_dB];
        TIME  = [TIME;TimeComps];

        figure('Position', [1000, 500, 500, 450]);

        if doPlotResiduals == 1
            sp1 = subplot(2,1,1);
            sp2 = subplot(2,1,2);
            sp1.Position = [0.1300    0.4381    0.7750    0.4869];
            sp2.Position = [0.1300    0.1024    0.7750    0.1833];
            axes(sp1)
        end

        plot(tTemplate,D1,'Color',[0.6 0.6 0.6],'LineWidth',0.5)
        hold on
        plot(tTemplate,d1,'r','LineWidth',2)

        plot(tM(1:halfN),d1M(1:halfN),'b','LineWidth',1)
        p = plot(tM(1:halfN),d1M(1:halfN),'o','Color',[1 1 1],'LineWidth',1);
        p.MarkerFaceColor = [1 1 1];
        p.MarkerSize = 6;
        p.MarkerEdgeColor = [0 0 1];

        plot(tM(halfN+1:end),d1M(halfN+1:end),'b','LineWidth',1)
        p = plot(tM(halfN+1:end),d1M(halfN+1:end),'o','Color',[1 1 1],'LineWidth',1);
        p.MarkerFaceColor = [1 1 1];
        p.MarkerSize = 6;
        p.MarkerEdgeColor = [0 0 1];

        xmin = 0;
        xmax = 8;
        ymin = 0;
        ymax = 4;

        if doPlotResiduals == 1
            for jj = 1:nRuns
                line([ThdOn(jj),ThdOn(jj)],[ymin,ymax],'LineWidth',0.5,'LineStyle',':','Color',[0 0 0])
                line([ThdOff(jj),ThdOff(jj)],[ymin,ymax],'LineWidth',0.5,'LineStyle',':','Color',[0 0 0])
            end
        end

        line([thdOn,thdOn],[ymin,ymax],'LineWidth',1,'LineStyle','--','Color',[0 0 0])
        line([thdOff,thdOff],[ymin,ymax],'LineWidth',1,'LineStyle','--','Color',[0 0 0])

        set(gca, 'FontSize', TS);
        grid off;
        xlim([xmin,xmax])
        ylim([ymin,ymax])
        xticks([0,1,2,3,4,5,6,7,8])
        xticklabels({'40','58','75','93','110','93','75','58','40'})
        xlabel('Elicitor Level (dB SPL)','FontSize', FS)
        ylabel('Total Change (dB)','FontSize', FS)
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

        saveName = ['Participant ',num2str(ii),'_Mepgod1.svg'];
        saveas(gcf,[savePath,saveName],'svg')

        Delta_dB  = [];
        TimeComps = [];
    end

    q1 = d1(:,1);
    t1 = tTemplate';
    [qmax,imax] = max(q1);
    criterion = 0.5;

    [~,indx1] = min(abs(q1(1:imax) - (qmax-criterion)));
    [~,indx2] = min(abs(q1(imax:end) - (qmax-criterion)));
    indx2 = indx2 + imax;

    qq = q1(indx1:indx2);
    tt = t1(indx1:indx2);

    pp = polyfit(tt,qq,2);
    yy = polyval(pp,tt);

    resid = sqrt(mean((yy-qq).^2)); % just a quick residual check

    dydt = gradient(yy)./gradient(tt);
    dydt2 = gradient(dydt)./gradient(tt);

    A2(:,ii) = median(dydt2);
    A1(:,ii) = polyder(polyder(pp));
end

gapStart = 3.50;
gapEnd   = 5.14;
keepIdx = (TIME < gapStart) | (TIME > gapEnd);

DELTAp = DELTA(keepIdx);
TIMEp  = TIME(keepIdx);

timeBins = unique(round(TIMEp*1e5)/1e5);
timeBins = sort(timeBins);

figure('Color','w','Position',[200 200 1100 700]);
ax = axes('Parent',gcf);
set(ax,'Color','w'); hold(ax,'on');

boxplot(DELTAp, TIMEp, ...
    'Positions', timeBins, ...
    'Widths', 0.22, ...
    'Symbol', '', ...
    'Whisker', 1);
hold on

set(findobj(ax,'Tag','Box'),    'Color',[0 0 1],'LineWidth',1.5);
set(findobj(ax,'Tag','Median'), 'Color',[1 0 0],'LineWidth',2);

set(findobj(ax,'Tag','Whisker'),'Visible','on');
set(findobj(ax,'Tag','Cap'),    'Visible','on');
set(findobj(ax,'Tag','Outliers'),'Visible','on');

set(ax,'XColor','k','YColor','k', 'FontSize', TS, 'Box','off', 'LineWidth',1.5);

tTicks = 0:0.5:8;
ax.XTick = tTicks;

Lticks = zeros(size(tTicks));
upIdx = tTicks <= 4;
dnIdx = tTicks > 4;

Lticks(upIdx) = 40 + 17.5 * tTicks(upIdx);
Lticks(dnIdx) = 110 - 17.5 * (tTicks(dnIdx) - 4);

ax.XTickLabel = compose('%.0f', Lticks);

ax.YTickMode = 'auto';
ax.YTickLabelMode = 'auto';

xlabel('Elicitor Level (dB SPL)','FontSize',FS);
ylabel('Δ (Swept - Mepani) (dB)','Interpreter','none');

ymin = -1;
ymax = 2;
xmin = 0;
xmax = 8;

xlim([xmin,xmax]);
ylim([ymin,ymax]);

line([xmin,xmax],[0,0],'LineWidth',1,'LineStyle','--','Color','k')
title('')

S = load(fullfile(savePath,'MEK_results_matched_discrete.mat'));
MEK = S.MEK_results;

elicAll = MEK.ELIC_ALL(:);
yMEK    = MEK.DIFF_ALL_THR;

nPts = numel(elicAll);
if mod(nPts,2) ~= 0
    warning('MEK overlay skipped: MEK.ELIC_ALL length is odd.');
else
    halfN = nPts/2;

    Lmin  = 40;
    Lmax  = 110;
    m     = 17.5;

    elicUp = elicAll(1:halfN);
    elicDn = elicAll(halfN+1:end);

    tUp = (elicUp - Lmin) / m;
    tDn = 4 + (Lmax - elicDn) / m;
    tMEK = [tUp; tDn];

    nPerPoint   = sum(~isnan(yMEK), 2);
    medPerPoint = median(yMEK, 2, 'omitnan');

    medMasked = medPerPoint;
    medMasked(nPerPoint < 5) = NaN; % don’t draw a line if too few people

    medUp = medMasked(1:halfN);
    medDn = medMasked(halfN+1:end);

    jitterAmp = 0.05;

    for s = 1:size(yMEK,2)
        ys = yMEK(:,s);
        ok = ~isnan(ys);

        xj = tMEK(ok) + jitterAmp*(rand(sum(ok),1)-0.5);
        plot(xj, ys(ok), 'o', ...
            'MarkerSize', 6.5, ...
            'MarkerFaceColor', 'none', ...
            'MarkerEdgeColor', 'k', ...
            'LineWidth', 1.2);
    end

    plot(tUp, medUp, 'k-', 'LineWidth', 3);
    plot(tDn, medDn, 'k-', 'LineWidth', 3);
end

saveName = 'FIG14_Grouping_Mepgod.svg';
saveas(gcf,[savePath,saveName],'svg');

muMEK = mean(yMEK(:),'omitnan');
sdMEK = std(yMEK(:),'omitnan');
fprintf('MEK (Swept - Mepani) diffs: mean = %.4f dB, std = %.4f dB (NaNs omitted)\n', muMEK, sdMEK);

end
