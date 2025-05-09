function [h1,h2] = ARLas_dpoaePlotting_fixedF1Continuous_v1(DPOAE)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [h] = ARLas_dpoaePlotting_fixedF1Continuous_v1(DPOAE);
%
% Plot swept DPOAE analyses; goes with ARLas_dpoae_fixedRatioContinuous_v1.m.
%
% Author: Shawn Goodman, PhD
% Auditory Research Lab, the University of Iowa
% Date: March 7, 2022
% Last Updated: March 7, 2022 -- ssg
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % unpack needed fields of the DPOAE structure ---
    subjectID = DPOAE.subjID;
    timeStamp = DPOAE.timeStamp;
    ear = DPOAE.ear;
    targetL1 = DPOAE.targetL1;
    nLevels = length(targetL1);
    targetL2 = DPOAE.targetL2;
    targetCalType = DPOAE.targetCalType;
    L1 = DPOAE.L1;
    N1 = DPOAE.N1;
    L2 = DPOAE.L2;
    N2 = DPOAE.N2;
    Ldp = DPOAE.Ldp;
    Ndp = DPOAE.Ndp;
    f1 = DPOAE.f1;
    f2 = DPOAE.f2;
    
    Ldp_3oct = DPOAE.Ldp_3oct;
    Ndp_3oct = DPOAE.Ndp_3oct;
    fc = DPOAE.fc;

    % setup for octave plots
    tickLocations = log2([1000,1500,2000,3000,4000,6000,8000,12000,16000]);
    tickLabels = {'1','1.5','2','3','4','6','8','12','16'};
    logf2 = log2(f2);
    logfc = log2(fc);
    logf1 = log2(f1);
    xmin = log2(min(f2));
    xmax = log2(max(f2));
    xmin1 = log2(min(f1));
    ymin = -25;
    ymax = 30;
    
    % plot primaries ------------------------------------------------------
    if strcmp(ear,'Left') || strcmp(ear,'left')
        h1 = figure(2);
    elseif strcmp(ear,'Right') || strcmp(ear,'right')
        h1 = figure(4);
    else
        h1 = figure(6);
    end
    co = get(gca,'colororder');    
    
        hold off
        for kk=1:nLevels
            plot(logf1(:,kk),L1(:,kk),'Color',co(kk,:),'LineWidth',1)
            hold on
        end
        for kk=1:nLevels
            plot(logf2(:,kk),L2(:,kk),'Color',co(kk,:),'LineWidth',1)
            hold on
        end
        for kk=1:nLevels
            plot(logf1(:,kk),N1(:,kk),'Color',co(kk,:),'LineStyle','--','LineWidth',0.5)
            hold on
        end
        for kk=1:nLevels
            plot(logf2(:,kk),N2(:,kk),'Color',co(kk,:),'LineStyle','-.','LineWidth',0.5)
            hold on
        end
        xticks(tickLocations)
        xticklabels(tickLabels)
        xlabel('Frequency (Hz)')
        ylabel('Primary F2 Level (dB SPL)')
        if xmax <= xmin
            xmax = xmin + 1;
        end
        xlim([xmin,xmax])
        title([subjectID,'  ',ear,' Ear  DPOAE Primary Levels  ',char(cellstr(timeStamp))])

        
    % plot dpoaes % -------------------------------------------------------
    if strcmp(ear,'Left') || strcmp(ear,'left')
        h2 = figure(3);
    elseif strcmp(ear,'Right') || strcmp(ear,'right')
        h2 = figure(5);
    else
        h2 = figure(7);
    end
        hold off
        for kk=1:nLevels
            plot(logf2(:,kk),Ldp(:,kk),'Color',co(kk,:),'LineStyle','-','LineWidth',0.5)
            hold on
        end
        for kk=1:nLevels
            plot(logf2(:,kk),Ndp(:,kk),'Color',co(kk,:),'LineStyle','--','LineWidth',0.5)
            hold on
            plot(logf2(:,kk),Ndp(:,kk)+6,'Color',[0 0 0],'LineStyle','--','LineWidth',0.5)
        end
        xticks(tickLocations)
        xticklabels(tickLabels)
        xlabel('F2 Frequency (Hz)')
        ylabel('DPOAE Level (dB SPL)')
        xlim([xmin,xmax])
        %ylim([ymin,ymax])
        txt1 = ['L1,L2 = ',round(num2str(targetL1(1))),',',round(num2str(targetL2(1))),' dB ',targetCalType];
        txt2 = 'Noise Floor';
        txt3 = 'noise Floor + 6 dB';
            legend(txt1,txt2,txt3,'Location','best')            
        title([subjectID,'  ',ear,' Ear  DPOAE Levels  ',char(cellstr(timeStamp))])

end

% OLD CODE ----------------------------------------------------------------
%         if nLevels > 1
%             txt2 = ['L1,L2 = ',round(num2str(targetL1(2))),',',round(num2str(targetL2(2))),' dB ',targetCalType];
%         end
%         if nLevels > 2
%             txt3 = ['L1,L2 = ',round(num2str(targetL1(3))),',',round(num2str(targetL2(3))),' dB ',targetCalType];
%         end
%         if nLevels == 1
%             legend(txt1)
%         elseif nLevels == 2
%             legend(txt1,txt2)
%         else
