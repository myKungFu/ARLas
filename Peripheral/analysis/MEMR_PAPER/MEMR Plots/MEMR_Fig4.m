function [] = MEMR_Fig4()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% series of figures to show how magnitude non-monotonicity can be avoided
% and variability significantly reduced by computing total change instead.
%
% Author: Shawn Goodman
% DAte: January 28, 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    load('C:\myWork\ARLas\MEMR_Analysis_Final_15\MEM06_Run1_Analysis1.mat')
    D = MEMR_mem.D(:,5:15);
    freq = MEMR_mem.freq(5:15);
    % ---------------------------------------
    figure
    plot(20*log10(abs(D)),'.-','Color',[.7 .7 .7],'LineWidth',0.5,'MarkerSize',10)
    hold on
    plot(20*log10(abs(D(:,7))),'.-','Color',[1 0 0],'LineWidth',0.5,'MarkerSize',10)
    plot(20*log10(abs(D(:,8))),'.-','Color',[0 0 1],'LineWidth',0.5,'MarkerSize',10)
    grid on
    line([0 160],[0 0],'Color',[0 0 0],'LineWidth',0.5,'LineStyle','-')
    xticks([0,20,40,60,80,100,120,140,160])
    xticklabels({'40','58','75','93','110','93','75','58','40'})
    %levels = 17.5*0.05([0,20,40,60,80,60,40,20,0]')+40
    xlabel('Elicitor Level (dB SPL)','FontSize',14)
    ylabel('\Delta Magnitude (dB)','FontSize',14)
    set(gca,'fontsize',11)
    xlim([0,160])
    ylim([-6 6])
    clear all

    load('C:\myWork\ARLas\MEMR_Analysis_Final_15\MEM06_Run1_Analysis1.mat')
    D = MEMR_mem.D(:,5:15);
    freq = MEMR_mem.freq(5:15);
    figure
    % ---------------------------------------
    % create a line segment for u1
    Theta = linspace(0,2*pi,2000)';
    u = exp(1i*Theta);
    plot(u,'k')
    hold on
    xlim([-1.25,1.25])
    ylim([-1.25,1.25])
    line([-1 1],[0 0],'Color',[0 0 0],'LineWidth',0.5,'LineStyle','-')
    line([0 0],[-1 1],'Color',[0 0 0],'LineWidth',0.5,'LineStyle','-')
    grid on
    set(gca,'fontsize',11)
    xlabel('Real','FontSize',14)
    ylabel('Imaginary','FontSize',14)
    axis square

    plot(D(:,11),'.-','Color',[.7 .7 .7],'LineWidth',0.5,'MarkerSize',8)
    xlim([-.5,.5]+.75)
    ylim([-.5,.5])

    plot(D,'.-','Color',[.7 .7 .7],'LineWidth',0.5,'MarkerSize',8)
    plot(D(:,7),'.-','Color',[1 0 0],'LineWidth',0.5,'MarkerSize',8)
    plot(D(:,8),'.-','Color',[0 0 1],'LineWidth',0.5,'MarkerSize',8)
    xlim([.4 1.4])
    ylim([-.7,.3])


    % --------------------------------------
    % --------------------------------------
    figure
    load('C:\myWork\ARLas\MEMR_Analysis_Final_15\MEM06_Run1_Analysis1.mat')
    D = MEMR_mem.D(:,5:15);
    % create a line segment for u1
    Theta = linspace(0,2*pi,2000)';
    u = exp(1i*Theta);
    plot(u,'k')
    hold on
    line([-1 1],[0 0],'Color',[0 0 0],'LineWidth',0.5,'LineStyle','-')
    line([0 0],[-1 1],'Color',[0 0 0],'LineWidth',0.5,'LineStyle','-')
    grid on
    set(gca,'fontsize',11)
    xlabel('Real','FontSize',14)
    ylabel('Imaginary','FontSize',14)
    axis square

    %plot(D,'.-','Color',[.7 .7 .7],'LineWidth',0.5,'MarkerSize',8)
    plot(D(:,7),'.-','Color',[1 0 0],'LineWidth',0.5,'MarkerSize',8)
    plot(D(:,8),'.-','Color',[0 0 1],'LineWidth',0.5,'MarkerSize',8)
    xlim([.4 1.4])
    ylim([-.7,.3])

    load('D:\MEMR_Analysis1\MEM06_Run2_Analysis1.mat')
    D = MEMR_mem.D(:,5:15);
    plot(D(:,7),'.-','Color',[1 0 0],'LineWidth',0.5,'MarkerSize',8)
    plot(D(:,8),'.-','Color',[0 0 1],'LineWidth',0.5,'MarkerSize',8)
    load('D:\MEMR_Analysis1\MEM06_Run3_Analysis1.mat')
    D = MEMR_mem.D(:,5:15);
    plot(D(:,7),'.-','Color',[1 0 0],'LineWidth',0.5,'MarkerSize',8)
    plot(D(:,8),'.-','Color',[0 0 1],'LineWidth',0.5,'MarkerSize',8)
    load('D:\MEMR_Analysis1\MEM06_Run4_Analysis1.mat')
    D = MEMR_mem.D(:,5:15);
    plot(D(:,7),'.-','Color',[1 0 0],'LineWidth',0.5,'MarkerSize',8)
    plot(D(:,8),'.-','Color',[0 0 1],'LineWidth',0.5,'MarkerSize',8)


   figure
    load('C:\myWork\ARLas\MEMR_Analysis_Final_15\MEM06_Run1_Analysis1.mat')
    D = MEMR_mem.D(:,5:15);
    plot(20*log10(abs(D(:,7))),'.-','Color',[1 0 0],'LineWidth',0.5,'MarkerSize',10)
    hold on
    plot(20*log10(abs(D(:,8))),'.-','Color',[0 0 1],'LineWidth',0.5,'MarkerSize',10)
    grid on
    line([0 160],[0 0],'Color',[0 0 0],'LineWidth',0.5,'LineStyle','-')
    xticks([0,20,40,60,80,100,120,140,160])
    xticklabels({'40','58','75','93','110','93','75','58','40'})
    xlabel('Elicitor Level (dB SPL)','FontSize',14)
    ylabel('\Delta Magnitude (dB)','FontSize',14)
    set(gca,'fontsize',11)
    xlim([0,160])
    ylim([-6 6])
    load('D:\MEMR_Analysis1\MEM06_Run2_Analysis1.mat')
    D = MEMR_mem.D(:,5:15);
    plot(20*log10(abs(D(:,7))),'.-','Color',[1 0 0],'LineWidth',0.5,'MarkerSize',10)
    plot(20*log10(abs(D(:,8))),'.-','Color',[0 0 1],'LineWidth',0.5,'MarkerSize',10)
    load('D:\MEMR_Analysis1\MEM06_Run3_Analysis1.mat')
    D = MEMR_mem.D(:,5:15);
    plot(20*log10(abs(D(:,7))),'.-','Color',[1 0 0],'LineWidth',0.5,'MarkerSize',10)
    plot(20*log10(abs(D(:,8))),'.-','Color',[0 0 1],'LineWidth',0.5,'MarkerSize',10)
    load('D:\MEMR_Analysis1\MEM06_Run4_Analysis1.mat')
    D = MEMR_mem.D(:,5:15);
    plot(20*log10(abs(D(:,7))),'.-','Color',[1 0 0],'LineWidth',0.5,'MarkerSize',10)
    plot(20*log10(abs(D(:,8))),'.-','Color',[0 0 1],'LineWidth',0.5,'MarkerSize',10)


 figure
    load('C:\myWork\ARLas\MEMR_Analysis_Final_15\MEM06_Run1_Analysis1.mat')
    D = MEMR_mem.D1(:,5:15);
    plot(20*log10((D(:,7))),'.-','Color',[1 0 0],'LineWidth',0.5,'MarkerSize',10)
    hold on
    plot(20*log10((D(:,8))),'.-','Color',[0 0 1],'LineWidth',0.5,'MarkerSize',10)
    grid on
    line([0 160],[0 0],'Color',[0 0 0],'LineWidth',0.5,'LineStyle','-')
    xticks([0,20,40,60,80,100,120,140,160])
    xticklabels({'40','58','75','93','110','93','75','58','40'})
    xlabel('Elicitor Level (dB SPL)','FontSize',14)
    ylabel('\Delta Total (dB)','FontSize',14)
    set(gca,'fontsize',11)
    xlim([0,160])
    ylim([-6 6])
    load('D:\MEMR_Analysis1\MEM06_Run2_Analysis1.mat')
    D = MEMR_mem.D1(:,5:15);
    plot(20*log10((D(:,7))),'.-','Color',[1 0 0],'LineWidth',0.5,'MarkerSize',10)
    plot(20*log10((D(:,8))),'.-','Color',[0 0 1],'LineWidth',0.5,'MarkerSize',10)
    load('D:\MEMR_Analysis1\MEM06_Run3_Analysis1.mat')
    D = MEMR_mem.D1(:,5:15);
    plot(20*log10((D(:,7))),'.-','Color',[1 0 0],'LineWidth',0.5,'MarkerSize',10)
    plot(20*log10((D(:,8))),'.-','Color',[0 0 1],'LineWidth',0.5,'MarkerSize',10)
    load('D:\MEMR_Analysis1\MEM06_Run4_Analysis1.mat')
    D = MEMR_mem.D1(:,5:15);
    plot(20*log10((D(:,7))),'.-','Color',[1 0 0],'LineWidth',0.5,'MarkerSize',10)
    plot(20*log10((D(:,8))),'.-','Color',[0 0 1],'LineWidth',0.5,'MarkerSize',10)
    

  % ---------------------------------------
    figure
    load('C:\myWork\ARLas\MEMR_Analysis_Final_15\MEM06_Run1_Analysis1.mat')
    D = MEMR_mem.D1(:,5:15);
    plot(20*log10(D),'.-','Color',[.7 .7 .7],'LineWidth',0.5,'MarkerSize',10)
    hold on
    plot(20*log10(abs(D(:,7))),'.-','Color',[1 0 0],'LineWidth',0.5,'MarkerSize',10)
    %plot(20*log10(abs(D(:,8))),'.-','Color',[0 0 1],'LineWidth',0.5,'MarkerSize',10)
    plot(20*log10(abs(D(:,5))),'.-','Color',[0 0 1],'LineWidth',0.5,'MarkerSize',10)
    % --> this is a cheat here, because 7 and 8 are right on top if each
    % other. In order to see that they are both plotted, I used freq 5
    % instead.
    grid on
    line([0 160],[0 0],'Color',[0 0 0],'LineWidth',0.5,'LineStyle','-')
    xticks([0,20,40,60,80,100,120,140,160])
    xticklabels({'40','58','75','93','110','93','75','58','40'})
    %levels = 17.5*0.05([0,20,40,60,80,60,40,20,0]')+40
    xlabel('Elicitor Level (dB SPL)','FontSize',14)
    ylabel('\Delta Total (dB)','FontSize',14)
    set(gca,'fontsize',11)
    xlim([0,160])
    ylim([-6 6])


   figure
    load('C:\myWork\ARLas\MEMR_Analysis_Final_15\MEM06_Run1_Analysis1.mat')
    D = MEMR_mem.D(:,5:15);   
    plot(abs(20*log10(abs(D))),'.-','Color',[.7 .7 .7],'LineWidth',0.5,'MarkerSize',10)
    hold on
    plot(abs(20*log10(abs(D(:,7)))),'.-','Color',[1 0 0],'LineWidth',0.5,'MarkerSize',10)
    plot(abs(20*log10(abs(D(:,8)))),'.-','Color',[0 0 1],'LineWidth',0.5,'MarkerSize',10)
    grid on
    line([0 160],[0 0],'Color',[0 0 0],'LineWidth',0.5,'LineStyle','-')
    xticks([0,20,40,60,80,100,120,140,160])
    xticklabels({'40','58','75','93','110','93','75','58','40'})
    %levels = 17.5*0.05([0,20,40,60,80,60,40,20,0]')+40
    xlabel('Elicitor Level (dB SPL)','FontSize',14)
    ylabel('\Delta Magnitude (dB)','FontSize',14)
    set(gca,'fontsize',11)
    xlim([0,160])
    ylim([-6 6])

% for the solution slide (how to do it)
 % ---------------------------------------
    load('C:\myWork\ARLas\MEMR_Analysis_Final_15\MEM06_Run1_Analysis1.mat')
    D = MEMR_mem.D(:,5:15);     
    Theta = linspace(0,2*pi,2000)';
    u = exp(1i*Theta);
    plot(u,'k')
    hold on
    xlim([-1.25,1.25])
    ylim([-1.25,1.25])
    line([-1 1],[0 0],'Color',[0 0 0],'LineWidth',0.5,'LineStyle','-')
    line([0 0],[-1 1],'Color',[0 0 0],'LineWidth',0.5,'LineStyle','-')
    grid on
    set(gca,'fontsize',11)
    xlabel('Real','FontSize',14)
    ylabel('Imaginary','FontSize',14)
    axis square
    plot(D(:,11),'.-','Color',[.7 .7 .7],'LineWidth',0.5,'MarkerSize',8)
    xlim([-.5,.5]+.75)
    ylim([-.5,.5])
    plot(D,'.-','Color',[.7 .7 .7],'LineWidth',0.5,'MarkerSize',8)
    plot(D(:,7),'.-','Color',[1 0 0],'LineWidth',0.5,'MarkerSize',8)
    plot(D(:,8),'.-','Color',[0 0 1],'LineWidth',0.5,'MarkerSize',8)
    ylim([-1.15 1.15])
    xlim([-1.15 1.15])

    load('C:\myWork\ARLas\MEMR_Analysis_Final_15\MEM06_Run1_Analysis1.mat')
    D = MEMR_mem.D(:,5:15);
    D = D - 1;
    figure
    Theta = linspace(0,2*pi,2000)';
    u = exp(1i*Theta);
    plot(u,'k')
    hold on
    xlim([-1.25,1.25])
    ylim([-1.25,1.25])
    line([-1 1],[0 0],'Color',[0 0 0],'LineWidth',0.5,'LineStyle','-')
    line([0 0],[-1 1],'Color',[0 0 0],'LineWidth',0.5,'LineStyle','-')
    grid on
    set(gca,'fontsize',11)
    xlabel('Real','FontSize',14)
    ylabel('Imaginary','FontSize',14)
    axis square
    plot(D(:,11),'.-','Color',[.7 .7 .7],'LineWidth',0.5,'MarkerSize',8)
    xlim([-.5,.5]+.75)
    ylim([-.5,.5])
    plot(D,'.-','Color',[.7 .7 .7],'LineWidth',0.5,'MarkerSize',8)
    plot(D(:,7),'.-','Color',[1 0 0],'LineWidth',0.5,'MarkerSize',8)
    plot(D(:,8),'.-','Color',[0 0 1],'LineWidth',0.5,'MarkerSize',8)
    ylim([-1.15 1.15])
    xlim([-1.15 1.15])




   figure
    load('C:\myWork\ARLas\MEMR_Analysis_Final_15\MEM06_Run1_Analysis1.mat')
    D = MEMR_mem.d1;   
    plot(20*log10(abs(D)),'.-','Color',[0 0 0],'LineWidth',0.5,'MarkerSize',10)
    xticks([0,20,40,60,80,100,120,140,160])
    xticklabels({'40','58','75','93','110','93','75','58','40'})
    %levels = 17.5*0.05([0,20,40,60,80,60,40,20,0]')+40
    xlabel('Elicitor Level (dB SPL)','FontSize',14)
    ylabel('\Delta Total (dB)','FontSize',14)
    set(gca,'fontsize',11)
    xlim([0,160])
    ylim([0 4.5])
    grid on






keyboard
    % ----------------------------------------
