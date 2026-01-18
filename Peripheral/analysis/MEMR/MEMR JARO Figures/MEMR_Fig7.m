 figure
    load('C:\myWork\ARLas\MEMR_Analysis_Final_15\MEM09_Run1_Analysis1.mat')
    D = MEMR_mem.D(:,5:15);
    % create a line segment for u1
    Theta = linspace(0,2*pi,2000)';
    u = exp(1i*Theta);
    plot(u,'k')
    hold on
    line([-1 1],[0 0],'Color',[0 0 0],'LineWidth',0.5,'LineStyle','-')
    line([0 0],[-1 1],'Color',[0 0 0],'LineWidth',0.5,'LineStyle','-')
    grid on
    set(gca,'fontsize',12)
    xlabel('Real','FontSize',14)
    ylabel('Imaginary','FontSize',14)
    axis square

    %plot(D,'.-','Color',[.7 .7 .7],'LineWidth',0.5,'MarkerSize',8)
    plot(D(:,7),'.-','Color',[1 0 0],'LineWidth',0.5,'MarkerSize',8)
    plot(D(:,8),'.-','Color',[0 0 1],'LineWidth',0.5,'MarkerSize',8)
    xlim([.4 1.4])
    ylim([-.7,.3])

    load('C:\myWork\ARLas\MEMR_Analysis_Final_15\MEM09_Run2_Analysis1.mat')
    D = MEMR_mem.D(:,5:15);
    plot(D(:,7),'.-','Color',[1 0 0],'LineWidth',0.5,'MarkerSize',8)
    plot(D(:,8),'.-','Color',[0 0 1],'LineWidth',0.5,'MarkerSize',8)
    load('C:\myWork\ARLas\MEMR_Analysis_Final_15\MEM09_Run3_Analysis1.mat')
    D = MEMR_mem.D(:,5:15);
    plot(D(:,7),'.-','Color',[1 0 0],'LineWidth',0.5,'MarkerSize',8)
    plot(D(:,8),'.-','Color',[0 0 1],'LineWidth',0.5,'MarkerSize',8)
    load('C:\myWork\ARLas\MEMR_Analysis_Final_15\MEM09_Run4_Analysis1.mat')
    D = MEMR_mem.D(:,5:15);
    plot(D(:,7),'.-','Color',[1 0 0],'LineWidth',0.5,'MarkerSize',8)
    plot(D(:,8),'.-','Color',[0 0 1],'LineWidth',0.5,'MarkerSize',8)