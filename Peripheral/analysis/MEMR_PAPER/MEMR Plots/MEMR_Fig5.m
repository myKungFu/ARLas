function [] = MEMR_Fig5(whichOne)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MEMR_Fig5(whichOne)
%
% series of figures to show results of specific measurements
% and variability related to ICC.
%
% Input Argument:
%   whichOne = 1, 2, 3, or 4
%              (for max total change, delay, thresholds, hysteresis)
% A common metric for quantifying variability is the ICC. The ICC is a 
% descriptive statistic that can be used to describe how strongly
% strongly measurements from the same individual resemble each other, 
% relative to measurements from other individuals. We computed ICC
% as an evaluation of test–retest reliability, based on a two-way
% mixed-effects model and absolute agreement between measurements 
% (Koo et al., 2016), denoted ICC(2,1) by Shrout and Fleiss [21]. 
% We used the first two measurements from each participant for computation. 
% The ICC was 0.80 (95% CI 0.64–0.89)for XXX. ICC is here computed using
% the function icc21.m, located in Analysis/SJ/ICC/.
%
% Author: Shawn Goodman
% DAte: February 2, 2024
% Last Updated: February 4, 2024 -- for ARO talk
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load('D:\MEMR_GroupData_Fix1\MEMR_groupData.mat')

FS = 18; % font size
TS = 16; % tick size

if nargin == 0
    whichOne = 1;
end

if whichOne == 1 % ------ Max Total Change ------------
    q = 20*log10(pAmp);
    myICC = ICC(3,'k',q); % this one is by Kevin Brownhill
    %q2 = q(:,1:2);
    %[r,varargout] = icc21(q);
    %ci(1,1) = varargout(1);
    %ci(1,2) = varargout(2);


    q1 = q(:);
    PD = fitdist(q1,'kernel');
    chain = PD.random(1000000,1);
    post = calculatePosterior_v2(chain);
    medq = median(q,2);
    [~,indx] = sort(medq);
    qs = q(indx,:);
    ss = (1:1:30)';
    
    qqq = median(qs,2);
    sm = 0.1; % was 0.0001; % was 0.00001;  smoothing factor (smaller numbers are more smooth)
    www = ones(size(ss));
    pp = csaps(ss,qqq,sm,[],www); % piecewise polynomial real coefficients
    qs_sm = ppval(pp,ss); % evaluate only at the original x values
       
    figure(1)
    for ii=1:30
    plot(ss(ii),qs(ii,:),'.','Color',[.7 .7 .7],'MarkerSize',14)
    hold on
    end
    for ii=1:30
    plot([ss(ii),ss(ii)],[min(qs(ii,:)),max(qs(ii,:))],'-','Color',[.7 .7 .7],'LineWidth',1)
    end
    plot(ss,qs_sm,'Color',[.5 0 0],'LineWidth',2)
    %line([0,30],[0 0],'Color',[0 0 0],'LineWidth',0.5,'LineStyle','--')
    grid off
    xmin = 0;
    xmax = 31;
    ymin = 0;
    ymax = 12
    ylim([ymin,ymax])
    xlim([xmin,xmax])
    set(gca,'fontsize',TS)
    xlabel('Participants (sorted)','FontSize',FS)
    ylabel('Max Total Change (dB)','FontSize',FS)
    
    figure(2)
    plot(post.xx,post.yy,'Color',[.5 0 0],'LineWidth',2)
    ymin = 0;
    ymax = max(post.yy);
    ymax = ymax * 1.1;
    xmin = min(post.xx);
    xmax = max(post.xx);
    hold on
    line([post.hdi95(1),post.hdi95(1)],[ymin ymax],'Color',[0 0 0],'LineWidth',0.5)
    line([post.hdi95(2),post.hdi95(2)],[ymin ymax],'Color',[0 0 0],'LineWidth',0.5)
    line([post.mod,post.mod],[ymin ymax],'Color',[0.5 0 0],'LineWidth',0.5,'LineStyle','--')
    grid on
    ylim([ymin,ymax])
    xlim([xmin,xmax])
    set(gca,'fontsize',TS)
    xlabel('Max Total Change (dB)','FontSize',FS)
    ylabel('Probability Density','FontSize',FS)
    
elseif whichOne == 2 % ------ Reflex Delay ------------
    q = 1000*delay; % put delay in ms
    myICC = ICC(3,'k',q); % this one is by Kevin Brownhill
    
    q1 = q(:);
    PD = fitdist(q1,'kernel');
    chain = PD.random(1000000,1);
    post = calculatePosterior_v2(chain);
    medq = median(q,2);
    [~,indx] = sort(medq);
    qs = q(indx,:);
    ss = (1:1:30)';
    
    qqq = median(qs,2);
    sm = 0.1; % was 0.0001; % was 0.00001;  smoothing factor (smaller numbers are more smooth)
    www = ones(size(ss));
    pp = csaps(ss,qqq,sm,[],www); % piecewise polynomial real coefficients
    qs_sm = ppval(pp,ss); % evaluate only at the original x values
       
    figure(1)
    for ii=1:30
    plot(ss(ii),qs(ii,:),'.','Color',[.7 .7 .7],'MarkerSize',14)
    hold on
    end
    for ii=1:30
    plot([ss(ii),ss(ii)],[min(qs(ii,:)),max(qs(ii,:))],'-','Color',[.7 .7 .7],'LineWidth',1)
    end
    plot(ss,qs_sm,'Color',[.5 0 0],'LineWidth',2)
    line([0,31],[0 0],'Color',[0 0 0],'LineWidth',0.5,'LineStyle','--')
    grid off
    xmin = 0;
    xmax = 31;
    ymin = -250; %min(post.xx);
    ymax = 500; %max(post.xx);
    ylim([ymin,ymax])
    xlim([xmin,xmax])
    set(gca,'fontsize',TS)
    yticks([-500,-250,0,250,500])
    yticklabels({'-500','-250','0','250','500'})
    xlabel('Participants (sorted)','FontSize',FS)
    ylabel('Reflex Delay (ms)','FontSize',FS)
    
    figure(2)
    plot(post.xx,post.yy,'Color',[.5 0 0],'LineWidth',2)
    ymin = 0;
    ymax = max(post.yy);
    ymax = ymax * 1.1;
    xmin = -200; %min(post.xx);
    xmax = 600; %max(post.xx);
    hold on
    line([post.hdi95(1),post.hdi95(1)],[ymin ymax],'Color',[0 0 0],'LineWidth',0.5)
    line([post.hdi95(2),post.hdi95(2)],[ymin ymax],'Color',[0 0 0],'LineWidth',0.5)
    line([post.mod,post.mod],[ymin ymax],'Color',[0.5 0 0],'LineWidth',0.5,'LineStyle','--')
    grid on
    ylim([ymin,ymax])
    xlim([xmin,xmax])
    set(gca,'fontsize',TS)
    xlabel('Reflex Delay (ms)','FontSize',FS)
    ylabel('Probability Density','FontSize',FS)
    
elseif whichOne == 3 % ------ Thresholds ------------
    q = thdOnsetLvl;
    q2 = thdOffsetLvl;
    myICC = ICC(3,'k',q); % this one is by Kevin Brownhill
    myICC2 = ICC(3,'k',q2);

    % for onset threshold --------------
    q1 = q(:); 
    PD = fitdist(q1,'kernel');
    chain = PD.random(1000000,1);
    post = calculatePosterior_v2(chain);
    medq = median(q,2);
    [~,indx] = sort(medq);
    qs = q(indx,:);
    ss = (1:1:30)';
    
    qqq = median(qs,2);
    sm = 0.1; % was 0.0001; % was 0.00001;  smoothing factor (smaller numbers are more smooth)
    www = ones(size(ss));
    pp = csaps(ss,qqq,sm,[],www); % piecewise polynomial real coefficients
    qs_sm = ppval(pp,ss); % evaluate only at the original x values
       
    figure(1)
    for ii=1:30
    plot(ss(ii),qs(ii,:),'.','Color',[1 .7 .7],'MarkerSize',14)
    hold on
    end
    for ii=1:30
    plot([ss(ii),ss(ii)],[min(qs(ii,:)),max(qs(ii,:))],'-','Color',[1 .7 .7],'LineWidth',1)
    end
    plot(ss,qs_sm,'Color',[.5 0 0],'LineWidth',2)

    % for offst threshold -------------
    q1 = q2(:); % for onset threshold
    PD2 = fitdist(q1,'kernel');
    chain = PD2.random(1000000,1);
    post2 = calculatePosterior_v2(chain);
    %medq = median(q2,2);
 %[~,indx] = sort(medq);
    qs = q2(indx,:);
    ss = (1:1:30)';
    
    qqq = median(qs,2);
    sm = 0.1; % was 0.0001; % was 0.00001;  smoothing factor (smaller numbers are more smooth)
    www = ones(size(ss));
    pp = csaps(ss,qqq,sm,[],www); % piecewise polynomial real coefficients
    qs_sm = ppval(pp,ss); % evaluate only at the original x values
       
    figure(1)
    hold on
    for ii=1:30
        plot(ss(ii),qs(ii,:),'.','Color',[.7 .7 1],'MarkerSize',14)
    end
    for ii=1:30
    plot([ss(ii),ss(ii)],[min(qs(ii,:)),max(qs(ii,:))],'-','Color',[.7 .7 1],'LineWidth',1)
    end
    plot(ss,qs_sm,'Color',[0 0 .5],'LineWidth',2)
    
    grid off
    xmin = 0;
    xmax = 31;
    %ymin = min(post.xx);
    %ymax = max(post.xx);
    ymax = 105;
    ymin = 55;
    ylim([ymin,ymax])
    xlim([xmin,xmax])
    set(gca,'fontsize',TS)
    xlabel('Participants (sorted *)','FontSize',FS)
    ylabel('Threshold (dB SPL)','FontSize',FS)
    
    figure(2)
    plot(post.xx,post.yy,'Color',[.5 0 0],'LineWidth',2)
    hold on
    plot(post2.xx,post2.yy,'Color',[0 0 .5],'LineWidth',2)
    ymin = 0;
    ymax = max(post.yy);
    ymax = ymax * 1.1;
    xmin = 40;
    xmax = 110;
    hold on
    line([post.hdi95(1),post.hdi95(1)],[ymin ymax],'Color',[0 0 0],'LineWidth',0.5)
    line([post.hdi95(2),post.hdi95(2)],[ymin ymax],'Color',[0 0 0],'LineWidth',0.5)
    line([post.mod,post.mod],[ymin ymax],'Color',[0.5 0 0],'LineWidth',0.5,'LineStyle','--')
    line([post2.hdi95(1),post2.hdi95(1)],[ymin ymax],'Color',[0 0 0],'LineWidth',0.5)
    line([post2.hdi95(2),post2.hdi95(2)],[ymin ymax],'Color',[0 0 0],'LineWidth',0.5)
    line([post2.mod,post2.mod],[ymin ymax],'Color',[0 0 .5],'LineWidth',0.5,'LineStyle','--')
    grid on
    ylim([ymin,ymax])
    xlim([xmin,xmax])
    set(gca,'fontsize',TS)
    xlabel('Threshold (dB SPL)','FontSize',FS)
    ylabel('Probability Density','FontSize',FS)
    
elseif whichOne == 4 % ------ Hysteresis ------------
    %q = 20*log10(hysteresis);
    q = thdOnsetLvl;
    q2 = thdOffsetLvl;
    %q = q2./q;
    q = q2 - q;
    myICC = ICC(3,'k',q); % this one is by Kevin Brownhill    
    
    q1 = q(:);
    PD = fitdist(q1,'kernel');
    chain = PD.random(1000000,1);
    post = calculatePosterior_v2(chain);
    medq = median(q,2);
    [~,indx] = sort(medq);
    qs = q(indx,:);
    ss = (1:1:30)';
    
    qqq = median(qs,2);
    sm = 0.1; % was 0.0001; % was 0.00001;  smoothing factor (smaller numbers are more smooth)
    www = ones(size(ss));
    pp = csaps(ss,qqq,sm,[],www); % piecewise polynomial real coefficients
    qs_sm = ppval(pp,ss); % evaluate only at the original x values
       
    figure(1)
    for ii=1:30
    plot(ss(ii),qs(ii,:),'.','Color',[.7 .7 .7],'MarkerSize',14)
    hold on
    end
    for ii=1:30
    plot([ss(ii),ss(ii)],[min(qs(ii,:)),max(qs(ii,:))],'-','Color',[.7 .7 .7],'LineWidth',1)
    end
    plot(ss,qs_sm,'Color',[.5 0 0],'LineWidth',2)
    line([0,31],[0 0],'Color',[0 0 0],'LineWidth',0.5,'LineStyle','--')
    grid off
    xmin = 0;
    xmax = 31;
    ymin = -30; %min(post.xx);
    ymax = 10; %max(post.xx);
    ylim([ymin,ymax])
    xlim([xmin,xmax])
    set(gca,'fontsize',TS)
    xlabel('Participants (sorted)','FontSize',FS)
    ylabel('Hysteresis','FontSize',FS)
    
    figure(2)
    plot(post.xx,post.yy,'Color',[.5 0 0],'LineWidth',2)
    ymin = 0;
    ymax = max(post.yy);
    ymax = ymax * 1.1;
    xmin = min(post.xx);
    xmax = max(post.xx);
    hold on
    line([post.hdi95(1),post.hdi95(1)],[ymin ymax],'Color',[0 0 0],'LineWidth',0.5)
    line([post.hdi95(2),post.hdi95(2)],[ymin ymax],'Color',[0 0 0],'LineWidth',0.5)
    line([post.mod,post.mod],[ymin ymax],'Color',[0.5 0 0],'LineWidth',0.5,'LineStyle','--')
    grid on
    ylim([ymin,ymax])
    xlim([xmin,xmax])
    set(gca,'fontsize',TS)
    xlabel('Hysteresis','FontSize',FS)
    ylabel('Probability Density','FontSize',FS)
elseif whichOne == 5     
    q = trend;
    q = 20*log10(q);
    myICC = ICC(3,'k',q); % this one is by Kevin Brownhill    
    
    q1 = q(:);
    PD = fitdist(q1,'kernel');
    chain = PD.random(1000000,1);
    post = calculatePosterior_v2(chain);
    medq = median(q,2);
    [~,indx] = sort(medq);
    qs = q(indx,:);
    ss = (1:1:30)';
    
    qqq = median(qs,2);
    sm = 0.1; % was 0.0001; % was 0.00001;  smoothing factor (smaller numbers are more smooth)
    www = ones(size(ss));
    pp = csaps(ss,qqq,sm,[],www); % piecewise polynomial real coefficients
    qs_sm = ppval(pp,ss); % evaluate only at the original x values
       
    figure(1)
    for ii=1:30
    plot(ss(ii),qs(ii,:),'.','Color',[.7 .7 .7],'MarkerSize',14)
    hold on
    end
    for ii=1:30
    plot([ss(ii),ss(ii)],[min(qs(ii,:)),max(qs(ii,:))],'-','Color',[.7 .7 .7],'LineWidth',1)
    end
    plot(ss,qs_sm,'Color',[.5 0 0],'LineWidth',2)
    %line([0,31],[0 0],'Color',[0 0 0],'LineWidth',0.5,'LineStyle','--')
    grid off
    xmin = 0;
    xmax = 31;
    ymin = 0; %min(post.xx);
    ymax = 12; %max(post.xx);
    ylim([ymin,ymax])
    xlim([xmin,xmax])
    set(gca,'fontsize',TS)
    xlabel('Participants (sorted)','FontSize',FS)
    ylabel('Trend','FontSize',FS)
    
    figure(2)
    plot(post.xx,post.yy,'Color',[.5 0 0],'LineWidth',2)
    ymin = 0;
    ymax = max(post.yy);
    ymax = ymax * 1.1;
    xmin = -0.5; %min(post.xx);
    xmax = 12; %max(post.xx);
    hold on
    line([post.hdi95(1),post.hdi95(1)],[ymin ymax],'Color',[0 0 0],'LineWidth',0.5)
    line([post.hdi95(2),post.hdi95(2)],[ymin ymax],'Color',[0 0 0],'LineWidth',0.5)
    line([post.mod,post.mod],[ymin ymax],'Color',[0.5 0 0],'LineWidth',0.5,'LineStyle','--')
    grid on
    ylim([ymin,ymax])
    xlim([xmin,xmax])
    set(gca,'fontsize',TS)
    xlabel('Trend','FontSize',FS)
    ylabel('Probability Density','FontSize',FS)
end

post.mod
post.hdi95
myICC
if whichOne == 3
    post2.mod
    post2.hdi95
    myICC2
end

keyboard

pc = getSpeechPC;
load('C:\MEMR_GroupData\MEMR_groupData.mat')
pAmp = 20*log10(pAmp);
delay = 1000*delay;
hysteresis = thdOffsetLvl - thdOnsetLvl;

figure
subplot(3,3,1)
plot(pAmp,delay,'.')
xlabel('pAmp')
ylabel('delay')
subplot(3,3,2)
plot(pAmp,hysteresis,'.')
xlabel('pAmp')
ylabel('hysteresis')
subplot(3,3,3)
plot(pAmp,thdOnsetLvl,'.')
xlabel('pAmp')
ylabel('Onset')
subplot(3,3,4)
plot(pAmp,thdOffsetLvl,'.')
xlabel('pAmp')
ylabel('Offset')
subplot(3,3,5)
plot(delay,hysteresis,'.')
xlabel('delay')
ylabel('hysteresis')
subplot(3,3,6)
plot(delay,thdOnsetLvl,'.')
xlabel('delay')
ylabel('Onset')
subplot(3,3,7)
plot(delay,thdOffsetLvl,'.')
xlabel('delay')
ylabel('Offset')


r = corrcoef(mean(pAmp,2),pc); % r = -0.1
r_pAmp = r(1,2);
R = corrcoef(mean(delay,2),pc); % -0.02
r_delay = r(1,2);
r = corrcoef(mean(hysteresis,2),pc); % r = 0.1796
r_hyst = r(1,2);
r = corrcoef(mean(thdOnsetLvl,2),pc); % r = -0.5337
r_thdOn = r(1,2);
r = corrcoef(mean(thdOffsetLvl,2),pc); % r = -0.4064
r_thdOff = r(1,2);

% --------------------------------------------------------
FS = 14;
TS = 12;
figure
subplot(2,2,1)
 plot(mean(pAmp,2),pc,'.','MarkerSize',12)
    % Linear model Poly1:
    %      f(x) = p1*x + p2
    % Coefficients (with 95% confidence bounds):
    %        p1 =       2.293  (-1.646, 6.231)
    %        p2 =       73.68  (63.32, 84.04)
    % 
    % Goodness of fit:
    %   SSE: 4245
    %   R-square: 0.05019
    %   Adjusted R-square: 0.01501
    %   RMSE: 12.54
 p1 =  2.293;
 p2 = 73.68;
 xx = linspace(0,10,100)';
 fx = p1 .* xx + p2;
 hold on
 plot(xx,fx,'k','LineWidth',1)
set(gca,'fontsize',TS)
xlabel('Max Total Change (dB)','FontSize',FS)
ylabel('QuickSin % Correct','FontSize',FS)

subplot(2,2,2)
 plot(mean(thdOnsetLvl,2),pc,'.','MarkerSize',12)
    % Linear model Poly1:
    %      f(x) = p1*x + p2
    % Coefficients (with 95% confidence bounds):
    %        p1 =     -0.8666  (-1.444, -0.2895)
    %        p2 =       153.5  (103.2, 203.7)
    % 
    % Goodness of fit:
    %   SSE: 4143
    %   R-square: 0.1837
    %   Adjusted R-square: 0.1545
    %   RMSE: 12.16   
 p1 =  -0.8666;
 p2 = 153.5;
 xx = linspace(50,110,100)';
 fx = p1 .* xx + p2;
 hold on
 plot(xx,fx,'k','LineWidth',1)
set(gca,'fontsize',TS)
xlabel('Onset Threshold (dB SPL)','FontSize',FS)
ylabel('QuickSin % Correct','FontSize',FS)

subplot(2,2,3)
 plot(mean(thdOffsetLvl,2),pc,'.','MarkerSize',12)
    % Linear model Poly1:
    %      f(x) = p1*x + p2
    % Coefficients (with 95% confidence bounds):
    %        p1 =      -0.577  (-1.172, 0.01779)
    %        p2 =       118.8  (75.47, 162.1)
    % 
    % Goodness of fit:
    %   SSE: 4298
    %   R-square: 0.153
    %   Adjusted R-square: 0.1228
    %   RMSE: 12.39
 p1 =  -0.577;
 p2 = 118.8;
 xx = linspace(50,110,100)';
 fx = p1 .* xx + p2;
 hold on
 plot(xx,fx,'k','LineWidth',1)
set(gca,'fontsize',TS)
xlabel('Offset Threshold (dB SPL)','FontSize',FS)
ylabel('QuickSin % Correct','FontSize',FS)

subplot(2,2,4)
 plot(mean(hysteresis,2),pc,'.','MarkerSize',12)
    % Linear model Poly1:
    %      f(x) = p1*x + p2
    % Coefficients (with 95% confidence bounds):
    %        p1 =    -0.05136  (-1.128, 1.025)
    %        p2 =       76.68  (59.87, 93.49)
    % 
    % Goodness of fit:
    %   SSE: 4582
    %   R-square: 0.000355
    %   Adjusted R-square: -0.03667
    %   RMSE: 13.03  
 p1 =  -0.05136;
 p2 =  76.68;
 xx = linspace(-30,10,100)';
 fx = p1 .* xx + p2;
 hold on
 plot(xx,fx,'k','LineWidth',1)
set(gca,'fontsize',TS)
xlabel('Hysteresis (dB)','FontSize',FS)
ylabel('QuickSin % Correct','FontSize',FS)














subplot(3,2,2)
 plot(mean(delay,2),pc,'.')
 xlabel('delay')
 ylabel('PC')
% Linear model Poly1:
%      f(x) = p1*x + p2
% Coefficients (with 95% confidence bounds):
%        p1 =     0.07066  (0.007656, 0.1337)
%        p2 =       60.87  (46.2, 75.53)
% 
% Goodness of fit:
%   SSE: 4128
%   R-square: 0.09925
%   Adjusted R-square: 0.06589
%   RMSE: 12.37
 
subplot(3,2,3)
 
% -->> here!
 plot(mean(thdOnsetLvl,2),pc,'.','MarkerSize',14)
 xlabel('Onset')
 ylabel('PC')
 p1 =  -0.8666;
 p2 = 153.5;
 xx = linspace(70,105,100)';
 fx = p1 .* xx + p2;
 hold on
 plot(xx,fx,'k','LineWidth',2)
 xlim([65,105])
 ylim([50,105])
set(gca,'fontsize',TS)
    xlabel('Onset Threshold (dB SPL)','FontSize',FS)
    ylabel('QuickSin % Correct','FontSize',FS)
% Linear model Poly1:
%      f(x) = p1*x + p2
% Coefficients (with 95% confidence bounds):
%        p1 =     -0.8666  (-1.444, -0.2895)
%        p2 =       153.5  (103.2, 203.7)
% 
% Goodness of fit:
%   SSE: 4143
%   R-square: 0.1837
%   Adjusted R-square: 0.1545
%   RMSE: 12.16

subplot(3,2,4)
 plot(mean(thdOffsetLvl,2),pc,'.')
 xlabel('Offset')
 ylabel('PC')
% Linear model Poly1:
%      f(x) = p1*x + p2
% Coefficients (with 95% confidence bounds):
%        p1 =      -0.577  (-1.172, 0.01779)
%        p2 =       118.8  (75.47, 162.1)
% 
% Goodness of fit:
%   SSE: 4298
%   R-square: 0.153
%   Adjusted R-square: 0.1228
%   RMSE: 12.39




subplot(3,2,5)
 plot(mean(hysteresis,2),pc,'.')
 xlabel('hysteresis')
 ylabel('PC')
% Linear model Poly1:
%      f(x) = p1*x + p2
% Coefficients (with 95% confidence bounds):
%        p1 =    -0.05136  (-1.128, 1.025)
%        p2 =       76.68  (59.87, 93.49)
% 
% Goodness of fit:
%   SSE: 4582
%   R-square: 0.000355
%   Adjusted R-square: -0.03667
%   RMSE: 13.03






end

% INTERNAL FUNCTIONS ------------------------------------------------------
function [pc] = getSpeechPC()
    % get speech (QuickSin) percent correct (average of 2 tests)
    pc = [86.0
92.0
100.0
68.0
84.0
100.0
80.0
72.0
82.0
86.0
82.0
90.0
100.0
52.0
64.0
76.0
68.0
54.0
72.0
58.0
74.0
98.0
64.0
86.0
64.0
76.0
80.0
74.0
78.0
86.0
];
end
function [post] = calculatePosterior_v2(chain)
    % calculate the posterior distribution using a kernel method
    dist = fitdist(chain,'kernel');
    minx = min(chain);
    maxx = max(chain);
    N = 1000;
    xx = linspace(minx,maxx,N);
    yy = dist.pdf(xx);
    alpha = 0.1; % this was 0.05 (95%). Here replaced with 0.1 (90%)
    [hdi95,HDIdensity,inOut95] = getHDI(alpha,xx,yy);
    alpha = 0.5;
    [hdi50,HDIdensity,inOut50] = getHDI(alpha,xx,yy);
    med = dist.median;
    [~,medIndx] = min(abs(med-xx));
    medY = yy(medIndx);
    modY = max(yy);
    indxMod = find(yy==modY);
    mod = xx(indxMod);
    post.dist = dist; % distribution object
    post.hdi95 = hdi95; % hdi interval for 95%
    post.hdi50 = hdi50; % hid interval for 50%
    post.med = med; % median of distribution
    post.medY = medY; % height of median of distribution
    post.xx = xx; % x-axis of distribution
    post.yy = yy; % y-axis of distribution
    post.mod = mod; % mode of distribution
    post.modY = modY; % height of mode of distribution
    post.inOut95 = inOut95;
    post.inOut50 = inOut50;
end
function [HDI,HDIdensity,inOut] = getHDI(alpha,xx,yy)
    % Calculate highest density interval on alpha
    [p,I] = sort(yy,'ascend'); % sort densities in ascending order
    cutValue = sum(p)*alpha; % find the alpha% cut value
    cutIndx = min(find(cumsum(p)>cutValue)); % find the location of the cut value
    waterline = p(cutIndx); % this is the cutoff density
    [goodValues,goodIndx] = find(yy >= waterline); % locate all values > cut
    HDI = [xx(min(goodIndx));xx(max(goodIndx))]; % determine the interval
    HDIdensity = waterline;
    inOut = (yy >= waterline);
end

% OLD CODE ----------------------------------------------------------------

% grid on
% yyy = PD.icdf(xxx);
% hold on
% plot(xxx*30,yyy,'r', 'LineWidth', 1.25)
% 
% d2 = pAmp;
% % d2(6,:) = [];
% d2 = d2(:);
% PD = fitdist(d2,'kernel');
% figure (02)
% xxx = linspace(-.5,1,1000)';
% dd = median(pAmp,2);
% [~,indx] = sort(dd);
% hold on 
% grid on
% plot(pAmp(indx,:),'.-','Color',[.7 .7 .7])
% grid on
% yyy = PD.icdf(xxx);
% hold on
% plot(xxx*30,yyy,'b', 'LineWidth', 1.25)
% 
% d2 = hysteresis;
% % d2(6,:) = [];
% d2 = d2(:);
% PD = fitdist(d2,'kernel');
% figure (03)
% xxx = linspace(-.5,1,1000)';
% dd = median(hysteresis,2);
% [~,indx] = sort(dd);
% hold on 
% grid on
% plot(hysteresis(indx,:),'.-','Color',[.7 .7 .7])
% grid on
% yyy = PD.icdf(xxx);
% hold on
% plot(xxx*30,yyy,'g', 'LineWidth', 1.25)

% figure(6)
% plot(ss,qs,'.-','Color',[.7 .7 .7])
% hold on
% plot(ss,qs_sm,'Color',[.5 0 0],'LineWidth',2)
% plot(ss,qqq,'.-','Color',[.5 0 0],'LineWidth',2)
% %plot(xs,ICDF)
% xxx = linspace(0,1,100)';
% xxx(end) = .999;
% xxxs = xxx*30;
% yyy = PD.icdf(xxx);
% plot(xxxs,yyy)
% 
%qqq = meanSmoother(qqq,5);


%plot(x,PDF,'Color',[.5 0 0],'LineWidth',2)
%hold on
%grid on
%plot(post.mod,post.modY,'b*')

%xlabel('\theta','FontSize',12)
%ylabel('Probability Density','FontSize',12)
%title(['Mode = ',num2str(post.mod),' HDI = ',num2str(post.hdi95(1)),' to ',num2str(post.hdi95(2))])

%x = linspace(xmin,xmax,1000)';
%PDF = PD.pdf(x);
%ICDF = PD.icdf(x);
%xs = (x*29)+1;
