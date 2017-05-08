function [] = ARLas_plotDPOAE(d,h,full)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [] = ARLas_plotDPOAE(d,h,full);
%
% Plot analyzed distortion product otoacoustic emissions using data obtained from
% ARLas and the experiment file ARLas_dpoae.m
%
% d = data structore
% h = figure handle
% full = 1 to make full figure, 0 to make subplot (2,1,2)
%
% Auditory Research Lab, The University of Iowa
% Deptartment of Communication Sciences & Disorders
% The University of Iowa
% Author: Shawn S. Goodman, PhD
% Date: April 5, 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% plot dp-gram
figure(h)
if full == 0
    subplot(2,1,2)
end
hold off
fmin = min(d.f2) - 10; % minumum frequency (Hz)
fmax = max(d.f2) + 10; % maximum frequency (Hz)
fff = round(linspace(fmin,fmax,10))';
fff = fff /1000; % put into kHz
logX = log2(fff); % convert x axis to octaves
xlim = [log2(fmin/1000),log2(fmax/1000)]; % set octave frequency limits
xtick = log2(d.f2/1000);
for ii=1:length(xtick)
    decimate = 3; % restrict to this number of digits
    xtickstr(ii) = cellstr(num2str(2.^xtick(ii),decimate)); % make tick labels
end
plot(log2(d.f2/1000),d.Ldp_cubic,'r-*')
hold on
plot(log2(d.f2/1000),d.Ldp_diff,'b-*')
plot(log2(d.f2/1000),d.Ndp_cubic,'--','Color',[0 0 0])
plot(log2(d.f2/1000),d.Ndp_diff,':','Color',[0 0 0])
plot(log2(d.f2/1000),d.L1,'c-v')
plot(log2(d.f2/1000),d.L2,'c-^')
if full == 1
    legend('2f1-f2','f2-f1','nf cubic','nf diff','L1','L2')
end
set(gca,'XLim',xlim)
set(gca,'XTick',xtick,'XTickLabel',xtickstr);    
xlabel('Frequency (kHz)','FontSize',12)
ylabel('Magnitude (dB SPL)','FontSize',12)
if full == 1
    title(['DP-GRAM:  ',d.subjName],'FontSize',12)    
else
    title('DP-GRAM','FontSize',12)    
end

