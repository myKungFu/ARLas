function [h] = DW_plotAudiogram(frequency,magnitude,h)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% h = DW_plotAudiogram(frequency,magnitude,h);
%
% function to plot audiogram data on a log2 (octave) axis
% frequency = frequency vector in kHz
% magnitude = corresponding magnitude vector in dB
%
% Author: Shawn Goodman
% Date: March 20, 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fmin = 1000; % minumum frequency (Hz)
fmax = 32000; % maximum frequency (Hz)
niceTicks = [1000, 1250, 1600, 2000, 2500, 3150, 4000, 5000,...
    6300, 8000, 10000, 12500, 16000, 20000, 25000, 32000];
octaveStepSize = 1/2; % fractional octave step size for tick marks
niceTicks = 2.^(log2(fmin):octaveStepSize:log2(fmax))';
niceTicks = niceTicks / 1000; % put into kHz

frequency = frequency(:); % force to column vector
magnitude = magnitude(:);

logX = log2(frequency); % convert x axis to octaves

xlim = [log2(fmin/1000),log2(fmax/1000)]; % set octave frequency limits
xtick = log2(niceTicks);
for ii=1:length(niceTicks)
    decimate = 1; % restrict to this number of digits
    xtickstr(ii) = cellstr(num2str(2.^xtick(ii))); % make tick labels
end

figure(h)
plot(logX,magnitude,'*-r');
% %
% frequency = frequency *1000;
% semilogx(frequency,magnitude,'*-r')
% %
%%%%%%%%%%%%%
hold on
%%%%%%%%%%%%%

ylim([0 100])

set(gca,'XLim',xlim)
set(gca,'XTick',xtick,'XTickLabel',xtickstr);    
xlabel('Frequency (kHz)','FontSize',12);
ylabel('Level (dB)','FontSize',12);
grid on

% %%%%%%%%%%%%%%%%%
d = capNorms();
normFreqs = d(:,1)/1000; % in kHz
normFreqsLog = log2(normFreqs);
figure(h)
plot(normFreqsLog,d(:,2),'b')
plot(normFreqsLog,d(:,3),'b:','LineWidth',2)
% %%%%%%%%%%%%%%%%%
legend('CAP Thresholds','+1 SEM Choongheon Norms','-1 SEM Choongheon Norms')
end

% internal functions ------------------------------------------------------

function [d] = capNorms() % CAP audiograms from all of Choongheon's experiments since being at Wash U
d = [
    22627	97.41792127	76.89116964
    19027	75.85339824	53.65248411
    16000	63.65664236	37.09987938
    13454	57.83444751	32.86990032
    11314	36.85921098	18.96687597
    9513	30.52930643	19.27069357
    8000	33.69814672	24.53663589
    6727	35.37776239	26.55267239
    5657	38.65958891	30.54910674
    4757	39.42099796	30.58769769
    4000	46.53483962	37.49124734
    3364	41.88737492	31.57349464
    2828	42.17228808	31.27988583
    2378	44.79304218	32.92869695
    2000	40.89587072	32.26934668
    1682	45.57039817	34.73394965
    1414	48.37633703	32.91931514
    1189	51.28485967	33.1934012
    1000	51.22312146	31.77687854
    ];
end