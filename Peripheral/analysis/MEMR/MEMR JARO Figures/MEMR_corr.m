```matlab
function [] = MEMR_corr(pickme)

load('C:\myWork\ARLas\MEMR_Data_Paper\MEMR_groupDataM.mat')

pAmp = 20 * log10(pAmp);
A = mean(pAmp, 2);
B = mean(thdOnsetLvl, 2);
C = mean(thdOffsetLvl, 2);
D = mean(hysteresis, 2);
E = mean(thdOnsetTime, 2);
F = mean(thdOffsetTime, 2);
G = mean(delay, 2);
H = mean(slopeUp, 2);

FS = 14; % font size
TS = 12; % tick size

% quick percent-correct scores (out of 30)
scores = [21.5, 23, 25, 17, 21, 26.5, 20, 18, 20.5, 21.5, ...
          20.5, 22.5, 25.5, 13, 16, 19, 17, 13.5, 18, 14.5, ...
          18.5, 24.5, 16, 21.5, 16, 19, 20, 18.5, 19.5, 21.5]';

if nargin == 0
    whichOne = 1;
end

if pickme == 1
    metric = A;

    percentage_correct = (scores / 30) * 100;
    mdl = fitlm(metric, percentage_correct, 'RobustOpts', 'bisquare');
    R2 = mdl.Rsquared.Ordinary;
    R = sqrt(R2);

    x_fit = linspace(min(metric), max(metric), 100)';
    [y_fit, y_CI] = predict(mdl, x_fit);

    figure('Position', [1000, 500, 500, 450]);
    set(gca, 'DataAspectRatio', [1 1 1]);

    scatter(metric, percentage_correct, 'filled', 'MarkerEdgeColor', 'k');
    hold on;
    plot(x_fit, y_fit, 'k-', 'LineWidth', 1.5); % best-fit line
    fill([x_fit; flipud(x_fit)], [y_CI(:,1); flipud(y_CI(:,2))], ...
         'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none'); % shaded CI band

    xlabel('Max Total Change (dB)', 'FontSize', FS);
    ylabel('Percentage Correct (%)', 'FontSize', FS);
    title(['R = ', num2str(R, '%.2f'), ' | R^2 = ', num2str(R2, '%.2f')], 'FontSize', FS);

    yticks([30,40,50,60,70,80,90,100])
    yticklabels({'30','40','50','60','70','80','90','100'})
    xticks([0,5,10,15,20,25,30])
    xticklabels({'0','5','10','15','20','25','30'})
    set(gca, 'FontSize', TS);
    grid off;

    CI = coefCI(mdl);
    ci_upper = CI(2, 2);
    ci_lower = CI(2, 1);

    if ci_lower <= 0 && ci_upper >= 0
        sig_text = 'Slope is NOT significantly different from zero.';
    else
        sig_text = 'Slope is significantly different from zero.';
    end
    text(mean(metric), max(percentage_correct), sig_text, 'FontSize', TS);

    % print stats so you can copy/paste if needed
    disp(['R: ', num2str(R)]);
    disp(['R^2: ', num2str(R2)]);
    disp(['Confidence intervals for the slope: [', num2str(ci_lower), ', ', num2str(ci_upper), ']']);
    disp(sig_text);

    hold off;

elseif pickme == 2

    metric = B;

    percentage_correct = (scores / 30) * 100;
    mdl = fitlm(metric, percentage_correct, 'RobustOpts', 'bisquare');
    R2 = mdl.Rsquared.Ordinary;
    R = sqrt(R2);

    x_fit = linspace(min(metric), max(metric), 100)';
    [y_fit, y_CI] = predict(mdl, x_fit);

    figure('Position', [1000, 500, 500, 450]);
    set(gca, 'DataAspectRatio', [1 1 1]);

    scatter(metric, percentage_correct, 'filled', 'MarkerEdgeColor', 'k');
    hold on;

    plot(x_fit, y_fit, 'k-', 'LineWidth', 1.5); % best-fit line

    fill([x_fit; flipud(x_fit)], [y_CI(:,1); flipud(y_CI(:,2))], ...
         'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none'); % shaded CI band

    xlabel('Onset Threshold (dB SPL)', 'FontSize', FS);
    ylabel('Percentage Correct (%)', 'FontSize', FS);
    title(['R = ', num2str(R, '%.2f'), ' | R^2 = ', num2str(R2, '%.2f')], 'FontSize', FS);

    set(gca, 'FontSize', TS);
    grid off;

    CI = coefCI(mdl);
    ci_upper = CI(2, 2);
    ci_lower = CI(2, 1);

    if ci_lower <= 0 && ci_upper >= 0
        sig_text = 'Slope is NOT significantly different from zero.';
    else
        sig_text = 'Slope is significantly different from zero.';
    end
    text(mean(metric), max(percentage_correct), sig_text, 'FontSize', TS);

    % print stats so you can grab them quickly
    disp(['R: ', num2str(R)]);
    disp(['R^2: ', num2str(R2)]);
    disp(['Confidence intervals for the slope: [', num2str(ci_lower), ', ', num2str(ci_upper), ']']);
    disp(sig_text);

    hold off;

elseif pickme == 3

    metric = C;

    percentage_correct = (scores / 30) * 100;
    mdl = fitlm(metric, percentage_correct, 'RobustOpts', 'bisquare');
    R2 = mdl.Rsquared.Ordinary;
    R = sqrt(R2);

    x_fit = linspace(min(metric), max(metric), 100)';
    [y_fit, y_CI] = predict(mdl, x_fit);

    figure('Position', [1000, 500, 500, 450]);
    set(gca, 'DataAspectRatio', [1 1 1]);

    scatter(metric, percentage_correct, 'filled', 'MarkerEdgeColor', 'k');
    hold on;

    plot(x_fit, y_fit, 'k-', 'LineWidth', 1.5); % best-fit line

    fill([x_fit; flipud(x_fit)], [y_CI(:,1); flipud(y_CI(:,2))], ...
         'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none'); % shaded CI band

    xlabel('Offset Threshold (dB SPL)', 'FontSize', FS);
    ylabel('Percentage Correct (%)', 'FontSize', FS);

    xlim([min(C) - 2.5, max(C) + 2.5])
    xticks([50,60,70,80,90,100])
    xticklabels({'50','60','70','80','90','100'})

    title(['R = ', num2str(R, '%.2f'), ' | R^2 = ', num2str(R2, '%.2f')], 'FontSize', FS);

    set(gca, 'FontSize', TS);
    grid off;

    CI = coefCI(mdl);
    ci_upper = CI(2, 2);
    ci_lower = CI(2, 1);

    if ci_lower <= 0 && ci_upper >= 0
        sig_text = 'Slope is NOT significantly different from zero.';
    else
        sig_text = 'Slope is significantly different from zero.';
    end
    text(mean(metric), max(percentage_correct), sig_text, 'FontSize', TS);

    disp(['R: ', num2str(R)]);
    disp(['R^2: ', num2str(R2)]);
    disp(['Confidence intervals for the slope: [', num2str(ci_lower), ', ', num2str(ci_upper), ']']);
    disp(sig_text);

    hold off;

elseif pickme == 4

    metric = D;

    percentage_correct = (scores / 30) * 100;
    mdl = fitlm(metric, percentage_correct, 'RobustOpts', 'bisquare');
    R2 = mdl.Rsquared.Ordinary;
    R = sqrt(R2);

    x_fit = linspace(min(metric), max(metric), 100)';
    [y_fit, y_CI] = predict(mdl, x_fit);

    figure('Position', [1000, 500, 500, 450]);
    set(gca, 'DataAspectRatio', [1 1 1]);

    scatter(metric, percentage_correct, 'filled', 'MarkerEdgeColor', 'k');
    hold on;

    plot(x_fit, y_fit, 'k-', 'LineWidth', 1.5); % best-fit line

    fill([x_fit; flipud(x_fit)], [y_CI(:,1); flipud(y_CI(:,2))], ...
         'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none'); % shaded CI band

    xlabel('Hysteresis (dB)', 'FontSize', FS);
    ylabel('Percentage Correct (%)', 'FontSize', FS);

    xlim([min(D) - 0.1, max(D) + 0.1])
    xticks([0.8,1,1.2,1.4,1.6,1.8,2.0])
    xticklabels({'0.8','1.0','1.2','1.4','1.6','1.8','2.0'})

    title(['R = ', num2str(R, '%.2f'), ' | R^2 = ', num2str(R2, '%.2f')], 'FontSize', FS);

    set(gca, 'FontSize', TS);
    grid off;

    CI = coefCI(mdl);
    ci_upper = CI(2, 2);
    ci_lower = CI(2, 1);

    if ci_lower <= 0 && ci_upper >= 0
        sig_text = 'Slope is NOT significantly different from zero.';
    else
        sig_text = 'Slope is significantly different from zero.';
    end
    text(mean(metric), max(percentage_correct), sig_text, 'FontSize', TS);

    disp(['R: ', num2str(R)]);
    disp(['R^2: ', num2str(R2)]);
    disp(['Confidence intervals for the slope: [', num2str(ci_lower), ', ', num2str(ci_upper), ']']);
    disp(sig_text);

    hold off;

elseif pickme == 5
    metric = G * 1000;

    percentage_correct = (scores / 30) * 100;
    mdl = fitlm(metric, percentage_correct, 'RobustOpts', 'bisquare');
    R2 = mdl.Rsquared.Ordinary;
    R = sqrt(R2);

    x_fit = linspace(min(metric), max(metric), 100)';
    [y_fit, y_CI] = predict(mdl, x_fit);

    figure('Position', [1000, 500, 500, 450]);
    set(gca, 'DataAspectRatio', [1 1 1]);

    scatter(metric, percentage_correct, 'filled', 'MarkerEdgeColor', 'k');
    hold on;

    plot(x_fit, y_fit, 'k-', 'LineWidth', 1.5); % best-fit line

    fill([x_fit; flipud(x_fit)], [y_CI(:,1); flipud(y_CI(:,2))], ...
         'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none'); % shaded CI band

    xlabel('Delay (ms)', 'FontSize', FS);
    ylabel('Percentage Correct (%)', 'FontSize', FS);

    xlim([min(metric)-50, max(metric)+50])

    title(['R = ', num2str(R, '%.2f'), ' | R^2 = ', num2str(R2, '%.2f')], 'FontSize', FS);

    set(gca, 'FontSize', TS);
    grid off;

    CI = coefCI(mdl);
    ci_upper = CI(2, 2);
    ci_lower = CI(2, 1);

    if ci_lower <= 0 && ci_upper >= 0
        sig_text = 'Slope is NOT significantly different from zero.';
    else
        sig_text = 'Slope is significantly different from zero.';
    end
    text(mean(metric), max(percentage_correct), sig_text, 'FontSize', TS);

    disp(['R: ', num2str(R)]);
    disp(['R^2: ', num2str(R2)]);
    disp(['Confidence intervals for the slope: [', num2str(ci_lower), ', ', num2str(ci_upper), ']']);
    disp(sig_text);

    hold off;

elseif pickme == 6

    metric = H;

    percentage_correct = (scores / 30) * 100;
    mdl = fitlm(metric, percentage_correct, 'RobustOpts', 'bisquare');
    R2 = mdl.Rsquared.Ordinary;
    R = sqrt(R2);

    x_fit = linspace(min(metric), max(metric), 100)';
    [y_fit, y_CI] = predict(mdl, x_fit);

    figure('Position', [1000, 500, 500, 450]);
    set(gca, 'DataAspectRatio', [1 1 1]);

    scatter(metric, percentage_correct, 'filled', 'MarkerEdgeColor', 'k');
    hold on;

    plot(x_fit, y_fit, 'k-', 'LineWidth', 1.5); % best-fit line

    fill([x_fit; flipud(x_fit)], [y_CI(:,1); flipud(y_CI(:,2))], ...
         'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none'); % shaded CI band

    xlabel('Hysteresis (dB)', 'FontSize', FS);
    ylabel('Percentage Correct (%)', 'FontSize', FS);

    title(['R = ', num2str(R, '%.2f'), ' | R^2 = ', num2str(R2, '%.2f')], 'FontSize', FS);

    set(gca, 'FontSize', TS);
    grid off;

    CI = coefCI(mdl);
    ci_upper = CI(2, 2);
    ci_lower = CI(2, 1);

    if ci_lower <= 0 && ci_upper >= 0
        sig_text = 'Slope is NOT significantly different from zero.';
    else
        sig_text = 'Slope is significantly different from zero.';
    end
    text(mean(metric), max(percentage_correct), sig_text, 'FontSize', TS);

    disp(['R: ', num2str(R)]);
    disp(['R^2: ', num2str(R2)]);
    disp(['Confidence intervals for the slope: [', num2str(ci_lower), ', ', num2str(ci_upper), ']']);
    disp(sig_text);

    hold off;

end
```
