load('C:\myWork\ARLas\MEMR_Data_Final\MEMR_groupData_15.mat') % Use one
% of these
% load('C:\myWork\ARLas\MEMR_Data_Final_nt\MEMR_groupData.mat') % 
pAmp = 20*log10(pAmp);
means = mean(pAmp);
std_devs = std(pAmp);

% Mean and Std
disp('Means of the four tests:');
disp(means);
disp('Standard deviations of the four tests:');
disp(std_devs);

% Boxplot
figure;
ylim([1, 3.5])
boxplot(pAmp, 'Labels', {'Run 1', 'Run 2', 'Run 3', 'Run 4'});
title('Boxplot of Peak Amplitudes Across Tests');
xlabel('Run');
ylabel('Peak Amplitude');

% Normality Test (Shapiro-Wilk test), download the extension from mathwork
% community website
normality_p_values = zeros(1, 4);
for i = 1:4
    [~, p] = swtest(pAmp(:, i));
    normality_p_values(i) = p;
end

disp('Shapiro-Wilk Test p-values:');
disp(normality_p_values);

% Multiple Paired t-tests
disp('Paired t-tests:');
for i = 1:3
    [h, p] = ttest(pAmp(:, i), pAmp(:, i+1));
    fprintf('Run %d vs Run %d: p-value = %.4f\n', i, i+1, p);
end

% Wilcoxon Signed-Rank Test (if pAmp is not normal)
disp('Wilcoxon Signed-Rank Tests:');
for i = 1:3
    p = signrank(pAmp(:, i), pAmp(:, i+1));
    fprintf('Run %d vs Run %d: p-value = %.4f\n', i, i+1, p);
end






load('D:\MEMR_GroupData\MEMR_groupData.mat') % Use one
% of these
% load('C:\myWork\ARLas\MEMR_Data_Final_nt\MEMR_groupData.mat') % 
pAmp = 20*log10(pAmp);
means = mean(pAmp);
std_devs = std(pAmp);

% Display means and standard deviations
disp('Means of the four runs:');
disp(means);
disp('Standard deviations of the four runs:');
disp(std_devs);

% Boxplot
figure;
ylim([1, 3.5])
hold on;
boxplot(pAmp, 'Labels', {'Run 1', 'Run 2', 'Run 3', 'Run 4'});

% Add participant lines
for i = 1:30
    plot(1:4, pAmp(i, :), '-o', 'Color', [0.7 0.7 0.7], 'LineWidth', 1);
end

title('Yes Detrending');
xlabel('Runs');
ylabel('Peak Amplitude');
hold off;

% Normality Test (Lilliefors test)
normality_p_values = zeros(1, 4);
for i = 1:4
    [~, p] = swtest(pAmp(:, i));
    normality_p_values(i) = p;
end

disp('Shapiro-Wilk Test p-values:');
disp(normality_p_values);

% Multiple Paired t-tests
disp('Paired t-tests:');
for i = 1:3
    [h, p] = ttest(pAmp(:, i), pAmp(:, i+1));
    fprintf('Run %d vs Run %d: p-value = %.4f\n', i, i+1, p);
end

% Wilcoxon Signed-Rank Test (if data is not normal)
disp('Wilcoxon Signed-Rank Tests:');
for i = 1:3
    p = signrank(pAmp(:, i), pAmp(:, i+1));
    fprintf('Run %d vs Run %d: p-value = %.4f\n', i, i+1, p);
end
