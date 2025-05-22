load('D:\MEMR_AnalysisM\MEM03_Mepani_Analysis1.mat')
S2 = MEMR_mem;
tt = S2.timeMepani(1:12);
load('D:\MEMR_AnalysisComp\MEM03_Run3_Analysis1.mat')
S1 = MEMR_mem;
% Extract
vec2 = S2.d1(1:12)';
vec1 = S1.d1 (1:58);
t1 = S2.timeMepani(2:12);
t2 = S2.timeMepani(13:23);
T1 = S1.t(1:80)';
T2 = S1.t(81:160)'';
mep = S2.d1(2:23)';
mep1 = S2.d1(2:12)';
mep2 = S2.d1(13:23)';

t = S2.timeMepani (2:end-1)';
T = S1.t(2:end-1)';

nearest_values = zeros(size(t1));
nearest_indices = zeros(size(t1));

% Iterate
for i = 1:length(t)

    min_diff = Inf;
    min_index = 0;

    % Iterate through T1
    for j = 1:length(T)
        % abs diff
        diff = abs(t(i) - T(j));

        % Update
        if diff < min_diff
            min_diff = diff;
            min_index = j;
        end
    end

    % Store
    nearest_values(i) = T(min_index);
    nearest_indices(i) = min_index;
end

% Display
disp("Nearest values in vector 2:");
disp(nearest_values);
disp("Corresponding indices in vector 2:");
disp(nearest_indices);
ssg1 = (1:1:22)';

for ii=1:length(t)
    ssg1(ii) = S1.d2(nearest_indices(ii));
end

x = ssg1; 
y = mep; 
% xq1 = tt;
% p = pchip(y,x,xq1);
% s = spline(y,x,xq1);
% 
% plot(x,y,'o',xq1,p)
% 
% id = interp1(tt, vec1, 1:12);
% correlation_coefficient = corrcoef(id, vec2);
% 
% disp(['Correlation Coefficient: ', num2str(correlation_coefficient(1,2))]);


% 
% function flooredVal = floorS(val, nS)
%     pw = ceil(log10(val)); % Find order of magnitude of val.
%     res = 10^(pw-nS); % Resolution to round to.
%     flooredVal = floor(val/res)*res; % < change floor() to ceil(), for ceiling equivalent.
% end