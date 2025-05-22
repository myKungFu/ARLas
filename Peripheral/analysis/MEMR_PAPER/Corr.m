load('D:\MEMR_AnalysisM\MEM03_Mepani_Analysis1.mat')
S2 = MEMR_mem;
tt = S2.timeMepani(1:12);
load('D:\MEMR_AnalysisComp\MEM03_Run1_Analysis1.mat')
S1 = MEMR_mem;
% Extract
vec2 = S2.d1(1:12)';
vec1 = S1.d1 (1:58);

x = vec2; 
y = vec1; 
xq1 = tt;
% p = pchip(y,x,xq1);
s = spline(y,x,xq1);
% m = makima(x,y,xq1);
plot(x,y,'o',xq1,p)

% Assuming you have the continuous data curve in matrix 'data_curve' and the discrete points in matrix 'discrete_points'
% 'tt' is the time vector

% Interpolate the continuous data curve at the discrete time points
id = interp1(tt, vec1, 1:12);

% Calculate the correlation coefficient between the interpolated data and discrete points
correlation_coefficient = corrcoef(id, vec2);

disp(['Correlation Coefficient: ', num2str(correlation_coefficient(1,2))]);



function flooredVal = floorS(val, nS)
    pw = ceil(log10(val)); % Find order of magnitude of val.
    res = 10^(pw-nS); % Resolution to round to.
    flooredVal = floor(val/res)*res; % < change floor() to ceil(), for ceiling equivalent.
end