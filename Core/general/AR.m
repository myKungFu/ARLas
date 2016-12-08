function [indx,nRejects] = AR(y,tolerance,doPlot)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [indx,nRejects] = AR(y,tolerance,doPlot);
%
% Performs artifact rejection based on quartiles (Tukey method).  
% This is a robust method because large outliers in the upper and lower 25% 
% of the data do not significantly affect the estimate. 
% This function does not actually delete artifacts; it simply returns an 
% array of indices showing the location of artifacts. Use AR_engine.m to 
% remove the outliers.
%
% y = input data vector
% tolerance = a string specifying the stingency of outlier detection:
%   'moderate' (default) uses 2.25 * the IQR as the criterion
%   'mild' = uses 1.5 * the IQR
%   'severe' uses 3 * the IQR (note that this one is a ssg invention)
% doPlot = a plotting switch.  If set = 1, will plot the sorted values of y
%   along with lines representing the IQRs and cutoff values used for
%   outlier detection. If not specified, will default to 0.
% indx = vector of indices indicating outliers
%
% Author: Shawn Goodman
% Date: 6/15/2005
% Modified: 7/15/2005 -- ssg cleaned up the code
%           7/28/2005 -- ssg added documentation
%           10/22/2013 -- ssg made to handle nan
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%doPlot = 0;

y = y(:); % ensure y is a column vector
%y = y(isfinite(y)); % get rid of NaN
warning('OFF')
if nargin <= 2,
    doPlot = 0; % set = 1 to look at outlier plot
end
if nargin == 1,
    tolerance = 'moderate';
end
switch tolerance
    case 'mild'
        multiplier = 1.5;
    case 'moderate'
        multiplier = 2.25;
    case 'severe'
        multiplier = 3;
    otherwise
        disp('Error: requested tolerance is not correctly specified.')
        return
end

xx = sort(y(:));
[dummy,indx] = max(xx);
xx = xx(1:indx); % get rid of any nan
N = length(xx); % number of observations
q = 100 *(0.5:N-0.5)./N;
xx = [min(xx); xx(:); max(xx)];
q = [0 q 100];
F1 = interp1(q,xx,25); % the first fourth, approx 25th precentile
F2 = interp1(q,xx,50); % the first half, approx 50th precentile
F3 = interp1(q,xx,75); % the third fourth, approx 75th percentile
IQR = F3 - F1; % the interquartile range

ArtifactVector = y >= (F1-multiplier*IQR) & y <= (F3 + multiplier*IQR);
[indx,val] = find(~ArtifactVector); % index the artifacts which should be rejected

nRejects = length(indx); % number of rejected buffers
percentRejected = (length(indx) / N) * 100;
% disp([num2str(percentRejected),' percent of the buffers rejected as artifacts.'])
% commented out previous line on 6/3/2010--ras and ssg

if doPlot == 1,
    figure(63)
    plot(q,xx,'*b-')
    hold on
    L1 = line([q(1) q(end)],[F1 F1]);
    L2 = line([q(1) q(end)],[F2 F2]);
    L3 = line([q(1) q(end)],[F3 F3]);
    L4 = line([q(1) q(end)],[F1-multiplier*IQR F1-multiplier*IQR]);
    L5 = line([q(1) q(end)],[F3 + multiplier*IQR F3 + multiplier*IQR]);
    set(L1,'Color',[0 1 0])
    set(L2,'Color',[0 0 0],'LineStyle',':')
    set(L3,'Color',[0 1 0])
    set(L4,'Color',[1 0 0])
    set(L5,'Color',[1 0 0])
    hold off
    %pause
end
warning('ON')