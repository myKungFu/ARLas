function [epsilon] = findOptimalLengths1(cavityLengths)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% epsilon = findOptimalLengths1(cavityLengths);
%
% Minimize the error based on ideal and calculated PRESSURE (Neely method)
% Based originally on code generously provided by Stephen T. Neely
% Deptartment of Communication Sciences & Disorders
% The University of Iowa
% Author: Shawn Goodman
% Date: July 21, 2015
% Updated: November 4, 2017 - ssg
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

sep = filesep;
q = load(['C:',sep,'myWork',sep,'ARLas',sep,'temp.mat']);
obj = q.obj;

ZLi = calculate_ZLi(cavityLengths,obj);
[PS,ZS] = calculate_source(ZLi,obj);
epsilon = calculate_error(ZLi,PS,ZS,obj);
% figure(113)
% plot(epsilon,'*')
% hold on
% pause(.01)

% subfunctions ------------------------------------------------------------
function [ZLi] = calculate_ZLi(cavityLengths,obj)
% calculate theoretical cavity impedance
R = exp(-2 * obj.wn * cavityLengths'); % wave equation exponent.  Volume velocity is 0 at the closed end of the tube.
%R = exp(-2 * obj.wn * cavityLengths(:)); 

% zc = z0 .* (1 + R) ./ (1 - R); % ideal cavity impedance; Alternative notation: Zc = -iZo*cot(kL), where k = wavenumber
% this formulation of zc assumes perfectly reflective cavities, which
% is not the case with real tubes. Reducing R in the numerator reduces
% the notch depth; reducing R in the denominator reduces the peak
% height. Both values of R can also be reduced. 
ZLi = obj.z0 .* (1 + (R * obj.reflNum)) ./ (1 - (R * obj.reflDenom)); % this reduces the notch depth

function [PS,ZS] = calculate_source(ZLi,obj)    
% calculate thevenin source pressure and impedance
nFreqs = size(obj.PL,1); % number of frequencies in analysis range
ZS = zeros(nFreqs,1); % initialize source impedance vector
PS = zeros(nFreqs,1); % initialize source pressure vector
for ii=1:nFreqs
   z = ZLi(ii,:).'; % cavity impedance (theoretical); looping across rows (frequency), taking all columns (cavities) simultaneously
   p = obj.PL(ii,:).'; % cavity pressure (measured); looping across rows (frequency), taking all columns (cavities) simultaneously
   A = [z -p]; % Zc - Pc = cavity impedance - cavity pressure (measured and defined as H)
   B = z .* p; % Zc * Pc = cavity impedance * cavity pressure
   x = A \ B; % matrix division
   PS(ii) = x(1); % Ps = source pressure
   ZS(ii) = x(2); % Zs = source impedance
end

function [epsilon] = calculate_error(ZLi,PS,ZS,obj)
% calculate the error in the fit
% An error (epsilon) value of 1 is considered to be acceptable with this method
doPlot = 0;
% calculate error term: the variance in the measured pressure NOT accounted for by theoretical (source) pressure
% calculate crosstalk; similar to "subtracting DC component of voltage when calculating the rms value"; allows another degree of freedom in the equations
s1 = 0; % initialize sums of squares variables
s2 = 0;
f1 = obj.fminIndx;
f2 = obj.fmaxIndx;
PL = obj.PL(f1:f2,:);
ZLi = ZLi(f1:f2,:);
PS = PS(f1:f2);
ZS = ZS(f1:f2);
[nFreqs,nCavities] = size(PL);
accumulatedPressure = zeros(nFreqs,1);
for ii=1:nCavities % loop across cavities
    PLi(:,ii) = PS .* ZLi(:,ii) ./ (ZS + ZLi(:,ii)); % idealized pressure
    %deltaP = PL(:,ii) - PS .* ZLi(:,ii) ./ (ZS + ZLi(:,ii)); % Numerator |Pc - Ps * Zc / (Zs + Zc)|
    deltaP = PL(:,ii) - PLi(:,ii);
    accumulatedPressure = accumulatedPressure + deltaP; % sum of pressure differences across cavities; only used if crosstalk is being considered
    s1 = s1 + sum(abs(deltaP).^2); % sum of squares:  square the numerator and sum across frequency;
    s2 = s2 + sum(abs(PL(:,ii)).^2); % sum of squares: square the measured cavity pressure and sum across frequency
end
s1 = s1 - sum(abs(accumulatedPressure).^2) / nCavities; % line not used if crosstalk isn't being considered
epsilon = 10000 * s1 / s2; % divide the sum of squares of the pressure difference with the sum of squares of the measured pressure and multiply by 10,000
if doPlot == 1
    % Plot to watch the fit converge
    figure(500)
    subplot(2,1,1)
    hold off
    plot(20*log10(abs(PLi)))
    hold on
    plot(20*log10(abs(PL)),':')
    subplot(2,1,2)
    plot(20*log10(abs(deltaP)),'r')
    pause
end

