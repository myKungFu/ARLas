function [epsilon] = findOptimalLengths2(cavityLengths)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% epsilon = findOptimalLengths2(cavityLengths);
%
% Minimize the error based on ideal and calculated IMPEDANCE (Goodman method)
%
% Deptartment of Communication Sciences & Disorders
% The University of Iowa
% Author: Shawn Goodman
% Date: July 21, 2015
% Updated: Jun2 13, 2017 - ssg; forced cavity lengths to be column vector
%                               on line 19
% Updated: November 4, 2017 - ssg
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

sep = filesep;
q = load(['C:',sep,'myWork',sep,'ARLas',sep,'temp.mat']);
obj = q.obj;

ZLi = calculate_ZLi(cavityLengths(:),obj);
[PS,ZS] = calculate_source(ZLi,obj);
ZL = calculate_ZL(PS,ZS,obj);
epsilon = calculate_error(ZL,ZLi,obj);

% subfunctions ------------------------------------------------------------
function [ZLi] = calculate_ZLi(cavityLengths,obj)
% calculate theoretical cavity impedance
R = exp(-2 * obj.wn * cavityLengths'); % wave equation exponent.  Volume velocity is 0 at the closed end of the tube.
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

function [ZL] = calculate_ZL(PS,ZS,obj)
% cavity impedance (calculated, not ideal)
nCavities = size(obj.PL,2); % number of cavities to analyze
ZL = zeros(size(obj.PL)); % initialize
for ii=1:nCavities  % calculate load impedance for each cavity
    ZL(:,ii) = ZS .* obj.PL(:,ii) ./ (PS - obj.PL(:,ii));
end    
    
function [epsilon] = calculate_error(ZL,ZLi,obj)
% calculate the error in the fit
doPlot = 0;
f1 = obj.fminIndx;
f2 = obj.fmaxIndx;
deltaZ = ZL(f1:f2,:) - ZLi(f1:f2,:);
N = size(deltaZ(:),1);
%mag = abs(deltaZ);
%phi = angle(deltaZ);
%eMag = sqrt(sum(sum(mag.^2)))/N;
%w = abs(mag).^2;
%ePhi = sum(sum(abs(phi).*w))./sum(sum(w));
eTotal = sum(sum(abs(deltaZ)))/N;
%epsilon = eMag;
%epsilon = ePhi;
epsilon = eTotal;
if doPlot == 1
    % Plot to watch the fit converge
    figure(500)
    subplot(2,1,1)
    hold off
    plot(20*log10(abs(ZLi(f1:f2,:))))
    hold on
    plot(20*log10(abs(ZL(f1:f2,:))),':')
    subplot(2,1,2)
    plot(20*log10(abs(deltaZ)),'r')
    pause
end
