function [stim,F] = ARLasChrp()

fs = 96000;
fmin = 100; % starting frequency (hz)
% find ending frequency:
%   2 samples = (half a cycle / f) * fs;
%   2 = (0.5/x)*fs
%   0.5*fs/2 = x
%   2400 = x;
fmax = 0.5*fs/2;

N = fs/8;
X = (0:1:N-1)';
X = X / max(X);

% formula = fmin*exp(-a*X)
%
% X(end) = 1, so the value of a is as follows:
% fmax = 100*exp(-aX)
% fmax / 100 = exp(-a)
% -log(fmax / 100) = a
a = -log(fmax / fmin);
F = 100*exp(-a*X);
lead = round(0.5/F(1)*fs); % lead samples for half a cycle
lag = round(0.25/F(end)*fs); % lead samples for a quarter cycle
F = [ones(lead,1)*F(1);F;ones(lag,1)*F(end)];

% Rotate the phase:
% At each time sample, time has advanced by the sampling period (1/fs).
% At each time sample, the phase has rotated by the frequency (f) times
% the sampling period: phase advance (in cycles) = f * (1/fs), or f/fs.
% To make this a radian change, multiply by 2pi: phi = (2*pi*f)/fs;
% Finally, add the newly accumulated phase to the previous phase value.
stim = zeros(length(F),1);
phi = 0; % starting phase is zero
for jj = 1:length(F)
    stim(jj,1) = sin(phi);
    phi = phi + ((2*pi*F(jj))/fs);
end

rampN = round(lead*2);
h = hann(rampN*2);
h = h(1:rampN);
stim(1:rampN) = stim(1:rampN) .* h;

rampN = round(lead*.5);
h = hann(rampN*2);
h = h(1:rampN);
stim = flipud(stim);
stim(1:rampN) = stim(1:rampN) .* h;
stim = flipud(stim);

padLen = 0.050; % zero-pad the stimulus with 5 ms
padN = round(padLen * fs);
pad = zeros(padN,1);
stim = [pad;stim;pad];
F = [pad;F;pad];
stim = stim / max(abs(stim));  % rescale so 1 is max out
stim = stim * .99;

return

b = flipud(p);
b = b(2:end);
stim = [p;b];
stim = stim(2:end);

plot([stim;stim],'b')
hold on
plot(stim,'r')
plot(p,'g')

N = length(stim);
N2 = N/2;

resp1 = stim(1:N2);
resp2 = flipud(stim(N2+1:end));
resp1 = resp1(1:end-1);
resp2 = resp2(2:end);

lead2 = round(0.25/F(1)*fs); % lead samples for a quarter cycle
resp1 = resp1(1:end-lead2+1);
resp2 = resp2(lead2:end);

% 
% 
% figure
% plot(resp1,'b')
% hold on
% plot(resp2,'r:')
% 
% 
% plot([p;b])
% hold on
% plot(p,'r')
% 

%----------------------------------------------------------------------------



% 
% 
% 
% 
% 
% 
% 
% 
% [off,indx] = min(abs(F-20001));
% F = F(1:indx);
%     
    
% General model Exp2:
%      f(x) = a*exp(b*x) + c*exp(d*x)
% Coefficients (with 95% confidence bounds):
%        a =       104.1  (103.9, 104.3)
%        b =       5.298  (5.295, 5.3)
%        c =   5.179e-11  (2.903e-11, 7.456e-11)
%        d =       30.93  (30.5, 31.37)
% Goodness of fit:
%   SSE: 3.323e+08
%   R-square: 0.9998
%   Adjusted R-square: 0.9998
%   RMSE: 76.57
% a =       104.1;
% b =       5.298;
% c =   5.179e-11;
% d =       30.93;
% F = a*exp(b*X) + c*exp(d*X);

%     
% lead = round(0.5/F(1)*fs); % lead samples for half a cycle
% lag = round(0.5/F(end)*fs); % lead samples for half a cycle
% 
% F = [ones(lead,1)*F(1);F;ones(lag,1)*F(end)];

% Rotate the phase:
% At each time sample, time has advanced by the sampling period (1/fs).
% At each time sample, the phase has rotated by the frequency (f) times
% the sampling period: phase advance (in cycles) = f * (1/fs), or f/fs.
% To make this a radian change, multiply by 2pi: phi = (2*pi*f)/fs;
% Finally, add the newly accumulated phase to the previous phase value.



% fs = 96000;
% fmin = 100;
% fmax = 22100;
% sweepRate = 21900; %200000; %30000; % Hz per second
% sweepLen = (fmax-fmin) / sweepRate;
% sweepN = round(sweepLen * fs); % number of samples in sweep
% if mod(sweepN,2)~=0
%     sweepN = sweepN + 1;
% end
% fStepSize = (fmax-fmin) / sweepN; % frequency step size (Hz);
% F = (fmin:fStepSize:fmax)'; % frequency vector
% F = F(1:sweepN);
% 
% 
%    % morph F to account for equal "on" time of periods
%     P = 1./F; % period in sec is inverse of frequency
%     P = P / (min(P)); % normalize so that the smallest period is 1 sample
%     P = round(P);
%     X = zeros(sum(P),1);
%     start = 1;
%     finish = start + P(1) - 1;
%     for ii=1:length(P)-1
%         X(start:finish) = F(ii);
%         start = finish + 1;
%         finish = start + P(ii+1)-1;
%     end
%     sweepN = length(X);
%     F = X;
%     F(end) = F(end-1);
%     
%     X = (0:length(F)-1)';
%     X = X / max(X);    
% 
%     
% % General model:
% %      f(x) = 100*exp(-a*x)
% % Coefficients (with 95% confidence bounds):
% %        a =      -5.369  (-5.369, -5.369)
% % Goodness of fit:
% %   SSE: 9.291e+08
% %   R-square: 0.9994
% %   Adjusted R-square: 0.9994
% %   RMSE: 127.7    
% %a = -5.369;

