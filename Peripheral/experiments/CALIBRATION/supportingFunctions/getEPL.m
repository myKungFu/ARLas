function [correction,canalLength,Fqw,RS,RL] = getEPL(z0,ZS,ZL,IPL,SPL,freq,c,doPlot)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [correction,canalLength,Fqw,RS,RL] = getEPL(z0,ZS,ZL,IPL,SPL,freq,c);
%
% Calculate emission pressure level (EPL) correction.
% Based on Karolina Charaziak and Chris Shera (2017).
%
% INPUT ARGUMENTS:
% z0 = estimated canal surge impedance
% ZS = complex source impedance
% ZL = complex load impedance
% freq = vector of frequencies (in Hz) associated with ZL and ZL 
% c = assumed speed of sound in air
% doPlot = turn on and off plotting (1 = on; 0 = off)
%
% OUTPUT ARGUMENTS:
% correction = mulitiply OAE (in complex freq domain) by this to get EPL
% canalLength = length of cavity (canal) (cm)
% Fqw = frequency of quarter-wavelength null (Hz)
% RS = source reflectance
% RL = load reflectance
%
% Auditory Research Lab, The University of Iowa
% Deptartment of Communication Sciences & Disorders
% The University of Iowa
% Author: Shawn Goodman
% Date: December 1, 2017
% Updated: October 16, 2018
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if nargin < 8
        doPlot = 0;
    end
    
    % average human ear canal length is 2.2 to 3.2 cm
    % Assuming the speed of sound is 34588 cm/s 
    % We are looking for quarter-wave nulls between 2700 and 5400 Hz
    fmin = 2700;
    fmax = 5400;

    % We are also looking for half-wave resonances between 5400 and 10800 Hz
    % But in this code, we concentrate our search on the quarter wave nulls.
    % The time domain does not provide fine enough sampling for accurate
    % detection. In the frequency domain, the null will be the frequency with
    % the largest difference between SPL (fpl and rpl sum out of phase) and IPL
    % (fpl and rpl sum in phase).
    
    [dummy,iMax] = min(abs(freq - fmax)); % the maximum frequency of the expected null
    [dummy,iMin] = min(abs(freq - fmin)); % the minimum frequency of the expected null
    ff = freq(iMin:iMax);
    delta = IPL - SPL;
    delta = delta(iMin:iMax);
    [dummy,iQW] = max(delta); % maximum of the ipl spl difference
    Fqw = ff(iQW); % frequency of the quarter-wave null
    Aqw = SPL(iQW+iMin); % amplitude of the quarter-wave null (for plotting)
    canalLength = c / (Fqw*4); % length of canal in cm
    canalLength = round(canalLength*1000)/1000; % control precision based on max detectable difference
    
    tau = canalLength / c; % estimated one-way ear canal delay
    omega = freq * 2 * pi; % angular frequency
    t = exp(1i*omega*tau); % the one-way ear canal delay 
    RS = (ZS - z0)./(ZS + z0); % source reflectance
    RL = (ZL - z0)./(ZL + z0); % load reflectance
    % EPL = OAE .* ((1-Rsr.*Rload)./(t.*(1+Rsr))); % Karolina syntax
    % EPL = OAE .* ((1-RS.*RL)./(t.*(1+RS))); % Goodman syntax
    correction = ((1-RS.*RL)./(t.*(1+RS))); % mulitiply OAE (in freq domain) by this to get EPL
    
    % figure showing the calculated diameter and length of ear canal
    if doPlot == 1
        flow = 500;
        fhigh = 20000;
        [dummy,iiMin] = min(abs(freq - flow));
        [dummy,iiMax] = min(abs(freq - fhigh));
        ymin = min(SPL(iiMin:iiMax)) * .95;
        ymax = max(IPL(iiMin:iiMax)) * 1.05;
        figure 
        plot(freq/1000,SPL,'b')
        hold on
        plot(freq/1000,IPL,'r')
        xlim([flow/1000,fhigh/1000])
        ylim([ymin,ymax])
        line([freq(iMin)/1000,freq(iMin)/1000],[ymin,ymax],'Color',[0 0 0],'LineStyle',':','LineWidth',0.5)
        line([freq(iMax)/1000,freq(iMax)/1000],[ymin,ymax],'Color',[0 0 0],'LineStyle',':','LineWidth',0.5)
        plot(freq(iQW+iMin)/1000,Aqw,'go')
        xlabel('Frequency (kHz)','FontSize',12)
        ylabel('Magnitude (Max Output)','FontSize',12)
        title(['Diameter = ',num2str(diam),'cm     Length = ',num2str(canalLength),' cm'])
        legend('SPL','IPL','fmin','fmax','quarter-wave null')    
    end
end

%--------------------------------------------------------------------------
% Calculations follow:
%  HALF-WAVE RESONANCES ----------------
% LONGEST
% half-wave resonsance is at 3.2*2 cm
% Rhalf = (3.2 * 2)/34588 = 1.850352723487915e-04 seconds (185 usec)
%         which corrsponds to Rhalf * fs = 8.16 samples
% Fhalf = 1 / Rhalf = 5404 Hz
% SHORTEST UN-OCCLUDED
% half-wave resonsance is at 2.2*2 cm
% Rhalf = (2.2 * 2)/34588 = 1.272117497397942e-04 seconds (127 usec)
%         which corrsponds to Rhalf * fs = 5.61 samples
% Fhalf = 1 / Rhalf = 7860 Hz
% SHORTEST OCCLUDED
% half-wave resonsance is at 2.2*2 cm - .6 cm (eartip) = 1.6*2 cm
% Rhalf = (1.6 * 2)/34588 = 9.251763617439575e-05 seconds (92.5 usec)
%         which corrsponds to Rhalf * fs = 4.08 samples
% Fhalf = 1 / Rhalf = 10808 Hz
%  QUARTER-WAVE NULLS ---------------
% LONGEST
% quarter-wave NULL is at 3.2*4 cm
% Rquart = (3.2 * 4)/34588 = 3.700705446975830e-04 seconds (370 usec)
%         which corrsponds to Rquart * fs = 16.32 samples
% Fhalf = 1 / Rquart = 2702 Hz
% SHORTEST UN-OCCLUDED
% quarter-wave null is at 2.2*4 cm
% Rquart = (2.2 * 4)/34588 = 2.544234994795883e-04 seconds (254 usec)
%         which corrsponds to Rquart * fs = 11.22 samples
% Fhalf = 1 / Rquart = 3930 Hz
% SHORTEST OCCLUDED
% quarter-wave null is at 2.2*4 cm - .6 cm (eartip) = 1.6*4 cm
% Rquart = (1.6 * 4)/34588 = 1.850352723487915e-04 seconds (185 usec)
%         which corrsponds to Rquart * fs = 8.16 samples
% Fhalf = 1 / Rquart = 5404 Hz

% This is for group delay estimate, which does not currently work (4/17/2018 ssg)
%     phi = unwrap(angle(PL(obj.fminIndx:end))); % unwrapped phase of recording
%     freq = obj.freq(obj.fminIndx:end);
%     phi = phi(1:1000);
%     freq = freq(1:1000);
%     omega = freq * 2 * pi;
%     tau = -diff(phi) ./ diff(omega); % group delay of recording
%     tau = tau'; % put into a row vector
%     tau = tau/2; % one-way delay
%     w = 20*log10(abs(PL(obj.fminIndx:end)));
%     w = w(1:1000);
%     w = w - min(w);
%     w = w(1:end-1);
%     w = w';
%     tau2 = sum(tau .* w) / sum(w);
%     [indx,nRejects] = AR(tau,'mild',0); % Tukey-style artifact rejection
%     tau = AR_engine(tau,indx);
%     w2 = AR_engine(w,indx);
%     tau3 = sum(tau .* w2) / sum(w2);
%     tau4 = mean(tau);
%     tau = median(tau); % time delay of chirp in seconds
%     % convert time delay to distance
%     % distance = speed of sound * time
%     [tau;tau2;tau3;tau4] * obj.c
%     % tau is the windsorized median
%     % tau2 is the weighted mean
%     % tau3 is the weighted mean, windsorized

% This is for a notch and peak based estimate of length, which is the best
% method currently (4/17/2018 ssg)
% average human ear canal length is 2.2 to 3.2 cm
% Assuming the speed of sound is 34588 cm/s 
% We are looking for quarter-wave nulls between 2700 and 5400 Hz
% We are looking for half-wave resonances between 5400 and 10800 Hz
% Here, we concentrate our search on the quarter wave nulls.
% freq = obj.freq;
% [~,iMax] = min(abs(freq - 5400));
% [~,iMin] = min(abs(freq - 2700));
% M = 20*log10(abs(PL(iMin:iMax)));
% F = freq(iMin:iMax);
% [~,iQW] = min(M);
% Fqw = F(iQW); % frequency of the quarter-wave null
% canalLength2 = obj.c / (Fqw*4); % length of canal, estimated directly from the null
