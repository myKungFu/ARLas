function [stimScaled,scaling,errors] = ARLas_applyLongTubeCal(LTC,targetLvl,f,stim)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [stimScaled,scaling,errors] = ARLas_applyLongTubeCal(LTC,targetLvl,f,stim);
%
% Apply in-situ calibration to a stimulus.
% 
% INPUT ARGUMENTS:
% LTC = Long Tube Calibration structure obtained using ARLas_longTubeCal.m
% targetLvl = desired stimulus level in dB SPL
% f = desired frequency bandwidth of click (Hz). Give two cutoff frequencies, e.g. [1000,8000].
%
% OUTPUT ARGUMENTS:
% stimScaled = stimulus scaled to result in achieving target level and zero phase.
% scaling = multiplier used to scale stimulus and achieve zero phase.
%       Scaling is the impulse response (FIR) used to filter stim to yield stimScaled.
% errors = report of whether target can be reached (1) or not (0).
%
% Auditory Research Lab, The University of Iowa
% Deptartment of Communication Sciences & Disorders
% The University of Iowa
% Author: Shawn Goodman
% This code based on ARLas_applyISC, Original Date: December 1, 2017
% Updated: October 16, 2018, July 9, 2019, November 5, 2019
% This new version for long tube calibration, Original Date: June 7, 2021
% Last Updated: June 9, 2021 -- ssg
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    fs = LTC.fs;
    if length(f)==1 % if the stimulus is narrowband
        % get the full-out amplitude from the in-situ cal measurement
        fo = interp1(LTC.freq,LTC.fo_spl,f,'pchip');
        phi = interp1(LTC.freq,LTC.phi,f,'pchip');
        d = targetLvl - fo; % calculate the dB difference between full out and desired (target) level
        k = 10^(d/20); % express dB difference as a linear multiplier
        if k >= 1 % if cannot achieve target output because desired > than full out
            k = 0; % set multiplier to zero
            errors = d; % error is the dB amount that exceeds full out
            warning('The requested dB value exceeds full out.')
        else
            errors = 0;
        end
        a = hilbert(stim); % create an analytic signal
        stimScaled = real(a .* k*exp(1i*phi)); % scale and shift stimulus
        scaling = k*cos(phi)+1i*k*sin(phi);
    else % the stimulus is broadband: create an inverse magnitude FIR filter 
        tau = 0.0015; % filter group delay (sec)
        L = round(tau*2 * fs); % calculate filter order
        if mod(L,2)~=1 % ensure filter order is even
            L = L + 1;
        end
        M = L-1; % filter order
        M2 = M/2; % filter group delay in samples
        N = fs; % number of samples in design (using FIR from the "window" method)
        F = (0:1:N-1)'*(fs/N); % make a frequency axis corresponding to X
        f1 = min(f); % minimum frequency to include (Hz)
        f2 = max(f); % maximum frequency to include (Hz)
        if f1 < LTC.fmin
            error('Requested fmin value is < long tube calibration fmin.')
        end
        if f2 > LTC.fmax
            error('Requested fmax value is > long tube calibration fmin.')
        end
        [~,indx1] = min(abs(F - f1)); % index of minimum frequency
        [~,indx2] = min(abs(F - f2)); % index of maximum frequency
        [~,indx4] = min(abs(F-(fs - F(indx1)))); % the locations of the alised portions of the spectrum
        [~,indx3] = min(abs(F-(fs - F(indx2))));
        Fcut = F(indx1:indx2); % cut frequency vector down to size

        fo =  interp1(LTC.freq,LTC.fo_spl,Fcut,'pchip'); % this is full out magnitude in dB
        phi = interp1(LTC.freq,LTC.phi,Fcut,'pchip'); % this is phase in radians

        % calculate the maximum possible output as an impulse
        Mag = zeros(fs,1); % initialize magnitude spectrum
        Phase = zeros(fs,1); % initialize phase spectrum
        bar = min(fo); % level of constant output (dB SPL)
        foInv = bar - fo; % inverse full output
        foInv = foInv + bar; % put back around min value
        %target = foInv; % (for plotting purposes only)
        foInv = 10.^(foInv/20)*0.00002; %foInv = 10.^(foInv/20);  % put into linear units (Pa)

        % create an impulse response from the presumed flat output spectrum
        Mag(indx1:indx2) = 10.^(bar/20)*0.00002;
        Mag(indx3:indx4) = flipud(10.^(bar/20)*0.00002);
        Phase(indx1:indx2) = 0;
        Phase(indx3:indx4) = flipud(0);
        X = Mag .* cos(Phase) + 1i*Mag.*sin(Phase); % spectrum, complex rectangular form
        x = real(ifft(X)); % time domain impulse
        x = [x(end-M2:end);x(1:M2)]; % make filter causal and cut down to size
        h = hann(length(x)); % window
        foImpulse = x .* h; % apply window to the impulse
        
        maxOut = bar; % maximum output in dB SPL
        d = targetLvl - maxOut; % calculate the dB difference between full out and desired (target) level
        k = 10^(d/20); % express dB difference as a linear multiplier
        
        % create an impulse response that will create a flat output
        Mag(indx1:indx2) = foInv;
        Mag(indx3:indx4) = flipud(foInv);
        Phase(indx1:indx2) = -phi;
        Phase(indx3:indx4) = flipud(phi);
        X = Mag .* cos(Phase) + 1i*Mag.*sin(Phase); % spectrum, complex rectangular form
        x = real(ifft(X)); % time domain impulse
        x = [x(end-M2:end);x(1:M2)]; % make filter causal and cut down to size
        h = hann(length(x)); % window
        a = x .* h; % apply window to the impulse
        
        
        a = a / max(abs(a));
        
        %scaling = a * k;
        scaling = a;
        
        errors = 0;
        stimScaled = fastFilter(scaling,stim);
    end
end

% INTERNAL FUNCTIONS ------------------------------------------------------
