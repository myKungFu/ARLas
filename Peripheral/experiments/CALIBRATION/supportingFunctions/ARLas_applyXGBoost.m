function [stimScaled,scaling,errors] = ARLas_applyXGBoost(isc,targetLvl,f,stim)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [stimScaled,scaling,errors] = ARLas_applyXGBoost(isc,targetLvl,type,f,stim);
%
% Apply XGBoost calibration to a stimulus.
% 
% INPUT ARGUMENTS:
% isc = structure obtained using XGBoost, containing the following fields
%    isc.freq = frequency in Hz
%    isc.spl = sound pressure level in dB SPL
%
% targetLvl = desired stimulus level in dB SPL
% f = dominant frequency of stimulus (Hz) 
%       For narrowband stimuli, give a single frequency in Hz.
%       For broadband stimuli, give two cutoff frequencies in Hz.
% stim = stimulus of unit amplitude (i.e. full output)
%
% OUTPUT ARGUMENTS:
% stimScaled = stimulus scaled to result in achieving target level.
% scaling = multiplier used to scale stimulus and achieve zero phase.
%       If stimulus is narrowband, scaling is a scalar multiplier (complex)
%       If stimulus is broadband, scaling is the impulse response (FIR) used to
%          filter stim to yield stimScaled.
% errors = report of whether target can be resached (1) or not (0).
%
% Auditory Research Lab, The University of Iowa
% Deptartment of Communication Sciences & Disorders
% The University of Iowa
% Author: Shawn Goodman & Daniel Tay
% Date: February 20, 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    correction = 76.97; % correction to get to full output
    isc.spl = isc.spl + correction;

    if length(f)>2 % if the stimulus is broadband
        % get the full-out amplitude from XGBoost
        fo = interp1(isc.freq,isc.spl,f,'pchip');
        phi = zeros(size(fo));

        d = targetLvl - fo; % calculate the dB difference between full out and desired (target) level
        k = 10.^(d/20); % express dB difference as a linear multiplier
        indx = find(f==0);
        k(indx) = 0;
        if any(k>1) % changed 8/12/2024 -- ssg --if k >= 1 % if cannot achieve target output because desired > than full out
            k = 0; % set multiplier to zero
            errors = d; % error is the dB amount that exceeds full out
            warning('The requested dB value exceeds full out.')
        else
            errors = 0;
        end
        a = hilbert(stim); % create an analytic signal
        stimScaled = real(a .* k.*exp(1i.*phi)); % scale and shift stimulus
        scaling = k.*cos(phi)+1i.*k.*sin(phi);
    elseif length(f)==1 % if the stimulus is narrowband ---------------------------------------------------
        % get the full-out amplitude from the in-situ cal measurement
        fo = interp1(isc.freq,isc.spl,f,'pchip');
        phi = zeros(size(fo));
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
        
    else % the stimulus is broadband impulsive create an inverse magnitude FIR filter  ---------------------
        tau = 0.0015; % filter group delay (sec)
        L = round(tau*2 * isc.fs); % calculate filter order
        if mod(L,2)~=1 % ensure filter order is even
            L = L + 1;
        end
        M = L-1; % filter order
        M2 = M/2; % filter group delay in samples
        N = isc.fs; % number of samples in design (using FIR from the "window" method)
        F = (0:1:N-1)'*(isc.fs/N); % make a frequency axis corresponding to X
        f1 = min(f); % minimum frequency to include (Hz)
        f2 = max(f); % maximum frequency to include (Hz)
        [~,indx1] = min(abs(F - f1)); % index of minimum frequency
        [~,indx2] = min(abs(F - f2)); % index of maximum frequency
        [~,indx4] = min(abs(F-(isc.fs - F(indx1)))); % the locations of the alised portions of the spectrum
        [~,indx3] = min(abs(F-(isc.fs - F(indx2))));
        Fcut = F(indx1:indx2); % cut frequency vector down to size
        
        fo = interp1(isc.freq,isc.spl,Fcut,'pchip'); % this is full out magnitude in dB
        phi = zeros(size(fo));
        
        % calculate the maximum possible output as an impulse
        Mag = zeros(isc.fs,1); % initialize magnitude spectrum
        Phase = zeros(isc.fs,1); % initialize phase spectrum
        bar = min(fo); % level of constant output (dB SPL)
        foInv = bar - fo; % inverse full output
        foInv = foInv + bar; % put back around min value
        target = foInv; % (for plotting purposes only)
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
        
        maxOut = bar; % maximum output in dB type (type = spl, rpl, fpl, or ipl)
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
        scaling = a * k;
        
        errors = 0;
        stimScaled = fastFilter(scaling,stim);
    end
end

% INTERNAL FUNCTIONS ------------------------------------------------------

% OLD CODE ----------------------------------------------------------------
%         A = abs(fft(a,isc.fs)); % magnitude of impulse, zero padded
%                 %A = A / (L/2); %A = A / (isc.fs/2); % scaled magnitude --> DO NOT SCALE FILTERS, DUMMY!
%         A = 20*log10(A/.00002); % fulloutput SPL of impulse
%         maxOut = bar; % maximum output in dB type (type = spl, rpl, fpl, or ipl)
%         d = targetLvl - maxOut; % calculate the dB difference between full out and desired (target) level
%         k = 10^(d/20); % express dB difference as a linear multiplier
%             scaling = a * k;
%            a = [a;0]; % make order ODD (don't want fast filter to apply group delay corection here!)


%         if k >= 1 % if cannot achieve target output because desired > than full out
%             errors = d; % error is the dB amount that exceeds full out
%             disp('Error in ARLas_applyISC.m! Line 161')
%             figure
%             plot(Fcut,fo,'r')
%             hold on
%             plt(F,A)
%             plt(Fcut,foInv,'g')
%             xlabel('Frequency (Hz)')
%             ylabel('Magnitude (dB type)')
%             title('Output Calibration')
%             legend('full out','inverse filter','target')
%             xlim(f)
%             %keyboard
%             %return
%         else

%         elseif strcmp(type,'pr') % make flat absorbance (still working throught this one!!)
%             fo = interp1(isc.freq,isc.fpl,Fcut,'pchip');
%             phi = interp1(isc.freq,isc.fplPhi,Fcut,'pchip');
%             pr = interp1(isc.freq,isc.pr,Fcut,'pchip');

                
        % NN = length(stim);
        % fs = isc.fs;
        % FF = (0:1:NN-1)'*(fs/NN);
        % [~,indx] = min(abs(FF-f));
        % X = fft(stim);
        % mag = abs(X);
        % mag = mag * k;
        % phase = angle(X);
        % phase(indx) = phase(indx) - phi;
        % phase = angle(cos(phase)+1i*sin(phase));
        % phase(NN-indx+2) = -phase(indx);
        % stimScaled = real(ifft(mag.*cos(phase)+1i*mag.*sin(phase)));
        % scaling = k*cos(phi)+1i*k*sin(phi);

%         if strcmp(type,'pr') % note: this is not correct!! just trying out ideas!!
%             foInv = foInv.*pr;
%         end


%         bar = mean(fo); % mean full output value
%         foInv = bar - fo; % inverse full output
%         foInv = foInv + bar; % put back around mean value
%         %foInv = foInv - min(foInv); % make minumum value = 0 dB -- why??
%         foInv = 10.^(foInv/20);     % put into linear units
%         [~,indx4] = min(abs(F-(isc.fs - F(indx1)))); % the locations of the alised portions of the spectrum
%         [~,indx3] = min(abs(F-(isc.fs - F(indx2))));
%         Mag(indx1:indx2) = foInv;
%         Mag(indx3:indx4) = flipud(foInv);
%         Phase = zeros(size(Mag));
%         X = Mag .* cos(Phase) + 1i*Mag.*sin(Phase);
%         x = real(ifft(X));
% 
%         L = round(tau*2 * isc.fs);
%         if mod(L,2)~=1
%             L = L + 1;
%         end
%         M = L-1;
%         M2 = M/2;
%         x = [x(end-M2:end);x(1:M2)];
%         h = hann(length(x));
%         a = x .* h;
%
%        stimScaled = fastFilter(a,stim);
%         A = A(indx1:indx2); % cut to to set of desired frequencies
%         reduction = 1./max(abs(a));
%         reduction_dB = 20*log10(reduction);
%         
%         a2 = a * reduction;        
%         A2 = abs(fft(a2,isc.fs));
%         A2 = A2 / (L/2); %A2 = A2 / (isc.fs/2);
%         A2 = 20*log10(A2/.00002);
%         A2 = A2(indx1:indx2); 
%         
%         maxOut = mean(A2);
%         ui = a2; % the unit impulse 
%         
%         fullCardOutput = 0.99;
%         maxOut = bar + 20*log10(fullCardOutput./max(abs(a)));


%         % double check the output to make sure no unwanted phase shifts are
%         % present
%         try
%             maxIterations = 3;
%             sampleOffset = 1;
%             counter = 1;
%             while sampleOffset ~= 0 && counter < maxIterations
%                 [stimScaled,sampleOffset] = adjustPhase(stim,stimScaled,isc);
%                 counter = counter + 1;
%             end
%             if sampleOffset ~= 0
%                 disp('Error: unable to correct group delay when applying ISC.')
%                 errors = 1;
%             end
%         catch
%             %keyboard
%         end
        
%         scaling = 1./max(abs(stimScaled));
%         dB_change = 20*log10(scaling);
%         stimScaled = stimScaled * scaling * 0.99;

%         keyboard
%         % make sure target output levels reached
%         d = targetLvl - fo; % calculate the dB difference between full out and desired (target) level
%         aa = 10.^(d/20); % express dB difference as a linear multiplier
%         if a >= 1 % if cannot achieve target output because desired > than full out
%             a = 0.99; % set multiplier to full out
%             errors = d; % error is the dB amount that exceeds full out
%         else
%             errors = 0;
%         end
%         stimScaled = stim * a;
% function [stimScaled,sampleOffset] = adjustPhase(stim,stimScaled,isc)
%     fs = isc.fs;
%     nfft = length(stim);
%     X1 = fft(stim,nfft);
%     X2 = fft(stimScaled,nfft);
%     F = (0:1:nfft-1)'*(fs/nfft); % make a frequency axis corresponding to X
%     [~,indx1] = min(abs(F-isc.fmin));
%     [~,indx2] = min(abs(F-isc.fmax));
%     X1 = X1(indx1:indx2);
%     X2 = X2(indx1:indx2);
%     F = F(indx1:indx2);
%     phi1 = angle(X1);
%     phi2 = angle(X2);
%     dphi = unwrap(phi2 - phi1);
%     omega = 2*pi*F;
%     order = 1;
%     warning off
%     p = polyfit(omega,dphi,order);
%     pd = polyder(p);
%     slope = -polyval(pd,omega);
%     warning on
%     slope = mean(slope);
%     sampleOffset = round(slope * fs);
%     if sampleOffset ~= 0 % if delay is changed, try to adjust time to give correct delay
%         if sampleOffset > 0 % if delay is positive
%             stimScaled = [stimScaled(sampleOffset+1:end);zeros(sampleOffset,1)]; % shift to left
%         else
%             stimScaled = [zeros(sampleOffset,1);stimScaled(sampleOffset+1:end-sampleOffset+1)]; % shift right
%         end
%     end
% end

        % if doSmoothing == 1 % added 8/12/2024 -- ssg -- can turn on smoothing to compensate for noisy calibration
        %     smHz = 100; % apply smoothing in 2 Hz steps
        %     smPts = round(smHz / median(gradient(f))); % smoothing points
        %     if smPts > 1
        %         fo = meanSmoother(fo,smPts); % smooth fo
        %     end
        % end
