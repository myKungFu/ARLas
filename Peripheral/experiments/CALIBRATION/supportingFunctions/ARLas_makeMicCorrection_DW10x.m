function [b] = ARLas_makeMicCorrection_DW10x(fplWaveform,refWaveform,fmin,fmax,fs)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% b = ARLas_makeMicCorrection_DW10x(fplWaveform,refWaveform,fmin,fmax,fs);
%
% Create microphone correction for the ER10X probe system. 
% This version is specifically for use with guinea pigs and an ear bar.
% Assuming frequency range from 200 - 32000 Hz. 
% Assuming measurements made in the custom calibration assembly:
%       ER10X probe -- ear bar -- brass coupler -- 1/8" condensor mic
%
% Input Arguments:
%  fplWaveform = calculated forward pressure level in the calibration assembly.
%  refWaveform = reference waveform measured by the condensor mic.
%  fmin = minimum frequency to include (200 Hz)
%  fmax = maximum frequency to include (32000 Hz)
%  fs = sampling rate (Hz; assuming 96000)
% 
% Authors: Shawn S Goodman, PhD & Jeffery T Lichtenhan PhD
% Auditory Research Lab, Dept. of Comm. Sciences & Disorders, The University of Iowa
% Auditory Physiology Lab, Dept. of Otolaryngology, Washington University School of Medicine, St. Louis
% Date: October 12, 2019
% Updated: October 12-22, 2019 -- ssg
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% --------------------- Adjustable Parameters -----------------------------
    doSmoothing = 0;
    doPhaseAdjust = 1;
    filterGroupDelay = 0.0015; % desired filter group delay (s)
                              % for DW, 1 ms group delay puts you within 1 dB of target values
    smoothing = 0.000000001;  % spling smoothing value
%--------------------------------------------------------------------------    
    
    N = length(fplWaveform); % length if input waveform
    N2 = length(refWaveform);
    if N ~= N2 % ensure the two waveforms are the same size
        error('Input waveforms must be the same size.')
    end
    if mod(N,2)~=0 % ensure an even number of samples
        fplWaveform = [fplWaveform;0];
        refWaveform = [refWaveform;0];
        N = length(fplWaveform); 
    end
    
    FPL = fft(fplWaveform,N); % put into the frequency domain
    REF = fft(refWaveform,N);
    Xfer = FPL ./ REF; % microphone transfer function (fpl re: reference)
    freq = (0:1:N-1)' * (fs/N); % frequency vector
    [~,indxCut1] = min(abs(freq-fmin)); % cutoff indices
    [~,indxCut2] = min(abs(freq-fmax));
    mag = abs(Xfer);
    ang = angle(Xfer);
    
    if doSmoothing == 1
        w = ones(size(freq)); % create vector of weights for the smoothing spline
        w(1:indxCut1) = 0; % weight outside of the fitting zone to zeros...
        w(indxCut2:end) = 0; % this keeps away crazy at the edges--important!
        % Xfer is complex; smooth on mag and phase individually
        mag = csaps(freq,mag,smoothing,freq,w); % the smoothed, weighted spline, densly-spaced estimate
        ang = csaps(freq,unwrap(ang),smoothing,freq,w); % the smoothed, weighted spline
        ang = angle(cos(ang) + 1i*sin(ang)); % wrap the phase back up
    end

    if doPhaseAdjust == 1
        % adjust the phase by subtracting off constant group delay that
        % represents the travel time in the tube
        angStub = unwrap(ang(indxCut1:indxCut2));
        freqStub = freq(indxCut1:indxCut2);
        p = polyfit(freqStub,angStub,1); % straight-line fit is constant group delay
        gd = p(1)*freqStub + p(2); % create group delay line
        angStub = angStub - gd; % subtract off group delay. 
        ang(indxCut1:indxCut2) = angStub; % remainder is phase that needs to be corrected
    else
        ang = ang * 0; % if phase will not be corrected, then set to zero delay
    end
    
    % Invert the mic transfer function and make into an FIR filter
    mag = 1./mag; % invert magnitude (multiplicative inverse)
    ang = -ang; % invert phase (additive inverse)
    % regardless of what came before, set all values outside fmin and fmax to zero
    ang(1:indxCut1-1) = 0;
    ang(indxCut2+1:end) = 0;
    mag(1:indxCut1-1) = 0;
    mag(indxCut2+1:end) = 0;
    % ------------ CREATE FIR FILTER USING WINDOW METHOD ------------------
    nyquist = (N/2)+1; % nyquist sample
    % create the aliased "half" of the spectrum
    mag(nyquist+1:end) = flipud(mag(2:nyquist-1));
    ang(nyquist+1:end) = -flipud(ang(2:nyquist-1));
    Y = mag.*cos(ang) +1i*mag.*sin(ang); % complex rectangular form
    y = real(ifft(Y)); % time domain waveform
    yy = [y(nyquist:end);y(1:nyquist-1)]; % make filter causal
    m = round(filterGroupDelay * fs); % filter order
    if mod(m,2)~= 0 % force group delay to be even so that can remove it after filtering
        m = m + 1;
    end
    m2 = m/2; % half group delay
    
    yy = yy(nyquist-m2:nyquist+m2); % cut down filter to desired size
    h = hann(length(yy)); % window to smooth edges
    b = yy .* h; % filter coefficients in the time domain.

end