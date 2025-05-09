function [stim,lpc,reduction,errors] = ARLas_applyISC_noise(isc,targetLvl,type,f,stimLen,fs)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [stim,lpc,reduction,errors] = ARLas_applyISC_noise(isc,targetLvl,type,f,stimLen,fs);
%
% Apply in-situ calibration to create bandpass noise.
% Here, f is a vector of frequencies to include.
% 
% INPUT ARGUMENTS:
% isc = in-situ calibration structure obtained using ARLas_inSituCal.m
% targetLvl = desired stimulus level in dB (dB SPL, dB FPL, etc.)
% type = type of calibration to apply ('spl','fpl','rpl',ipl')
% f = dominant frequency of stimulus (Hz) 
%       For broadband stimuli, give all of the frequencies to include.
%       Amplitude will be calibrated; phase will be assigned randomly
% stim = stimulus of unit amplitude (i.e. full output)
%
% OUTPUT ARGUMENTS:
% stimScaled = stimulus scaled to result in achieving target level and zero phase.
% scaling = multiplier used to scale stimulus and achieve zero phase.
%       If stimulus is narrowband, scaling is a scalar multiplier (complex)
%       If stimulus is broadband, scaling is the impulse response (FIR) used to
%          filter stim to yield stimScaled.
% errors = report of whether target can be resached (1) or not (0).
%
% Auditory Research Lab, The University of Iowa
% Deptartment of Communication Sciences & Disorders
% The University of Iowa
% Author: Shawn Goodman
% Date: July 7, 2021
% Last Updated: July 7, 2021 -- ssg -- added capability for swept tones
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



    % get the full-out amplitude from the in-situ cal measurement
    if strcmp(type,'spl')
        fo = interp1(isc.freq,isc.spl,f,'pchip');
        phi = interp1(isc.freq,isc.splPhi,f,'pchip');
    elseif strcmp(type,'fpl')
        fo = interp1(isc.freq,isc.fpl,f,'pchip');
        phi = interp1(isc.freq,isc.fplPhi,f,'pchip');
    elseif strcmp(type,'rpl')
        fo = interp1(isc.freq,isc.rpl,f,'pchip');
        phi = interp1(isc.freq,isc.rplPhi,f,'pchip');
    elseif strcmp(type,'ipl')
        fo = interp1(isc.freq,isc.ipl,f,'pchip');
        phi = interp1(isc.freq,isc.fplPhi,f,'pchip');
        % ipl technically has no phase; here, we use fpl phase.
    else
        error('Unrecognized calibration type. Must be a string: spl, fpl, rpl, or ipl.')
    end        
    
    
    if length(targetLvl)==1
        origTargetLvl = targetLvl;
        targetLvl = ones(size(f))*targetLvl;
    else
        origTargetLvl = [];
    end
    d = targetLvl - fo; % calculate the dB difference between full out and desired (target) level
    k = 10.^(d/20); % express dB difference as a linear multiplier
    if k >= 1 % if cannot achieve target output because desired > than full out
        k = 0; % set multiplier to zero
        errors = d; % error is the dB amount that exceeds full out
        warning('The requested dB value exceeds full out.')
    else
        errors = 0;
    end
    
    nSamples = round(stimLen * fs);
    t = (0:1:nSamples-1)'/fs;
    nFreqs = length(f);
    p = rand(nFreqs,1)*(2*pi);
    x = zeros(nSamples,1);
    for ii=1:nFreqs
        dummy = k(ii)*cos(2*pi*f(ii)*t + p(ii));
        x = x + dummy;
    end
    %reduction = 1 ./ (3*std(x));
    reduction = 10^(-40/20);
    x = x * reduction;
    [stim,errors] = lowNoiseNoise(x);
    
    if ~isempty(origTargetLvl)
        targetLvl = origTargetLvl;
    end
    lpc = targetLvl + 20*log10(reduction);
    reduction = 20*log10(reduction);
    
    % In shawn's ear, this full-scaled output gave a noise amplitude of 82
    % dB SPL rms. Shawn has a small ear, so this level would generally be
    % smaller, say up to 10 dB smaller? 
    % Say then that maximum output of the system is ~75 dB SPL rms.
    
    
    
end

% INTERNAL FUNCTIONS ------------------------------------------------------
function [noise,errors] = lowNoiseNoise(noise)
    errors = 0;
    [Rows,Cols] = size(noise);
    noise = noise(:);
    N  = length(noise);
    noiseMax = 1;
    indx = find(abs(noise) >= noiseMax); % find the samples that exceed the max
    
    criterion = 0.1;
    if (length(indx)/N) > criterion
        disp(['WARNING: more than ',num2str(criterion*100),' % of noise samples exceed amplitude of 1.'])
        errors = 1;
    end
    
    replacement = (rand(length(indx),1)*2)-1;
    noise(indx) = replacement;
    noise = noise * .995; 
    
    noise = reshape(noise,Rows,Cols);
end


% 
% % NOTES:
% dBp2p - dBbw = dBlpc
% p2p / bw = lpc
% 
% p2p = 1.4572;
% dBp2p = 20*log10(p2p/.00002);
% bw = 16000 - 100;
% dBbw = 10*log10(bw); % 42
% dBlpc = dBp2p - dBbw; % 55.2358
% 
% ft = abs(fft(clicksImpulse));
% ft = ft / (length(ft)/2);
% ft = mean(ft(5:50));
% dBft = 20*log10(ft/.00002); % 60
% 
% rms = mean(sqrt(clicksImpulse.^2));
% dBrms = 20*log10(rms/.00002); % 68
% 


% So if fo was 55 dB SPL, then
% working backwards,
% peak level would be roughly 
% dBp2p = dBlpc + dBbw
% dBp2p = 55    +  42 = 97 dB

% if fo was 110, then p2p would be 110 + 42 = 152.
% peak level is half that (assuming vertical symmetry, which isn't quite
% true), but say -6 dB, so peak is 146 dB.
% Say that yields a click with amplitude of 18. Then we would have to
% divide by 18 to make a full output, which would reduce