function [isc] = ARLas_inSituCal(obj,recording,stimulus)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [isc] = ARLas_inSituCal(obj,recording,stimulus);
%
% For use with ARLas (Auditory Research Laboratory auditory software)
% Used in applying Thevenin source calibration to ear canal (or other cavity)recordings.
%
% obj = object of class thevCal (this contains the Thevenin source calibration)
% recording = sound pressure level recording made in the ear canal
% stimulus = the electrical stimulus used to make the recording
%
% isc is a structure containing the following fields:
% isc.fs = sampling rate (Hz)
% isc.fmin = minimum frequency of calibration
% isc.fmax = maximum frequency of calibration
% isc.Amax = maximum output level of system
% isc.nfft = number of samples in the fft
% isc.PS = source pressure
% isc.ZS = source impedance
% isc.PL = measured load pressure (Pa)
% isc.ZL = calculated load impedance (ohms)
% isc.z0 = calculated characteristic (surge) impedance of load (ohms)
% isc.diam = cavity diameter estimated from surge impedance (cm)
% isc.pr = power relectance
% isc.length = cavity length, in cm, estimated from rpl/fpl calculations
% isc.length2 = cavity length, in cm, estimated directly from null
% isc.FPL = forward pressure level (complex) using full scaled output stimulus
% isc.RPL = reverse pressure level (complex) using full-scaled output stimulus
% isc.PRR = pressure reflectance (complex)
% isc.spl = maximum sound pressure level magnitude (dB re: 20 uPa)
% isc.fpl = maximum forward pressure level magnitude (dB re: 20 uPa)
% isc.rpl = maximum reverse pressure level magnitude (dB re: 20 uPa)
% isc.ipl = maximum integrated pressure level magnitude (dB re: 20 uPa)
% isc.freq = frequency vector associated with the above values (in Hz)
% isc.splPhi = unwrapped phase for spl
% isc.fplPhi = unwrapped phase for fpl
% isc.rplPhi = % unwrapped phase for rpl
% isc.eplCorrection = correction; % epl correction (multiply SPL by this)
%
% Auditory Research Lab, The University of Iowa
% Deptartment of Communication Sciences & Disorders
% The University of Iowa
% Author: Shawn Goodman
% Date: December 1, 2017
% Updated: October 16, 2018
% Updated: July 9, 2019 -- ssg & sb; added phase outputs
% Updated: May 21, 2025 -- ssg -- added epl correction as an output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N = length(stimulus);
nfft = obj.nfft;
if N > nfft
    error('Calibration stimulus is longer than Thevenin calibration stimulus. Cannot compute values.')
end

S = fft(stimulus,nfft); % stimulus vector
R = fft(recording,nfft); % recordings vector
PL = R ./ S; % load pressure
PL = PL(1:obj.fmaxIndx); % cut to highest calibrated frequency
PL = PL(obj.fminIndx:end);
ZL = (obj.ZS .* PL) ./ (obj.PS - PL); % calculate load impedance
[z0,diam] = getZ0(obj,ZL); % Use "surge" method to estimate surge impedance and cavity diameter:
freq = obj.freq; % get frequency vector
FPL = (PL/2) .* (1 + (z0 ./ ZL)); % complex FPL
RPL = (PL/2) .* (1 - (z0 ./ ZL)); % complex RPL
PRR = RPL ./ FPL;                 % complex pressure reflectance PR
pr = abs(RPL).^2 ./ abs(FPL).^2; % power reflectance

Amax = 1; % system maximum output (usually 1)
pRef = 0.00002; % pressure reference (20 micropascals)
spl = 20*log10(Amax*abs(PL)/pRef);  % load sound pressure magnitude (dB spl)
fpl = 20*log10(Amax*abs(FPL)/pRef); % load forward pressure magnitude (dB fpl)
rpl = 20*log10(Amax*abs(RPL)/pRef); % load reverse pressure magnitude (dB rpl)
ipl = 20*log10((Amax*abs(FPL) + Amax*abs(RPL))/pRef); % load integrated pressure level (dB ipl)
splPhi = unwrap(angle(PL)); % phase; added 7/9/2019 Shawn & Sriram
fplPhi = unwrap(angle(FPL));
rplPhi = unwrap(angle(RPL));

c =  3.4933386657e+04; % speed of sound in air
[correction,canalLength,Fqw,RS,RL] = getEPL(z0,obj.ZS,ZL,ipl,spl,freq,c,diam); % canal length here is estimated from FPL and RPL

isc.fs = obj.fs;     % sampling rate (Hz)
isc.fmin = obj.fmin; % minimum frequency of calibration
isc.fmax = obj.fmax; % maximum frequency of calibration
isc.Amax = Amax;     % maximum output level of system
isc.nfft = obj.nfft; % number of samples in the fft
isc.PS = obj.PS;     % source pressure
isc.ZS = obj.ZS;     % source impedance
isc.PL = PL;         % measured load pressure (Pa)
isc.ZL = ZL;         % calculated load impedance (ohms)
isc.z0 = z0;         % calculated characteristic (surge) impedance of load (ohms)
isc.diam = diam;     % cavity diameter estimated from surge impedance (cm)
isc.length = canalLength; % cavity length, in cm, estimated from rpl/fpl calculations
isc.pr = pr;         % power relectance
isc.FPL = FPL;       % forward pressure level (complex) using full scaled output stimulus
isc.RPL = RPL;       % reverse pressure level (complex) using full-scaled output stimulus
isc.PRR = PRR;       % pressure reflectance (complex)
isc.spl = spl;       % maximum sound pressure level magnitude (dB re: 20 uPa)
isc.fpl = fpl;       % maximum forward pressure level magnitude (dB re: 20 uPa)
isc.rpl = rpl;       % maximum reverse pressure level magnitude (dB re: 20 uPa)
isc.ipl = ipl;       % maximum integrated pressure level magnitude (dB re: 20 uPa)
isc.freq = obj.freq; % frequency vector associated with the above values (in Hz)
isc.Fqw = Fqw;       % location of quarter wavelength null (Hz)
isc.RS = RS;         % Source reflectance
isc.RL = RL;         % Load reflectance
isc.splPhi = splPhi; % unwrapped phase for spl
isc.fplPhi = fplPhi; % unwrapped phase for fpl
isc.rplPhi = rplPhi; % unwrapped phase for rpl
isc.eplCorrection = correction; % epl correction (multiply SPL by this)

end

% internal functions ------------------------------------------------------
function [z0,diam] = getZ0(obj,ZL)
    % Use "surge" method to estimate surge impedance and cavity diameter:
    minFreq = obj.fmin * 1.25; % use the lowest possible frequency, plus 25%
    [~,indxMin] = min(abs(obj.freq-minFreq));
    maxFreq = obj.fmax;
    [~,indxMax] = min(abs(obj.freq-maxFreq));
    winN = (indxMax-indxMin)+1;
    w = blackman(winN*2);
    w = w(winN+1:end);
    pad = zeros(indxMin-1,1);
    w = [pad;w];
    % z0 = mean(real(ZL(indxMin:indxMax)));
    z0 = abs(sum(w.*(real(ZL)))/sum(w));
    nIterations = 8;
    %figure(13); hold on
    %nn = (1:nIterations)';
    %zz = zeros(nIterations,1);
    for kk=1:nIterations % originally 4; changed to 8 per Steve and Daniel
        wr = sum(w.*(real(((ZL - z0)./(ZL + z0)))))/sum(w);
        z0 = z0 * (1 + wr);
        % zz(kk,1) = z0;
        % plot(nn(1:kk),zz(1:kk),'g-o')
        % pause(.1)
    end
    diam = getDiameter(z0);
    % figure(113)
    % plot(real(ZL),'b')
    % hold on
    % plot(imag(ZL),'r')
    % xlim([54,400])
    % pause(.5)
end
function diam = getDiameter(z0)
    % get cavity diameter, based on calcualted surge impedance (z0)
    a = 1.219;  
    b = -0.03033 ; 
    c = 0.8667 ;
    d = -0.002683 ; 
    diam = a*exp(b*z0) + c*exp(d*z0);
end

% 
% [isc.splmax,isc.SPLmax,isc.maxFreq,isc.splFlatmax] = makeImpulseResp(PL,Amax,fmin,fmax,fs,nfft);
% [isc.fplmax,isc.FPLmax,~,isc.fplFlatmax] = makeImpulseResp(FPL,Amax,fmin,fmax,fs,nfft);
% [isc.rplmax,isc.RPLmax,~,isc.rplFlatmax] = makeImpulseResp(RPL,Amax,fmin,fmax,fs,nfft);
% [isc.iplmax,isc.IPLmax,~,isc.iplFlatmax] = makeImpulseResp(abs(FPL) + abs(RPL),Amax,fmin,fmax,fs,nfft);

% function [xmax,Xmax,freq,flatMax] = makeImpulseResp(X,Amax,fmin,fmax,fs,nfft)
%     % X is the magnitude of the obtained value, re: stimulus, in linear units
%     %X = 10.^(X/20); % put into into linear units
%     %Amax = 1; 			% maximum amplitude of digital output
%     Xmax = Amax * X;	% maximum output
%     Xmax = abs(Xmax);  	% magnitude
%     H = -Xmax;          % transfer function
%     freq = (0:1:nfft-1)'*(fs/nfft); % frequency vector (Hz)
%     [~,indxMin] = min(abs(freq-fmin)); % location of frequency min and max
%     [~,indxMax] = min(abs(freq-fmax));
%     [~,indxNyquist] = min(abs(freq-(fs/2))); % Nyquist frequency location
%     %H = H(1:indxNyquist); % cut to include only un-aliased parts
%     freq = freq(1:indxNyquist);
%     H(1:indxMin) = 0; % zero frequencies outside of fmin and fmax
%     %H(indxMax:end) = 0;
%     H = [H;zeros(indxNyquist-length(H),1)];
%     M = round(.02 * fs); % filter order is 2% of the sampling rate
%     if mod(M,2)~= 0 % filter order must be an even number
%         M = M + 1;
%     end
%     pRef = 0.0002; % 20 uPa is pressure reference
%     n = length(H);
%     X = [H;flipud(H(2:end-1))]; % add aliased portion of the DFT
%     x = ifft(X); % create in impulse response
%     N = length(x);
%     x = [x(N/2+2:end);x(1:N/2+1)]; % make filter causal
%     x = x((N/2)-(M/2):(N/2)+(M/2)); % take only the desired center portion
%     w = hann(M+1); % create hann window
%     xmax = x .* w; % window to smooth edges; This is the impulse response
%     Xmax = fft(xmax,N); % put back into the frequency domain
%     pRef = 0.00002;  % pressure reference: 20 uPa
%     Xmax = 20*log10(abs(Xmax)/pRef); % convert to dB SPL
%     Xmax = Xmax(indxMin:indxMax);
%     freq = (0:1:N-1)'*(fs/N); % frequency vector (Hz)
%     freq = freq(indxMin:indxMax);
%     flatMax = min(Xmax); % maximum output if the stimulus is broadband flattened
% end
% 
% function [H] = makeImpulseResponse(refGain,N,fmin,fmax,M,samplingRate)
%     refGain = 10.^(refGain/20); % put refGain into linear units
%     ff = (0:1:N-1)'*(samplingRate/N); % create a frequency vector (Hz)
%     [dummy,indxStart] = min(abs(ff-fmin)); % find the index corresonding to start frequency
%     [dummy,indxFinish] = min(abs(ff-fmax)); % find the index corresonding to end frequency
%     H = zeros(M+1,size(refGain,2));
%     for ii=1:size(refGain,2)
%         RefGain = zeros(size(ff));
%         RefGain(indxStart:indxFinish) = refGain(:,ii);
%         RefGain = RefGain(1:floor(N/2)+1);
%         Mag = [RefGain;flipud(RefGain(2:end-1))];
%         Phi = zeros(size(Mag));
%         RefGain = (Mag.*cos(Phi) + j*Mag.*sin(Phi));
%         h = (real(ifft(RefGain))); % convert to time domain
%         designN = length(h);
%         h = [h(designN/2+1:end);h(1:designN/2)]; % shift the impulse
%         %N = M+1; % number of samples in the filter
%         n = (0:(1/M):1)'; % ... make a hann window
%         w = 0.5 * (1-cos(2*pi*n)); % ...hann window
%         middle = (designN/2) + 1; 
%         halfw = (length(w)-1)/2;
%         h = h(middle-halfw:middle+halfw); % truncate
%         H(:,ii) = h .* w; % apply the window to create the new correction filter;
%     end
% end
% 


