function [isc,z0] = calculateFPL(obj,recording,stimulus)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [isc,z0] = calculateFPL(obj,recording,stimulus);
%
% For use with ARLas (Auditory Research Laboratory auditory software)
% Used in applying Thevenin source calibration to ear canal (or other cavity)recordings.
%
% obj = object of class thevCal (this contains the Thevenin source calibration)
% recording = sound pressure level recording made in the ear canal
% stimulus = the electrical stimulus used to make the recording
% isc = a structur containing spl, fpl, rpl, ipl, pr, etc.
% z0 = the estimated surge (characteristic) impedance
%
% Auditory Research Lab, The University of Iowa
% Deptartment of Communication Sciences & Disorders
% The University of Iowa
% Author: Shawn Goodman
% Date: July 7, 2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% still need to check for zero padding here------------------------------------------------------
N = length(stimulus);
nfft = obj.nfft;
if N > nfft
    error('Experiment stimulus is longer than Thevenin calibration stimulus. Cannot compute FPL.')
end

S = fft(stimulus,nfft); % stimulus matrix
R = fft(recording,nfft); % recordings matrix
PL = R ./ S; % load pressure
PL = PL(1:obj.fmaxIndx);
    
ZL = (obj.ZS .* PL) ./ (obj.PS - PL); % load impedance

% Use "surge" method to estimate cavity diameter:
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

FPL = (PL/2) .* (1 + (z0 ./ ZL)); % complex FPL
RPL = (PL/2) .* (1 - (z0 ./ ZL)); % complex RPL
PR = RPL ./ FPL; % complex pressure reflectance        
spl = 20*log10(abs(PL)); % load pressure magnitude (dB) (SPL)
fpl = 20*log10(abs(FPL)); % magnitude (dB) fpl
rpl = 20*log10(abs(RPL)); % magnitude (dB) rpl
ipl = fpl + rpl; % integrated pressure level
pr = abs(RPL).^2 ./ abs(FPL).^2; % power reflectance

% make transfer functions
filterGroupDelay = .0025; % desired filter group delay in seconds
M = round(filterGroupDelay * obj.fs); % filter order
if mod(M,2) ~= 0 % force filter order to be even
    M = M + 1;
end
H_spl = spl(obj.fminIndx:end);
H_fpl = fpl(obj.fminIndx:end);
H_rpl = rpl(obj.fminIndx:end);
H_ipl = ipl(obj.fminIndx:end);
H_pr = pr(obj.fminIndx:end);
try
    h_spl = makeImpulseResponse(-H_spl,N,obj.fmin,obj.fmax,M,obj.fs);
    h_fpl = makeImpulseResponse(-H_fpl,N,obj.fmin,obj.fmax,M,obj.fs);
    h_rpl = makeImpulseResponse(-H_rpl,N,obj.fmin,obj.fmax,M,obj.fs);
    h_ipl = makeImpulseResponse(-H_ipl,N,obj.fmin,obj.fmax,M,obj.fs);
end
isc.H_spl = -H_spl;  
isc.H_fpl = -H_fpl;  
isc.H_rpl = -H_rpl;  
isc.H_ipl = -H_ipl;
try
    isc.h_spl = h_spl;  
    isc.h_fpl = h_fpl;  
    isc.h_rpl = h_rpl;  
    isc.h_ipl = h_ipl;
end
isc.H_freq = obj.freq(obj.fminIndx:end);
isc.freq = obj.freq;
isc.H_pr = H_pr;
isc.SPL = PL;
isc.FPL = FPL; % forward pressure level
isc.RPL = RPL; % reverse pressure level
isc.PR = PR; % pressure reflectance
isc.zl = ZL;
isc.spl = spl;
isc.fpl = fpl; % forward pressure level (mag dB)
isc.rpl = rpl; % reverse pressure level (mag dB)
isc.ipl = ipl;
isc.prc = pr; % power reflectance
isc.z0 = z0; % calculated characteristic (surge) impedance
isc.diam = diam; % cavity diamter estimated from surge impedance

% internal functions ------------------------------------------------------

  function [H] = makeImpulseResponse(refGain,N,fmin,fmax,M,samplingRate)
    refGain = 10.^(refGain/20); % put refGain into linear units
    ff = (0:1:N-1)'*(samplingRate/N); % create a frequency vector (Hz)
    [dummy,indxStart] = min(abs(ff-fmin)); % find the index corresonding to start frequency
    [dummy,indxFinish] = min(abs(ff-fmax)); % find the index corresonding to end frequency
    H = zeros(M+1,size(refGain,2));
    for ii=1:size(refGain,2)
        RefGain = zeros(size(ff));
        RefGain(indxStart:indxFinish) = refGain(:,ii);
        RefGain = RefGain(1:floor(N/2)+1);
        Mag = [RefGain;flipud(RefGain(2:end-1))];
        Phi = zeros(size(Mag));
        RefGain = (Mag.*cos(Phi) + j*Mag.*sin(Phi));
        h = (real(ifft(RefGain))); % convert to time domain
        designN = length(h);
        h = [h(designN/2+1:end);h(1:designN/2)]; % shift the impulse
        %N = M+1; % number of samples in the filter
        n = (0:(1/M):1)'; % ... make a hann window
        w = 0.5 * (1-cos(2*pi*n)); % ...hann window
        middle = (designN/2) + 1; 
        halfw = (length(w)-1)/2;
        h = h(middle-halfw:middle+halfw); % truncate
        H(:,ii) = h .* w; % apply the window to create the new correction filter;
    end
    
function diam = getDiameter(z0)
a = 1.219;  
b = -0.03033 ; 
c = 0.8667 ;
d = -0.002683 ; 
diam = a*exp(b*z0) + c*exp(d*z0);