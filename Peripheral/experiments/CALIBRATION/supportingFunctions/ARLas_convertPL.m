function [pl,Pl,phi,other,wf] = ARLas_convertPL(X,isc)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [pl,PL,phi,other,wf] = ARLas_convertPL(x,isc);
%
% Convert pressure level to various other forms (fpl, etc.)
%
% x = input sound pressure waveform (column vector of amplitudes in Pa)
% isc = in-situ calibration structure; obtained from ARLas_inSituCal.m
%
% function returns four structures:
%   1) pl = pressure magnitudes in dB. Includes the following:
%           spl, fpl, rpl, ipl,epl, pr (pressure reflectance), pr (power
%           reflectance), associated noise floors, and frequency.
%   2) PL = complex pressure values. Includes the following:
%           SPL, FPL, RPL, EPL, PRR (pressure reflectance, and frequency.
%   3) phi = associated phases in radians
%   4) other = other related values: ZL, ZS, z0, freq, RS, RL, etc.
%   5) wf = waveforms of original spl and epl, along with time vector
%
% Auditory Research Lab, The University of Iowa
% Deptartment of Communication Sciences & Disorders
% The University of Iowa
% Author: Shawn Goodman
% Date: December 1, 2017
% Updated: October 16, 2018
% Updated: March 15, 2019 -- ssg; added time waveforms
% Updated: July 9, 2019 -- ssg & sb; fixed error relating to long input signals
% Updated: March 30, 2020 -- ssg; changed smoothing from a weighted spline
%                            to a simple mean smoother.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ZL = isc.ZL; % load impedance
    ZS = isc.ZS; % source impedance
    z0 = isc.z0; % characteristic (surge) impedance
    freq = isc.freq; % associated frequency vector
    fmax = isc.fmax; % highest calibrated frequency
    fmin = isc.fmin; % lowest calibrated frequency
    if fmax> freq(end)
        fmax = freq(end);
    end
    if fmin< freq(1)
        fmin = freq(1);
    end
    
    fs = isc.fs; % sampling rate
    % adjust fft size so that the largest of x and isc is used
    if size(X,1) > isc.nfft % if x is longer than the chirp cal
        nfft = size(X,1);   % no zero padding will be added to x, but cal must be interpolated
        zpadN = 0; % number of samples in the zero padding
        freqNew = (0:1:nfft-1)'*(fs/nfft);
        % [~,fminIndx] = min(abs(freqNew-freq(1))); removed and replaced
        % with following line 9/7/2022 ssg
        [~,fminIndx] = min(abs(freqNew-fmin));
        [~,fmaxIndx] = min(abs(freqNew-fmax));
        freqNew = freqNew(fminIndx:fmaxIndx);
        ZL = interp1(freq,ZL,freqNew,'pchip');
        ZS = interp1(freq,ZS,freqNew,'pchip');
        freq = freqNew;
    elseif size(X,1) <= isc.nfft % if x is shorter than the chirp cal
        nfft = isc.nfft; % zero padding will be added to x, but not to cal
        zpadN = nfft - size(X,1); % number of samples of zero padding added
    end
    nyquist = (nfft/2)+1;

    SPL = fft(X,nfft); % measured sound presure level
    frequency = (0:1:nfft-1)'*(isc.fs/nfft);
    [~,fminIndx] = min(abs(frequency-fmin));
    [~,fmaxIndx] = min(abs(frequency-fmax));
    
    spl_orig = X; % orignal time waveform
    time = (0:1:length(X)-1)'/isc.fs; % associated time vector (s)

    scaling = ((nfft-zpadN)/2); % fft scaling factor
    SPL = SPL / scaling; % scale to give correct output
    SPL = SPL(fminIndx:fmaxIndx,:); % cut to highest calibrated frequency
    SPLraw = SPL; 
    
    ZLmatrix = repmat(ZL,1,size(SPL,2));
    FPL = (SPL/2) .* (1 + (z0 ./ ZLmatrix)); % complex FPL
    RPL = (SPL/2) .* (1 - (z0 ./ ZLmatrix)); % complex RPL
    PRR = RPL ./ FPL;                 % complex pressure reflectance PR
    IPL = abs(FPL) + abs(RPL);        % integrated pressure level
    PR = abs(RPL).^2 ./ abs(FPL).^2;  % power reflectance

    [SPL,splnf] = getSEM(SPL);
    [FPL,fplnf] = getSEM(FPL);
    [RPL,rplnf] = getSEM(RPL);
    [IPL,iplnf] = getSEM(IPL);
    iplnf = iplnf * 2;
    [PRR,PRRnf] = getSEM(PRR);
    [PR,PRnf] = getSEM(PR);

% -----    
    FPL = smoother(freq,FPL,SPL);
    RPL = smoother(freq,RPL,SPL);
    IPL = abs(FPL) + abs(RPL);        % integrated pressure level
    
    IPL2 = IPL.*cos(angle(FPL)) + 1i*IPL.*sin(angle(FPL)); % ipl using fpl phase
    
    PR = abs(RPL).^2 ./ abs(FPL).^2;  % power reflectance
    ZL = smoother(freq,ZL,SPL);
% -----    
    
    pRef = 0.00002; % pressure reference (uPa)
    spl = 20*log10(abs(SPL)/pRef); % express as magnitude re: 20 uPa
    ipl = 20*log10(IPL/pRef);
    cc =  3.493338665700000e+04; % speed of sound in air
    [correction,canalLength,Fqw,RS,RL] = getEPL(z0,ZS,ZL,ipl,spl,freq,cc);
    EPL = mean(SPLraw,2) .* correction; % emission pressure level
     
% -----    
    EPL = smoother(freq,EPL,SPL);
% -----    
    
    % calculate waveforms ------------------------------------------------=
    pad1 = zeros(fminIndx,1);
    pad2 = zeros(nyquist-fmaxIndx,1);
    EPL_full = [pad1;EPL*scaling;pad2]; % emission pressure level
    EPL_full = [EPL_full(1:nyquist);flipud(conj(EPL_full(2:nyquist-1)))];
    epl_full = real(ifft(EPL_full));
    epl_full = epl_full(1:length(spl_orig));
    epl_full = ARLas_ramp(epl_full,fs,0.001);
    
    FPL_full = [pad1;FPL*scaling;pad2]; % forward pressure level
    FPL_full = [FPL_full(1:nyquist);flipud(conj(FPL_full(2:nyquist-1)))];
    fpl_full = real(ifft(FPL_full));
    fpl_full = fpl_full(1:length(spl_orig));
    fpl_full = ARLas_ramp(fpl_full,fs,0.001);

    IPL2_full = [pad1;IPL2*scaling;pad2]; % intergrated pressure level using forward phase
    IPL2_full = [IPL2_full(1:nyquist);flipud(conj(IPL2_full(2:nyquist-1)))];
    ipl2_full = real(ifft(IPL2_full));
    ipl2_full = ipl2_full(1:length(spl_orig));
    ipl2_full = ARLas_ramp(ipl2_full,fs,0.001);    

    IPL_full = [pad1;IPL*scaling;pad2]; % integrated pressure level
    IPL_full = [IPL_full(1:nyquist);flipud(conj(IPL_full(2:nyquist-1)))];
    ipl_full = real(ifft(IPL_full));
    nn = length(ipl_full);
    ipl_full = [ipl_full(nn/2:end);ipl_full(1:nn/2+1)];
    m2 = length(spl_orig)/2;
    center = length(ipl_full)/2;
    ipl_full = ipl_full(center-m2+1:center+m2);
    ipl_full = ARLas_ramp(ipl_full,fs,0.001);
    
    RPL_full = [pad1;RPL*scaling;pad2]; % reverse pressure level
    RPL_full = [RPL_full(1:nyquist);flipud(conj(RPL_full(2:nyquist-1)))];
    rpl_full = real(ifft(RPL_full));
    rpl_full = rpl_full(1:length(spl_orig));
    rpl_full = ARLas_ramp(rpl_full,fs,0.001);
    
    epl = 20*log10(abs(EPL)/pRef); % expressed as magnitude re: 20 uPa
    fpl = 20*log10(abs(FPL)/pRef);
    rpl = 20*log10(abs(RPL)/pRef);

    splnf = 20*log10(splnf/pRef); % express noise floors as magnitude re: 20 uPa
    iplnf = 20*log10(iplnf/pRef);
    fplnf = 20*log10(fplnf/pRef);
    rplnf = 20*log10(rplnf/pRef);
    % noise floor correction cannot be applied for epl

    % function returns five structures:
    % 1) complex values (PL)
    % 2) pressure magnitudes in dB (pl)
    % 3) phases in radians (phi)
    % 4) other values related to the calculations (other)
    % 5) time waveforms for selected calculations

    pl.spl = spl(fminIndx:end); % sound pressure level (dB: re 20 uPa)
    pl.fpl = fpl(fminIndx:end); % forward pressure level (dB: re 20 uPa)
    pl.rpl = rpl(fminIndx:end); % reverse pressure level (dB: re 20 uPa)
    pl.ipl = ipl(fminIndx:end); % integrated pressure level (dB: re 20 uPa)
    pl.epl = epl(fminIndx:end); % emission pressure level (dB: re 20 uPa)
    pl.pr = PR(fminIndx:end);   % power reflectance
    pl.splnf = splnf(fminIndx:end);
    pl.fplnf = fplnf(fminIndx:end);
    pl.rplnf = rplnf(fminIndx:end);
    pl.iplnf = iplnf(fminIndx:end);
    pl.f = freq(fminIndx:end);

    Pl.SPL = SPL(fminIndx:end);
    Pl.FPL = FPL(fminIndx:end);
    Pl.RPL = RPL(fminIndx:end);
    Pl.EPL = EPL(fminIndx:end);
    Pl.PRR = PRR(fminIndx:end);
    Pl.f = freq(fminIndx:end);

    phi.spl = angle(SPL(fminIndx:end));     % phase of spl (rad)
    phi.fpl = angle(FPL(fminIndx:end));     % phase of fpl
    phi.rpl = angle(RPL(fminIndx:end));     % phase of rpl
    phi.prr = angle(PRR(fminIndx:end));     % phse of pressure reflectance
    phi.epl = angle(EPL(fminIndx:end));     % phase of epl
    phi.f = freq(fminIndx:end);
    % power reflectance and ipl have no phase; being composed of sum or ratio of magnitudes

    other.ZL = ZL; % load impedance
    other.ZS = ZS; % cource impedance
    other.z0 = z0; % characteristic (surge) impedance
    other.freq = freq; % frequency vector (Hz)
    other.fmax = fmax; % highest calibrated frequency
    other.fmin = fmin; % lowest calibrated frequency
    other.correction = correction;
    other.canalLength = canalLength; 
    other.diam = isc.diam; % ear canal diameter
    other.Fqw = Fqw; % quarter-wavelength frequency (from EPL calibration)
    other.RS = RS; % source reflectance
    other.RL = RL; % load reflectance

    wf.time = time; % time waveform (seconds)
    wf.spl = spl_orig; % original sound pressure level recording
    wf.epl = epl_full; % epl corrected recording (use this for emissions only)
    wf.fpl = fpl_full; % fpl corrected recording (use this for stimulus only)
    wf.rpl = rpl_full; % rpl corrected recording (use this for stimulus only)
    wf.ipl = ipl_full; % this has a constant, but arbitrary group delay.
    wf.ipl2 = ipl2_full; % ipl corrected recording, using fpl phase
end

% internal functions ------------------------------------------------------
function [x,sem] = getSEM(X)
% Calculate the noise floor, which in the case of of multiple time-locked
% responses, is the standard error of the mean (sem)
    K = size(X,2);
    x = mean(X,2);
    if K <=1
        sem = [];
        return
    end
    XBAR = repmat(x,1,K);
    S2 = (1/(K-1)) * sum((X - XBAR) .* conj(X - XBAR),2); % variance
    sem = sqrt((1/K) * S2); % standard error of mean
end
function [y] = smoother(XX,YY,SPL)
    y = abs(YY);
    mag = meanSmoother(y,6);
    y = unwrap(angle(YY));
    ang = meanSmoother(y,6);
    ang = angle(cos(ang) + 1i*sin(ang));
    y = mag.*cos(ang) +1i*mag.*sin(ang);
end
function [y] = smoother2(XX,YY,SPL)
try
    smoothing = .000001;
    w = ones(size(XX));
    xx = XX;
    x = XX;

    q = 20*log10(abs(mean(SPL,2)));
    mu = mean(q);
    q = q - mu;
    w = ones(size(q));
    criterion = 0;
    w(find(q<criterion)) = 0;
    dw = diff(w);
    indxStart = find(dw<-.1) + 1;
    indxFinish = find(dw>.1);
    if length(indxFinish) > 1
        if indxFinish(1) < indxStart(1)
            indxFinish(1) = [];
        end
        if indxStart(end) > indxFinish(end)
            indxStart(end) = [];
        end
    else
        dummy = indxStart;
        indxStart = indxFinish;
        indxFinish = dummy;
    end
    
    w = ones(size(q));
    for ii=1:length(indxStart)
        n = indxFinish(ii) - indxStart(ii) + 1;
        h = -blackman(n) + 1;
        w(indxStart(ii):indxFinish(ii)) = h;
    end
    
    y = abs(YY);
    mag = csaps(x,y,smoothing,xx,w); % the smoothed, weighted spline, densly-spaced estimate
    y = unwrap(angle(YY));
    ang = csaps(x,y,smoothing,xx,w); % the smoothed, weighted spline, densly-spaced estimate
    ang = angle(cos(ang) + 1i*sin(ang));
    y = mag.*cos(ang) +1i*mag.*sin(ang);
    
    %plt(XX,abs(YY),'b')
    %hold on
    %plt(XX,abs(y),'r')
catch ME
    keyboard
end
end

%------------------- OLD CODE
% PL = fft(x,nfft); % measured sound presure level
% scaling = ((nfft-zpadN)/2); % fft scaling factor
% PL = PL / scaling; % scale to give correct output
% 
% [~,fminIndx] = min(abs(freq-fmin));
% [~,fmaxIndx] = min(abs(freq-fmax));
% PL = PL(1:fmaxIndx); % cut to highest calibrated frequency
% 
% FPL = (PL/2) .* (1 + (z0 ./ ZL)); % complex FPL
% RPL = (PL/2) .* (1 - (z0 ./ ZL)); % complex RPL
% PRR = RPL ./ FPL;                 % complex pressure reflectance PR
% 
% pr = abs(RPL).^2 ./ abs(FPL).^2;  % power reflectance
% ipl = abs(FPL) + abs(RPL);        % integrated pressure level
% 
% pRef = 0.00002; % pressure reference (uPa)
% spl = 20*log10(abs(PL)/pRef); % express as magnitude re: 20 uPa
% ipl = 20*log10(ipl/pRef);
% c =  3.493338665700000e+04; % speed of sound in air
% [correction,canalLength,Fqw,RS,RL] = getEPL(z0,ZS,ZL,ipl,spl,freq,c);
% 
% EPL = PL .* correction; % emission pressure level
% epl = 20*log10(abs(EPL)/pRef); % expressed as magnitude re: 20 uPa
% fpl = 20*log10(abs(FPL)/pRef);
% rpl = 20*log10(abs(RPL)/pRef);
% 
% % complex forms are in upper case; magnitudes are in lower case
% % cut to lowest calibrated frequency
% 
% freq = freq(fminIndx:end); % associated frequency vector
% ZS = ZS(fminIndx:end);   % source impedance
% ZL = ZL(fminIndx:end);   % load (ear canal) impedance
% SPL = PL(fminIndx:end);  % sound pressure level
% FPL = FPL(fminIndx:end); % forward pressure level
% RPL = RPL(fminIndx:end); % reverse pressure level
% PRR = PRR(fminIndx:end); % pressure reflectance
% EPL = EPL(fminIndx:end); % emission pressure level
% 
% 
% % function returns four structures:
% % 1) complex values (PL)
% % 2) pressure magnitudes in dB (pl)
% % 3) phases in radians (phi)
% % 4) other values reltated to the calculations (other)
% 
% pl.spl = spl(fminIndx:end); % sound pressure level (dB: re 20 uPa)
% pl.fpl = fpl(fminIndx:end); % forward pressure level (dB: re 20 uPa)
% pl.rpl = rpl(fminIndx:end); % reverse pressure level (dB: re 20 uPa)
% pl.ipl = ipl(fminIndx:end); % integrated pressure level (dB: re 20 uPa)
% pl.epl = epl(fminIndx:end); % emission pressure level (dB: re 20 uPa)
% pl.pr = pr(fminIndx:end);   % power reflectance 
% pl.f = freq;
% 
% Pl.SPL = SPL;
% Pl.FPL = FPL;
% Pl.RPL = RPL;
% Pl.EPL = EPL;    
% Pl.PRR = PRR;
% Pl.f = freq;
% 
% phi.spl = angle(SPL);     % phase of spl (rad)
% phi.fpl = angle(FPL);     % phase of fpl
% phi.rpl = angle(RPL);     % phase of rpl
% phi.prr = angle(PRR);     % phse of pressure reflectance
% phi.epl = angle(EPL);     % phase of epl
% phi.f = freq;
% % power reflectance and ipl have no phase; being composed of sum or ratio of magnitudes
% 
% other.ZL = ZL; % load impedance
% other.ZS = ZS; % cource impedance
% other.z0 = z0; % characteristic (surge) impedance
% other.freq = freq; % frequency vector (Hz)
% other.fmax = fmax; % highest calibrated frequency
% other.fmin = fmin; % lowest calibrated frequency
% other.correction = correction;
% other.canalLength = canalLength; 
% other.diam = isc.diam; % ear canal diameter
% other.Fqw = Fqw; % quarter-wavelength frequency (from EPL calibration)
% other.RS = RS; % source reflectance
% other.RL = RL; % load reflectance
% 
% 
% 
% 
% 
% if ~isempty(f) % if specific frequencies are requested, assume narrowband
%     nFreqs = length(f); % number of frequencies to extract
%     pl.spl = zeros(nFreqs,1); % pre-allocate space
%     pl.fpl = zeros(nFreqs,1);
%     pl.rpl = zeros(nFreqs,1);
%     pl.ipl = zeros(nFreqs,1);
%     pl.prr = zeros(nFreqs,1);
%     pl.pr = zeros(nFreqs,1);
%     for ii=1:nFreqs % loop over frequencies
%         [~,indx] = min(abs(f(ii)-freq)); % find the index of the frequency
%         pl.spl(ii,1) = spl(indx);
%         pl.fpl(ii,1) = fpl(indx);
%         pl.rpl(ii,1) = rpl(indx);
%         pl.ipl(ii,1) = ipl(indx);
%         pl.prr(ii,1) = PRR(indx);
%         pl.pr(ii,1) = pr(indx);
%         pl.epl(ii,1) = epl(indx);
%     end
% else % otherwise assume broadband stimulus
%     pl.spl = spl;
%     pl.fpl = fpl;
%     pl.rpl = rpl;
%     pl.ipl = ipl;
%     pl.prr = PRR;
%     pl.pr =  pr;
%     pl.epl = epl;    
%     f = freq;
% end

    %[~,fminIndx] = min(abs(freq-fmin));
    %[~,fmaxIndx] = min(abs(freq-fmax));
    %SPL = SPL(1:fmaxIndx,:); % cut to highest calibrated frequency
%     SPL_orig = SPL; % original frequency spectrum
%     SPL_full = SPL; % full frequency spectrum, with low and high frequencies set to zero
%     SPL_full(1:fminIndx) = 0;
%     SPL_full(fmaxIndx:nyquist) = 0;
%     SPL_full = [SPL_full(1:nyquist);flipud(conj(SPL_full(2:nyquist-1)))];
%     spl_full = real(ifft(SPL_full));

