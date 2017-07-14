function [d] = ARLas_analyzeDPOAE(header,Data,d,h)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% d = ARLas_analyzeDPOAE(header,Data,d,h);
%
% Analyze distortion product otoacoustic emissions using data obtained from
% ARLas and the experiment file ARLas_dpoae.m
% d = a data structure. Input can be empty, if no data have been collected yet.
% h = a figure handle. Optional input; if not an argument, will not plot
% analysis.
%
% Auditory Research Lab, The University of Iowa
% Deptartment of Communication Sciences & Disorders
% The University of Iowa
% Author: Shawn S. Goodman, PhD
% Date: April 5, 2017
% Updated: Jun2 14, 2017 - ssg; fixed error; filtering was occurring after
%                               FFT; now occurs before.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

onBandNoise = 1; % if =1, will estimate noise in same bin as signal; else an average of 10 bins either side

ref = 0.00002; % reference is 20 uPa
fs = header.fs; % sampling rate (Hz)
f1 = header.userInfo.f1; % first primary frequency
f2 = header.userInfo.f2; % second primary frequency
fdp_cubic = header.userInfo.fdp_cubic; % cubic distortion frequency
fdp_diff = header.userInfo.fdp_diff; % difference distortion frequency

% lowest frequency is the difference frequency; highest is f2. Make an fir
% filter that encompases these frequencies, and apply it.
edgeHz = f2 * .05;
lowCut = fdp_diff - edgeHz;
highCut = f2 + edgeHz;
Data = filterMe(lowCut,highCut,Data,fs);

[frequency,signal,noiseFloor,phase] = ARLas_fda(Data,fs,ref); % frequency domain analysis (fda)

[~,indx_L1] = min(abs(frequency - f1)); % find L1 (actual)
L1 = signal(indx_L1); % level of f1
P1 = phase(indx_L1); % phase of f1
[~,indx_L2] = min(abs(frequency - f2)); % find L2 (actual)
L2 = signal(indx_L2); % level of f2
P2 = phase(indx_L2); % phase of f2

[~,indx_cubic] = min(abs(frequency - fdp_cubic)); % find the cubic distortion tone
Ldp_cubic = signal(indx_cubic); % level of cubic distortion (dB SPL)
if onBandNoise == 1
    Ndp_cubic = noiseFloor(indx_cubic); % level of noise floor (dB SPL) -- ON BAND   
else
    Ndp_cubic = mean([signal(indx_cubic-10),signal(indx_cubic+10)]); % level of noise floor (dB SPL) -- OFF BAND
end
Pdp_cubic = phase(indx_cubic); % phase of cubic distortion (rad)
Pdp_cubic = Pdp_cubic - (2*P1-P2); % subtract off the stimulus phase

[~,indx_diff] = min(abs(frequency - fdp_diff)); % find the difference tone
Ldp_diff = signal(indx_diff); % level of difference tone (dB SPL)
if onBandNoise == 1
    Ndp_diff = noiseFloor(indx_diff); % level of noise floor (dB SPL) -- ON BAND   
else
    Ndp_diff = mean([signal(indx_diff-10),signal(indx_diff+10)]); % level of noise floor -- OFF BAND
end
Pdp_diff = phase(indx_diff); % phase of difference distortion (rad)
Pdp_diff = Pdp_diff - (P2-P1); % subtract off the stimulus phase

if isempty(d)
    d.ref = ref;
    d.fs = fs;
    indx = 1;
else
    n = length(d.f1);
    indx = n+1;
end

d.f1(1,indx) = f1; % f1 frequency (Hz)
d.f2(1,indx) = f2; % f2 frequency (Hz)
d.fdp_cubic(1,indx) = fdp_cubic; % cubic distortion frqeuency (Hz)
d.fdp_diff(1,indx) = fdp_diff; % difference frequency (Hz)
d.L1(1,indx) = L1; % level of f1 primary (dB SPL)
d.L2(1,indx) = L2; % level of f2 primary (dB SPL)
d.Ldp_cubic(1,indx) = Ldp_cubic; % level of cubic distortion (dB SPL)
d.Ldp_diff(1,indx) = Ldp_diff; % level of difference distortion (dB SPL)
d.Ndp_cubic(1,indx) = Ndp_cubic; % noise floor of cubic distortion (dB SPL)
d.Ndp_diff(1,indx) = Ndp_diff; % noise floor of difference distortion (dB SPL)
d.P1(1,indx) = P1;
d.P2(1,indx) = P2;
d.Pdp_cubic(1,indx) = Pdp_cubic;
d.Pdp_diff(1,indx) = Pdp_diff;
d.frequency = frequency; % this is always the same; no need to indx.
d.signal(:,indx) = signal;
d.noiseFloor(:,indx) = noiseFloor;
d.waveform = mean(Data,2); % the DPOAE waveform (averaged)

if nargin == 4 % plot analysis results, if figure handle was passed as an argument
    figure(h)
    subplot(2,1,1)
    hold off
    plot(frequency(indx_cubic)/1000,Ldp_cubic,'*r') % plot first for legend
    hold on
    plot(frequency(indx_diff)/1000,Ldp_diff,'*b')
    plot(frequency(indx_L1)/1000,L1,'*c')
    plot(frequency(indx_L2)/1000,L2,'*c')
    plot(frequency/1000,signal,'Color',[.7 .7 .7],'LineWidth',1)
    plot(frequency/1000,noiseFloor,'Color',[0 0 0],'LineWidth',1)
    plot(frequency(indx_cubic)/1000,Ldp_cubic,'*r') % replot to put on top
    plot(frequency(indx_diff)/1000,Ldp_diff,'*b')
    plot(frequency(indx_L1)/1000,L1,'*c')
    plot(frequency(indx_L2)/1000,L2,'*c')
    xmin = fdp_diff - 100;
    xmax = f2 + 100;
    [~,indxMin] = min(abs(frequency - xmin));
    [~,indxMax] = min(abs(frequency - xmax));
    ymin = min([min(noiseFloor(indxMin:indxMax)),L1,L2,Ldp_cubic,Ldp_diff]);
    ymax = max([max(noiseFloor(indxMin:indxMax)),max(signal(indxMin:indxMax))]);
    ymin = ymin - 5;
    ymax = ymax + 5;
    xlim([xmin/1000 xmax/1000])
    ylim([ymin ymax])
    xlabel('Frequency (kHz)','FontSize',12)
    ylabel('Magnitude (dB SPL)','FontSize',12)
    title(['DPOAE SPECTRUM: ',num2str(d.f2(end)),' Hz'],'FontSize',12)
    legend('2f1-f2','f2-f1','f1,f2','Location','NorthWest')
    pause(0.25)
end

function [X] = filterMe(lowCut,highCut,X,fs)
% ---- bandpass fir filtering (for matrices) -----
coeff = nbf(lowCut,highCut,fs);
[rows,cols] = size(X);
X = fastFilter(coeff,X(:)); 
X = reshape(X,rows,cols);


function b = nbf(cf1,cf2,fs)
% ---- narrow-band filter -----
Fstop1 = cf1-20;         % First Stopband Frequency
Fpass1 = cf1;            % First Passband Frequency
Fpass2 = cf2;            % Second Passband Frequency
Fstop2 = cf2+20;         % Second Stopband Frequency
Dstop1 = 0.001;          % First Stopband Attenuation
Dpass  = 0.057501127785; % Passband Ripple
Dstop2 = 0.001;          % Second Stopband Attenuation
flag   = 'scale';        % Sampling Flag
[N,Wn,BETA,TYPE] = kaiserord([Fstop1 Fpass1 Fpass2 Fstop2]/(fs/2), [0 1 0], [Dstop1 Dpass Dstop2]);
if mod(N,2)~=0
    N = N + 1;
end
b  = fir1(N, Wn, TYPE, kaiser(N+1, BETA), flag)';
