% solution to the discrete versus swept levels problem

% load the discrete data here
cd('C:\myWork\ARLas\Data\julia pair 8\julia pair 8_12SEP2024')
load('Ch3_ER10xA_dpoae_0007.mat')

nSweeps = header.userInfo.nSweeps;
Data = data(:);
Data = reshape(Data,length(Data)/nSweeps,nSweeps);
[Rows,Cols] = size(Data);

f1 = header.userInfo.F1(1);
f2 = header.userInfo.F2(1);
fRatio = header.userInfo.fRatio;
fdp = header.userInfo.Fdp(1);
fs = header.userInfo.fs;
ref = 0.00002;
[frequency,mag,nf,phase] = ARLas_fda(Data(:,2:end-1),fs,ref);

[~,indx1] = min(abs(frequency - f1));
[~,indx2] = min(abs(frequency - f2));
[~,indxdp] = min(abs(frequency - fdp));

L1 = mag(indx1);
L2 = mag(indx2);
Ldp = mag(indxdp);
N1 = nf(indx1);
N2 = nf(indx2);
Ndp = nf(indxdp);
P1 = phase(indx1);
P2 = phase(indx2);
Pdp = phase(indxdp) - (2*P1-P2);

% double check analysis codes for fixed L1
targetL2 = [10 20 30 40 50 60 70]; % L2 in dB FPL
targetL1 = ones(size(targetL2)) * 65; % target levels for primaries (dB FPL); 
% 
% target L1 was 65, achieved was 61.32 (3.68 dB low). Because of RMS???
% target L2 was 70, achieved 64.779 (5.221 dB low) ???

% noise floor is quite a bit elevated on-band for the primaries. Suggests
% that the stimulus is changing over time.

% try looking at each second (there are 12 1-second sweeps) for systematic
% changes
for ii=1:nSweeps
    [~,magx,nfx,phasex] = ARLas_fda(Data(:,ii),fs,ref);
    L1x(ii,1) = magx(indx1);
    L2x(ii,1) = magx(indx2);
    Ldpx(ii,1) = magx(indxdp);
end

figure
subplot(3,1,1)
plot(L1x,'^b-')
ylabel('L1 (dB SPL')
xlabel('Sweep Number')
subplot(3,1,2)
plot(L2x,'vb-')
ylabel('L2 (dB SPL')
xlabel('Sweep Number')
subplot(3,1,3)
plot(Ldpx,'or-')
ylabel('Ldp (dB SPL')
xlabel('Sweep Number')

% compare with a LSF
N = size(Data,1);
t = (0:1:N-1)'/fs;
yc1 = cos(2*pi*f1*t);
ys1 = -sin(2*pi*f1*t);
yc2 = cos(2*pi*f2*t);
ys2 = -sin(2*pi*f2*t);
ycdp = cos(2*pi*fdp*t);
ysdp = -sin(2*pi*fdp*t);

data = mean(Data,2);
w = ones(N,1);
SM = [yc1,ys1,yc2,ys2,ycdp,ysdp,w];
alfa = 0.05;
[coeff,anova,modelSummary] = OLSfit(SM,data,alfa,w);

mag1 = abs(coeff.b(1) + 1i*coeff.b(2));
mag2 = abs(coeff.b(3) + 1i*coeff.b(4));
mag3 = abs(coeff.b(5) + 1i*coeff.b(6));
mag1 = 20*log10(mag1/ref);
mag2 = 20*log10(mag2/ref);
mag3 = 20*log10(mag3/ref);

% note that these LSF magnitudes are 3 dB different from the FFT version.
% the 3 dB difference comes from computing rms (fft version) instead of
% peak version (LSF). LSF is 3 dB higher primaries (but not in actual fact,
% just in the numbers computed
mag1 - L1 % delta = 3.0123 dB
mag2 - L2 % delta = 3.0054 dB
mag3 - Ldp % delta = 3.1150 dB

% target L1 was 65, achieved LSF was 64.3 (0.6671 dB low).
% target L2 was 70, achieved 67.7844 (2.2 dB low)

% now compare with heterodyne filter
cedp = ycdp - 1i*ysdp;
model = coeff.b(1)*yc1 + coeff.b(2)*ys1 + coeff.b(3)*yc2 +coeff.b(4)*ys2;
residual = data - model;
h = hann(1000);
h = h(1:500);
residual(1:500) = residual(1:500).*h;
residual = flipud(residual);
residual(1:500) = residual(1:500).*h;
residual = flipud(residual);
z = residual .* cedp;

SOS = [1.000000000000000  -1.999999498135409   1.000000000000000   1.000000000000000  -1.999917173175506   0.999917185170538
1.000000000000000  -1.999997074920576   1.000000000000000   1.000000000000000  -1.999796028594313   0.999796040830916];        
G = [0.023900934554411
0.004183340469401
1.000000000000000];

re = real(z);
im = imag(z);
re = filtfilt(SOS, G, re);
im = filtfilt(SOS, G, im);

hatdp = 2*abs(re+1i*im);
hatdp = mean(hatdp);
mag3_het = 20*log10(hatdp/.00002); % 8.6218 dB SPL

% FFT = Ldp = 5.7920
% LSF = mag3 = 8.9070
% het = mag3_het = 8.6218 






