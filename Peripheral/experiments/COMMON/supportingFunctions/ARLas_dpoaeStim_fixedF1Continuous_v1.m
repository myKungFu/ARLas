function [Y1,Y2,F1,F2,Fdp,Time,Y1c,Y2c,Ydp,Ydpc] = ARLas_dpoaeStim_fixedF1Continuous_v1(fs,f1,ratioMin,ratioMax,sweepRate,targetL1,targetL2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [Y1,Y2,F1,F2,Fdp,Time,Y1c,Y2c,Ydp,Ydpc] = ARLas_dpoaeStim_fixedF1Continuous_v1(fs,f1,ratioMin,ratioMax,sweepRate,targetL1,targetL2);
%
% Create two vectors of DPOAE stimuli, f1 and f2.
% The vectors are then folded into matrices for efficiency presenting
% through ARLas; goes with ARLas_dpoae_fixedF1Continuous_v1.m.
%
% Author: Shawn Goodman, PhD
% Auditory Research Lab, the University of Iowa
% Date: March 14, 2022
% Last Updated: March 14, 2022 -- ssg
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    octPerSec = sweepRate;

    foldLength = 0.1; % fold into a matrix with columns this length (sec)
    foldN = round(foldLength * fs); % number of samples in each column
    
    fmin = f1 * ratioMin;
    fmax = f1 * ratioMax;

    nOctaves = log2(fmax)-log2(fmin); % numbmer of octaves to sweep over
    len = nOctaves / octPerSec; % length of total stimulus (sec)
    N = round(fs*len); % number of samples in total stimulus

    oct = linspace(0,nOctaves,N)'; % octave spacing
    f2 = 1000 * 2.^(oct); % put into linear space.
    k = fmin / 1000; % k is a multiplier to transform to the desired frequency range
    f2 = f2 * k; % f2 primary frequencies

    rampLen1 = 0.003; %5/fmin; % need some ramp on and off to avoid frequency splatter
    rampN1 = round(fs * rampLen1);
    rampLen2 = 0.003; %5/fmax;
    rampN2 = round(fs * rampLen2);
    pad1 = ones(rampN1,1)*fmin;
    pad2 = ones(rampN2,1)*fmax;
    f2 = [pad1;f2;pad2];

    f1 = ones(size(f2))*f1; % get f1 vector of constant frequencies
    fdp = 2*f1 - f2;
    fRatio = f2 ./ f1;
    
    %[f1,fdp,fRatio] = getRatios(targetL1,targetL2,f2);
    N = length(f1); % recalculate new number of samples

    phi1 = 0; % starting phase
    phi2 = 0;
    phiDP = 0;
    y1 = zeros(N,1); % initialize output
    y2 = zeros(N,1); % initialize output
    ydp = zeros(N,1); % initialize output
    y1c = zeros(N,1); % initialize output
    y2c = zeros(N,1); % initialize output
    ydpc = zeros(N,1); % initialize output

    deltaT = 1/fs; % sampling period
    for ii=1:N
        phaseChange1 = 2*pi* (deltaT * f1(ii));
        phi1 = phi1 + phaseChange1;

        phaseChange2 = 2*pi* (deltaT * f2(ii));
        phi2 = phi2 + phaseChange2;

        phaseChangeDP = 2*pi* (deltaT * fdp(ii));
        phiDP = phiDP + phaseChangeDP;

        y1(ii,1) = sin(phi1);
        y2(ii,1) = sin(phi2);
        ydp(ii,1) = sin(phiDP);

        y1c(ii,1) = cos(phi1);
        y2c(ii,1) = cos(phi2);    
        ydpc(ii,1) = cos(phiDP);

    end

    h = hann(rampN1*2);
    h = h(1:rampN1);
    y1(1:rampN1) = y1(1:rampN1) .* h;
    y2(1:rampN1) = y2(1:rampN1) .* h;
    ydp(1:rampN1) = ydp(1:rampN1) .* h;
    y1c(1:rampN1) = y1c(1:rampN1) .* h;
    y2c(1:rampN1) = y2c(1:rampN1) .* h;
    ydpc(1:rampN1) = ydpc(1:rampN1) .* h;
    h = hann(rampN2*2);
    h = h(1:rampN2);
    h = flipud(h);
    y1(end-rampN2+1:end) = y1(end-rampN2+1:end) .* h;
    y2(end-rampN2+1:end) = y2(end-rampN2+1:end) .* h;
    ydp(end-rampN2+1:end) = ydp(end-rampN2+1:end) .* h;
    y1c(end-rampN2+1:end) = y1c(end-rampN2+1:end) .* h;
    y2c(end-rampN2+1:end) = y2c(end-rampN2+1:end) .* h;
    ydpc(end-rampN2+1:end) = ydpc(end-rampN2+1:end) .* h;


    NN = length(y1);
    nFolds = ceil(NN / foldN);
    extra = mod(NN,foldN);
    if extra ~= 0
        pad = zeros(foldN - extra,1);
        f1 = [f1;pad];
        f2 = [f2;pad];
        fdp = [fdp;pad];
        y1 = [y1;pad];
        y2 = [y2;pad];
        ydp = [ydp;pad];
        y1c = [y1c;pad];
        y2c = [y2c;pad];
        ydpc = [ydpc;pad];
    end
    Y1 = reshape(y1,foldN,nFolds);
    Y2 = reshape(y2,foldN,nFolds);
    Ydp = reshape(ydp,foldN,nFolds);
    Y1c = reshape(y1c,foldN,nFolds);
    Y2c = reshape(y2c,foldN,nFolds);
    Ydpc = reshape(ydpc,foldN,nFolds);
    F1 = reshape(f1,foldN,nFolds);
    F2 = reshape(f2,foldN,nFolds);
    Fdp = reshape(fdp,foldN,nFolds);
    NN = length(y1);
    time = (0:1:NN-1)'/fs;
    Time = reshape(time,foldN,nFolds);
end

% % INTERNAL FUNCTIONS ----------------------------------------------------

function [f1,fdp,fRatio] = getRatios(L1,L2,f2)
% Values taken from Samantha Ginter dissertaion with Sumit Dhar, 2021
    %if L1==65 && L2==55
    if L1<75 && L2<75
        f = [0.5   0.75  1     1.5   2     3     4     6     8     10    11.2  12.5  14    15    16    17    18    19];
        r = [1.26  1.26  1.24  1.24  1.24  1.20  1.20  1.18  1.16  1.16  1.16  1.16  1.16  1.16  1.16  1.16  1.14  1.14];
    elseif L1>=75 && L2>=75
        f = [0.5   0.75  1     1.5   2     3     4     6     8     10    11.2  12.5  14    15    16    17    18    19];
        r = [1.28  1.28  1.28  1.28  1.28  1.24  1.22  1.20  1.18  1.18  1.18  1.18  1.18  1.18  1.18  1.18  1.18  1.18];
    else
        f = [0.5   0.75  1     1.5   2     3     4     6     8     10    11.2  12.5  14    15    16    17    18    19];
        r = [1.22  1.22  1.22  1.22  1.22  1.22  1.22  1.22  1.22  1.22  1.22  1.22  1.22  1.22  1.22  1.22  1.22  1.22];
    end
    f = f * 1000;
    fRatio = interp1(f,r,f2,'pchip');
    f1 = f2 ./ fRatio;
    fdp = 2*f1 - f2;
end


% OLD CODE ----------------------------------------------------------------
