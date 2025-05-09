function [DPOAE,DPOAE_ref] = ARLas_dpoaeAnalysis_sweepL2Continuous_v1(header,Data,DPOAE,RefData)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [DPOAE,DPOAE_ref] = ARLas_dpoaeAnalysis_sweepL2Continuous_v1(header,Data,DPOAE,RefData);
%
% Analyze swept DPOAE recordings; goes with ARLas_dpoae_fixedRatioContinuous_v1.m.
%
% Author: Shawn Goodman, PhD
% Auditory Research Lab, the University of Iowa
% Date: March 14, 2022
% Last Updated: March 14, 2022 -- ssg
% Last Updated: July 18, 2023 -- ssg -- updating for RB experiments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    [fmin,fmax,Y1,Y2,Y1c,Y2c,Ydpc,Ydp,Time,F1,F2,Fdp,fs,sweepRate,nReps,nSweeps] = unpack(header);
    Q = Data(:);
    QQ = reshape(Q,length(Q)/header.userInfo.nSweeps,header.userInfo.nSweeps);
    Data = QQ;
    
    b = getbpf; % bandpass filter from 1-16 kHz
    [Rows,Cols] = size(Data);
    Data = Data(:);
    Data = fastFilter(b,Data);
    Data = reshape(Data,Rows,Cols);
    
%     if ~isempty(RefData)
%         [Rows,Cols] = size(RefData);
%         RefData = RefData(:);
%         RefData = fastFilter(b,RefData);
%         RefData = reshape(RefData,Rows,Cols);
%         Y1c_r = Y1c;
%         Y1_r = Y1;
%         Y2c_r = Y2c;
%         Y2_r = Y2;
%         Ydpc_r = Ydpc;
%         Ydp_r = Ydp;
%         F2_r = F2;
%         F1_r = F1;
%         Fdp_r = Fdp;
%         Time_r = Time;
%     end
%     
    
    [X,Y1c,Y1,Y2c,Y2,Ydpc,Ydp,F2,F1,Fdp,Time] = clean(Data,nSweeps,Y1c,Y1,Y2c,Y2,Ydpc,Ydp,F2,F1,Fdp,Time);
    if ~isempty(RefData)
        Xref = clean(RefData,nSweeps,Y1c_r,Y1_r,Y2c_r,Y2_r,Ydpc_r,Ydp_r,F2_r,F1_r,Fdp_r,Time_r);
    end    
    [Delta_t,BW] = getAnalysisWindows(F2,sweepRate);
    wn = round(Delta_t * fs); % number of samples in each analysis window

    % find the number of analysis steps -------------
    stepSize = 12; % 12 max, 6 min in a good range; original at 8 for fixed ratio
    jj = 1;
    start = 1; % analysis start sample
    finish = start + wn(jj); % analysis finish sample
    NN = size(X,1);
    while finish < NN 
        step = round(wn(jj) / stepSize);
        start = start + step; 
        finish = finish + step;
        jj = jj + 1;
    end
    nIterations = jj-1;
    jj = 1;
    start = 1; % analysis start sample
    finish = start + wn(jj); % analysis finish sample
    
    
% ->> process here
    [L1,N1,P1,L2,N2,P2,Ldp,Ndp,Pdp,f1,f2,fdp] = doAnalysis(nIterations,X,F1,F2,start,finish,Y1c,Y1,Y2c,Y2,Fdp,fs,wn,stepSize);
    if ~isempty(RefData)
        [L1_r,N1_r,P1_r,L2_r,N2_r,P2_r,Ldp_r,Ndp_r,Pdp_r] = doAnalysis(nIterations,Xref,F1,F2,start,finish,Y1c,Y1,Y2c,Y2,Fdp,fs,wn,stepSize);
        DPOAE_ref.L1_r = L1_r;
        DPOAE_ref.N1_r = N1_r;
        DPOAE_ref.P1_r = P1_r;
        DPOAE_ref.L2_r = L2_r;
        DPOAE_ref.N2_r = N2_r;
        DPOAE_ref.P2_r = P2_r;
        DPOAE_ref.Ldp_r = Ldp_r;
        DPOAE_ref.Ndp_r = Ndp_r;
        DPOAE_ref.Pdp_r = Pdp_r;
        DPOAE_ref.f1 = f1;
        DPOAE_ref.f2 = f2;
        DPOAE_ref.fdp = fdp;
    else
        DPOAE_ref = [];
    end
    
% put into octave averages

    % partition the signal energy into 1/3 octaves.  
    % frequency = frequency vector (Hz)
    % Center frequencies are derived from ANSI S1.11-1986, 
    % "Specification for Octave-Band and Fractional Octave Band Analog and Digital Filters", 
    % and from ANSI S1.6-1984, "Preferred Frequencies, Frequency Levels, and Band Numbers for
    % Acoustical Measurements".  The defining formulas yield some values that
    % are not integral numbers.
    fc = [1000, 1250, 1600, 2000, 2500, 3150, 4000, 5000,...
        6300, 8000, 10000, 12500, 16000]; % filter center frequencies
    fL = [900, 1120, 1400, 1800, 2240, 2800, 3550, 4500,...
        5600, 7100, 9000, 11200, 14000]; % lower filter cutoff
    fu = [1120, 1400, 1800, 2240, 2800, 3550, 4500,...
        5600, 7100, 9000, 11200, 14000, 18000]; % upper filter cutoff

    nn = length(fc);
    Ldp_3oct = zeros(nn,1);
    Ndp_3oct = zeros(nn,1);
    for ii=1:nn
        [~,indxL] = min(abs(f2-fL(ii)));
        [~,indxH] = min(abs(f2-fu(ii)));
        chunkS = Ldp(indxL:indxH);
        chunkN = Ndp(indxL:indxH);
        snr = chunkS - chunkN;
        snr(find(snr<=0)) = 0.01; % set all values of snr < 0 to a small number greater than zero
        chunkS = 10.^(chunkS/20);
        chunkN = 10.^(chunkN/20);
        q = sum(chunkS .* snr) / sum(snr);
        if ~isreal(q)
            keyboard
        end
        Ldp_3oct(ii,1) = 20*log10(q);
        q = sum(chunkN .* snr) / sum(snr);
        if ~isreal(q)
            keyboard
        end
        Ndp_3oct(ii,1) = 20*log10(q);
    end
    
    counter = size(DPOAE.L1,2);
    counter = counter + 1;
    DPOAE.targetL1(:,counter) = header.userInfo.targetL1;
    DPOAE.targetL2(:,counter) = header.userInfo.targetL2;
    DPOAE.L1(:,counter) = L1;
    DPOAE.N1(:,counter) = N1;
    DPOAE.P1(:,counter) = P1;
    DPOAE.L2(:,counter) = L2;
    DPOAE.N2(:,counter) = N2;
    DPOAE.P2(:,counter) = P2;
    DPOAE.Ldp(:,counter) = Ldp;
    DPOAE.Ndp(:,counter) = Ndp;
    DPOAE.Pdp(:,counter) = Pdp;
    DPOAE.f1(:,counter) = f1;
    DPOAE.f2(:,counter) = f2;
    DPOAE.fdp(:,counter) = fdp;
    DPOAE.Ldp_3oct(:,counter) = Ldp_3oct;
    DPOAE.Ndp_3oct(:,counter) = Ndp_3oct;
    DPOAE.fc(:,counter) = fc;
    
    %toc/60
end

% INTERNAL FUNCTIONS  -----------------------------------------------------
function [L1,N1,P1,L2,N2,P2,Ldp,Ndp,Pdp,f1,f2,fdp] = doAnalysis(nIterations,X,F1,F2,start,finish,Y1c,Y1,Y2c,Y2,Fdp,fs,wn,stepSize)
    doConstantFit = 1; % fit for a constant dp frequency in each window
                       % otherwise use the swept target.
                       
    f1 = zeros(nIterations,1);
    f2 = zeros(nIterations,1);
    fdp = zeros(nIterations,1);
    Ldp = zeros(nIterations,1);
    Ndp = zeros(nIterations,1);
    Pdp = zeros(nIterations,1);
    L1 = zeros(nIterations,1);
    N1 = zeros(nIterations,1);
    P1 = zeros(nIterations,1);
    L2 = zeros(nIterations,1);
    N2 = zeros(nIterations,1);
    P2 = zeros(nIterations,1);
                       
    for jj=1:nIterations
        if mod(jj,100)==0
            disp(['   analyzing ',num2str(jj),' of ',num2str(nIterations)])
        end

        Signal = X(start:finish,:);
        [rows,cols] = size(Signal);

        % frequency of current analysis
        f1(jj,1) = median(F1(start:finish));
        f2(jj,1) = median(F2(start:finish));
        fdp(jj,1) = median(Fdp(start:finish));

        y1c = Y1c(start:finish);
        y1s = Y1(start:finish);
        y2c = Y2c(start:finish);
        y2s = Y2(start:finish);

        if doConstantFit == 1
            t = (0:1:rows-1)'/fs;
            ydpc = cos(2*pi*fdp(jj,1)*t);
            ydps = -sin(2*pi*fdp(jj,1)*t);
        else
            ydpc = Ydpc(start:finish);
            ydps = Ydp(start:finish);
        end
        
        h = blackman(rows); % method 2.5: blackman window -- has slightly better snr than hann window

        % first, look to reject individual samples
        v = var(Signal,[],2); % variance of the signal agross time
        n = round(size(Signal,1)); % use a max smoother to figure out time regions over wich to downweight
        n = n * .03;
        v2 = v;
        [indx,~] = AR(v2,'mild',0); % this identifies which time samples have noise
        if ~isempty(indx)
            nBad = length(indx);
        else
            nBad = 0;
        end
        W = ones(size(Signal)); % matrix of weights; initialize to ones

        counter = 0;
        for kk=1:nBad % loop across bad time samples
            q = Signal(indx(kk),:);
            [indx2,nRejects] = AR(q,'mild',0);
            if ~isempty(indx2)
                %keyboard
                W(indx(kk),indx2) = 0;
                counter = counter + nRejects;
            end
        end

        % look to reject any full buffers, after removing noisy samples
        wrms = sqrt(sum((Signal .* W).^2,1) ./ sum(W,1)); % weighted rms energy in each sweep
        [indxFull,~] = AR(wrms,'mild',0);

        if ~isempty(indxFull)
            W(:,indxFull) = [];
            Signal(:,indxFull) = [];
        end
        newCols = size(Signal,2);    
        clear B1 B2 Bdp ww    
        for ii=1:newCols
            w = W(:,ii);
            %w = meanSmoother(w,n);
            ww(ii,1) = sum(w);
            [B1(ii,1),B2(ii,1),Bdp(ii,1)] = ARLas_dpoae_OLSfit_v1(y1c,y1s,y2c,y2s,ydpc,ydps,Signal(:,ii),h.*w);
        end
        
        [indx,~] = AR(abs(Bdp),'mild',0);
        www = ww;
        www(indx) = 0;
        [Ldp(jj,1),Ndp(jj,1),Pdp(jj,1)] = getSN(Bdp,www);
        
        [indx,~] = AR(abs(B1),'mild',0);
        www = ww;
        www(indx) = 0;
        [L1(jj,1),N1(jj,1),P1(jj,1)] = getSN(B1,www);

        [indx,~] = AR(abs(B2),'mild',0);
        www = ww;
        www(indx) = 0;
        [L2(jj,1),N2(jj,1),P2(jj,1)] = getSN(B2,www);
        
% if L2(jj,1) >95
%     keyboard
% end
        
        step = round(wn(jj) / stepSize);
        start = start + step; 
        finish = finish + step;
    end
end

function [signal,noise,phase] = getSN(B,W)
    if nargin == 2
        % make sure the sum of the weights equals 1
        W = W ./ sum(W);
        B = B(:)';
        W = W(:)';
        K = length(B);
        Xk = B; % take the DFT of each buffer in the matrix (down columns)
        Xbar = sum(W .* Xk,2) / sum(W); % signal is the coherent weighted mean
        Xbar2 = abs(Xbar) .^2; % signal energy
        phase = angle(Xbar); % signal phase
        XBAR = repmat(Xbar,1,K); % replicate to matrix (to avoid running a loop)
        S2 = sum(W .* ((Xk - XBAR) .* conj(Xk - XBAR))) / (sum(W)-(sum(W)/K)); % weighted variance (this is the HA noise floor)
        Se2 = (1/K) * S2; % energy of the standard error
        pRef = 0.00002;
        signal = 10*log10(Xbar2/(pRef^2));
        noise = 10*log10(Se2/(pRef^2));
    else
        B = B(:)';
        K = length(B);
        Xk = B; % take the DFT of each buffer in the matrix (down columns)
        Xbar = (1/K) * sum(Xk,2); % signal is the coherent mean
        Xbar2 = abs(Xbar) .^2; % signal energy
        phase = angle(Xbar); % signal phase
        XBAR = repmat(Xbar,1,K); % replicate to matrix (to avoid running a loop)
        S2 = (1/(K-1)) * sum((Xk - XBAR) .* conj(Xk - XBAR),2); % variance (this is the HA noise floor)
        % Se = sqrt((1/K) * S2); % standard error of mean (calculation not necessary)
        Se2 = (1/K) * S2; % energy of the standard error
        pRef = 0.00002;
        signal = 10*log10(Xbar2/(pRef^2));
        noise = 10*log10(Se2/(pRef^2));
    end
end

function [Delta_t,BW] = getAnalysisWindows(F2,sweepRate)
    % Delta_t = BW / sweepRate;
    %Delta_t = alpha / (NR(foae(t)) * d*log(foae(t)) / dt); % Abdala eq. 3
    % foae(t) describes how the oae frequency changes with time
    % NR is the reflecton component delay measured in perods
    % alpha is a constant chosen to yield a bandwidth of 0.0625 oct at 1 kHz.
    % For logarithmically swept tones at rate r, the function
    % d*ln(foae(t))/dt becomes r * ln(2).
    % alpha / (NR * r * log(2))
    % to solve this, once we know freq(t) frequency as a function of time,
    % we can get the expected reflection delay, express as a number of
    % cycles of f. 
    
    [tau,tau5,tau95] = sfoaeLatency(F2); % tau is the reflection delay in seconds.
    % here, use tau5, which is shorter delay, because we are using higher stimulus levels. 
    NR = tau ./ (1./F2);

    % Calculate Delta_t, the analysis length.
    %Delta_t = alpha / (NR(foae(t)) * d*log(foae(t)) / dt); % Abdala eq. 3
    %Delta_t = alpha / (NR * r * log(2));
    %Delta_t = alpha ./ (NR * sweepRate * log(2));
    %
    % In this code, Delta_t = BW / r
    % We want BW to be 0.0625 at 1 kHz
    % Re-arranging the first equation,
    alpha = 1;
    Delta_t = alpha ./ (NR * sweepRate * log(2));
    [~,indx1k] = min(abs(F2 - 1000));
    BW1k = Delta_t(indx1k) * sweepRate;
    target = 0.0625;
    alpha = target / BW1k;
    Delta_t = alpha ./ (NR * sweepRate * log(2));
    BW1k = Delta_t(indx1k) * sweepRate;
    if abs(BW1k - target) > 0.0001
        disp('Warning: target bandwidth not achieved!')
        keyboard
    end
    BW = Delta_t * sweepRate;
end

function [fmin,fmax,Y1,Y2,Y1c,Y2c,Ydpc,Ydp,Time,F1,F2,Fdp,fs,sweepRate,nReps,nSweeps] = unpack(header)
    fmin = header.userInfo.fmin;
    fmax = header.userInfo.fmax;
    Y1 = header.userInfo.Y1;
    Y2 = header.userInfo.Y2;
    Y1c = header.userInfo.Y1c;
    Y2c = header.userInfo.Y2c;
    Ydpc = header.userInfo.Ydpc;
    Ydp = header.userInfo.Ydp;
    Time = header.userInfo.Time;
    F1 = header.userInfo.F1;
    F2 = header.userInfo.F2;
    Fdp = header.userInfo.Fdp;
    %LTC1 = header.userInfo.LTC1;
    %LTC2 = header.userInfo.LTC2;
    fs = header.userInfo.fs;
    sweepRate = header.userInfo.sweepRate;
    nReps = header.userInfo.nReps; % number of reps
    nSweeps = header.userInfo.nSweeps; % number of sweeps
end

function [X,Y1c,Y1,Y2c,Y2,Ydpc,Ydp,F2,F1,Fdp,Time] = clean(Data,nSweeps,Y1c,Y1,Y2c,Y2,Ydpc,Ydp,F2,F1,Fdp,Time)
    Rows = size(Y1(:),1);
    Cols = nSweeps;
    X = reshape(Data,Rows,Cols);
    Y1c = Y1c(:);
    Y1 = Y1(:);
    Y2c = Y2c(:);
    Y2 = Y2(:);
    Ydpc = Ydpc(:);
    Ydp = Ydp(:);
    F2 = F2(:);
    F1 = F1(:);
    Fdp = Fdp(:);
    Time = Time(:);

    indx = find(F2>0);
    X = X(indx,:);
    Y1c = Y1c(indx);
    Y1 = Y1(indx);
    Y2c = Y2c(indx);
    Y2 = Y2(indx);
    Ydpc = Ydpc(indx);
    Ydp = Ydp(indx);
    F2 = F2(indx);
    F1 = F1(indx);
    Fdp = Fdp(indx);
    Time = Time(indx);

%     indx1 = min(find(F2>1000));
%     indx16 = min(find(F2>=16000));
%     X = X(indx1:indx16,:);
%     Y1c = Y1c(indx1:indx16);
%     Y1 = Y1(indx1:indx16);
%     Y2c = Y2c(indx1:indx16);
%     Y2 = Y2(indx1:indx16);
%     Ydpc = Ydpc(indx1:indx16);
%     Ydp = Ydp(indx1:indx16);
%     F2 = F2(indx1:indx16);
%     F1 = F1(indx1:indx16);
%     Fdp = Fdp(indx1:indx16);
%     Time = Time(indx1:indx16);
     
end

function b = getbpf()
    Fs = 96000;
    Fstop1 = 100;             % First Stopband Frequency 
    Fpass1 = 400;            % First Passband Frequency -->for a low freq of 750, the dp is 450 Hz
    Fpass2 = 20000;           % Second Passband Frequency
    Fstop2 = 24000;           % Second Stopband Frequency
    Dstop1 = 0.001;           % First Stopband Attenuation
    Dpass  = 0.057501127785;  % Passband Ripple
    Dstop2 = 0.0001;          % Second Stopband Attenuation
    flag   = 'scale';         % Sampling Flag
    % Calculate the order from the parameters using KAISERORD.
    [N,Wn,BETA,TYPE] = kaiserord([Fstop1 Fpass1 Fpass2 Fstop2]/(Fs/2), [0 ...
                                 1 0], [Dstop1 Dpass Dstop2]);
    if mod(N,2)~= 0
        N = N + 1;
    end
    b  = fir1(N, Wn, TYPE, kaiser(N+1, BETA), flag)';
end

% OLD CODE ----------------------------------------------------------------
