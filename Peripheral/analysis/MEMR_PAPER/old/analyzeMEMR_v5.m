function [MEMR] = analyzeMEMR_v5(pathName,fileNameL,fileNameR,subjectName,runName)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [h10,h11] = analyzeMEMR_v5()
% (Data,time,fs,indx)
%
% Author: Shawn Goodman
% Date: August 6, 2022
% Last Updated: August 9, 2022
% Last Updated: October 10, 2022 -- ssg -- version 2; goes with
%                       goodmanMEMR_v2.m.
% Last Updated: October 11, 2022 -- ssg -- added figure handles to output.
% Last Updated: February 27, 2023 -- ssg --
% Last Updated: March 9, 2023 -- ssg -- moved on to v4. Calculating the arc
%               lengths instead of magnitude
% Last Updated: March 13, 2023 -- ssg -- still working on definding the
%               noise, so can find threshold
% Last Updated: March 14-17, 2023 -- ssg --
%
% NOTES: Still to do: 
% choose peak frequencies by picking the larges values
%   over a range. Currently we just use fixed values, which are not always
%   the best choice (i.e., our current frequencies will not usually yield
%   the largest value). This won't match with the Mepani method.
%
% Add a data structre to save the analysis. Include template, templateSNR, etc.
% analyize actual noise (rather than intended noise and report that too
% make sure that just taking abs and turned delta do give nonmonotonic
% results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
    % ---------------------------------------------------------------------
    % type 1 = curve length
    % type 2 = abs(delta) -- straight magnitude
    % type 3 = abs(delta-1) -- turned delta

    octStepSize = 1/2; % octave step size for analysis
    fmin = 250; %250; % minimum frequency to examine
    fmax = 4000; % maximum frequency to examine
    % frequencies of interest:
    peakFreqs = round(2.^(log2(fmin):octStepSize:log2(fmax))');

    % extract the raw recordings -----
    load([pathName,fileNameL]) % the clicks
    micCorrection = header.userInfo.micCorrectionL;
    Data = applyMicCorrection(data,micCorrection);
     clear data header
    load([pathName,fileNameR]) % the noise
    micCorrection = header.userInfo.micCorrectionL;
    Noise = applyMicCorrection(data,micCorrection);
    
    chunkSize = header.userInfo.chunkSize; % 9600
    nChunks = header.userInfo.nChunks; % 80
    nSweeps = header.userInfo.nSweeps; % 15
    Data = reshape(Data,chunkSize*nChunks,nSweeps);
    Noise = reshape(Noise,chunkSize*nChunks,nSweeps);
    fs = header.fs; % sampling rate (96 kHz)
    modLen = 4; % modulation time (4 seconds increasing, 4 seconds decreasing)
    time = (0:1:size(Data,1)-1)'/fs; % time vector
    indx = header.userInfo.clickIndx; % location of the clicks
     clear data header
    time = time / time(end); % express time from 0 to 1

    % When considering which portions of each time window to analyze:
    % Within each click window, we want to include several round trip
    % travel times of the click in the canal, but no low-frequency OAEs. We
    % should consider only the first 2 ms of the click, which is 
    % round(0.002*fs) = 192 samples. 
    % We used a noise activator with a rise time of 4 seconds and a fall
    % time of 4 seconds. This gives 8 seconds / .05 = 160 clicks.
    nSamples = 193; % number of time samples in each click analysis window
    nClicks = length(indx); % number of total clicks to analyze

% location of noise levels
    modLen = 4; % modulation duration in sec (from onset of noise to offset)
    [nLoops,clickTrain,noiseSamples,nReps,H] = getStimuli(fs,modLen);
    % Here, H is the (linear) amplitude vector for the noise. It is 8
    % seconds long total--4 seconds rise and 4 seconds fall.
    % The noise floor in the ear canal is roughly 60 dB SPL.
    % The noise itself goes from 50 dB SPL to 120, covering a change of 70
    % dB in 4 seconds, so a rate of 17.5 dB/second
    h = H(:,1)*(10^(120/20)*.00002); % h is the vector of noise levels (intended)
    %plot(20*log10(abs(Noise(:,1))/.00002))
    %hold on
    %plot(20*log10(h/.00002),'r','LineWidth',2)
    noiseLvl = 20*log10(h/.00002); % noise level
    noiseRate = (120 - 70) / 4; % noise changes at 12.5 dB per second
            % multiply this value by the timeShift value to see how many dB
            % to subtract from the noise vector when calculating thresholds

    % run the data analysis program
    MEMR = runme(Data,nSamples,nClicks,indx,time,fs,peakFreqs,noiseLvl,noiseRate,subjectName,runName);
keyboard

end
% INTERNAL FUNCTIONS ------------------------------------------------------
function [MEMR] = runme(Data,nSamples,nClicks,indx,time,fs,peakFreqs,noiseLvl,noiseRate,subjectName,runName)
    %T = fft(template,fs); % put template into the frequency domain
    % note that we are not scaling, because we will take a ratio later, so it doesn't matter
    F = (0:1:fs-1)'*(fs/fs); % full frequency vector (Hz)
    % we're not examining the entire frequency range, so find the cut points
    [~,fIndx1] = min(abs(min(peakFreqs)-F));
    [~,fIndx2] = min(abs(max(peakFreqs)-F));
    F = F(fIndx1:fIndx2);
    %T = T(fIndx1:fIndx2);

    for jj=1:nClicks % loop across each click position in the sweep
        timeChunk(1,jj) = time(indx(jj)); % the temporal postion of each click position
    end
    % time went from 0 to 1, but actually, the length is 8 seconds
    timeChunk = timeChunk * 8; % this gives the correct 50 ms interval between clicks
    
    % find the location of the desired analysis frequencies
    for ii=1:length(peakFreqs)
        [~,peakIndx(ii,1)] = min(abs(F-peakFreqs(ii)));
    end
    nFreqs = length(peakIndx); % number of frequencies to analyze
    nSweeps = size(Data,2); % number of sweeps
    
    % we will make a multidimensional matrix to store the analysis 
    % nClicks x nSweeps x nFreqs 
    % so each column is a single sweep, each row contains the click change for a single frequency
    % M is size 160 clicks x 15 sweeps x n frequencies
    M = zeros(nClicks,nSweeps,nFreqs);  % unnormalized analysis
    Mnorm = zeros(nClicks,nSweeps,nFreqs); % normalized analysis

    for ii=1:nClicks % loop across each click position
        Chunk = Data(indx(ii):indx(ii)+nSamples,:); % contains all the sweeps
        FT = fft(Chunk,fs); % contains all sweeps at all frequencies for a given click position
        FT = FT(fIndx1:fIndx2,:); % cut down to size
        R = FT;
        for jj=1:nFreqs % extract each desired frequency
            r = R(peakIndx(jj),:);
            M(ii,:,jj) = r.'; % save in the appropriate position in the matrix
        end
    end

    % there is slow drift across sweeps (not related to the MEMR)
    % subtract this off
    for jj=1:nFreqs
        m = M(:,:,jj);
        [rows,cols] = size(m);
        mr = real(m(:));
        mi = imag(m(:));
        t = (1:1:length(mr))';
        smoothing = 0.0000000009;
        ppr = csaps(t,mr,smoothing);
        ppi = csaps(t,mi,smoothing);
        mr_sm = ppval(ppr,t);
        mi_sm = ppval(ppi,t);
        % plot(mr)
        % hold on
        % plot(mr_sm)
        %pause
        mr = mr - mr_sm;
        mi = mi - mi_sm;
        m = reshape(mr + 1i*mi,rows,cols);
        M(:,:,jj) = m;
    end

    % Normalize each sweep by the start of that sweep
    % use the average of baselineN clicks to do the normalization
    baselineN = 3;
    for jj=1:nFreqs
        m = M(:,:,jj);
        baseline = m(1:baselineN,:);
        baseline2 = m(end-baselineN+1:end,:);
        baseline = mean([baseline;baseline2],1); % take from the end, too
        mnorm = m ./ baseline; % normalized m 
        %mnorm = m; % use this to NOT normalize
        Mnorm(:,:,jj) = mnorm;
    end

    % find the time shift of the MEMR
    for jj=1:nFreqs % loop across frequencies
        m = Mnorm(:,:,jj);
        Swirls(:,jj) = mean(m,2);
        type = 1;
        [X(:,jj),snr(jj,1)] = getTimeShift(m,type);
    end
    % obviously, there can only be 1 time shift, and each frequency gives
    % an estimate. But they estimates are not all equal. Get an overall
    % estimate by taking the weighted least-squraes fit, weighting by snr values
    % Also, don't fit across entire time window:
    i1 = 70; % start of fitting time window
    i2 = 100; % end of fitting time window
    XX = X(i1:i2,:); % just the parts that will be fit
    tt = timeChunk(i1:i2)'; % 
    TT = repmat(tt,1,size(XX,2));
    TTT = TT(:); % time vecgtor
    YYY = XX(:); % observed data points
    X = [TTT.^2,TTT,ones(size(TTT))]; % solution matrix
    WWW = ones(size(XX)); % set the weights to the SNRs of each frequency
    for kk=1:nFreqs
        WWW(:,kk) = WWW(:,kk) * snr(kk);
    end
    WWW = WWW(:);
    [B,yhat,residuals] = ARLas_wlsf(YYY,X,WWW); % solve the weighted least-squares problem
    timeShift = roots(polyder(B)) - 4; % this is the peak time shift in seconds
    % plotting
    %     plot(timeChunk,x,'k.-')
    %     hold on
    %     plot(TT,XX,'.-')
    %     yy = polyval(p,tt);
    %     plot(tt,yy,'r','LineWidth',2)
    %     yy = polyval(B,tt);
    %     plot(tt,yy,'g','LineWidth',2)
    

    % get the cumulative arc lengths
    for jj=1:nFreqs
    %jj=9;        
        m = Mnorm(:,:,jj);
        type = 1; % this is for cumulative (complex) arc lenghts (CAL)
        [CAL(:,jj),CALn(:,jj),Swirls_sm(:,jj)] = crawl(timeChunk,m,type,timeShift,subjectName,runName);

        type = 2; % this is for magnitude only (MAG)
        [MAG(:,jj),MAGn(:,jj)] = crawl(timeChunk,m,type,timeShift,subjectName,runName);
    end
    CAL = abs(CAL); % signal must be 0 or greater

    % analyze each frequency band: 
    %   threshold up
    %   threshold down
    %   peak activation level
    %   amount of hysteresis (area under the curve, or threshold difference
    pw1 = 78; % window where the peak should be
    pw2 = 88;
    snrCut = 6; % cutoff for saying signal is present above the noise
    TS = 0.05; %8 / size(CAL,1); % here, TS (sampling period) is 50 ms (between each click)
    time2 = (0:1:size(CAL,1)-1)'* TS;
    for ii=1:length(time2)
        [~,indx] = min(abs(time2(ii)-time*8));
        noiseLvl2(ii,1) = noiseLvl(indx);
    end
    
    present = zeros(size(CAL,2),1); % MEMR present (1) or not present (0)
    peakAmp = zeros(size(CAL,2),1); % MEMR amplitude (linear ratio change)
    peakTime = zeros(size(CAL,2),1); % time at which peak occurs
    peakIndx = zeros(size(CAL,2),1); % index at which peak occurs
    thdUp_time = zeros(size(CAL,2),1); % time at which MEMR starts
    thdDn_time = zeros(size(CAL,2),1); % time at which MEMR finishes
    thdUp_lvl = zeros(size(CAL,2),1); % noise level at which MEMR starts
    thdDn_lvl = zeros(size(CAL,2),1); % noise level at which MEMR finishes
    thdUp_indx = zeros(size(CAL,2),1); % index at which MEMR starts
    thdDn_indx = zeros(size(CAL,2),1); % index at which MEMR finishes
    hysteresis = zeros(size(CAL,2),1); % mean absolute difference between up and down
    hysteresis2 = zeros(size(CAL,2),1); % threshold differences in seconds

    for kk=1:nFreqs
        Signal = CAL(:,kk);
        Noise = CALn(:,kk);
        Snr = 20*log10(Signal ./ Noise);
        
        % if the SNR is acceptable in the peak window, then MEMR is present
        snr2 = Snr(pw1:pw2);
        if max(snr2) >= snrCut
            present(kk,1) = 1;
        end
        % if MEMR is present, look for peak and thresholds
        if present(kk,1) == 1
            signal = Signal(pw1:pw2);
            [peakAmp(kk,1),peakIndx(kk,1)] = max(signal);
            peakIndx(kk,1) = peakIndx(kk,1) + pw1 - 1;
            peakTime(kk,1) = time2(peakIndx(kk,1));

            % look for the low-side (up) threshold
            done = 0;
            counter = peakIndx(kk,1);
            while done == 0
                next = Snr(counter-1);
                if next > snrCut
                    counter = counter - 1;
                else
                    done = 1;
                end
                if counter < 2
                    done = 1;
                end
            end
            thdUp_indx(kk,1) = counter;
            thdUp_time(kk,1) = time2(counter);
            thdUp_lvl(kk,1) = noiseLvl2(counter);

            % look for the high-side (down) threshold
            NNN = length(Snr);
            done = 0;
            counter = peakIndx(kk,1);
            while done == 0
                %current = Snr(counter);
                next = Snr(counter+1);
                if next > snrCut
                    counter = counter + 1;
                else
                    done = 1;
                end
                if counter > NNN-1
                    done = 1;
                end
            end
            thdDn_indx(kk,1) = counter;
            thdDn_time(kk,1) = time2(counter);
            thdDn_lvl(kk,1) = noiseLvl2(counter);
            
            % calculate hysteresis
            try
            leftSide = CAL(thdUp_indx(kk):peakIndx(kk),kk);
            rightSide = CAL(peakIndx(kk):thdDn_indx(kk),kk);
            leftN = length(leftSide);
            rightN = length(rightSide);
            bigN = max([leftN,rightN]);
            leftSide = CAL(peakIndx(kk)-bigN+1:peakIndx(kk),kk);
            rightSide = CAL(peakIndx(kk):peakIndx(kk)+bigN-1,kk);
            leftSide = flipud(leftSide);
            hysteresis(kk,1) = mean(abs(rightSide - leftSide));
            hysteresis2(kk,1) = (abs(thdUp_indx(kk) - peakIndx(kk)) - abs(thdDn_indx(kk) - peakIndx(kk)))*0.05;

            %figure(111)
            %plot(leftSide)
            %hold on
            %plot(rightSide)
            %hold off
            %pause
            %close(111)
            catch
            end

        end

        figure % plotting thresholds and max
        mult = 10^(snrCut/20);
        plot(time2,CAL(:,kk),'r')
        hold on
        plot(time2,CALn(:,kk)*mult,'k')
        if present(kk,1) == 1
            plot(thdUp_time(kk,1),CAL(thdUp_indx(kk,1),kk),'+b')
            plot(thdDn_time(kk,1),CAL(thdDn_indx(kk,1),kk),'+b')
            plot(peakTime(kk,1),peakAmp(kk,1),'*b')
            line([peakTime(kk,1),peakTime(kk,1)],[0,peakAmp(kk,1)],'Color',[0 0 1],'LineStyle','-','LineWidth',0.5)
            line([thdUp_time(kk,1),thdUp_time(kk,1)],[0,CAL(thdUp_indx(kk,1),kk)],'Color',[0 0 1],'LineStyle','--','LineWidth',0.5)
            line([thdDn_time(kk,1),thdDn_time(kk,1)],[0,CAL(thdDn_indx(kk,1),kk)],'Color',[0 0 1],'LineStyle','--','LineWidth',0.5)
        end
        xlabel('Time (s)')
        ylabel('Relative Change')
        title([subjectName(1:end-1),'     ',runName(1:end-1),'     ',num2str(peakFreqs(kk)),' kHz'])
        %keyboard

    end

    figure % comparing CAL to MAG size -------
    subplot(1,2,1)
    plot(peakFreqs,20*log10(max(CAL)),'r')
    hold on
    plot(peakFreqs,20*log10(max(MAG)),'b')
    xlabel('Frequency (Hz)')
    ylabel('Peak Magnitude (dB)')
    legend('CAL','MAG')
    subplot(1,2,2)
    plot(20*log10(snr),20*log10(max(CAL)) - 20*log10(max(MAG)),'k*')
    xlabel('SNR (dB)')
    ylabel('CAL - MAG (dB)')

    MEMR.present = present;
    MEMR.peakFreqs = peakFreqs;
    MEMR.nFreqs = length(peakFreqs);
    MEMR.timeShift = timeShift;
    MEMR.snr = snr;
    MEMR.peakAmp = peakAmp;
    MEMR.peakTime = peakTime;
    MEMR.peakIndx = peakIndx;
    MEMR.thdUp_time = thdUp_time;
    MEMR.thdDn_time = thdDn_time;
    MEMR.thdUp_lvl = thdUp_lvl - (timeShift * noiseRate);
    MEMR.thdDn_lvl = thdDn_lvl - (timeShift * noiseRate);
    MEMR.thdUp_indx = thdUp_indx;
    MEMR.thdDn_indx = thdDn_indx;
    MEMR.hysteresis = hysteresis;
    MEMR.hysteresis2 = hysteresis2;
    MEMR.CAL = CAL; % complex arc length -- signal
    MEMR.CALn = CALn; % complex arc length -- noise
    MEMR.MAG = MAG; % complex arc length -- signal
    MEMR.MAGn = MAGn; % complex arc length -- noise
    MEMR.time = time2; % time across each sweep (8 seconds)
    MEMR.noiseLvl = noiseLvl2; % noise level as a function of time
    MEMR.noiseLvlCorrected = noiseLvl2 - (timeShift * noiseRate); % noise level corrected for shift
    MEMR.Swirls = Swirls;
    MEMR.Swirls_sm = Swirls_sm;

    keyboard

end
function [xr] = newAR(xr,multiplier)
    q25 = prctile(xr,25);
    q75 = prctile(xr,75);
    qIQR = iqr(xr);
    indxBadU = find(xr>q75+multiplier*qIQR);
    indxBadL = find(xr<q25-multiplier*qIQR);
    if isempty(indxBadU) & isempty(indxBadL)
        indxBad = [];
    elseif ~isempty(indxBadU) & ~isempty(indxBadL)
        indxBad = unique([indxBadU(:);indxBadL(:)]);
    elseif ~isempty(indxBadU)
        indxBad = indxBadU;
    elseif ~isempty(indxBadL)
        indxBad = indxBadL;
    end
    if ~isempty(indxBad) % if there are rejections to make
        nRejects = length(indxBad);
        nn = length(xr);
        for rindx = 1:nRejects
            if indxBad(rindx) == 1
                xr(indxBad(rindx)) = mean([xr(indxBad(rindx)+1),xr(indxBad(rindx)+2)]);
            elseif indxBad(rindx) == nn
                xr(indxBad(rindx)) = mean([xr(indxBad(rindx)-1),xr(indxBad(rindx)-2)]);
            else
                xr(indxBad(rindx)) = mean([xr(indxBad(rindx)-1),xr(indxBad(rindx)+1)]);
            end
        end
    end
end
function [x,snr] = getTimeShift(x,type)
    % MEMR is delayed by some amount (~200 ms?); find the delay ("timeShift")
    medSm = 5;
    meanSm = 15;
    [xrs,xis,xrn,xin,sn_r,sn_i] = doSmoothing(x,type,medSm,meanSm);
    x = abs(mean(xrs,2) + 1i*mean(xis,2));
    snr = max([sn_r,sn_i]);
end
function [xrs,xis,xrn,xin,sn_r,sn_i] = doSmoothing(x,type,medSm,meanSm)
    % mean and median smoothing points can be specified.
    % if not specified, defaults are 5 (median) and 2 (mean).
    if nargin < 3
        medSm = 5; % number of samples (to both left and right) over which to smooth for median
        meanSm = 2; % same, for mean
    end
    if type == 1
        x = x-1;
    end
    xr = real(x);
    xi = imag(x);
    % We will take the sum to estimate the signal and the difference to
    % estimate the noise. Therefore, we need an even number. Only use 14.
    % this will leave averaging and estimating across 7 buffers
    xr = xr(:,1:14);
    xi = xi(:,1:14);
    xrs = (xr(:,[1,3,5,7,9,11,13]) + xr(:,[2,4,6,8,10,12,14]))/2; % signal
    xrn = (xr(:,[1,3,5,7,9,11,13]) - xr(:,[2,4,6,8,10,12,14]))/2; % noise
    xis = (xi(:,[1,3,5,7,9,11,13]) + xi(:,[2,4,6,8,10,12,14]))/2; % signal
    xin = (xi(:,[1,3,5,7,9,11,13]) - xi(:,[2,4,6,8,10,12,14]))/2; % noise

    % put matrices into single column vectors for artifact rejection and smoothing
    [rows,cols] = size(xrs); % this is so you can put back into matrices later
    xrs = xrs(:); % signal ---
    xis = xis(:);
    xrn = xrn(:); % noise ---
    xin = xin(:);

    % artifact rejection for exteme values
    multiplier = 3; % outliers will be considered this value time the IQR.
                    % outliers will be replaced by the average of surround values
    xrs = newAR(xrs',multiplier)'; % signal ---
    xis = newAR(xis',multiplier)';
    xrn = newAR(xrn',multiplier)'; % noise ---
    xin = newAR(xin',multiplier)';

    % smoothing is done by a median smoother, followed by a mean smoother
    xrs = medianSmoother(xrs,medSm); % signal ---
    xrs = meanSmoother(xrs,meanSm);
    xis = medianSmoother(xis,medSm);
    xis = meanSmoother(xis,meanSm);
    xrn = medianSmoother(xrn,medSm); % noise ---
    xrn = meanSmoother(xrn,meanSm);
    xin = medianSmoother(xin,medSm);
    xin = meanSmoother(xin,meanSm);

    sn_r = sqrt(mean(xrs.^2)) / sqrt(mean(xrn.^2));
    sn_i = sqrt(mean(xis.^2)) / sqrt(mean(xin.^2));

    % put back into matrices (7 columns)
    xrs = reshape(xrs,rows,cols); % signal ---
    xis = reshape(xis,rows,cols);
    xrn = reshape(xrn,rows,cols); % noise ---
    xin = reshape(xin,rows,cols);
end
function [delta,nf,swirl] = crawl(t,x,type,timeShift,subjectName,runName)
    % calculate the amount of change by finding the length of the smoothed curve
    % t = timeChunk (1 x 160 vector)
    % x = complex vector of changes across time
    % arcLen_sm = the smoothed cumulative arc length
    [xrs,xis,xrn,xin] = doSmoothing(x,type);
    swirl = mean(xrs,2)+1i*mean(xis,2);
    Q = [];
    Qn = [];

    %ttt = t;
    %tt = (t(1):0.001:t(end));
    %t = tt;
    dt = median(gradient(t));
    for iii=1:size(xrs,2) % loop over the 7 buffers
        sr = xrs(:,iii)';
        si = xis(:,iii)';
        nr = xrn(:,iii)';
        ni = xin(:,iii)';

        %sr = spline(ttt,sr,tt);
        %si = spline(ttt,si,tt);
        %nr = spline(ttt,nr,tt);
        %ni = spline(ttt,ni,tt);
        %[~,midPoint] = min(abs(tt-4));
        %foldIndx = midPoint + round(timeShift);
        [~,midPoint] = min(abs(t-4));
        foldIndx = midPoint + round(timeShift/0.05);
        foldIndx2 = (80 +  80 - foldIndx) + 1;
         nn = length(sr);

        % CALCULATE FOR SIGNAL --------------------------------------------
        % integrate from the left
        if type == 1
            grad = sqrt((  (gradient(sr,t)).^2   + (gradient(si,t)).^2  )  );
            k = ones(size(grad));
            k(foldIndx:end) = -k(foldIndx:end);
            q1 = cumsum(grad.*k) * dt;
            sr = fliplr(sr);
            si = fliplr(si);
            grad = sqrt((  (gradient(sr,t)).^2   + (gradient(si,t)).^2  )  );
            k = ones(size(grad));
            k(foldIndx2:end) = -k(foldIndx2:end);
            q2 = fliplr(cumsum(grad.*k) * dt);
            wf1 = linspace(1,0,nn); % weighting function
            wf2 = linspace(0,1,nn);
            arcLen = ((q1.*wf1)+(q2.*wf2)) ./ (wf1+wf2); % weighted sum 

%             cc = 0;
%             nn = length(sr);
%             for ii=2:foldIndx
%                 cc(ii,1) = cc(ii-1) + (sqrt(    ((sr(ii+1)-sr(ii))/dt).^2    +      ((si(ii+1)-si(ii))/dt).^2    )       *dt);
%             end
%             for ii=foldIndx+1:nn-1
%                 cc(ii,1) = cc(ii-1) - (sqrt(    ((sr(ii+1)-sr(ii))/dt).^2    +      ((si(ii+1)-si(ii))/dt).^2    )       *dt);
%             end
%             arcLen = cc';
%             
%             sr = fliplr(sr);
%             si = fliplr(si);
%             foldIndx2 = (80 +  80 - foldIndx) + 1;
%             cc = 0;
%             for ii=2:foldIndx2
%                 cc(ii,1) = cc(ii-1) + (sqrt(sr(ii).^2 + si(ii).^2)*dt);
%             end
%             for ii=foldIndx2+1:nn
%                 cc(ii,1) = cc(ii-1) - (sqrt(sr(ii).^2 + si(ii).^2)*dt);
%             end
%             arcLen2 = flipud(cc)';
%             
%             wf = linspace(1,0,nn);
%             wf2 = linspace(0,1,nn);
%             arcLen = ((arcLen.*wf)+(arcLen2.*wf2)) ./ (wf+wf2); 
%         
        %plot(arcLen)
        %hold on
        %plot(arcLen2)            
            %sr = sr - sr(1); % make sure starts integrating from zero
            %si = si - si(1) + 1;
            %arcLenL = cumsum(sqrt((  (gradient(sr,t)).^2   + (gradient(si,t)).^2  )  )) * dt;
        elseif type == 2
            arcLenL = abs(complex(sr+1i*si));
        elseif type == 3
            arcLenL = abs(complex((sr+1i*si)-1));
        end
        % flip over a the center point to make the descending I/O function
        %midPoint = 81; % center point of the I/O function
         %q1 = arcLenL(1:midPoint);
         %q2 = arcLenL(midPoint:end);
         %q2 = q2 - q2(1);
         %q2 = - q2;
         %q2 = q2 + q1(end);
         %arcLenL = [q1,q2];
        % intergrate from the right
        if type == 1
            sr = fliplr(sr);
            si = fliplr(si);
            sr = sr - sr(1); % make sure starts integrating from zero
            si = si - si(1);
            arcLenR = cumsum(sqrt((  (gradient(sr,t)).^2   + (gradient(si,t)).^2  )  )) * dt;
            %arcLenR = fliplr(cumsum(fliplr(sqrt((  (gradient(sr,t)).^2   + (gradient(si,t)).^2  )  ))) * dt);
            % flip over a the center point to make the descending I/O function
             arcLenR = fliplr(arcLenR);
             %q1 = arcLenR(1:midPoint);
             %q2 = arcLenR(midPoint:end);
             %q2 = q2 - q2(1);
             %q2 = - q2;
             %q2 = q2 + q1(end);
             %arcLenR = [q1,q2];
             %arcLenR = fliplr(arcLenR);
            %arcLen = mean([arcLenL;arcLenR],1);
            
       %arcLen = [arcLenL(1:foldIndx),arcLenR(foldIndx:1:end)];

            %hat = arcLenL(foldIndx);
            %r = arcLenL(foldIndx:end);
            %r = (-(r - hat)) + hat;
            %arcLen = [arcLenL(1:foldIndx),r];

        else
            arcLen = arcLenL;
        end
        meanSm = 2; 
        arcLen = meanSmoother(arcLen',meanSm)';


        % CALCULATE FOR NOISE ---------------------------------------------
        % integrate from the left
        if type == 1
            grad = sqrt((  (gradient(nr,t)).^2   + (gradient(ni,t)).^2  )  );
            k = ones(size(grad));
            k(foldIndx:end) = -k(foldIndx:end);
            q1 = cumsum(grad.*k) * dt;
            nr = fliplr(nr);
            ni = fliplr(ni);
            grad = sqrt((  (gradient(nr,t)).^2   + (gradient(ni,t)).^2  )  );
            k = ones(size(grad));
            k(foldIndx2:end) = -k(foldIndx2:end);
            q2 = fliplr(cumsum(grad.*k) * dt);
            wf1 = linspace(1,0,nn); % weighting function
            wf2 = linspace(0,1,nn);
            arcLen_n = ((q1.*wf1)+(q2.*wf2)) ./ (wf1+wf2); % weighted sum 

%              cc = 0;
%             nn = length(nr);
%             for ii=2:foldIndx
%                 cc(ii,1) = cc(ii-1) + (sqrt(nr(ii).^2 + ni(ii).^2)*dt);
%             end
%             for ii=foldIndx+1:nn
%                 cc(ii,1) = cc(ii-1) - (sqrt(nr(ii).^2 + ni(ii).^2)*dt);
%             end
%             arcLen_n = cc';
% 
%             nr = fliplr(nr);
%             ni = fliplr(ni);
%             cc = 0;
%             for ii=2:foldIndx2
%                 cc(ii,1) = cc(ii-1) + (sqrt(nr(ii).^2 + ni(ii).^2)*dt);
%             end
%             for ii=foldIndx2+1:nn
%                 cc(ii,1) = cc(ii-1) - (sqrt(nr(ii).^2 + ni(ii).^2)*dt);
%             end
%             arcLen2_n = flipud(cc)';
% 
%             arcLen_n = ((arcLen_n.*wf)+(arcLen2_n.*wf2)) ./ (wf+wf2); 
%             
            %nr = nr - nr(1); % make sure starts integrating from zero
            %ni = ni - ni(1);        
            %arcLenLn = cumsum(sqrt((  (gradient(nr,t)).^2   + (gradient(ni,t)).^2  )  )) * dt;
        elseif type == 2
            arcLenLn = abs(complex(nr+1i*ni));
        elseif type == 3
            arcLenLn = abs(complex(nr+1i*ni)-1);
        end
        % flip over a the center point to make the descending I/O function
         %q1 = arcLenLn(1:midPoint);
         %q2 = arcLenLn(midPoint:end);
         %q2 = q2 - q2(1);
         %q2 = - q2;
         %q2 = q2 + q1(end);
         %arcLenLn = [q1,q2];
        % intergrate from the right 
        if type == 1
%             nr = fliplr(nr);
%             ni = fliplr(ni);
%             nr = nr - nr(1); % make sure starts integrating from zero
%             ni = ni - ni(1);
%             arcLenRn = cumsum(sqrt((  (gradient(nr,t)).^2   + (gradient(ni,t)).^2  )  )) * dt;
%             %arcLenRn = fliplr(cumsum(fliplr(sqrt((  (gradient(nr,t)).^2   + (gradient(ni,t)).^2  )  ))) * dt);
%             % flip over a the center point to make the descending I/O function
%              arcLenRn = fliplr(arcLenRn);
             %q1 = arcLenRn(1:midPoint);
             %q2 = arcLenRn(midPoint:end);
             %q2 = q2 - q2(1);
             %q2 = - q2;
             %q2 = q2 + q1(end);
             %arcLenRn = [q1,q2];
             %arcLenRn = fliplr(arcLenRn);
            %arcLen_n = mean([arcLenLn;arcLenRn],1);
     %arcLen_n = [arcLenLn(1:foldIndx),arcLenRn(foldIndx:1:end)];

            %hat = arcLenLn(foldIndx);
            %r = arcLenLn(foldIndx:end);
            %r = (-(r - hat)) + hat;
            %arcLen_n = [arcLenLn(1:foldIndx),r];

        else
            arcLen_n = arcLenLn;
        end
        arcLen_n = meanSmoother(arcLen_n',meanSm)';

        Q = [Q;arcLen];
        Qn = [Qn;arcLen_n];
    end
    if type == 1
        Delta = (Q-Qn)';
        %nf = std(Delta,[],2) / sqrt(size(Delta,2));
        delta = mean(Delta,2);
        nf1 = std(Q',[],2) / sqrt(size(Q',2));
        nf2 = std(Qn',[],2) / sqrt(size(Qn',2));
        nf = mean([nf1,nf2],2);
    else
        Delta = Q';
        DeltaN = Qn';
        delta = mean(Delta,2);
        %nf = std(DeltaN,[],2) / sqrt(size(DeltaN,2));
        nf = mean(DeltaN,2);
    end

    %     plot(Delta,'Color',[.7 .7 .7])
%     hold on
%     plot(delta,'b','LineWidth',2)
%     plot(nf*2,'k','LineWidth',2)
%     pause
%     c

    delta = delta(:);
    nf = nf(:);
end

function [nLoops,clickTrain,noiseSamples,nReps,H] = getStimuli(fs,modLen)
    % create click stimuli ---------------------------------------------------------
    %fs = obj.fs; % get the system sampling rate
    clickN = 1; % number of samples in the click
    clickLen = 0.05; % desired stimulus length in seconds
    nSamples = round(clickLen * fs); % total number of samples
    if mod(nSamples,2) ~= 0 % make total number of samples even
        nSamples = nSamples + 1;
    end
    click = zeros(nSamples,1);
    startDelay = 0.003; % silence before the click onset (s)
    startSample = round(startDelay * fs);
    click(startSample:startSample + (clickN-1)) = 1;
     dBatten = 23; % dB attenuation re: full out
     multiplier = 10^(-dBatten/20);
     click = click * multiplier;
    % Paige's threshold is 58 dB attenuation re: full output
    % Presentation level = 35 dB SL = 58 - 35 = 23 dB atten
    % IEC711 coupler results: 
    %   0.03 Pa peak
    %   0.03 + 0.024 = 0.054 peak to peak
    %   0.054 / 2 = 0.027 peak equivalent
    %   63.5 dB pSPL
    %   68.6 dB peak to peak SPL
    %   62.6 dB peSPL


    % create noise ------------------------------------------------------------
    clicksPerCondition = 1000; % total number of clicks recorded per condition (noise or no noise)
    %modLen = 8; % modulation duration in sec (from onset of noise to offset)
    loopDuration = modLen * 30; %240; %60; % seconds per loop
    clicksPerLoop = loopDuration; % because one second gives you 1 click at each level
    nLoops = ceil(clicksPerCondition / clicksPerLoop); % number of loops
    noiseSamples = round(modLen * fs); % total 
    nReps = loopDuration / (modLen * 2); % number of repeitions per loop
    %nReps = loopDuration / (modLen * 4); % number of repeitions per loop
    %h = linspace(0,70,noiseSamples/2)';
    h = linspace(0,70,noiseSamples)';
    h(1) = eps;
    h = [h;flipud(h)];
    h = 10.^(h/20); % noise amplitude in linear units
    h = h / max(abs(h)) * 1; % rescale to unit amplitude
    H = repmat(h,1,nReps);
    dBattenN = 1; % dB attenuation re: full out
    multiplierNoise = 10^(-dBattenN/20);
    % Paige's threshold is 65 dB attenuation re: full output
    % Presentation level = full output - 1 dB = 64 dB SL
    % IEC711 coupler results: 
    %   1.3 Pa peak
    %   96.3 dB SPL at peak intensity

    % finish making click train
    clicksPerTrain = (noiseSamples / nSamples)*2;
    clickTrain = repmat(click,1,clicksPerTrain);
    clickTrain = clickTrain(:);
end

function [maxTime,maxIndx,dx,fitresult,gof] = createFit(xx,qq,doPlot)
    % fit with a 3-gaussian mixture
    if nargin < 3
        doPlot = 0;
    end
    [xData, yData] = prepareCurveData(xx,qq);
    % Set up fittype and options.
    ft = fittype( 'gauss3' );
    opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
    opts.Display = 'Off';
    opts.Lower = [-Inf -Inf 0 -Inf -Inf 0 -Inf -Inf 0];
    opts.Robust = 'Bisquare';
    opts.StartPoint = [1 164 9.60201287675777 0.877600381176352 183 16.355870659702 0.867957144362498 145 8.69431682989013];
    % Fit model to data.
    [fitresult, gof] = fit(xData,yData,ft,opts);
    
    a1 = fitresult.a1;
    b1 = fitresult.b1;
    c1 = fitresult.c1;
    a2 = fitresult.a2;
    b2 = fitresult.b2;
    c2 = fitresult.c2;
    a3 = fitresult.a3;
    b3 = fitresult.b3;
    c3 = fitresult.c3;
    step = xx(2)-xx(1); % original sampling (50 ms steps)
    %nn = length(xx);
    % we want 1 ms steps sizes, so increase sampling by 50 times
    newStep = step / 50; % sample
    xxx = (xx(1):newStep:xx(end))';
    dx = xxx(2)-xxx(1);
    yyy = a1*exp(-((xxx-b1)./c1).^2) + a2*exp(-((xxx-b2)./c2).^2) + a3*exp(-((xxx-b3)./c3).^2);
    [~,maxIndx] = max(yyy);
    maxTime = xxx(maxIndx);

    if doPlot == 1
        % Plot fit with data.   
        %figure
        plot(xx,qq,'.b')
        hold on
        plot(xxx,yyy,'r-')
        grid on
        %h = plot( fitresult, xData, yData );
        %legend( h, 'qq vs. xx', 'untitled fit 1', 'Location', 'NorthEast', 'Interpreter', 'none' );
        % Label axes
        %xlabel( 'xx', 'Interpreter', 'none' );
        %ylabel( 'qq', 'Interpreter', 'none' );
        %grid on
    end
end
function [B,yhat,residuals] = ARLas_wlsf(y,X,w)
    % weighted least-squares fit
    % Ordinary Least-Squares Regression for the model
    % y = b1*X1 + b2*X2 + ... + bk*Xk.
    % X = matrix of independent variables (your "solution matrix" of signals of interest)
    % y = vector of dependent variables (your recorded data)
    % B = coefficients
    % yhat = predicted values
    % residuals = residuals

    [rows,columns] = size(X); % size of the solution matrix
    sqrtw = sqrt(w); % and take the square root of the weighting vector
    X = repmat(sqrtw,1,columns).*X; % replicate sqrtw to matrix the size of X; J is the weighted version of X
    y = sqrtw.*y; % Dy is the weighted version of y
    [Q,R] = qr(X,0); % orthogonal-triangular decomposition; R is the Cholesky factor of the X matrix
                     % the input argument 0 makes an "economy sized" decomposition, so that
                     % [nSamples,nIVs] = size(Q), and 
                     % [nIvs,nIVs] = size(R).
    B = full(R\(Q'*y)); % Same as p = D*X\(D*y); b is a vector of coefficients; also same as b = (X'*X)\(X'*Y).
    yhat = X*B;      % Predicted responses at each data point.
    residuals = y - yhat;
end

% OLD CODE ----------------------------------------------------------------
% function [m] = ARfix(m)
%     % fix artifact problems
%     % check for and remove single-sample artifacts
%     mr = real(m);
%     mi = imag(m);
%     [rows,cols] = size(mr);
%     mr = mr(:).';
%     mi = mi(:).';
%     nn = length(mr);
%     doPlot = 0;
%     [indxAR,nRejects] = AR(mr','severe',doPlot);
%     if ~isempty(indxAR) % if there are rejections to make
%         for rindx = 1:nRejects
%             if indxAR(rindx) == 1
%                 mr(indxAR) = mean([mr(indxAR+1),mr(indxAR+2)]);
%             elseif indxAR(rindx) == nn
%                 mr(indxAR) = mean([mr(indxAR-1),mr(indxAR-2)]);
%             else
%                 mr(indxAR) = mean([mr(indxAR-1),mr(indxAR+1)]);
%             end
%         end
%     end
%     [indxAR,nRejects] = AR(mi','mild',doPlot);
%     if ~isempty(indxAR) % if there are rejections to make
%         for rindx = 1:nRejects
%             if indxAR(rindx) == 1
%                 mi(indxAR) = mean([mi(indxAR+1),mi(indxAR+2)]);
%             elseif indxAR(rindx) == nn
%                 mi(indxAR) = mean([mi(indxAR-1),mi(indxAR-2)]);
%             else
%                 mi(indxAR) = mean([mi(indxAR-1),mi(indxAR+1)]);
%             end
%         end
%     end
%     m = mr + 1i*mi;
%     m = reshape(m,rows,cols);
% end
%     arcLen_sm = Q-Qn;
%     arcLen_nf = std(arcLen_sm,[],1) / sqrt(size(arcLen_sm,1));
%     arcLen_sm = mean(arcLen_sm,1);

% 
    % A two second rise fall has four seconds, so 80 clicks, of which 6
    % (3.75 from each edge) form the template.
    
%     sigStart = 1; % signal start index
%     sigFinish = 15; % signal finish index
%     windowN = round(fs*0.01);
% 
%     startOffset = 1; %13;
%     finishOffset = 60; % 60
%     whichOne = 1; % which fft bin
%     counter = 1;

%     for ii=1:size(Data,2)
%         offset = ii*.001;
%         runme(template,Data(:,ii),nSamples,nClicks,indx,time,fs,iscS1L,offset)
%     end

% we were going to convert to FPL and RPL, but doesn't really make any
% difference? And may increase likelihood of errors from calibration?
%     tic
%     for kk=1:size(Data,2)
%         disp(['Converting ',num2str(kk),' of ',num2str(size(Data,2))])
%         [pl,Pl,phi,other,wf] = ARLas_convertPL(Data(:,kk),iscS1L);
%         FPL(:,kk) = wf.fpl;
%         RPL(:,kk) = wf.rpl;
%     end
%     PRR = RPL ./ FPL;
%     Data = PRR;
%     toc / 60

%delta(:,jj) = mean(Chunk,2) - template;
%[frequency,signal(:,jj),noiseFloor(:,jj)] = ARLas_fda(chunk,fs,0.00002,fs);
%S(:,jj) = fft(mean(chunk,2));
%rms = sqrt(mean(delta.^2,1));
%     if ~exist('iscS1L','var') % in-situ calibration data
%         iscS1L = [];
%     end

% function [] = crawl()
%     t = timeChunk;
%     dt = median(gradient(t));
%     for jj=1:length(peakIndx)
%         x = (complex(Delta(peakIndx(jj),:)));
%         
%         xr = real(x);
%         xi = imag(x);
%         smoothing = 0.999999;
%         ppr = csaps(t,xr,smoothing);
%         ppi = csaps(t,xi,smoothing);
%         xr_sm = ppval(ppr,t);
%         xi_sm = ppval(ppi,t);
%         resid_xr = xr - xr_sm; % residuals after subtracting the smoothing spline
%         resid_xi = xi - xi_sm;
%         ppr = csaps(t,resid_xr,smoothing);
%         ppi = csaps(t,resid_xi,smoothing);
%         xr_noise = ppval(ppr,t);
%         xi_noise = ppval(ppi,t);
%     
%     
%     
%     
%         int = cumsum(( (gradient(xr,t))   + (gradient(xi,t))  )) * dt;
%         int_smL = cumsum(( (gradient(xr_sm,t))   + (gradient(xi_sm,t))  )) * dt;
%         int_smR = fliplr(cumsum(fliplr(( (gradient(xr_sm,t))   + (gradient(xi_sm,t))  ))) * dt);
%     
%     
%     
%         %arcLen = cumsum(sqrt((  (gradient(xr,t)).^2   + (gradient(xi,t)).^2  )  )) * dt;
%         % integrate from the left
%             arcLen_smL = cumsum(sqrt((  (gradient(xr_sm,t)).^2   + (gradient(xi_sm,t)).^2  )  )) * dt;
%             arcLen_noiseL = cumsum(sqrt((  (gradient(xr_noise,t)).^2   + (gradient(xi_noise,t)).^2  )  )) * dt;
%         % intergrate from the right    
%             arcLen_smR = fliplr(cumsum(fliplr(sqrt((  (gradient(xr_sm,t)).^2   + (gradient(xi_sm,t)).^2  )  ))) * dt);
%         % find where the two functions cross: this is the peak
%             indxLT = find(arcLen_smL <= arcLen_smR);
%             [~,peakIndx2] = max(indxLT);
%             peak = arcLen_smL(peakIndx2);
%         % flip the descending side
%             arcLen_sm = [arcLen_smL(1:peakIndx2),arcLen_smR(peakIndx2+1:end)];
%             arcLen_sm(peakIndx2) = (arcLen_sm(peakIndx2-1) + arcLen_sm(peakIndx2+1))/2;
%         
%         plot(arcLen_smL,'b*-')
%         hold on
%         plot(arcLen_smR,'r*-')
%         plot(peakIndx2,peak,'go')
%     
%     end
% end



% % plotting -------------------------------------------------------
%     h10 = figure(10); hold on % plot complex data at desired frequencies
%     for jj=1:length(peakIndx)
%         plot(complex(Delta(peakIndx(jj),:)),'.-')
%     end
%     NN = 1000; % number of samples in the ellipses
%     t = linspace(0,2*pi,NN); 
%     eMinus0 = [cos(t);sin(t)]'; % create a unit circle
%     hold on; box on;
%     plot(eMinus0(:,1),eMinus0(:,2),'LineWidth',1,'Color',[0 0 0])
%     % plot the 1 to 6 dB reduction lines
%     eMinus1 = eMinus0 * (10^(-1/20));
%     eMinus2 = eMinus0 * (10^(-2/20));
%     eMinus3 = eMinus0 * (10^(-3/20));
%     eMinus4 = eMinus0 * (10^(-4/20));
%     eMinus5 = eMinus0 * (10^(-5/20));
%     eMinus6 = eMinus0 * (10^(-6/20));
%     innerColor = [0.5 0.5 0.5];
%     LW = 1; % line width
%     plot(eMinus1(:,1),eMinus1(:,2),'LineStyle','--','LineWidth',LW,'Color',innerColor)
%     plot(eMinus2(:,1),eMinus2(:,2),'LineStyle','--','LineWidth',LW,'Color',innerColor)
%     plot(eMinus3(:,1),eMinus3(:,2),'LineStyle','--','LineWidth',LW,'Color',innerColor)
%     plot(eMinus4(:,1),eMinus4(:,2),'LineStyle','--','LineWidth',LW,'Color',innerColor)
%     plot(eMinus5(:,1),eMinus5(:,2),'LineStyle','--','LineWidth',LW,'Color',innerColor)
%     plot(eMinus6(:,1),eMinus6(:,2),'LineStyle','--','LineWidth',LW,'Color',innerColor)
%     % Format the plot
%     Lh = line([-1 1],[0 0]); 
%     set(Lh,'Color',[0 0 0],'LineStyle','-','LineWidth',1)
%     Lv = line([0 0],[-1 1]);
%     set(Lv,'Color',[0 0 0],'LineStyle','-','LineWidth',1)
%     axis square
%     legend(['f=',num2str(round(freq(2))),' Hz'],['f=',num2str(round(freq(3))),' Hz'],['f=',num2str(round(freq(4))),' Hz'],['f=',num2str(round(freq(5))),' Hz'],['f=',num2str(round(freq(6))),' Hz'])
%     xlim([0.5 1.3])
%     ylim([-0.4 0.4])
%     xlabel('Real')
%     ylabel('Imaginary')
% 
%     %-----------------------------------
%     CCC = [0 0 0
%         1 0 0
%         .5 .5 0
%         0 .5 .5
%         0 1 0
%         0 0 .5
%         0 0 1];
%     h11 = figure(11); hold on
%     for jj=1:length(peakIndx)
%         delta = (complex(Delta(peakIndx(jj),:)));
%         if max(delta) >=1 %max(abs(Delta(jj,:)))>=1
%             multiplier = 1;
%         else
%             multiplier = -1;
%         end
%         %plot(timeChunk,multiplier*20*log10(abs(complex(Delta(peakIndx(jj),:))-1)),'.-') %,'Color',CCC(jj,:)
%         %plot(timeChunk,multiplier*(abs(complex(Delta(peakIndx(jj),:))  )),'.-')
%         plot(timeChunk,multiplier*(abs(complex(delta-1))),'.-')
% 
%         multiplier
%     end
%     %legend(['f=',num2str(round(freq(2))),' Hz'],['f=',num2str(round(freq(3))),' Hz'],['f=',num2str(round(freq(4))),' Hz'],['f=',num2str(round(freq(5))),' Hz'],['f=',num2str(round(freq(6))),' Hz'])
%     legend(num2str(peakIndx))
% 
% q = 20*log10(abs(complex(Delta(peakIndx(5),:))));
% A = 1;
% B = 1;
% %obj = fitIBeta(timeChunk,q,A,B);
% obj = [];

%         for ii=1:size(m,2)
%             mm = m(:,ii);
%             [CAL(:,ii),CALn(:,ii)] = crawl(timeChunk,mm);
%         end
%         CALd = CAL - CALn;
%         cald = mean(CALd,2);
%         %caln = mean(CALn,2);
%         cal_nf = std(CAL,[],2) / sqrt(size(CAL,2));
%         cald_nf = std(CALd,[],2) / sqrt(size(CALd,2));
%         caln_nf = std(CALn,[],2) / sqrt(size(CALn,2));
        

%         %plot(caln+2*caln_nf,'b')
%         %hold on
%         %plot(caln-2*caln_nf,'b')
%         %plot(caln,'b--')
%         %plot(cal,'r')
%         plot(cald,'b')
%         hold on
%         plot(cal_nf*2,'k--')
%         %plot(caln,'r')
%         plot(caln_nf,'k')
%         plot(caln_nf*2,'g')
%         hold off
% 
%         plot(CAL,'Color',[.7 .7 .7])
%         hold on
%         plot(cal,'r','LineWidth',2)
%         plot(cal_nf,'k')
%         keyboard
%  [frequency,signal,noiseFloor] = ARLas_fda(x,96000,0.00002);
% plot(frequency,signal)
% hold on
% plot(frequency,noiseFloor,'k')
% plot(frequency,noiseFloor+6,'k')
%     mask = ones(size(xr));
%     for ii=2:2:size(mask,2)
%         mask(:,ii) = -1;
%     end
%     %mask = mask(:,1:14);

    %bar_xr = mean(xr,2);
    %noi_xr = mean(xr.*mask,2);

    %xr = medianSmoother(xr,2);
    %xi = medianSmoother(xi,2);
    %smoothing = 0.999; % 0.999999;
    %ppr = csaps(t,xr,smoothing);
    %ppi = csaps(t,xi,smoothing);
    %xr_sm = ppval(ppr,t);
    %xi_sm = ppval(ppi,t);
% pdr = makedist('normal','mu',0,'sigma',std(xr_n));
% pdi = makedist('normal','mu',0,'sigma',std(xr_n));
% xr_sm = reshape(xr_sm,rows,cols);
% xi_sm = reshape(xi_sm,rows,cols);
% xr_sm = mean(xr_sm,2)';
% xi_sm = mean(xi_sm,2)';
% xr_sm = xr_sm - xr_sm(1);
% xi_sm = xi_sm - xi_sm(1);

%     nnn = length(xr_n);
%     indxR = rand(n,1);
%     [~,indx] = sort(indxR);
%     xr_n = xr_n(indx);
%     indxR = rand(n,1);
%     [~,indx] = sort(indxR);
%     xi_n = xi_n(indx);
%     
%     qr_n = pdr.random(1,nnn);
%     qi_n = pdi.random(1,nnn);
% 
%     qr_n = medianSmoother(qr_n',5);
%     qr_n = meanSmoother(qr_n,2)';
%     qi_n = medianSmoother(qi_n',5);
%     qi_n = meanSmoother(qi_n,2)';
%     qr_n = reshape(qr_n,rows,cols);
%     qi_n = reshape(qi_n,rows,cols);
%     qr_n = mean(qr_n,2)';
%     qi_n = mean(qi_n,2)';

%         % find where the two functions cross: this is the peak
%         indxLT = find(arcLen_smL <= arcLen_smR);
%         %[~,peakIndx2] = max(indxLT);
%         peakIndx2 = 82;
%         peak = arcLen_smL(peakIndx2);
%     % flip the descending side
%         arcLen_sm = [arcLen_smL(1:peakIndx2),arcLen_smR(peakIndx2+1:end)];
%         %[~,peakIndx3] = max(arcLen_sm);
%         peakIndx3 = 81;
%         arcLen_sm(peakIndx3) = mean([arcLen_sm(peakIndx3-1),arcLen_sm(peakIndx3),arcLen_sm(peakIndx2+1)]);

% 
%     % now same thing for the noise
%     % integrate from the left
%         arcLen_smL = cumsum(sqrt((  (gradient(nr,t)).^2   + (gradient(ni,t)).^2  )  )) * dt;
%         %arcLen_smL = cumsum(((  (gradient(xr.'-xr_sm,t))   + (gradient(xi.'-xi_sm,t))  )  )) * dt;
%         q1 = arcLen_smL(1:81);
%         q2 = arcLen_smL(81:end);
%         q2 = q2 - q2(1);
%         q2 = - q2;
%         q2 = q2 + q1(end);
%         arcLen_smL = [q1,q2];
% 
% %         % intergrate from the right
% %         arcLen_smR = fliplr(cumsum(fliplr(sqrt(( (gradient(nr,t)).^2   + (gradient(ni,t)).^2  )  ))) * dt);
% %         %arcLen_smR = fliplr(cumsum(fliplr(((  (gradient(xr.'-xr_sm,t))   + (gradient(xi.'-xi_sm,t))  )  ))) * dt);
% %     % find where the two functions cross: this is the peak
% %         %indxLT = find(arcLen_smL <= arcLen_smR);
% %     % flip the descending side
% %         arcLen_n = [arcLen_smL(1:peakIndx2),arcLen_smR(peakIndx2+1:end)];
% %         %arcLen_n = arcLen_smL;
% %         arcLen_n(peakIndx3) = mean([arcLen_n(peakIndx3-1),arcLen_n(peakIndx3),arcLen_n(peakIndx2+1)]);
    %     plot(arcLen_smL,'b*-')
    %     hold on
    %     plot(arcLen_smR,'r*-')
    %     plot(peakIndx2,peak,'go')
    %    plot(arcLen_sm,'k.-')
    %    keyboard
    % plot(arcLen_sm,'b')
    % hold on
    % plot(arcLen_n,'k')
    % hold off


%     for jj=1:nClicks % loop across each click position in the sweep
%         Chunk = Data(indx(jj):indx(jj)+nSamples,:); % current chunk to analyze (contains the click)
%         FT = fft(Chunk,fs);
%         %S(:,jj) = fft(mean(chunk,2),fs);
%     end
%     
%     Template = [S(:,1:3),S(:,nClicks-2:nClicks)]; % template in the frequency domain (template = first 3 and last 3 clicks)
%     Template = mean(Template,2);
%     for jj=1:length(indx) 
%         Delta(:,jj) = S(:,jj)./Template;
%     end
%     freq = (0:1:size(chunk,1)-1)*(fs/size(chunk,1));
% 
% fff = (1:1:fs)-1'; % full frequency vector
% % find the locatino of the desired center frequencies
% for ii=1:length(peakFreqs)
%     [~,peakIndx(ii)] = min(abs(fff-peakFreqs(ii)));
% end
%peakIndx = [349,933,1619,3312, 4039,7294,14109];    

%-------- keyboard
% for jj=1:length(peakIndx)-1
%         D = complex(Delta(peakIndx(jj):peakIndx(jj+1),:));
%         f = fff(peakIndx(jj):peakIndx(jj+1));
%         [~,fhatIndx] = max(mean(abs(D-1),2));
%         peakFreqs2(jj) = f(fhatIndx);
% end

%     % reshape the click stimuli into 100 ms chunks
%     [rows,cols] = size(clickTrain);
%     chunkLen = 0.1; % chunk size (sec)
%     chunkSize = round(chunkLen*fs); % chunk size (samples)
%     nChunks = floor(rows/chunkSize);
%     time = (0:1:length(clickTrain)-1)'/fs;
%     ClickTrain = reshape(clickTrain,chunkSize,nChunks);
%     clickIndx = find(clickTrain>.01); % location of clicks
%     S1L = ClickTrain;
%     S1R = S1L * 0;
%   % generate noise -----
%     multiplierNoise = 0.25;
%     [rows,cols] = size(S1L);
%     noise = randn(rows,cols*2);
%     noise = noise(:)';
%     [indx,nRejects] = AR(noise,'moderate',0);
%     noise = AR_engine(noise,indx);
%     noise = noise(1:rows*cols);
%     noise = noise(:);
%     h = H(:,1);
%     Noise = noise .* h;
%     Noise = Noise * multiplierNoise;
%     Noise = reshape(Noise,chunkSize,nChunks);
%     S2R = Noise;
%     S2L = Noise * 0;

%     dt = median(gradient(t));
% 
%     i1 = 79;
%     i2 = 89;
%     x = abs(mean(xrs,2) + 1i*mean(xis,2));
%     p = polyfit(t(i1:i2),x(i1:i2),2);
%     timeShift = roots(polyder(p)) - 4;
%     %gof = [];
%     snr = max([sn_r,sn_i]);
%     return
% 
%     % find the overall delay
%     tri = [(1:1:80)';(80:-1:1)']; % triangle shape representing growth of the noise
%     q = mean(xrs,2)+1i*mean(xis,2);
%     qqq = abs(xcorr(q,tri));
%     rrr = abs(xcorr(tri,tri));
%     qqq = qqq / max(qqq);
%     rrr = rrr / max(rrr);
%     % the max indx when you autocorrelate tri is 160.
%     xxx = (0:1:length(qqq)-1)'*dt;
%     % look across a limited time window (2.5 ms equivalent) from
%     i1 = 140; %110; % starting index of time window 
%     i2 = 180; %210; % ending index of time window
%     doPlot = 0;
%     [maxTimeQ,maxIndxQ,dx,fitresultQ,gofQ] = createFit(xxx(i1:i2),qqq(i1:i2),doPlot);
%     [maxTimeR,maxIndxR,dx,fitresultR,gofR] = createFit(xxx(i1:i2),rrr(i1:i2),doPlot);
%     indxShift = maxIndxQ - maxIndxR;
%     timeShift = maxTimeQ - maxTimeR;
%     % index shift can only be positive (a time delay) -- you can't go the othe way
%     if indxShift < 0
%         %indxShift = 0;
%         timeShift = 0;
%     end
%     % time shift should not be > 300 ms. If it is, play conservatively and
%     % assign a shift of 150 ms
% %     if timeShift > .300
% %         timeShift = .150;
% %     end
%     gof = min([gofQ.adjrsquare,gofR.adjrsquare]); % take the minimum R^2 as a goodness of fit measure for later weighting

            
%     % compare changes in clicks by making a "Template".
%     % Get the template by averaging the first 3 clicks and last 3 clicks.
%     % NOTE: The MEMR has hysteresis, so that you can't also average across
%     % multiple sweeps! Sweeps can show drift. Also, you aren't guaranteed
%     % to return to baseline by the end of the sweep, so you should probably
%     % avoid just averaging the start and end also?
%     % Each click window is 50 ms long, so the average of 3 is 150 ms.
%     templateN = 3; % number of clicks (from both beginning and end) to include in the template
%     for jj=1:templateN 
%         chunk = Data(indx(jj):indx(jj)+nSamples,:);
%         % artifact rejection:
%         hat = max(chunk,[],1); % get the maximum value in each sweep
%         doPlot = 0;
%         [indxAR,nRejects] = AR(hat(:)','mild',doPlot);
%         if ~isempty(indxAR) % if there are rejections to make
%             chunk(:,indxAR) = [];
%         end
%         template(:,jj) = mean(chunk,2);
%     end
%     counter = 4;
%     for jj=nClicks-2:nClicks % loop across click
%         chunk = Data(indx(jj):indx(jj)+nSamples,:);
%         % artifact rejection:
%         hat = max(chunk,[],1); % get the maximum value in each sweep
%         doPlot = 0;
%         [indxAR,nRejects] = AR(hat(:)','mild',doPlot);
%         if ~isempty(indxAR) % if there are rejections to make
%             chunk(:,indxAR) = [];
%         end
%         template(:,counter) = mean(chunk,2);
%         counter = counter + 1;
%     end
%     templateNoise = std(template,[],2)/sqrt(size(template,2));
%     template = mean(template,2);
%     % average SNR in the template window
%     templateSNR = 20*log10(mean(abs(template)./templateNoise));
