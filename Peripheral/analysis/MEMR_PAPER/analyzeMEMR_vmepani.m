function [MEMR_inc,MEMR_mem] = analyzeMEMR_vmepani(header,Clicks,headerN,Noise,subjectName)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [analysisPart1_inc,analysisPart1_mem] = analyzeMEMR_vmepani(pathName,subjectName,runNumber)
%
% Analyze MEMR data, mepani test.
%
% Author: Shawn Goodman
% Date: August 6, 2022
% Last Updated: January 26, 2024 -- ssg
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
    % ---------------------------------------------------------------------
    doJitterFix = 1;
    % ---------------------------------------------------------------------

    % extract the raw recordings -----
    micCorrection = header.userInfo.micCorrectionL;
    Clicks = applyMicCorrection(Clicks,micCorrection);
    micCorrectionN = headerN.userInfo.micCorrectionR;
    Noise = applyMicCorrection(Noise,micCorrectionN);
    %[pl,Pl,phi,other,wf] = ARLas_convertPL(Clicks(:,1),header.userInfo.iscS1L);

    
    % Jitter Fix! ---------------------------------------------------------
    % if strcmp(subjectName,'MEM10') & runNumber == 1
    %     doJitterFix = 0;
    % end
    % if strcmp(subjectName,'MEM37') & runNumber == 1
    %     doJitterFix = 0;
    % end
    % if strcmp(subjectName,'MEM12') & runNumber == 1
    %     doJitterFix = 0;
    % end
  



    fs = header.fs;
    nSamples = header.nSamples;
    nReps = header.userInfo.nReps;
    nSweeps = header.userInfo.nSweeps;
    noiseTarget = header.userInfo.targetLevelNoise;
    probeEar = header.userInfo.probeEar;
    elicitorEar = header.userInfo.activatorEar;
    templateClick = header.userInfo.stimClick;
    templateNoise = header.userInfo.stimNoise;

    totalSamples = nSamples * nReps; % total samples in 12 sweeps
    nSamples = totalSamples / nSweeps; % number of samples in each sweep
    sweepLength = nSamples / fs; % length of each sweep (s)
    Clicks = reshape(Clicks,nSamples,nSweeps);
    Noise = reshape(Noise,nSamples,nSweeps);

    Clicks = Clicks(1:51840*2,:);
    Noise = Noise(1:51840*2,:);

    if doJitterFix == 1
        Clicks = unjitter(Clicks,subjectName);
    end

    time = (0:1:size(Clicks,1)-1)'/fs; % time vector
    %indx = header.userInfo.clickIndx; % location of the clicks
    indx = [482,52322];
    clicks = mean(Clicks,2);
    clicks = clicks(1:500);
    [~,maxIndx] = max(clicks);
    if abs((maxIndx - indx(1))) > 2
        keyboard
    end
    clicks = mean(Clicks,2);
    clicks = clicks(52200:52400);
    [~,maxIndx] = max(clicks);
    maxIndx = maxIndx + 52200 -1;
    if abs((maxIndx - indx(2))) > 2
        keyboard
    end
    
    clear data %header
    elicitor = Noise(:,4); 
    % ----> what is special about this? Nothing, we just chose an exemplar for plotting


    % When considering which portions of each time window to analyze:
    % Within each click window, we want to include several round trip
    % travel times of the click in the canal, but no low-frequency OAEs. We
    % should consider only the first 2 ms of the click, which is 
    % round(0.002*fs) = 192 samples. 
    nClicks = length(indx); % number of total clicks to analyze (2)

    % % location of noise levels
    % modLen = 4; % modulation duration in sec (from onset of noise to offset)
    % [nLoops,clickTrain,noiseSamples,nReps,H] = getStimuli(fs,modLen);
    % % Here, H is the (linear) amplitude vector for the noise. It is 8
    % % seconds long total--4 seconds rise and 4 seconds fall.
    % % The noise floor in the ear canal is roughly 60 dB SPL.
    % % The noise itself goes from 50 dB SPL to 120, covering a change of 70
    % % dB in 4 seconds, so a rate of 17.5 dB/second
    % h = H(:,1)*(10^(120/20)*.00002); % h is the vector of noise levels (intended, in Pa)
    % %plot(20*log10(abs(Noise(:,1))/.00002))
    % %hold on
    % %plot(20*log10(h/.00002),'r','LineWidth',2)
    %noiseLvl = noiseTarget; % noise level (intended, dB SPL)
    Clicks = ARLas_hpFilter(Clicks,fs,100);


    % get the stimulus levels ---------------------------------------------
    t = (0:1:length(elicitor)-1)'/fs ;
    E = elicitor(4450:51000); % elicitor waveform
    C = Clicks(indx(1)-100:indx(1)+100,:); % click waveform
    C = median(C,2);
    RMS = sqrt(mean(E.^2)); % elicitor RMS
    RMS = 20*log10(RMS/.00002); % elicitor RMS in dB SPL

    RMSC = sqrt(mean(C.^2)); % rms of the click
    RMSC = 20*log10(RMSC / .00002); % rms of the click in dB SPL
    pSPL = 20*log10(max(abs(C))/.00002); % level of the click in pSPL
    

    %----------------------------------------------------------------------

    % analyze MEMR --------------------------------------------------------
    stabilityCheck = 0;
    MEMR_mem = runme(Clicks,nClicks,indx,time,fs,stabilityCheck,subjectName);
    % analyze ear canal stability
    stabilityCheck = 1;
    MEMR_inc = runme(Clicks,nClicks,indx,time,fs,stabilityCheck,subjectName);

    MEMR_mem.elicitor = elicitor;
    MEMR_inc.elicitor = elicitor;
    MEMR_mem.RMS = RMS; % rms of the elicitor Pa
    MEMR_mem.RMSC = RMSC; % rms of the click Pa
    MEMR_mem.pSPL = pSPL; % click peak in pSPL
    MEMR_mem.noiseTarget = noiseTarget; % intended elicitor level from header file
   

    %  plotting -----------------------------------------
    % try
    %     h1 = figure;
    %     sp1 = subplot(2,1,1);
    %     sp2 = subplot(2,1,2);
    %     sp1.Position = [0.1300    0.7230    0.7750    0.2020];
    %     sp2.Position = [0.1300    0.1100    0.7750    0.4924];
    %   axes(sp1)
    %     ax = gca;
    %     ax.FontSize = 11; 
    %     plot(MEMR_mem.timeTrend,20*log10(MEMR_mem.trend),'b')
    %     hold on
    %      plot(MEMR_mem.timeTrend,20*log10(MEMR_inc.trend),'b--')
    %     %Q = repmat(20*log10(MEMR_mem.d1),1,15);
    %     %Q = Q(:);
    %     %plot(MEMR_mem.timeTrend,20*log10(MEMR_mem.trend)+Q,'k')
    %     xlabel('Time (s)','FontSize',11)
    %     ylabel('Change (dB)','FontSize',11)
    %     title([subjectName,'    Elicitor Level ',num2str(noiseLvl),' dB SPL rms'])
    %     grid on
    %   axes(sp2)
    %     ax = gca;
    %     ax.FontSize = 11; 
    %     plot(MEMR_mem.t,20*log10(MEMR_mem.d1),'k')
    %     hold on
    %     plot(MEMR_mem.t,20*log10(MEMR_inc.d1),'k--')
    %     ymax = max(20*log10(MEMR_mem.d1))*1.1;
    %     line([4 4],[0 ymax],'Color',[0 0 0],'LineWidth',0.5,'LineStyle',':')
    %     plot(MEMR_mem.peakTime,20*log10(MEMR_mem.peakAmp),'r.','MarkerSize',14)
    %     line([MEMR_mem.peakTime,MEMR_mem.peakTime],[0 20*log10(MEMR_mem.peakAmp)],'Color',[1 0 0],'LineWidth',0.5,'LineStyle','-')
    %     line([MEMR_mem.thdOnsetTime,MEMR_mem.thdOnsetTime],[0 20*log10(MEMR_mem.thdAmp)],'Color',[1 0 0],'LineWidth',0.5,'LineStyle','-')
    %     line([MEMR_mem.thdOffsetTime,MEMR_mem.thdOffsetTime],[0 20*log10(MEMR_mem.thdAmp)],'Color',[1 0 0],'LineWidth',0.5,'LineStyle','-')
    %     line([0 8],[20*log10(MEMR_mem.thdAmp),20*log10(MEMR_mem.thdAmp)],'Color',[0 0 0],'LineWidth',0.5,'LineStyle',':')
    %     plot(MEMR_mem.thdOnsetTime,20*log10(MEMR_mem.thdAmp),'r.','MarkerSize',14)
    %     plot(MEMR_mem.thdOffsetTime,20*log10(MEMR_mem.thdAmp),'r.','MarkerSize',14)
    %     ylim([0,ymax])
    %     xlabel('Time (s)','FontSize',11)
    %     ylabel('Change (dB)','FontSize',11)
    %     keyboard
    % catch ME
    %     keyboard
    % end

    pause(0.01)

end
% INTERNAL FUNCTIONS ------------------------------------------------------
function [MEMR] = runme(Clicks,nClicks,indx,time,fs,stabilityCheck,subjectName)
    timeChunk = zeros(1,nClicks);
    for jj=1:nClicks % loop across each click position in the sweep
        timeChunk(1,jj) = time(indx(jj)); % the temporal postion of each click position
    end
    nSweeps = size(Clicks,2); % number of sweeps

    % Window all clicks and save into a matrix C --------------------------
        % round trip travel time between ear canal probe and eardrum:
        % depends on depth of insertion, but probably 1.5-2.5 cm
        % assume speed of sound is 34400 cm/s
        % round((1/((34400/1.5)))*fs)*2
        % 8-14 samples
    winN = 14; % number of samples to reduce (before reflection comes back)
    if stabilityCheck == 1 % this is for incident part of the click -----
        clickN = winN*2+1; % number of analysis samples in each (windowed click)
        h = hann(clickN); % make a hann window to window the edges 
        H = repmat(h,1,nSweeps);
        C = zeros(clickN,nSweeps,nClicks); % initialize matrix of windowed clicks (C for clicks)
        Cn = zeros(clickN,nSweeps,nClicks); % initialize noise matrix
        for ii=1:nClicks % loop across each click position (1-160)
            q = Clicks(indx(ii)-winN:indx(ii)+winN,:);
            q = q .* H;
            C(:,:,ii) = q;

            qn = Clicks(indx(ii)+1000-winN:indx(ii)+1000+winN,:); % window for noise estimate
            Cn(:,:,ii) = qn; % matrix of clicks
        end
    else % this is for the reflected (MEMR) part of the click -----------
        clickN = 100 + winN; % number of analysis samples in each (windowed click)
        h = hann(winN*2); % make a hann window to window the edges 
        h = h(1:winN);
        H = repmat(h,1,nSweeps);
        C = zeros(clickN,nSweeps,nClicks); % initialize matrix of windowed clicks (C for clicks)
        Cn = zeros(clickN+winN,nSweeps,nClicks);
        for ii=1:nClicks % loop across each click position (1-160)
            q = Clicks(indx(ii):indx(ii)+clickN-1,:); % reflected part of the sweep starting 14 samples after incident click
            qn = Clicks(indx(ii)-winN:indx(ii)+clickN-1,:); % *** include the incident part! ***
            q(1:winN,:) = q(1:winN,:) .* H;
            qn(1:winN,:) = qn(1:winN,:) .* H;
            q = flipud(q);
            qn = flipud(qn);
            q(1:winN,:) = q(1:winN,:) .* H;
            qn(1:winN,:) = qn(1:winN,:) .* H;
            q = flipud(q);
            qn = flipud(qn);

            C(:,:,ii) = q; % matrix of clicks
            Cn(:,:,ii) = qn; % matrix of clicks

            % The following is to calculate load impedance -----------
            %[pl,Pl,phi,other,wf] = ARLas_convertPL(mean(q,2),header.userInfo.iscS1L);
            %zl = other.ZL;
            % %S = fft(stimulus,nfft); % stimulus vector
            % nfft = header.userInfo.iscS1R.nfft;
            % R = fft(mean(q,2),nfft); % recordings vector
            % %PL = R ./ S; % load pressure
            % ff = (0:1:nfft-1)'*(fs/nfft);
            % [~,indx1] = min(abs(ff-header.userInfo.iscS1R.fmin));
            % [~,indx2] = min(abs(ff-header.userInfo.iscS1R.fmax));
            % PL = R(indx1:indx2);
            % %PL = PL(1:obj.fmaxIndx); % cut to highest calibrated frequency
            % %PL = PL(obj.fminIndx:end);
            % ZS = header.userInfo.iscS1R.ZS;
            % PS = header.userInfo.iscS1R.PS;
            % ZL(:,ii) = (ZS .* PL) ./ (PS - PL); % calculate load impedance
            %zl = pl.fpl;
            %ZL = repmat(q,1,size(Clicks,2));
            %q(1:winN,:) = q(1:winN,:) .* H;
            %q = flipud(q);
            %q(1:winN,:) = q(1:winN,:) .* H;
            %q = flipud(q);
            % % % q2(1:winN,:) = q2(1:winN,:) .* H;
            % % % q2 = flipud(q2);
            % % % q2(1:winN,:) = q2(1:winN,:) .* H;
            % % % q2 = flipud(q2);
            %ZL(:,ii) = zl;
        end
        % % % pr = PR(:,1);
        % % % pr = repmat(pr,1,size(PR,2));
        % % % PR = PR ./ pr;
        % % % plot(pl.f,20*log10(PR))

    end

    fmin = 100; % minimum frequency to analyze (Hz)
    fmax = 4000; % maximym frequency to analyze (Hz)


     % this is for load impedance -----------
     % freq = header.userInfo.iscS1R.freq;
     %    [~,indx1] = min(abs(fmin-freq)); % index of minimum frequency
     %    [~,indx2] = min(abs(fmax-freq)); % index of maximum frequency   
     %    ZL = ZL(indx1:indx2,:);
     %    freq = freq(indx1:indx2,:);
     %    zl = ZL(:,1);
     %    zl = repmat(zl,1,size(ZL,2));
     %    Q = ZL ./ zl;
     %    figure
     %    plot(freq,abs(Q))
     %    xlim([0 5000])
     %    pause(0.01)


    % Convert clicks into the frequency domain ----------------------------
    nfft = 960; % size of fft--chosen to give 100-Hz bin width
    frequency = (0:1:nfft-1)'*(fs/nfft); % frequency vector
    [~,indx1] = min(abs(fmin-frequency)); % index of minimum frequency
    [~,indx2] = min(abs(fmax-frequency)); % index of maximum frequency
    freq = frequency(indx1:indx2);
    nFreqs = length(freq); % number of frequencies to analyze
    originalN = size(squeeze(C(:,1,:)),1);
    scale = originalN / 2;
    
    CFT = zeros(nFreqs,nSweeps,nClicks); % initialize matrix of Fourier transformed clicks
    CFTn = zeros(nFreqs,nSweeps,nClicks); % initialize noise matrix
    TREND = zeros(nFreqs,nSweeps,nClicks); % initialize saved matrix of trend line
    TRENDn = zeros(nFreqs,nSweeps,nClicks);
    W = zeros(nFreqs,nSweeps,nClicks); % initialize weighting matrix (downweight bad samples)
    Wn = zeros(nFreqs,nSweeps,nClicks);
    
    warning off
    for ii=1:nSweeps % loop across 12 sweeps
        % for signal
        cc = squeeze(C(:,ii,:));

% figure(5)
% hold off
% plot(cc)
% ylabel('Pa')
% xlabel('Time (samples')


        FT = fft(cc,nfft); % contains all sweeps at all frequencies for a given click position
        FT = FT(indx1:indx2,:); % keep only the frequencies of interest
        FT  = FT ./ scale; % scale the magnitude appropriately
        CFT(:,ii,:) = complex(FT); % Fourier transform of C, cut to frequenc6ies of interest

        % for noise 
        ccn = squeeze(Cn(:,ii,:));
        FT = fft(ccn,nfft); % contains all sweeps at all frequencies for a given click position
        FT = FT(indx1:indx2,:); % keep only the frequencies of interest
        CFTn(:,ii,:) = complex(FT); % Fourier transform of C, cut to frequenc6ies of interest
    end

    if stabilityCheck ~=5 % this is hack to make this always happen

        %q = squeeze(CFT(5,:,:));
        %q = q.';
        %plot(abs(q(:)))
        %hold on
        
        % find long-term ipsilateral trends and extract them ---------
        [CFT,TREND,W] = memrDetrend(CFT);
        [CFTn,TRENDn,Wn] = memrDetrend(CFTn);
       for jj=1:nFreqs
            dummy = squeeze(TREND(jj,:,:));
            dummy = dummy.';
            dummy = dummy(:);
            dummy = dummy ./ dummy(1);
            dummy = abs(dummy - 1) + 1;
            %dummy = 20*log10(dummy);
            Trend(:,jj) = dummy;

            dummy = squeeze(TRENDn(jj,:,:));
            dummy = dummy.';
            dummy = dummy(:);
            dummy = dummy ./ dummy(1);
            dummy = abs(dummy - 1) + 1;
            %dummy = 20*log10(dummy);
            Trendn(:,jj) = dummy;
        end

        % Magnitude and phase analysis together --------------------------
        x = (0:1:nClicks-1);
        Z = zeros(nClicks,nFreqs);
        Z_sm = zeros(nClicks,nFreqs);
        % want to begin and end at the same place, so need to make a "knot"
        % condition at the ends
        xxx = (1:1:nClicks*3)-nClicks;
        for jj=1:nFreqs
            m = squeeze(CFT(jj,:,:)); % size m = 160 x 15
            w = squeeze(W(jj,:,:));
            Mr = real(m); % matrix of 15 sweeps at one click
            Mi = imag(m);
            sumw = sum(w,1); % take the weighted mean at one click (and one frequency)

% figure(6)
% subplot(2,1,1)
% hold off
% plot(Mr)
% ylabel('Real')
% subplot(2,1,2)
% hold off
% plot(Mi)
% ylabel('Imag')
% %pause

            maxIterations = 50;
            [signal,noise,weightsR1] = bisquareWeights(Mr(:,1).',maxIterations);
            [signal,noise,weightsR2] = bisquareWeights(Mr(:,2).',maxIterations);
            [signal,noise,weightsI1] = bisquareWeights(Mi(:,1).',maxIterations);
            [signal,noise,weightsI2] = bisquareWeights(Mi(:,2).',maxIterations);
            wr = [weightsR1',weightsR2'];
            wi = [weightsI1',weightsI2'];
            sumwr = sum(wr,1);
            sumwi = sum(wi,1);

            %mr = sum(Mr(:,:).*w,1)./ sumw;
            %mi = sum(Mi(:,:).*w,1)./ sumw;
            mr = sum(Mr(:,:).*wr,1)./ sumwr;
            mi = sum(Mi(:,:).*wi,1)./ sumwi;
            % smooth the response (effectively lowpass filter) -- this
            % worked when we had 160 clicks, but here, we have only 2, so
            % don't use it
            % mrrr = [mr,mr,mr];
            % miii = [mi,mi,mi];
            % www = [sumw,sumw,sumw];
            % sm = 0.0001; % was 0.00001;  smoothing factor (smaller numbers are more smooth)
            % ppr = csaps(xxx,mrrr,sm,[],www); % piecewise polynomial real coefficients
            % ppi = csaps(xxx,miii,sm,[],www); % piecewise polynomial imaginary coefficients
            % mr_sm = ppval(ppr,x); % evaluate only at the original x values
            % mi_sm = ppval(ppi,x);
            mr_sm = mr;
            mi_sm = mi;

            z = mr + 1i*mi;
            Z(:,jj) = z; % raw Z
            z_sm = mr_sm + 1i*mi_sm; % smoothed raw Z

            %figure
            %plot(z)
            %hold on
            %plot(z_sm,'r')
            %plot(z_sm(end),'r.')
            %plot(z_sm(1),'r*')


            % calculate length of curve -- this is problematic because it
            % doesn't return to where started. Calculating the tip is also
            % problematic. Probably don't use.
             t = x * 0.05; % time vector
                        %             dy = gradient(imag(z_sm)./gradient(t));
                        %             dx = gradient(real(z_sm)./gradient(t));
                        %             d = sqrt(dy.^2 + dx.^2);
                        %             [dmax,dmaxIndx] = max(d);
                        %             k = ones(size(d));
                        %             k(dmaxIndx+1:end) = -1;
                        %             D = cumsum(d.*k);

            baseline = z_sm(1);
            d = z_sm ./ baseline; % normalize z_sm
            d1 = abs(d-1)+1; % the combined mag+phase change
            d2 = abs(abs(d)-1)+1; % magnitude only change

            Z_sm(:,jj) = z_sm;
            D(:,jj) = d;
            D1(:,jj) = d1; % use this!
            D2(:,jj) = d2;
        end
        
    end

    % average over these frequencies:
    lowCutoffHz = 500;
    highCutoffHz = 1500;
    [~,lowIndx] = min(abs(freq -lowCutoffHz));
    [~,highIndx] = min(abs(freq -highCutoffHz));

    % rather than a straight average, zero out extreme values
    DD1 = D1(:,lowIndx:highIndx);
    peaks = max(DD1,[],1);
    doPlotAR = 0;
    multiplier = 1.5;
    [rejectIndx,nRejects] = newAR2(peaks,multiplier,doPlotAR); 
    ww1 = ones(size(peaks));
    ww1(rejectIndx) = 0;
    d = sum(D(:,lowIndx:highIndx).*ww1,2) / sum(ww1);
    d1 = sum(D1(:,lowIndx:highIndx).*ww1,2) / sum(ww1);
    d2 = sum(D2(:,lowIndx:highIndx).*ww1,2) / sum(ww1);
    if ~isempty(rejectIndx)
        %keyboard
    end

    trend =  sum(Trend(:,lowIndx:highIndx).*ww1,2) / sum(ww1);
    trendn =  sum(Trendn(:,lowIndx:highIndx).*ww1,2) / sum(ww1);
    timeTrend = (0:1:length(trend)-1)' * 0.05;

    % cut down to size (11 frequencies)
    D = D(:,lowIndx:highIndx);
    D1 = D1(:,lowIndx:highIndx);
    D2 = D2(:,lowIndx:highIndx);
    Z = Z(:,lowIndx:highIndx);
    Z_sm = Z_sm(:,lowIndx:highIndx);
    freq = freq(lowIndx:highIndx);


    MEMR.Trend = Trend;
    MEMR.Trendn = Trendn;
    MEMR.trend = trend;
    MEMR.trendn = trendn;
    MEMR.timeTrend = timeTrend;
    MEMR.D = D;
    MEMR.D1 = D1;
    MEMR.D2 = D2;
    MEMR.d = d;
    MEMR.d1 = d1;
    MEMR.d2 = d2;
    MEMR.x = x;
    MEMR.t = t;
    MEMR.freq = freq;
    MEMR.Z = Z;
    MEMR.Z_sm = Z_sm;
   % MEMR.elicitorLevel = elicitorLevel;
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
function [indx,nRejects] = newAR2(y,multiplier,doPlot)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Riff on AR.
% This version takes tolerance as any real number multiplier instead of a
% string (moderate, mild, etc.)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%doPlot = 0;

y = y(:); % ensure y is a column vector
%y = y(isfinite(y)); % get rid of NaN
warning('OFF')

xx = sort(y(:));
[dummy,indx] = max(xx);
xx = xx(1:indx); % get rid of any nan
N = length(xx); % number of observations
q = 100 *(0.5:N-0.5)./N;
xx = [min(xx); xx(:); max(xx)];
q = [0 q 100];
F1 = interp1(q,xx,25); % the first fourth, approx 25th precentile
F2 = interp1(q,xx,50); % the first half, approx 50th precentile
F3 = interp1(q,xx,75); % the third fourth, approx 75th percentile
IQR = F3 - F1; % the interquartile range

ArtifactVector = y >= (F1-multiplier*IQR) & y <= (F3 + multiplier*IQR);
[indx,val] = find(~ArtifactVector); % index the artifacts which should be rejected

nRejects = length(indx); % number of rejected buffers
percentRejected = (length(indx) / N) * 100;
% disp([num2str(percentRejected),' percent of the buffers rejected as artifacts.'])
% commented out previous line on 6/3/2010--ras and ssg

if doPlot == 1,
    figure(63)
    plot(q,xx,'*b-')
    hold on
    L1 = line([q(1) q(end)],[F1 F1]);
    L2 = line([q(1) q(end)],[F2 F2]);
    L3 = line([q(1) q(end)],[F3 F3]);
    L4 = line([q(1) q(end)],[F1-multiplier*IQR F1-multiplier*IQR]);
    L5 = line([q(1) q(end)],[F3 + multiplier*IQR F3 + multiplier*IQR]);
    set(L1,'Color',[0 1 0])
    set(L2,'Color',[0 0 0],'LineStyle',':')
    set(L3,'Color',[0 1 0])
    set(L4,'Color',[1 0 0])
    set(L5,'Color',[1 0 0])
    hold off
    %pause
end
warning('ON')
end
function [Clicks] = unjitter(Clicks,subjectName)
    Q = Clicks;
    Q = Q(1:2000,:); % isolate the first click (of two in each column)
    badIndx = []; % location of the "bad" (i.e., jittered) columns
    for ii=1:size(Q,2)
        [~,Indx(ii,1)] = max(Q(:,ii));
    end
    badIndx = find(diff(Indx)~=0);
    badIndx = badIndx + 1; % location of glitch in the matrix
    indxRef = Indx(1);
    offset = Indx(badIndx) - indxRef; % how many samples of glitch
    % return to baseline isn't wrong, and is coded as zero. Get ridDPOAE_L of these
    indx = find(offset~=0);
    if ~isempty(indx)
        offset = offset(indx);
        badIndx = badIndx(indx);
    end
    % if true peak is between two samples, noise jitters peak between the two; Get rid of these
    indx = find(abs(offset)>1);
    if ~isempty(indx)
        offset = offset(indx);
        badIndx = badIndx(indx);
    end
    % the following subjects are not actually problems. leave them alone
    if ~isempty(badIndx)
        disp('Found Jitter!')
        % line everything up with peaks at the mode
        winner = 482; %290;
        D2 = reshape(Clicks,size(Clicks,1)/2,size(Clicks,2)*2);
        for jj=1:size(D2,2)
            [~,maxIndx] = max(D2(1:4000,jj));
            shift = maxIndx - winner;
            if shift ~= 0
                dummy = D2(:,jj);
                if shift < 0
                    d2 = [zeros(abs(shift),1);dummy(1:end-abs(shift))];
                    D2(:,jj) = d2;
                elseif shift > 0
                    d2 = [dummy(abs(shift)+1:end);zeros(abs(shift),1)];
                    D2(:,jj) = d2;
                end
                [~,maxIndx] = max(D2(1:4000,jj));
                if maxIndx ~= winner
                    keyboard
                end
            end
        end
        Clicks = reshape(D2,size(D2,1)*2,size(D2,2)/2);
    end
end
function [CFT,TREND,W] = memrDetrend(CFT)
    nFreqs = size(CFT,1);
    for ii=1:nFreqs % do this individually for each frequency
        cc = squeeze(CFT(ii,:,:)); % pick the current frequency
        mr = real(cc); % the real part
        mr = mr'; % transpose
        [rows,cols] = size(mr); % use this later to reshape ack
        mr = mr(:); % force to a single column
        mi = imag(cc); % now do the same for the imaginary part
        mi = mi';
        mi = mi(:);
        x = (1:1:length(mr))'; % x-axis vector for spline fitting
        w = ones(size(mr)); % weighting vector (set to ones)
        
        doPlotAR = 0;
        multiplier = 1.5;
        [rejectIndx,nRejects] = newAR2(mr',multiplier,doPlotAR); 
        w(rejectIndx) = 0;
        
        smoothing = 0.01; % smoothing factor--0 is a straight line, 1 is cubic spline
        ppr = csaps(x,mr,smoothing,[],w); % piecewise polynomial real coefficients
        ppi = csaps(x,mi,smoothing,[],w); % piecewise polynomial imaginary coefficients
        mr_sm = ppval(ppr,x);
        mi_sm = ppval(ppi,x);
        
        %         figure
        %         plot(mr.*w,'b.')
        %         hold on
        %         plot(mr_sm,'r')
        %         keyboard

        trend = mr_sm+1i*mi_sm;
        offset = mean(trend); % mean location before detrending
        mrfixed = mr - mr_sm; % "fixed", i.e., detrended version of mr
        mifixed = mi - mi_sm;
        mrfixed = reshape(mrfixed,rows,cols);
        mifixed = reshape(mifixed,rows,cols);
        trend = reshape(trend,rows,cols);
        w = reshape(w,rows,cols);
        z = mrfixed + 1i*mifixed + offset; % add the offset back in to the complex form
        CFT(ii,:,:) = z.';
        TREND(ii,:,:) = trend.';
        W(ii,:,:) = w.';
    end
end

% OLD CODE ------------------------------    %
