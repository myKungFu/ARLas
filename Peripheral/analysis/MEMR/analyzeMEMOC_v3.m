function [MEMR_inc,MEMR_mem,h1,h2] = analyzeMEMOC_v3(header,Clicks,headerN,Noise,subjectName,runNumber)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [MEMR_inc,MEMR_mem,h1] = analyzeMEMOC_v3(header,Clicks,headerN,Noise,subjectName,runNumber)
%
% Originally was taken from analyzeMEMOC_v16
%
%
% Author: Shawn Goodman
% Date: August 19, 2024
% Last Updated: August 19, 2024 -- ssg
% Last Updated: August 30, 2024 -- ssg
% Last Updated: October 3, 2024 -- ssg
%   To fix:
%           ensure switching between ipsi, contra, bilateral
%           What should we hold constant?
%               elicitor range?  -- can't do this, because of MOC vs MEMR
%               rate of elicitor change? -- if do this, will make MOC sweep much shorter
%               length of sweep (seconds)? -- if do this, will have same
%                                           averaging time, but different sweep rates
% 
% 50 to 120 in 4 seconds = (120-50)/4 = 17.5 dB/s
% -20  to 50 in 4 seconds = (50--20)/4 = 17.5 dB/s
% 
% 45 to 115 in 4 seconds = (115-40)/4 = 17.5 dB/s
% -20  to 50 in 4 seconds = (50--20)/4 = 17.5 dB/s
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp(' ')
disp('Loading saved data....')
tic
%basePath = 'C:\myWork\ARLas\Data\MEMOC\ssgTrial_22AUG2024_moc_run1\';
basePath =  'C:\myWork\ARLas\Data\MEMOC\ssgTrial_22AUG2024_dummy_run2\';
cd(basePath)
dummy = load('Ch3_ER10xA_memoc_0001.mat');
header = dummy.header;
Clicks = dummy.data;
load('Ch4_ER10xB_memoc_0001.mat')
headerN = dummy.header;
Noise = dummy.data;
clear dummy
runNumber = 1;
subjectName = header.subjID;
cd( 'C:\myWork\ARLas\Peripheral\analysis\MEMOC')
disp('Done loading data....')
toc

    % ---------------------------------------------------------------------
    doJitterFix = 0;
    % ---------------------------------------------------------------------

    % extract the raw recordings -----
    micCorrection = header.userInfo.micCorrectionL;
    Clicks = applyMicCorrection(Clicks,micCorrection);
    micCorrectionN = headerN.userInfo.micCorrectionR;
    Noise = applyMicCorrection(Noise,micCorrectionN);
    %[pl,Pl,phi,other,wf] = ARLas_convertPL(Clicks(:,1),header.userInfo.iscS1L);

    % reshape to correct sizes
    chunkSize = header.userInfo.chunkSize; % 9600
    nChunks = header.userInfo.nChunks; % 80
    nSweeps = header.userInfo.nSweeps; % variable
    Clicks = reshape(Clicks,chunkSize*nChunks,nSweeps);
    Noise = reshape(Noise,chunkSize*nChunks,nSweeps);
    fs = header.fs; % sampling rate (96 kHz)
    time = (0:1:size(Clicks,1)-1)'/fs; % time vector
    indx = header.userInfo.clickIndx; % location of the clicks
    time = time / time(end); % express time from 0 to 1
    %elicitor = Noise(:,4); 
    % ----> what is special about this? Nothing, we just chose an exemplar for plotting
    h2 = figure;
    y = 20*log10(abs(Clicks(:,4))/.00002);
    plot(time,y,'k')
    ymax = max(y)+(0.05*max(y));
    ymin = 40;
    ylim([ymin,ymax])
    xlim([time(1),time(end)])
    xlabel('Time (s)','FontSize',12)
    ylabel('Magnitude (dB SPL)','FontSize',12)
    title(['Ipsilateral ',subjectName])

    % When considering which portions of each time window to analyze for MEMR:
    % Within each click window, we want to include several round trip
    % travel times of the click in the canal, but no low-frequency OAEs. We
    % should consider only the first 2 ms of the click, which is 
    % round(0.002*fs) = 192 samples. 
    % We used a noise activator with a rise time of 4 seconds and a fall
    % time of 4 seconds. This gives 8 seconds / .05 = 160 clicks.
    nClicks = length(indx); % number of total clicks to analyze

    % location of noise levels
    modLen = 4; % modulation duration in sec (from onset of noise to offset)
    [nLoops,clickTrain,noiseSamples,nReps,H] = getStimuli(fs,modLen);
    % Here, H is the (linear) amplitude vector for the noise. It is 8
    % seconds long total--4 seconds rise and 4 seconds fall.
    % The noise floor in the ear canal depends on subject but is usually roughly 40-60 dB SPL.
    % The elicitor noise goes from 
    % itself goes from 50 dB SPL to 120, covering a change of 70
    % dB in 4 seconds, so a rate of 17.5 dB/second
    h = H(:,1)*(10^(header.userInfo.targetNoise/20)*.00002); % h is the vector of noise levels (intended, in Pa)
    %plot(20*log10(abs(Noise(:,1))/.00002))
    %hold on
    %plot(20*log10(h/.00002),'r','LineWidth',2)
    noiseLvlTarget = 20*log10(h/.00002); % noise level (intended, dB SPL)
    % --> what are these numbers???    
    %noiseRate = (115 - 65) / 4; % noise changes at 12.5 dB per second
            % multiply this value by the timeShift value to see how many dB
            % to subtract from the noise vector when calculating thresholds

    % run the data analysis program
    %clicks = Clicks(1:2000,:);
    %Clicks = ARLas_hpFilter(Clicks,fs,100); % this doesn't add anything
    %here. Must have already been done.


    % get the stimulus levels ---------------------------------------------
    % note: elicitor here is the 4th (arbitrarily chosen) column of Noise.
    % This no longer makes sense here, given the new paradigm.
    % to get the click and noise levels, use median to avoid noise
    % contamination, then separate out the clicks from the noise portions.
    % for noise, should include the clicks and not for comparison.
    Q = median(Clicks,2);
    t = (0:1:length(Q)-1)'/fs ;
    % the peaks of the clicks are located very close to indx. 
    % go 192 samples on either side (total length 384 samples = 4 ms)
    offset = 192;
    ref = 0.00002;
    % get overall click sound pressures
    for jj=1:nClicks
        dummy = Q(indx(jj)-offset:indx(jj)+offset);
        noise_target(jj,1) = mean(noiseLvlTarget(indx(jj)-offset:indx(jj)+offset));
        click_pSPL(jj,1) = 20*log10(max(abs(dummy))/ref);
        click_peSPL(jj,1) = 20*log10(((max(dummy) + abs(min(dummy))) / 2 * sqrt(.5))/ref);
        click_rmsSPL(jj,1) = 20*log10((sqrt(mean(abs(dummy.^2))))/ref);
        time_levelEstimates(jj,1) = mean(t(indx(jj)-offset:indx(jj)+offset)); % mean time of each 50 ms chunk
    end
    % get overall noise sound pressures with clicks included
    offset2 = round(median(diff(indx)) / 2) - 2;
    for jj=1:nClicks
        if jj==1
            dummy = Q(1:indx(jj)+offset2);
        elseif jj == nClicks
            dummy = Q(indx(jj)-offset2:end);
        else
            dummy = Q(indx(jj)-offset2:indx(jj)+offset2);
        end
        noise_rmsSPL_withClick(jj,1) = 20*log10((sqrt(mean(abs(dummy.^2))))/ref);
    end
    % get overall noise sound pressures without clicks included
    QQ = Q;
    for jj=1:nClicks
        QQ(indx(jj)-offset:indx(jj)+offset) = 0;
        if jj==1
            dummy = QQ(1:indx(jj)+offset2);
        elseif jj == nClicks
            dummy = QQ(indx(jj)-offset2:end);
        else
            dummy = QQ(indx(jj)-offset2:indx(jj)+offset2);
        end
        noise_rmsSPL_withoutClick(jj,1) = 20*log10((sqrt(mean(abs(dummy.^2))))/ref);
    end
    figure
    plot(time_levelEstimates,click_peSPL,'b')
    hold on
    plot(time_levelEstimates,click_pSPL,'c')
    plot(time_levelEstimates,click_rmsSPL,'g')
    plot(time_levelEstimates,noise_rmsSPL_withoutClick,'r')
    plot(time_levelEstimates,noise_rmsSPL_withClick,'m')
    plot(time_levelEstimates,noise_target,'k')
    xlabel('Time (s)')
    ylabel('Level (dB SPL)')
    legend('Click peSPL','Click pSPL','Click rmsSPL','noise SPL with click','noiseSPL without click','noise target')
    %----------------------------------------------------------------------

    analyze MEMR --------------------------------------------------------
    stabilityCheck = 0;
    MEMR_mem = runme(Clicks,nClicks,indx,time,fs,stabilityCheck);
    % analyze ear canal stability
    stabilityCheck = 1; % inc stands for incident
    MEMR_inc = runme(Clicks,nClicks,indx,time,fs,stabilityCheck);

    % additionally save the stimulus levels
    MEMR_mem.time_levelEstimates = time_levelEstimates; % (average) time at which stimulus level estimates are made
    MEMR_mem.click_peSPL = click_peSPL; % peSPL of the clicks
    MEMR_mem.click_pSPL = click_pSPL; % pSPL of the clicks
    MEMR_mem.click_rmsSPL = click_rmsSPL; % rms SPL of the clicks
    MEMR_mem.noise_rmsSPL_withoutClick = noise_rmsSPL_withoutClick; % noise level (rms) without clicks (clicks set to zero)
    MEMR_mem.noise_rmsSPL_withClick = noise_rmsSPL_withClick; % noise level (rms) with clicks included
    MEMR_mem.noise_target = noise_target;
    MEMR_inc.time_levelEstimates = time_levelEstimates; % (average) time at which stimulus level estimates are made
    MEMR_inc.click_peSPL = click_peSPL; % peSPL of the clicks
    MEMR_inc.click_pSPL = click_pSPL; % pSPL of the clicks
    MEMR_inc.click_rmsSPL = click_rmsSPL; % rms SPL of the clicks
    MEMR_inc.noise_rmsSPL_withoutClick = noise_rmsSPL_withoutClick; % noise level (rms) without clicks (clicks set to zero)
    MEMR_inc.noise_rmsSPL_withClick = noise_rmsSPL_withClick; % noise level (rms) with clicks included
    MEMR_inc.noise_target = noise_target;

    %  plotting -----------------------------------------------------------
    try
        h1 = figure;
        sp1 = subplot(2,1,1);
        sp2 = subplot(2,1,2);
        sp1.Position = [0.1300    0.7230    0.7750    0.2020];
        sp2.Position = [0.1300    0.1100    0.7750    0.4924];
      axes(sp1)
        ax = gca;
        ax.FontSize = 11; 
        plot(MEMR_mem.timeTrend,20*log10(MEMR_mem.trend),'b')
        hold on
        plot(MEMR_mem.timeTrend,20*log10(MEMR_inc.trend),'b--')
        Q = repmat(20*log10(MEMR_mem.d1),1,nSweeps);
        Q = Q(:);
        plot(MEMR_mem.timeTrend,20*log10(MEMR_mem.trend)+Q,'k')
        xlabel('Time (s)','FontSize',11)
        ylabel('Change (dB)','FontSize',11)
        title([subjectName,'    Run# ',num2str(runNumber)])
        grid on
      axes(sp2)
        ax = gca;
        ax.FontSize = 11; 
        plot(MEMR_mem.t,20*log10(MEMR_mem.d1),'k')
        hold on
        plot(MEMR_mem.t,20*log10(MEMR_inc.d1),'k--')
        ymax = max(20*log10(MEMR_mem.d1))*1.1;
        line([4 4],[0 ymax],'Color',[0 0 0],'LineWidth',0.5,'LineStyle',':')
        plot(MEMR_mem.peakTime,20*log10(MEMR_mem.peakAmp),'r.','MarkerSize',14)
        line([MEMR_mem.peakTime,MEMR_mem.peakTime],[0 20*log10(MEMR_mem.peakAmp)],'Color',[1 0 0],'LineWidth',0.5,'LineStyle','-')
        line([MEMR_mem.thdOnsetTime,MEMR_mem.thdOnsetTime],[0 20*log10(MEMR_mem.thdAmp)],'Color',[1 0 0],'LineWidth',0.5,'LineStyle','-')
        line([MEMR_mem.thdOffsetTime,MEMR_mem.thdOffsetTime],[0 20*log10(MEMR_mem.thdAmp)],'Color',[1 0 0],'LineWidth',0.5,'LineStyle','-')
        line([0 8],[20*log10(MEMR_mem.thdAmp),20*log10(MEMR_mem.thdAmp)],'Color',[0 0 0],'LineWidth',0.5,'LineStyle',':')
        plot(MEMR_mem.thdOnsetTime,20*log10(MEMR_mem.thdAmp),'r.','MarkerSize',14)
        plot(MEMR_mem.thdOffsetTime,20*log10(MEMR_mem.thdAmp),'r.','MarkerSize',14)
        ylim([0,ymax])
        xlabel('Time (s)','FontSize',11)
        ylabel('Change (dB)','FontSize',11)
        %keyboard
    catch ME
        %keyboard
    end
    pause(0.01)

%     % analyze MOCR --------------------------------------------------------
%     stabilityCheck = 0;
%     MEMOC_moc = runmeMOC(Clicks,nClicks,indx,time,fs,stabilityCheck,header);
% 
%     MEMOC_moc.time_levelEstimates = time_levelEstimates; % (average) time at which stimulus level estimates are made
%     MEMOC_moc.click_peSPL = click_peSPL; % peSPL of the clicks
%     MEMOC_moc.click_pSPL = click_pSPL; % pSPL of the clicks
%     MEMOC_moc.click_rmsSPL = click_rmsSPL; % rms SPL of the clicks
%     MEMOC_moc.noise_rmsSPL_withoutClick = noise_rmsSPL_withoutClick; % noise level (rms) without clicks (clicks set to zero)
%     MEMOC_moc.noise_rmsSPL_withClick = noise_rmsSPL_withClick; % noise level (rms) with clicks included
%     MEMOC_moc.noise_target = noise_target;
% 
%     %  plotting -----------------------------------------------------------
%     try
%         h2 = figure;
%         ax = gca;
%         ax.FontSize = 11; 
%         plot(MEMOC_moc.t,20*log10(MEMOC_moc.d1),'k')
%         hold on
%         ymax = max(20*log10(MEMOC_moc.d1))*1.1;
%         line([4 4],[0 ymax],'Color',[0 0 0],'LineWidth',0.5,'LineStyle',':')
%         plot(MEMOC_moc.peakTime,20*log10(MEMOC_moc.peakAmp),'r.','MarkerSize',14)
%         line([MEMOC_moc.peakTime,MEMOC_moc.peakTime],[0 20*log10(MEMOC_moc.peakAmp)],'Color',[1 0 0],'LineWidth',0.5,'LineStyle','-')
%         line([MEMOC_moc.thdOnsetTime,MEMOC_moc.thdOnsetTime],[0 20*log10(MEMOC_moc.thdAmp)],'Color',[1 0 0],'LineWidth',0.5,'LineStyle','-')
%         line([MEMOC_moc.thdOffsetTime,MEMOC_moc.thdOffsetTime],[0 20*log10(MEMOC_moc.thdAmp)],'Color',[1 0 0],'LineWidth',0.5,'LineStyle','-')
%         line([0 8],[20*log10(MEMOC_moc.thdAmp),20*log10(MEMOC_moc.thdAmp)],'Color',[0 0 0],'LineWidth',0.5,'LineStyle',':')
%         plot(MEMOC_moc.thdOnsetTime,20*log10(MEMOC_moc.thdAmp),'r.','MarkerSize',14)
%         plot(MEMOC_moc.thdOffsetTime,20*log10(MEMOC_moc.thdAmp),'r.','MarkerSize',14)
%         ylim([0,ymax])
%         xlabel('Time (s)','FontSize',11)
%         ylabel('Change (dB)','FontSize',11)
%         title('MOCR')
%         %keyboard
%     catch ME
%         keyboard
%     end
%     pause(0.01)
% keyboard
end


% INTERNAL FUNCTIONS ------------------------------------------------------
function [MEMR] = runme(Clicks,nClicks,indx,time,fs,stabilityCheck)
    timeChunk = zeros(1,nClicks);
    for jj=1:nClicks % loop across each click position in the sweep
        timeChunk(1,jj) = time(indx(jj)); % the temporal postion of each click position
    end
    timeChunk = timeChunk * 8; % this gives the correct 50 ms interval between clicks
                               % time went from 0 to 1, but actually, the length is 8 seconds
    nSweeps = size(Clicks,2); % number of sweeps

    % Window all clicks and save into a matrix C --------------------------
        % round trip travel time between ear canal probe and eardrum:
        % depends on depth of insertion, but probably 1.5-2.5 cm
        % assume speed of sound is 34400 cm/s
        % round((1/((34400/1.5)))*fs)*2
        % 8-12 samples
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
        end

    end

    fmin = 100; % minimum frequency to analyze (Hz)
    fmax = 4000; % maximym frequency to analyze (Hz)

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
    for ii=1:nSweeps % loop across 15 sweeps
        % for signal
        cc = squeeze(C(:,ii,:));

        FT = fft(cc,nfft); % contains all sweeps at all frequencies for a given click position
        FT = FT(indx1:indx2,:); % keep only the frequencies of interest
        FT = FT ./ scale; % scale the magnitude appropriately
        CFT(:,ii,:) = complex(FT); % Fourier transform of C, cut to frequencies of interest

        % for noise 
        ccn = squeeze(Cn(:,ii,:));
        FT = fft(ccn,nfft); % contains all sweeps at all frequencies for a given click position
        FT = FT(indx1:indx2,:); % keep only the frequencies of interest
        CFTn(:,ii,:) = complex(FT); % Fourier transform of C, cut to frequenc6ies of interest
    end

    if stabilityCheck ~=5 % this is hack to make this always happen
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
            mr = sum(Mr(:,:).*w,1)./ sumw;
            mi = sum(Mi(:,:).*w,1)./ sumw;
            % smooth the response (effectively lowpass filter)
            mrrr = [mr,mr,mr];
            miii = [mi,mi,mi];
            www = [sumw,sumw,sumw];
            sm = 0.001; % was 0.0001; % was 0.00001;  smoothing factor (smaller numbers are more smooth)
            ppr = csaps(xxx,mrrr,sm,[],www); % piecewise polynomial real coefficients
            ppi = csaps(xxx,miii,sm,[],www); % piecewise polynomial imaginary coefficients
            mr_sm = ppval(ppr,x); % evaluate only at the original x values
            mi_sm = ppval(ppi,x);

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

            baseline = mean([z_sm(1:3),z_sm(end-2:end)]);
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


    % EXTRACT THE NEEDED METRICS ------------------------------------------
    if stabilityCheck ~= 1
        % peak delay -------------------
        sm = 1; % smoothing factor (smaller numbers are more smooth)
        n = 160; % number of clicks
        x = linspace(0,8,n)'; % x-axis for smoothing (click number)
        w = ones(size(x)); % weighting factor
        pp = csaps(x,d1,sm,[],w); % piecewise polynomial object
        dfdx = fnder(pp); % take derivative and solve for zero slope
        peakXX = fnzeros(dfdx); % peak location in seconds
        peakXX = peakXX(1,:);
        
        try
            
        peakYY = ppval(pp,peakXX);
        [~,peakIndx] = max(peakYY); 
        peakX = peakXX(peakIndx);
        peakY = peakYY(peakIndx);
      
        % if isempty(peakN)
        %     peakN = 1;
        % end
        % peakX = peakX(peakN); % function returns two outputs (both the same) take the first
        % peakY = ppval(pp,peakNf); % max average activation
        delay = peakX-4; % reflex delay in seconds. Peak of the noise occurred at t=4 seconds
        catch ME
            keyboard
        end
        %dY = ppval(dfdx,x);
        
        % thresholds ------------------
        %   use the Q concept: reduce 12 from peak
        d1_scaled = (d1-1)./(peakY-1);
        thd = 10.^(-12/20); % threshold is "Q12", or 12 dB down from peak
        d1_scaled = d1_scaled - thd;
        sm = 1; % smoothing factor (smaller numbers are more smooth)
        pp = csaps(x,d1_scaled,sm,[],w); % piecewise polynomial object    
        z2 = fnzeros(pp);
        Thd = z2(1,:); % threshold times
        try
            ons = Thd(find(Thd<peakX));
            thdOnsetTime = max(ons);
            offs = Thd(find(Thd>peakX));
            thdOffsetTime = min(offs);
            %thdOnsetTime = Thd(1); % onset threshold
            %thdOffsetTime = Thd(2); % offset threshold
        catch ME
            thdOnsetTime = NaN;
            thdOffsetTime = NaN;
        end
        % convert threshold times to thresholds re: nominal elicitor level
        % noise went from 45 to 115 in 4 seconds (17.5 dB/s)
        % y = mx + b
        % stimLevel = 17.5x + 45
        % but also account for reflex delay
        try
            % issue is whether or not to include delay. Theoretically it
            % seems right, but practically it causes problems, including
            % extra noise and issues with occasional negative delays.
            elicitorLevel = 17.5*x + 40; % convert time in seconds to elicitor level in dB SPL rms
            elicitorLevel = [elicitorLevel(1:80);flipud(elicitorLevel(1:80))];
            ThdLvl(1) = 17.5*thdOnsetTime + 40;
            ThdLvl(2) = -17.5*(thdOffsetTime-4) + 110;
            thdOnsetLvl = ThdLvl(1); % onset threshold re: stim level
            thdOffsetLvl = ThdLvl(2); % offset threshold re: stim level

            % if thdOffsetLvl > 95
            %     keyboard
            % end

            %ThdLvl(1) = 17.5*(Thd(1)-delay) + 45;
            %thdOnsetLvl = ThdLvl(1); % onset threshold re: stim level
            %ThdLvl(2) = -17.5*(Thd(2)-4-delay) + 115;
            %thdOffsetLvl = ThdLvl(2); % offset threshold re: stim level
            %elicitorLevel = 17.5*(x-delay) + 45;
        catch ME
            ThdLvl(1) = NaN;
            thdOnsetLvl = NaN; % onset threshold re: stim level
            ThdLvl(2) = NaN;
            thdOffsetLvl = NaN; % offset threshold re: stim level
            %elicitorLevel = 17.5*(x-delay) + 45;
            elicitorLevel = 17.5*x + 40; % convert time in seconds to elicitor level in dB SPL rms
            elicitorLevel = [elicitorLevel(1:80);flipud(elicitorLevel(1:80))];
        end
        % go back and get threshold amplitudes from non-scaled d1
        sm = 1; % smoothing factor (smaller numbers are more smooth)
        n = 160; % number of clicks
        x = linspace(0,8,n)'; % x-axis for smoothing (click number)
        w = ones(size(x)); % weighting factor
        pp = csaps(x,d1,sm,[],w); % piecewise polynomial object
        thdAmp = ppval(pp,thdOnsetTime);
        
        % hysteresis ---------------------
        try
            pp = csaps(x,d1,sm,[],w); % piecewise polynomial object
            hh = fnint(pp); % integrate the smoothed spline
            A = ppval(hh,thdOnsetTime);
            B = ppval(hh,peakX);
            C = ppval(hh,thdOffsetTime);
            aucLeft = B-A; % area under the curve left
            aucRight = C-B; % area under the curve right
            hyst = aucRight / aucLeft; % hysteresis as a ratio of area under the curves
            if hyst < 0
                keyboard
            end
            hysteresis = hyst;
        catch
            hysteresis = NaN;
        end
        % Calculate the slopes ------------------
        try
            peak_index = round(peakX /.05);
            index_xa = round(thdOnsetTime / .05);
            index_xb = round(thdOffsetTime / .05);
            part1_x = x(index_xa:peak_index);
            part1_y = d1(index_xa:peak_index);
            part2_x = x(peak_index+1:index_xb);
            part2_y = d1(peak_index+1:index_xb);
            slope_ascending = polyfit(part1_x,part1_y,1);
            slope_descending = polyfit(part2_x,part2_y,1);
            slopeUp = slope_ascending(1);
            slopeDn = slope_descending(1);
        catch
            slopeUp = NaN;
            slopeDn = NaN;
        end
    end
%---------------------------------------------------------------------------

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
    if stabilityCheck ~=1
        MEMR.peakTime = peakX;
        MEMR.peakAmp = peakY;
        MEMR.delay = delay;
        MEMR.thdOnsetTime = thdOnsetTime;
        MEMR.thdOffsetTime = thdOffsetTime;
        MEMR.thdOnsetLvl = thdOnsetLvl;
        MEMR.thdOffsetLvl = thdOffsetLvl;
        MEMR.hysteresis = hysteresis;
        MEMR.slopeUp = slopeUp;
        MEMR.slopeDn = slopeDn;
        MEMR.thd = thd;
        MEMR.thdAmp = thdAmp;
        MEMR.elicitorLevel = elicitorLevel;
    end

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
     % dBatten = 23; % dB attenuation re: full out
     % multiplier = 10^(-dBatten/20);
     % click = click * multiplier;
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
    %dBattenN = 1; % dB attenuation re: full out
    %multiplierNoise = 10^(-dBattenN/20);
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
        winner = 290;
        D2 = reshape(Clicks,size(Clicks,1)/2,size(Clicks,2)*2);
        for jj=1:size(D2,2)
            [~,maxIndx] = max(D2(:,jj));
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
                [~,maxIndx] = max(D2(:,jj));
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

        w(rejectIndx) = 0; % downweight noisy samples to zero
        
        % was 0.00000000001;
        smoothing = 0.000000001; % smoothing factor--0 is a straight line, 1 is cubic spline
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
function b = bpf()
% FIR Window Bandpass filter designed using the FIR1 function.
    Fs = 96000;  % Sampling Frequency
    Fstop1 = 700;              % First Stopband Frequency
    Fpass1 = 1000;             % First Passband Frequency
    Fpass2 = 3500;            % Second Passband Frequency
    Fstop2 = 4500;            % Second Stopband Frequency
    Dstop1 = 0.001;           % First Stopband Attenuation
    Dpass  = 0.057501127785;  % Passband Ripple
    Dstop2 = 0.001;           % Second Stopband Attenuation
    flag   = 'scale';         % Sampling Flag
    % Calculate the order from the parameters using KAISERORD.
    [N,Wn,BETA,TYPE] = kaiserord([Fstop1 Fpass1 Fpass2 Fstop2]/(Fs/2), [0 ...
                                 1 0], [Dstop1 Dpass Dstop2]);
    if mod(N,2)~=0
        N = N + 1;
    end
    % Calculate the coefficients using the FIR1 function.
    b  = fir1(N, Wn, TYPE, kaiser(N+1, BETA), flag);
    b = b(:);
end
function [MEMR] = runmeMOC(Clicks,nClicks,indx,time,fs,stabilityCheck,header)
    %clickIndx = header.userInfo.clickIndx;
    holeN1 = header.userInfo.holeN1; %288 beginning of the noise notch
    holeN2 = header.userInfo.holeN2; % 1152 ending of the noise notch

    timeChunk = zeros(1,nClicks);
    for jj=1:nClicks % loop across each click position in the sweep
        timeChunk(1,jj) = time(indx(jj)); % the temporal postion of each click position
    end
    timeChunk = timeChunk * 8; % this gives the correct 50 ms interval between clicks
                               % time went from 0 to 1, but actually, the length is 8 seconds
    nSweeps = size(Clicks,2); % number of sweeps

    % Window all clicks and save into a matrix C --------------------------
        % round trip travel time between ear canal probe and eardrum:
        % depends on depth of insertion, but probably 1.5-2.5 cm
        % assume speed of sound is 34400 cm/s
        % round((1/((34400/1.5)))*fs)*2
        % 8-12 samples
    if stabilityCheck == 1 % this is for incident part of the click -----
        % an option for memr but not for moc
        MEMR = [];
        return
    else % this is for the oae (MOCR) part of the click -----------
        %clickN = 100 + winN; % number of analysis samples in each (windowed click)
        winN = round(fs*.003); % 288 samples.
        h = hann(winN*2); % make a hann window to window the edges 
        h = h(1:winN);
        H = repmat(h,1,nSweeps);
        clickN = holeN2 - holeN1;
        C = zeros(clickN,nSweeps,nClicks); % initialize matrix of windowed clicks (C for clicks)
        Cn = C;
        Accumulator2 = zeros(nSamples,nClicks);
        Accumulator2n = zeros(size(qn,1),nClicks);
        for ii=1:nClicks % loop across each click position (1-160)
            disp(['Analyzing click ',num2str(ii),' of ',num2str(nClicks)])
            %q = Clicks(indx(ii):indx(ii)+clickN-1,:); % reflected part of the sweep starting 14 samples after incident click
            
            %q0 = Clicks(indx(ii)-holeN1+1:indx(ii)+holeN2,:); % this is the click waveform in the entire notch
            %q = Clicks(indx(ii)+holeN1+1:indx(ii)+holeN2,:);  % this is the waveform, exlcuding the stimulus click
                % q runs from 3 ms after the peak of the click to 12 ms
                % after the peak. Note: might be able to push this 1 ms
                % earlier, so holeN1-96?

            QQ = Clicks;
            [rows,cols] = size(QQ);
            b = bpf;
            QQ = fastFilter(b,QQ(:));
            QQ = reshape(QQ,rows,cols);
            %q = QQ(indx(ii)+holeN1-96:indx(ii)+holeN2,:);
            % was only 12 ms post peak of click. 1 kHz emission could be
            % longer (up to 17.1 ms), so add an additional 8 ms (768
            % samples), so including 20 ms post-click
            % NO!!! Must do only 11 ms, otherwise run into the noise
            q = QQ(indx(ii)+holeN1-96:indx(ii)+holeN2-96,:);
            qn = QQ(indx(ii)-96+1:indx(ii)+96,:);
            h = hann(96*2);
            h = h(1:96);
            H = repmat(h,1,nSweeps);
            q(1:96,:) = q(1:96,:) .* H;
            q = flipud(q);
            q(1:96,:) = q(1:96,:) .* H;
            q = flipud(q);
            h = hann(size(qn,1));
            H = repmat(h,1,nSweeps);
            qn = qn .* H;
            
            N = round(fs*0.025);
            t = (0:1:N-1)'/fs;
            t = t(1:length(q));

            fmin = 1000;
            fmax = 4000;
            CF = round(2.^(log2(fmin):1/6:log2(fmax)))';
            nFilters = length(CF);
            nSamples = size(q,1);
            Accumulator1 = zeros(nSamples,nFilters);
            for kk=1:nFilters
                disp(['   Analyzing filter ',num2str(ii),' of ',num2str(nFilters)])
                cf = CF(kk);
                [G,cf,tauHat,filtSamp] = GTFB(cf,cf,fs);
                clear Mag Wf
                for kk2=1:nSweeps
                    [Mag(:,kk2),frequency,time,Wf(:,kk2)] = applyGTfilterBank(q(:,kk2),G,cf,tauHat,filtSamp,fs);
                    [Magn(:,kk2),~,~,Wfn(:,kk2)] = applyGTfilterBank(qn(:,kk2),G,cf,tauHat,filtSamp,fs);
                end
                % specify time regions to include (from SS_analyze_getCeoaeCuts_v9.m)
                aa =       418.4;
                bb =      -0.4574;
                cc =      -4.247;
                tau = aa.*cf.^bb+cc;
                aa =        109.2;
                bb =    -0.3291;
                cc =     -3.854;
                tau05 = aa.*cf.^bb+cc;
                aa =      -1.22;
                bb =      0.3751;
                cc =       33.39;
                tau95 = aa.*cf.^bb+cc;
                [~,cut1] = min(abs(t-(tau05/1000)));
                [~,cut2] = min(abs(t-(tau95/1000)));
                winN = cut2-cut1 + 1;
                h = blackman(winN);
                H = repmat(h,1,nSweeps);
                Wf(cut1:cut2,:) = Wf(cut1:cut2,:) .* H;
                Wf(1:cut1,:) = 0;
                Wf(cut2:end,:) = 0;
                Accumulator1(:,kk) = mean(Wf,2);
                Accumulator2(:,ii) = Accumulator2(:,ii) + mean(Wf,2);
                Accumulator2n(:,ii) = Accumulator2n(:,ii) + mean(Wfn,2);
            
                nfft = size(Wf,1);
                originalN = winN;
                ref = 0.00002;
                [frequency,signal,noiseFloor,phase] = ARLas_fda(Wf,fs,0.00002,nfft,originalN);
                [~,peakIndx] = min(abs(frequency - cf));
                mag = signal(peakIndx);
                nf = noiseFloor(peakIndx); 
                phi = phase(peakIndx); 
                snr(kk,ii) = mag - nf;
                mag = 10^(mag/20)*ref;
                z(kk,ii) = mag*exp(1i*phi); % save results as complex coefficients for each filter and each click location

                [frequency,signal,noiseFloor,phase] = ARLas_fda(Wfn,fs,0.00002,nfft,originalN);
                %[~,peakIndx] = min(abs(frequency - cf));
                mag = signal(peakIndx);
                nf = noiseFloor(peakIndx); 
                phi = phase(peakIndx); 
                snrn(kk,ii) = mag - nf;
                mag = 10^(mag/20)*ref;
                zn(kk,ii) = mag*exp(1i*phi); % save results as complex coefficients for each filter and each click location
                
            end

%[frequency,signal,noiseFloor] = ARLas_fda(Wf,fs,0.00002);
            % 
            % 
            %     % winN is 288 samples, or 3 ms. This is quite a long for
            %     % the onset??
            %     % given the expected emission delay, and that you want to
            %     % include emissions from 1-3 kHz, you should  
            % 
            % qn = q; % the "noise" version ??? what is this here, just a place holder?
            % q(1:winN,:) = q(1:winN,:) .* H;
            % qn(1:winN,:) = qn(1:winN,:) .* H;
            % q = flipud(q);
            % qn = flipud(qn);
            % q(1:winN,:) = q(1:winN,:) .* H;
            % qn(1:winN,:) = qn(1:winN,:) .* H;
            % q = flipud(q);
            % qn = flipud(qn);
            % 
            % b = bpf;
            % q = fastFilter(b,q);
            % 
            % C(:,:,ii) = q; % matrix of clicks
            % Cn(:,:,ii) = qn; % matrix of clicks
        end

    end

keyboard    
% -----------------------------------------

Base = Accumulator2(:,1);
Base = repmat(Base,1,size(Accumulator2,2));
AC2 = Accumulator2 - Base;
RMS2 = sqrt(mean(AC2.^2,1));
RMS2(1) = RMS2(2);
base2 = mean(RMS2(1:3));
RMS21 = RMS2 ./ base2;
c1 = cumsum(RMS21);
w = ones(size(c1));
sm = 0.01; % was 0.0001; % was 0.00001;  smoothing factor (smaller numbers are more smooth)
t = (1:1:size(z,2))';
pp = csaps(t,c1,sm,[],w); % piecewise polynomial real coefficients
%c1sm = ppval(pp,t); % evaluate only at the original x values
pd = fnder(pp);
qsm = ppval(pd,t); % evaluate only at the original x values
baseline = mean([qsm(1:10);qsm(end-10:end)]);
qsm21 = qsm ./ baseline;


Base = Accumulator2n(:,1);
Base = repmat(Base,1,size(Accumulator2n,2));
AC2 = Accumulator2n - Base;
RMS2 = sqrt(mean(AC2.^2,1));
RMS2(1) = RMS2(2);
base2 = mean(RMS2(1:3));
RMS2 = RMS2 ./ base2;
c1 = cumsum(RMS2);
w = ones(size(c1));
sm = 0.05; % was 0.0001; % was 0.00001;  smoothing factor (smaller numbers are more smooth)
t = (1:1:size(z,2))';
pp = csaps(t,c1,sm,[],w); % piecewise polynomial real coefficients
%c1sm = ppval(pp,t); % evaluate only at the original x values
pd = fnder(pp);
qsm = ppval(pd,t); % evaluate only at the original x values
baseline = mean([qsm(1:10);qsm(end-10:end)]);
qsm = qsm ./ baseline;

plot(20*log10(RMS21),'m')
hold on
plot(20*log10(RMS2),'c')
plot(20*log10(qsm21),'r')
plot(20*log10(qsm),'b')

% ----------------------------------------
for ii=1:size(z,1)
    %m1 = abs(z(ii,:))';

    m1 = (zn(ii,:))';
    m1 = m1 ./ m1(1);
    m1 = abs(m1-1);


    t = (1:1:size(z,2))';
    c1 = cumsum(m1);
    w1 = snr(ii,:);
    w1 = 10.^(w1/10);
    
    sm = 0.001; % was 0.0001; % was 0.00001;  smoothing factor (smaller numbers are more smooth)
    pp = csaps(t,c1,sm,[],w1); % piecewise polynomial real coefficients
    %c1sm = ppval(pp,t); % evaluate only at the original x values
    pd = fnder(pp);
    m1sm = ppval(pd,t); % evaluate only at the original x values
    baseline = mean([m1sm(1:10);m1sm(end-10:end)]);
    m1sm = m1sm ./ baseline;

    figure(10)
    hold on
    plot(20*log10(m1sm))
    pause(0.01)

end
ABSZ = abs(z);
W = 10.^(snr./10);
Q = sum(ABSZ.*W,1) ./ sum(W,1);
c1 = cumsum(Q);
w = mean(W,1);
pp = csaps(t,c1,sm,[],w); % piecewise polynomial real coefficients
%c1sm = ppval(pp,t); % evaluate only at the original x values
pd = fnder(pp);
qsm = ppval(pd,t); % evaluate only at the original x values
baseline = mean([qsm(1:10);qsm(end-10:end)]);
qsm = qsm ./ baseline;
plot(20*log10(qsm),'r','LineWidth',2)




for ii=1:size(zn,1)
    m1 = (zn(ii,:))';
    m1 = m1 ./ m1(1);
    m1 = abs(m1-1);
    t = (1:1:size(zn,2))';
    c1 = cumsum(m1);
    w1 = snrn(ii,:);
    w1 = 10.^(w1/10);
    
    sm = 0.001; % was 0.0001; % was 0.00001;  smoothing factor (smaller numbers are more smooth)
    pp = csaps(t,c1,sm,[],w1); % piecewise polynomial real coefficients
    c1sm = ppval(pp,t); % evaluate only at the original x values
    pd = fnder(pp);
    m1sm = ppval(pd,t); % evaluate only at the original x values
    baseline = mean([m1sm(1:10);m1sm(end-10:end)]);
    m1sm = m1sm ./ baseline;

    figure(11)
    hold on
    plot(20*log10(m1sm))
    pause(0.01)
end

ABSZ = abs(zn);
W = 10.^(snrn./10);
Q = sum(ABSZ.*W,1) ./ sum(W,1);
c1 = cumsum(Q);
w = mean(W,1);
pp = csaps(t,c1,sm,[],w); % piecewise polynomial real coefficients
%c1sm = ppval(pp,t); % evaluate only at the original x values
pd = fnder(pp);
qsm = ppval(pd,t); % evaluate only at the original x values
baseline = mean([qsm(1:10);qsm(end-10:end)]);
qsm = qsm ./ baseline;
plot(20*log10(qsm),'r','LineWidth',2)



% -----------------------------------------
% ----------------------------------------
    CEOAE = C;
    fmin = 1000; % minimum frequency to analyze (Hz)
    fmax = 4000; % maximym frequency to analyze (Hz)

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
    for ii=1:nSweeps % loop across 15 sweeps
        % for signal
        cc = squeeze(C(:,ii,:));

        FT = fft(cc,nfft); % contains all sweeps at all frequencies for a given click position
        FT = FT(indx1:indx2,:); % keep only the frequencies of interest
        FT = FT ./ scale; % scale the magnitude appropriately
        CFT(:,ii,:) = complex(FT); % Fourier transform of C, cut to frequencies of interest

        % for noise 
        ccn = squeeze(Cn(:,ii,:));
        FT = fft(ccn,nfft); % contains all sweeps at all frequencies for a given click position
        FT = FT(indx1:indx2,:); % keep only the frequencies of interest
        CFTn(:,ii,:) = complex(FT); % Fourier transform of C, cut to frequenc6ies of interest
    end

    if stabilityCheck ~=5 % this is hack to make this always happen
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
            mr = sum(Mr(:,:).*w,1)./ sumw;
            mi = sum(Mi(:,:).*w,1)./ sumw;
            % smooth the response (effectively lowpass filter)
            mrrr = [mr,mr,mr];
            miii = [mi,mi,mi];
            www = [sumw,sumw,sumw];
            sm = 0.001; % was 0.0001; % was 0.00001;  smoothing factor (smaller numbers are more smooth)
            ppr = csaps(xxx,mrrr,sm,[],www); % piecewise polynomial real coefficients
            ppi = csaps(xxx,miii,sm,[],www); % piecewise polynomial imaginary coefficients
            mr_sm = ppval(ppr,x); % evaluate only at the original x values
            mi_sm = ppval(ppi,x);

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

            baseline = mean([z_sm(1:3),z_sm(end-2:end)]);
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


    % EXTRACT THE NEEDED METRICS ------------------------------------------
    if stabilityCheck ~= 1
        % peak delay -------------------
        sm = 1; % smoothing factor (smaller numbers are more smooth)
        n = 160; % number of clicks
        x = linspace(0,8,n)'; % x-axis for smoothing (click number)
        w = ones(size(x)); % weighting factor
        pp = csaps(x,d1,sm,[],w); % piecewise polynomial object
        dfdx = fnder(pp); % take derivative and solve for zero slope
        peakXX = fnzeros(dfdx); % peak location in seconds
        peakXX = peakXX(1,:);
        
        try
            
        peakYY = ppval(pp,peakXX);
        [~,peakIndx] = max(peakYY); 
        peakX = peakXX(peakIndx);
        peakY = peakYY(peakIndx);
      
        % if isempty(peakN)
        %     peakN = 1;
        % end
        % peakX = peakX(peakN); % function returns two outputs (both the same) take the first
        % peakY = ppval(pp,peakNf); % max average activation
        delay = peakX-4; % reflex delay in seconds. Peak of the noise occurred at t=4 seconds
        catch ME
            keyboard
        end
        %dY = ppval(dfdx,x);
        
        % thresholds ------------------
        %   use the Q concept: reduce 12 from peak
        d1_scaled = (d1-1)./(peakY-1);
        thd = 10.^(-12/20); % threshold is "Q12", or 12 dB down from peak
        d1_scaled = d1_scaled - thd;
        sm = 1; % smoothing factor (smaller numbers are more smooth)
        pp = csaps(x,d1_scaled,sm,[],w); % piecewise polynomial object    
        z2 = fnzeros(pp);
        Thd = z2(1,:); % threshold times
        try
            ons = Thd(find(Thd<peakX));
            thdOnsetTime = max(ons);
            offs = Thd(find(Thd>peakX));
            thdOffsetTime = min(offs);
            %thdOnsetTime = Thd(1); % onset threshold
            %thdOffsetTime = Thd(2); % offset threshold
        catch ME
            thdOnsetTime = NaN;
            thdOffsetTime = NaN;
        end
        % convert threshold times to thresholds re: nominal elicitor level
        % noise went from 45 to 115 in 4 seconds (17.5 dB/s)
        % y = mx + b
        % stimLevel = 17.5x + 45
        % but also account for reflex delay
        try
            % issue is whether or not to include delay. Theoretically it
            % seems right, but practically it causes problems, including
            % extra noise and issues with occasional negative delays.
            elicitorLevel = 17.5*x + 40; % convert time in seconds to elicitor level in dB SPL rms
            elicitorLevel = [elicitorLevel(1:80);flipud(elicitorLevel(1:80))];
            ThdLvl(1) = 17.5*thdOnsetTime + 40;
            ThdLvl(2) = -17.5*(thdOffsetTime-4) + 110;
            thdOnsetLvl = ThdLvl(1); % onset threshold re: stim level
            thdOffsetLvl = ThdLvl(2); % offset threshold re: stim level

            if thdOffsetLvl > 95
                keyboard
            end

            %ThdLvl(1) = 17.5*(Thd(1)-delay) + 45;
            %thdOnsetLvl = ThdLvl(1); % onset threshold re: stim level
            %ThdLvl(2) = -17.5*(Thd(2)-4-delay) + 115;
            %thdOffsetLvl = ThdLvl(2); % offset threshold re: stim level
            %elicitorLevel = 17.5*(x-delay) + 45;
        catch ME
            ThdLvl(1) = NaN;
            thdOnsetLvl = NaN; % onset threshold re: stim level
            ThdLvl(2) = NaN;
            thdOffsetLvl = NaN; % offset threshold re: stim level
            %elicitorLevel = 17.5*(x-delay) + 45;
            elicitorLevel = 17.5*x + 40; % convert time in seconds to elicitor level in dB SPL rms
            elicitorLevel = [elicitorLevel(1:80);flipud(elicitorLevel(1:80))];
        end
        % go back and get threshold amplitudes from non-scaled d1
        sm = 1; % smoothing factor (smaller numbers are more smooth)
        n = 160; % number of clicks
        x = linspace(0,8,n)'; % x-axis for smoothing (click number)
        w = ones(size(x)); % weighting factor
        pp = csaps(x,d1,sm,[],w); % piecewise polynomial object
        thdAmp = ppval(pp,thdOnsetTime);
        
        % hysteresis ---------------------
        try
            pp = csaps(x,d1,sm,[],w); % piecewise polynomial object
            hh = fnint(pp); % integrate the smoothed spline
            A = ppval(hh,thdOnsetTime);
            B = ppval(hh,peakX);
            C = ppval(hh,thdOffsetTime);
            aucLeft = B-A; % area under the curve left
            aucRight = C-B; % area under the curve right
            hyst = aucRight / aucLeft; % hysteresis as a ratio of area under the curves
            if hyst < 0
                keyboard
            end
            hysteresis = hyst;
        catch
            hysteresis = NaN;
        end
        % Calculate the slopes ------------------
        try
            peak_index = round(peakX /.05);
            index_xa = round(thdOnsetTime / .05);
            index_xb = round(thdOffsetTime / .05);
            part1_x = x(index_xa:peak_index);
            part1_y = d1(index_xa:peak_index);
            part2_x = x(peak_index+1:index_xb);
            part2_y = d1(peak_index+1:index_xb);
            slope_ascending = polyfit(part1_x,part1_y,1);
            slope_descending = polyfit(part2_x,part2_y,1);
            slopeUp = slope_ascending(1);
            slopeDn = slope_descending(1);
        catch
            slopeUp = NaN;
            slopeDn = NaN;
        end
    end
%---------------------------------------------------------------------------

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
    MEMR.CEOAE = CEOAE;

    MEMR.peakTime = peakX;
    MEMR.peakAmp = peakY;
    MEMR.delay = delay;
    MEMR.thdOnsetTime = thdOnsetTime;
    MEMR.thdOffsetTime = thdOffsetTime;
    MEMR.thdOnsetLvl = thdOnsetLvl;
    MEMR.thdOffsetLvl = thdOffsetLvl;
    MEMR.hysteresis = hysteresis;
    MEMR.slopeUp = slopeUp;
    MEMR.slopeDn = slopeDn;
    MEMR.thd = thd;
    MEMR.thdAmp = thdAmp;
    MEMR.elicitorLevel = elicitorLevel;
    

end

% [EOF]

% OLD CODE ------------------------------    %
