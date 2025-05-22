function [MEMR_inc,MEMR_mem,h1] = analyzeMEMR_ntrendv16(header,Clicks,headerN,Noise,subjectName,runNumber)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [analysisPart1_inc,analysisPart1_mem] = analyzeMEMR_v16(pathName,subjectName,runNumber)
%
% Analyze MEMR data, new goodman test.
% Call this function using goodmanMEMR_analysis_v2.m.
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
% Last Updated: March 18, 2023 -- ssg -- Version 6
% Last Updated: April 17, 2023 -- ssg -- version 7
% Last Updated: September 15, 2023 -- version 8
% Last Updated: September 25, 2023 -- version 9
% Last Updated: October 2, 2023 -- version 11 -- solving the drift problem
% Last Updated: October 3, 2023 -- ssg -- version 12. S''till to do: Analyze
%               click peak to show that it is stable.
% Last Updated: October 20, 2023 -- shows that change in ZL is the same
% result as change in SPL.
% Last Updated: October 26, 2023 -- getting noise windows to estimate snr
%               of integration
% Last Updated: November 13, 2023 -- getting ready to batch file.
%               Updated to version 15
% Last Updated: January 16, 2024 -- ssg
% Last Updated: January 22, 2024 -- ssg
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
    if strcmp(subjectName,'MEM10') & runNumber == 1
        doJitterFix = 0;
    end
    if strcmp(subjectName,'MEM37') & runNumber == 1
        doJitterFix = 0;
    end
    if strcmp(subjectName,'MEM12') & runNumber == 1
        doJitterFix = 0;
    end
    if doJitterFix == 1
        Clicks = unjitter(Clicks,subjectName);
    end

    chunkSize = header.userInfo.chunkSize; % 9600
    nChunks = header.userInfo.nChunks; % 80
    nSweeps = header.userInfo.nSweeps; % 15
    Clicks = reshape(Clicks,chunkSize*nChunks,nSweeps);
    Noise = reshape(Noise,chunkSize*nChunks,nSweeps);
    fs = header.fs; % sampling rate (96 kHz)
    %modLen = 4; % modulation time (4 seconds increasing, 4 seconds decreasing)
    time = (0:1:size(Clicks,1)-1)'/fs; % time vector
    indx = header.userInfo.clickIndx; % location of the clicks
    clear data %header
    time = time / time(end); % express time from 0 to 1
    %plot(20*log10(abs(Noise(:,4))/.00002),'k')
    elicitor = Noise(:,4); 
    % ----> what is special about this? Nothing, we just chose an exemplar for plotting


    % When considering which portions of each time window to analyze:
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
    % The noise floor in the ear canal is roughly 60 dB SPL.
    % The noise itself goes from 50 dB SPL to 120, covering a change of 70
    % dB in 4 seconds, so a rate of 17.5 dB/second
    h = H(:,1)*(10^(120/20)*.00002); % h is the vector of noise levels (intended, in Pa)
    %plot(20*log10(abs(Noise(:,1))/.00002))
    %hold on
    %plot(20*log10(h/.00002),'r','LineWidth',2)
    noiseLvl = 20*log10(h/.00002); % noise level (intended, dB SPL)
    % --> what are these numbers???    
    noiseRate = (115 - 65) / 4; % noise changes at 12.5 dB per second
            % multiply this value by the timeShift value to see how many dB
            % to subtract from the noise vector when calculating thresholds

    % run the data analysis program
    %clicks = Clicks(1:2000,:);
    Clicks = ARLas_hpFilter(Clicks,fs,100);


    % get the stimulus levels ---------------------------------------------
    t = (0:1:length(elicitor)-1)'/fs ;
    [rows,cols] = size(elicitor);
    chunkSize = round(fs*0.05);
    E = reshape(elicitor,chunkSize,rows/chunkSize);
    C = reshape(Clicks(:,1),chunkSize,rows/chunkSize);
    T = reshape(t,chunkSize,rows/chunkSize);
    [rows,cols] = size(E);
    for ii=1:cols
        RMS(ii,1) = sqrt(mean(E(:,ii).^2)); % elicitor RMS
        RMSC(ii,1) = sqrt(mean(C(:,ii).^2)); % rms of the click
        RMST(ii,1) = mean(T(:,ii)); % mean time of each 50 ms chunk
    end
    pSPL = 20*log10(max(C(:,ii)/.00002));
    rmsSPL = 20*log10(RMSC/.00002);

    %figure
    %plot(RMST,20*log10(RMS/.00002))
    %hold on
    %line([0,8],[pSPL,pSPL],'Color',[1 0 0])
    %line([0,8],[rmsSPL,rmsSPL],'Color',[0 1 0])
    %----------------------------------------------------------------------

    % analyze MEMR --------------------------------------------------------
    stabilityCheck = 0;
    MEMR_mem = runme(Clicks,nClicks,indx,time,fs,stabilityCheck,subjectName,runNumber);
    % analyze ear canal stability
    stabilityCheck = 1;
    MEMR_inc = runme(Clicks,nClicks,indx,time,fs,stabilityCheck,subjectName,runNumber);

    MEMR_mem.elicitor = elicitor;
    MEMR_inc.elicitor = elicitor;
    MEMR_mem.RMS = RMS; % rms of the elicitor Pa
    MEMR_mem.RMSC = RMSC; % rms of the click Pa
    MEMR_mem.RMST = RMST; % rms time vector Pa
    MEMR_mem.pSPL = pSPL; % click peak in pSPL
    MEMR_mem.rmsSPL = rmsSPL; % click rms in SPL

    %  plotting -----------------------------------------
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
        Q = repmat(20*log10(MEMR_mem.d1),1,15);
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
        keyboard
    end

    pause(0.01)

end
% INTERNAL FUNCTIONS ------------------------------------------------------
function [MEMR] = runme(Clicks,nClicks,indx,time,fs,stabilityCheck,subjectName,runNumber)
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
    for ii=1:nSweeps % loop across 15 sweeps
        % for signal
        cc = squeeze(C(:,ii,:));

% % show that same thing can be done in rms
% cc80 = squeeze(C(:,:,80));
% ccb = squeeze(C(:,:,1));
% cc80 = mean(cc80,2);
% ccb = mean(ccb,2);
% b = bpf;
% cc80 = fastFilter(b,cc80);
% ccb = fastFilter(b,ccb);
% 
% CC80 = fft(cc80);
% CCB = fft(ccb);
% DD = (CC80 ./ CCB)-1;
% dd = (cc80 - ccb);
% 
% z = sqrt(mean(dd.^2));
% z = z / sqrt(mean(ccb.^2));
% 
% Z = sqrt(mean(abs(DD).^2,1));

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

        %q = squeeze(CFT(5,:,:));
        %q = q.';
        %plot(abs(q(:)))
        %hold on
        
        % find long-term ipsilateral trends and extract them ---------

% % this is for the ARO talk, fig 3 ----
%     % call using MEMR_fig3.m
%     [CFT2,TREND,W] = memrDetrend(CFT);
% 
%     jj = 15;
%     m = squeeze(CFT(jj,:,:)); % no detrending
%     q = abs(m.');
%     ttt = squeeze(TREND(jj,:,:));
%     ttt = abs(ttt.');
%     m2 = squeeze(CFT2(jj,:,:)); % with detrending
%     q2 = abs(m2.');
%     longTime = (0:1:length(q2(:))-1)'*0.05;
%     shortTime = (0:1:size(q2,1)-1)'*0.05;
% 
%     ymax = .032; % for mem
%     ymin = .026;
%     %ymax = .23; % for incident
%     %ymin = .17;
% 
%    figure(3)
%    subplot(2,1,1)
%     plot(longTime,q(:),'Color',[.7 .7 .7])
%     hold on
%     plot(longTime,ttt(:),'b','LineWidth',2)
%     xlim([0,120])
%     ylim([ymin,ymax])
%     xlabel('Time (s)','FontSize',14)
%     hy1 = ylabel('Magnitude','FontSize',14);
%     %hy1.Position = [ -0.3991    0.0000   -1.0000];
%     set(gca,'fontsize',11)
%     xticks((8:16:120))
%     grid on
%     box on
%   subplot(2,1,2)
%     plot(longTime,q2(:),'Color',[.5 .5 .5])
%     xlim([0,120])
%     ylim([ymin,ymax])
%     xlabel('Time (s)','FontSize',14)
%     hy1 = ylabel('Magnitude','FontSize',14);
%     %hy1.Position = [ -0.3991    0.0000   -1.0000];
%     set(gca,'fontsize',11)
%     xticks((8:16:120))
%     grid on
%     box on
%    tbA = annotation('textbox',...
%     [0.15 0.65 0.3 0.15],...
%     'String',{'A'},...
%     'FontSize',18,...
%     'FontName','Arial',...
%     'LineStyle','none',...
%     'EdgeColor',[0 0 0],...
%     'Position',[0.8554    0.8571    0.0580    0.0810]);
%     tbB = annotation('textbox',...
%     [0.15 0.65 0.3 0.15],...
%     'String',{'B'},...
%     'FontSize',18,...
%     'FontName','Arial',...
%     'LineStyle','none',...
%     'EdgeColor',[0 0 0],...
%     'Position',[0.8554    0.3809    0.0580    0.0810]);
% 
%   % ------ % show smoothing
%     w = squeeze(W(jj,:,:));
%     x = (0:1:nClicks-1);
%     xxx = (1:1:nClicks*3)-nClicks;
%     Mr = real(m2); % matrix of 15 sweeps at one click
%     Mi = imag(m2);
%     sumw = sum(w,1); % take the weighted mean at one click (and one frequency)
%     mr = sum(Mr(:,:).*w,1)./ sumw;
%     mi = sum(Mi(:,:).*w,1)./ sumw;
%     % smooth the response (effectively lowpass filter)
%     mrrr = [mr,mr,mr];
%     miii = [mi,mi,mi];
%     www = [sumw,sumw,sumw];
%     sm = 0.001; % was 0.0001; % was 0.00001;  smoothing factor (smaller numbers are more smooth)
%     ppr = csaps(xxx,mrrr,sm,[],www); % piecewise polynomial real coefficients
%     ppi = csaps(xxx,miii,sm,[],www); % piecewise polynomial imaginary coefficients
%     mr_sm = ppval(ppr,x); % evaluate only at the original x values
%     mi_sm = ppval(ppi,x);
%     z = mr + 1i*mi;
%     z_sm = mr_sm + 1i*mi_sm; % smoothed raw Z
% 
%     baseline = mean([z_sm(1:3),z_sm(end-2:end)]);
%     d = z_sm ./ baseline; % normalize z_sm
%    
%    figure(4)
%    subplot(2,1,1)
%     plot(shortTime,q,'Color',[.7 .7 .7])
%     hold on    
%     plot(shortTime,abs(z_sm),'Color',[0 0 1],'LineWidth',2)
%     %plot(shortTime,mean(q,2),'r')
%     xlim([0,8])
%     ylim([ymin,ymax])
%     xlabel('Time (s)','FontSize',14)
%     hy1 = ylabel('Magnitude','FontSize',14);
%     %hy1.Position = [ -0.3991    0.0000   -1.0000];
%     set(gca,'fontsize',11)
%     xticks((0:1:8))
%     grid on
%     box on
%   subplot(2,1,2)
%     plot(shortTime,abs(d),'Color',[0 0 1],'LineWidth',2)
%     %plot(shortTime,mean(q,2),'r')
%     xlim([0,8])
%     ylim([0.9 1.01])
%     xlabel('Time (s)','FontSize',14)
%     hy1 = ylabel('Normalized Mag.','FontSize',14);
%     hy1.Position = [-0.8897    0.9550   -1.0000];
%     set(gca,'fontsize',11)
%     xticks((0:1:8))
%     grid on
%     box on    
% 
%   tbA = annotation('textbox',...
%     [0.15 0.65 0.3 0.15],...
%     'String',{'A'},...
%     'FontSize',18,...
%     'FontName','Arial',...
%     'LineStyle','none',...
%     'EdgeColor',[0 0 0],...
%     'Position',[0.8554    0.8571    0.0580    0.0810]);
% 
%   tbB = annotation('textbox',...
%     [0.15 0.65 0.3 0.15],...
%     'String',{'B'},...
%     'FontSize',18,...
%     'FontName','Arial',...
%     'LineStyle','none',...
%     'EdgeColor',[0 0 0],...
%     'Position',[0.8554    0.3452    0.0580    0.0810]);

%     keyboard
% % end for ARO talk, fig 3 ------------




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
    Fstop1 = 10;              % First Stopband Frequency
    Fpass1 = 500;             % First Passband Frequency
    Fpass2 = 1500;            % Second Passband Frequency
    Fstop2 = 2000;            % Second Stopband Frequency
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

% [EOF]

% OLD CODE ------------------------------    %
% calculate low-frequency drift across memr recordings ----------------
%     for ii=1:nFreqs % do this individually for each frequency
%         cc = squeeze(CFT(ii,:,:)); % pick the current frequency
%         mr = real(cc); % the real part
%         mr = mr'; % transpose
%         [rows,cols] = size(mr); % use this later to reshape ack
%         mr = mr(:); % force to a single column
%         mi = imag(cc); % now do the same for the imaginary part
%         mi = mi';
%         mi = mi(:);
%         x = (1:1:length(mr))'; % x-axis vector for spline fitting
%         w = ones(size(mr)); % weighting vector (set to ones)
%         
%         doPlotAR = 0;
%         multiplier = 1.5;
%         [rejectIndx,nRejects] = newAR2(mr',multiplier,doPlotAR); 
%         w(rejectIndx) = 0;
%         
%         smoothing = 0.00000000001; % smoothing factor--0 is a straight line, 1 is cubic spline
%         ppr = csaps(x,mr,smoothing,[],w); % piecewise polynomial real coefficients
%         ppi = csaps(x,mi,smoothing,[],w); % piecewise polynomial imaginary coefficients
%         mr_sm = ppval(ppr,x);
%         mi_sm = ppval(ppi,x);
%         
%         %         figure
%         %         plot(mr.*w,'b.')
%         %         hold on
%         %         plot(mr_sm,'r')
%         %         keyboard
% 
%         trend = mr_sm+1i*mi_sm;
%         offset = mean(trend); % mean location before detrending
%         mrfixed = mr - mr_sm; % "fixed", i.e., detrended version of mr
%         mifixed = mi - mi_sm;
%         mrfixed = reshape(mrfixed,rows,cols);
%         mifixed = reshape(mifixed,rows,cols);
%         trend = reshape(trend,rows,cols);
%         w = reshape(w,rows,cols);
%         z = mrfixed + 1i*mifixed + offset; % add the offset back in to the complex form
%         CFT(ii,:,:) = z.';
%         TREND(ii,:,:) = trend.';
%         W(ii,:,:) = w.';
%     end
% 
%     % calculate drift across the noise window ----------------------------
%     for ii=1:nFreqs % do this individually for each frequency
%         cc = squeeze(CFTn(ii,:,:)); % pick the current frequency
%         mr = real(cc); % the real part
%         mr = mr'; % transpose
%         [rows,cols] = size(mr); % use this later to reshape ack
%         mr = mr(:); % force to a single column
%         mi = imag(cc); % now do the same for the imaginary part
%         mi = mi';
%         mi = mi(:);
%         x = (1:1:length(mr))'; % x-axis vector for spline fitting
%         w = ones(size(mr)); % weighting vector (set to ones)
%         
%         doPlotAR = 0;
%         multiplier = 1.5;
%         [rejectIndx,nRejects] = newAR2(mr',multiplier,doPlotAR); 
%         %if ~isempty(rejectIndx)
%         %    keyboard
%             w(rejectIndx) = 0;
%         %end
%         
%         smoothing = 0.00000000001; % smoothing factor--0 is a straight line, 1 is cubic spline
%         ppr = csaps(x,mr,smoothing,[],w); % piecewise polynomial real coefficients
%         ppi = csaps(x,mi,smoothing,[],w); % piecewise polynomial imaginary coefficients
%         mr_sm = ppval(ppr,x);
%         mi_sm = ppval(ppi,x);
%         
%         trend = mr_sm+1i*mi_sm;
%         offset = mean(trend); % mean location before detrending
%         mrfixed = mr - mr_sm; % "fixed", i.e., detrended version of mr
%         mifixed = mi - mi_sm;
%         mrfixed = reshape(mrfixed,rows,cols);
%         mifixed = reshape(mifixed,rows,cols);
%         trend = reshape(trend,rows,cols);
%         w = reshape(w,rows,cols);
%         z = mrfixed + 1i*mifixed + offset; % add the offset back in to the complex form
%         CFTn(ii,:,:) = z.';
%         TRENDn(ii,:,:) = trend.';
%         Wn(ii,:,:) = w.';
%     end
% 
% 
%     % this just to double check we did it right
%     % for ii=1:nFreqs
%     %     for jj=1:nClicks
%     %     cc = squeeze(CFT(jj,:,jj));
%     %     ccend = squeeze(CFT(ii,:,jj));
%     %     cc1 = squeeze(CFT(ii,:,1));
%     %     plot(cc1)
%     %     hold on
%     %     plot(ccend)        
%     % end
% 
%     % calculate weighted means --------------------------------------------
%     x = (1:1:nClicks);
%     %Mr_sm = zeros(nClicks,nFreqs); % initializing the outputs
%     %Mi_sm = zeros(nClicks,nFreqs);
%     %Mrn_sm = Mr_sm;
%     %Min_sm = Mi_sm;
%     for jj=1:nFreqs
%         m = squeeze(CFT(jj,:,:)); % size m = 160 x 15
%         w = squeeze(W(jj,:,:));
%         Mr = real(m); % matrix of 15 sweeps at one click
%         Mi = imag(m);
%         sumw = sum(w,1); % take the weighted mean at one click (and one frequency)
%         mr = sum(Mr(:,:).*w,1)./ sumw;
%         mi = sum(Mi(:,:).*w,1)./ sumw;
%         % smooth the response (effectively lowpass filter)
%         sm = 0.00001;
%         ppr = csaps(x,mr,sm,[],sumw); % piecewise polynomial real coefficients
%         ppi = csaps(x,mi,sm,[],sumw); % piecewise polynomial imaginary coefficients
%         mr_sm = ppval(ppr,x);
%         mi_sm = ppval(ppi,x);
% 
%         %z = mr_sm + 1i*mi_sm;
%         z = mr + 1i*mi;
%         %z = z ./ z(1);
%         Z(:,jj) = z;
% 
% 
% 
%         m = squeeze(CFTn(jj,:,:)); % size m = 160 x 15
%         w = squeeze(Wn(jj,:,:));
%         Mr = real(m); % matrix of 15 sweeps at one click
%         Mi = imag(m);
%         sumw = sum(w,1); % take the weighted mean at one click (and one frequency)
%         mr = sum(Mr(:,:).*w,1)./ sumw;
%         mi = sum(Mi(:,:).*w,1)./ sumw;
%         % smooth the response (effectively lowpass filter)
%         %if any(sumw == 0)
%         %    keyboard
%         %end
%         ppr = csaps(x,mr,sm,[],sumw); % piecewise polynomial real coefficients
%         ppi = csaps(x,mi,sm,[],sumw); % piecewise polynomial imaginary coefficients
%         mrn_sm = ppval(ppr,x);
%         min_sm = ppval(ppi,x);
% 
%         %zn = mrn_sm + 1i*min_sm;
%         zn = mr + 1i*mi;
%         %zn = zn ./ z(1);
%         Zn(:,jj) = zn;
%     end
% 
%     analysisPart1.Z = Z;
%     analysisPart1.Zn = Zn;
%     analysisPart1.TREND = TREND;
%     analysisPart1.TRENDn = TRENDn;
%     analysisPart1.nFreqs = nFreqs;
%----------------------------------
        %q = squeeze(CFT(5,:,:));
        %q = q.';
        %plot(abs(q(:)))
        %t = squeeze(TREND(5,:,:));
        %t = t.';
        %plot(abs(t(:)),'g')
        
        % Xk = squeeze(CFT(:,:,1));
        % K = 15;
        % Xbar = (1/K) * sum(Xk,2); % signal is the coherent mean
        % Xbar2 = abs(Xbar) .^2; % signal energy
        % Xbar = (1/K) * sum(Xk,2); % signal is the coherent mean
        % Xbar2 = abs(Xbar) .^2; % signal energy
        % XBAR = repmat(Xbar,1,K); % replicate to matrix (to avoid running a loop)
        % S2 = (1/(K-1)) * sum((Xk - XBAR) .* conj(Xk - XBAR),2); % variance (this is the noise floor)
        % Se2 = (1/K) * S2; % energy of the standard error
        % ref = 0.00002;
        % signal1 = 10*log10(Xbar2/(ref^2));
        % %noise = 10*log10(S2/(ref^2));    
        % noiseFloor1 = 10*log10(Se2/(ref^2)); 
        % 
        % Xk = squeeze(CFT(:,:,80));
        % K = 15;
        % Xbar = (1/K) * sum(Xk,2); % signal is the coherent mean
        % Xbar2 = abs(Xbar) .^2; % signal energy
        % Xbar = (1/K) * sum(Xk,2); % signal is the coherent mean
        % Xbar2 = abs(Xbar) .^2; % signal energy
        % XBAR = repmat(Xbar,1,K); % replicate to matrix (to avoid running a loop)
        % S2 = (1/(K-1)) * sum((Xk - XBAR) .* conj(Xk - XBAR),2); % variance (this is the noise floor)
        % Se2 = (1/K) * S2; % energy of the standard error
        % ref = 0.00002;
        % signal80 = 10*log10(Xbar2/(ref^2));
        % %noise = 10*log10(S2/(ref^2));    
        % noiseFloor80 = 10*log10(Se2/(ref^2));
        % 
        % d = noiseFloor80 - noiseFloor1;
        % [nfMax,nfMaxIndx] = max(d);
        % d = signal80 - signal1;
        % [sMax,sMaxIndx] = max(abs(d));
        % nfF = frequency(nfMaxIndx);
        % sF = frequency(sMaxIndx);
    
    %octStepSize = 1/2; % octave step size for analysis
    %fmin = 125; % minimum frequency to examine
    %fmax = 8000; % maximum frequency to examine
    % frequencies of interest:
    %peakFreqs = round(2.^(log2(fmin):octStepSize:log2(fmax))'); % analysis center frequencies (Hz)


            % Standard magnitude analysis only -----------------------------------------
%         for ii=1:160 % loop across time
%             XX = squeeze(CFT(:,:,ii));
%             WW = squeeze(W(:,:,ii));
%             for jj=1:nFreqs % loop across frequencies
%                 Xk = XX(jj,:);
%                 Wk = WW(jj,:);
%                 badIndx = find(Wk==0);
%                 if ~isempty(badIndx)
%                     %keyboard
%                     Xk(badIndx) = [];
%                 end
%                 K = length(Xk);
%                 Xbar = (1/K) * sum(Xk,2); % signal is the coherent mean
%                 Xbar2 = abs(Xbar) .^2; % signal energy
%                 Xbar = (1/K) * sum(Xk,2); % signal is the coherent mean
%                 Xbar2 = abs(Xbar) .^2; % signal energy
%                 XBAR = repmat(Xbar,1,K); % replicate to matrix (to avoid running a loop)
%                 S2 = (1/(K-1)) * sum((Xk - XBAR) .* conj(Xk - XBAR),2); % variance (this is the noise floor)
%                 Se2 = (1/K) * S2; % energy of the standard error
%                 ref = 0.00002;
%                 signal(ii,jj) = 10*log10(Xbar2/(ref^2));
%                 noiseFloor(ii,jj) = 10*log10(Se2/(ref^2));
%             end
%         end
% 
% 
%         % apply smoothing -------------------------------------------------
%         sm = 0.95; % smoothing factor (smaller numbers are more smooth)
%         n = 160; % number of clicks
%         x = linspace(0,8,n)'; % x-axis for smoothing (click number)
%         xxx = linspace(0,8*3,n*3)';
%         for ii=1:size(signal,2) % loop across frequencies
% 
%             q = signal(:,ii); % get the mean (across frequencies) growth function
%             w = ones(size(x)); % weighting factor
%                 doPlotAR = 0;
%                 multiplier = 1.5;
%                 [rejectIndx,~] = newAR2(q',multiplier,doPlotAR); 
%                 w(rejectIndx) = 0;
%             qqq = [q;q;q];
%             www = [w;w;w];
%             pp = csaps(xxx,qqq,sm,[],www); % piecewise polynomial object
%             signal_sm(:,ii) = ppval(pp,x); 
% 
%             q = noiseFloor(:,ii); % get the mean (across frequencies) growth function
%             w = ones(size(q));
%                 doPlotAR = 0;
%                 multiplier = 1.5;
%                 [rejectIndx,~] = newAR2(q',multiplier,doPlotAR); 
%                 w(rejectIndx) = 0;
%             qqq = [flipud(q);q;q];
%             www = [flipud(w);w;w];
%             pp = csaps(xxx,qqq,sm,[],www); % piecewise polynomial object
%             noiseFloor_sm(:,ii) = ppval(pp,x);             
%         end
% 
%         % shift starting point to zero for all frequencies
%         aveNum = 10;
%         for ii=1:size(signal,2) % loop across frequencies
%             q = signal_sm(:,ii);
%             offset = mean([q(1:aveNum);q(end-aveNum+1:end)]);
%             q = q - offset;
%             signal_sm(:,ii) = q;
%             q = signal(:,ii);
%             q = q - offset;
%             signal(:,ii) = q;
% 
%             q = noiseFloor_sm(:,ii);
%             offset = mean([q(1:aveNum);q(end-aveNum+1:end)]);
%             q = q - offset;
%             noiseFloor_sm(:,ii) = q;
%             q = noiseFloor(:,ii);
%             q = q - offset;
%             noiseFloor(:,ii) = q;
%         end
% 
%         % sort smallest to largest
%         [~,sortIndx] = sort(nanmean(abs(signal_sm)));
%         signal_sorted = signal_sm(:,sortIndx);
%         freqSig_sorted = frequency(sortIndx);
%         sF = mean(freqSig_sorted(end));
% %         trend = squeeze(TREND(sortIndx(end),:,:));
% %         trend = trend.';
% %         trend = trend(:);
% %         trend = abs(trend);
% %         trend = trend ./ trend(1);
% %         trend = 20*log10(trend);
% % 
% %         trendn = squeeze(TRENDn(sortIndx(end),:,:));
% %         trendn = trendn.';
% %         trendn = trendn(:);
% %         trendn = abs(trendn);
% %         trendn = trendn ./ trendn(1);
% %         trendn = 20*log10(trendn);
% 
%         [~,sortIndx] = sort(nanmean(abs(noiseFloor_sm)));
%         noiseFloor_sorted = noiseFloor_sm(:,sortIndx);
%         freqNf_sorted = frequency(sortIndx);
%         nF = mean(freqNf_sorted(end));
%         sigHat = abs(signal_sorted(:,end));
%         %sigHat = mean(signal_sorted(:,24:26),2);
%         %sigHat = mean(abs(signal_sorted),2);
%         sigHatNF = mean(signal_sorted(:,1:3),2);
%         %nfHat = mean(noiseFloor_sorted(:,24:26),2);
%         nfHat = mean(abs(noiseFloor_sorted),2);
%         nfHatNF = mean(noiseFloor_sorted(:,1:3),2);
% 
%         % [~,maxSigIndx] = max(nanmean(abs(signal_sm)));
%         % [~,maxNfIndx] = max(nanmean(abs(noiseFloor_sm)));
%         % sF = frequency(maxSigIndx);
%         % nF = frequency(maxNfIndx);
%         % [~,minSigIndx] = min(nanmean(abs(signal_sm)));
%         % [~,minNfIndx] = min(nanmean(abs(noiseFloor_sm)));
% 
% 
%         figure
%         subplot(2,1,1)
%         %plot(signal,'Color',[.7 .7 .7])
%         hold on
%         plot(signal_sm,'Color',[0 0 1])
%         %plot(abs(signal_sm(:,maxSigIndx)),'g','LineWidth',2)
%         %plot(abs(signal_sm(:,minSigIndx)),'k','LineWidth',1)
%         plot(abs(sigHat),'g','LineWidth',2)
%         plot(abs(sigHatNF),'k','LineWidth',1)
%         title([num2str(sF),' Hz'])
%         xlim([1 160])
%         ymax = max(sigHat);
%         ymax = ymax * 1.1;
%         ymin = -ymax;
%         ylim([ymin,ymax])
%         line([thdOnIndx,thdOnIndx],[ymin,ymax],'Color',[1 0 0],'LineWidth',1)
%         line([thdOffIndx,thdOffIndx],[ymin,ymax],'Color',[1 0 0],'LineWidth',1)
% 
%         subplot(2,1,2)
% 
%         %plot(noiseFloor,'Color',[.7 .7 .7])
%         hold on
%         plot(abs(noiseFloor_sm),'Color',[1 0 0])
%         %plot(abs(noiseFloor_sm(:,maxNfIndx)),'g','LineWidth',2)
%         %plot(abs(noiseFloor_sm(:,minNfIndx)),'k','LineWidth',1)
%         plot(abs(nfHat),'g','LineWidth',2)
%         plot(abs(nfHatNF),'k','LineWidth',1)        
%         title([num2str(nF),' Hz'])
%         xlim([1 160])
% 
% %         figure
% %         Q = repmat(sigHat,1,15);
% %         Q = Q(:);
% %         plot(Q+trend)
% %         hold on
% %         plot(trend,'r')
% 
% %         keyboard
% 
% %         analysisPart1.Z = Z;
% %         analysisPart1.Zn = Zn;
% %         analysisPart1.TREND = TREND;
% %         analysisPart1.TRENDn = TRENDn;
% %         analysisPart1.nFreqs = nFreqs;
% 
% 

    % analysisPart1.signal_sm = signal_sm;
    % analysisPart1.sigHat = sigHat;
    % analysisPart1.sigHatNF = sigHatNF;
    % analysisPart1.sF = sF;
    % analysisPart1.thdOnIndx = thdOnIndx;
    % analysisPart1.thdOffIndx = thdOffIndx;
    % analysisPart1.noiseFloor_sm = noiseFloor_sm;
    % analysisPart1.nfHat = nfHat;
    % analysisPart1.nfHatNF = nfHatNF;
    % analysisPart1.nF = nF;
    % stuff we intend to use


%     % find thresholds, max, and hystersis ---------------------------------
% 
%         % find threshold by finding where largest response crosses the
%         % noise, defined by the largest value in the first 20 samples
%         d1dB = 20*log10(d1);
%         cutoff = 0.25; % threshold cutoff in dB. Determined emperically
%         [hat,hatIndx] = max(d1dB);
%         done = 0;
%         counter = 1;
%         indx = hatIndx;
%         while done == 0
%             q = d1dB(indx);
%             if q > cutoff
%                 indx = indx - 1;
%             else
%                 done = 1;
%             end
%             if indx < 1
%                 done = 1;
%             else
%                 counter = counter + 1;
%             end
%         end
%         thdOnIndx = indx; % location of onset threshold
%         indx = hatIndx;
%         done = 0;
%         counter = 1;
%         while done == 0
%             q = d1dB(indx);
%             if q > cutoff
%                 indx = indx + 1;
%             else
%                 done = 1;
%             end
%             if indx > 160
%                 done = 1;
%             else
%                 counter = counter + 1;
%             end
%         end        
%         thdOffIndx = indx;
