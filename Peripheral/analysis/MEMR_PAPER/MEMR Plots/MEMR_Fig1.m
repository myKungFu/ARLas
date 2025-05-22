function [] = MEMR_Fig1()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [] = MEMR_Fig1()
%
% Figure 1 for paper and ARO talk
%
% Author: Shawn Goodman
% Date: January 15, 2024
% Last Updated: January 15, 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
subjectName = 'MEM35'; %'MEM10'; % 'MEM12', 'MEM09', 'MEM17'
runNumber = 1;

    % 20 

    % ---------------------------------------------------------------------
    doJitterFix = 1;
    parentDrive = 'D'; % use 'C' or 'D'
    dataPathName = [parentDrive,':\MEMR_DATA\']; % location of raw data
    fileNameL = 'Ch3_ER10xA_memr_0001.mat';
    fileNameR = 'Ch4_ER10xB_memr_0001.mat';
    savePath = [parentDrive,':\MEMR_Analysis1\']; % where to save analyzed data
    
    % ---------------------------------------------------------------------

    d = dir([dataPathName,subjectName]); % read in file directory information
    dummy = load([dataPathName,subjectName,'\',d(runNumber+2).name,'\',fileNameL]); % the clicks
    header = dummy.header;
    Clicks = dummy.data;
    clear dummy
    dummy = load([dataPathName,subjectName,'\',d(runNumber+2).name,'\',fileNameR]); % the noise
    headerN = dummy.header;
    Noise = dummy.data;
    clear dummy 


    micCorrection = header.userInfo.micCorrectionL;
    Clicks = applyMicCorrection(Clicks,micCorrection);
    micCorrectionN = headerN.userInfo.micCorrectionR;
    Noise = applyMicCorrection(Noise,micCorrectionN);

    
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
    %nSamples = 193; % number of time samples in each click analysis window
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

% -------------------------------------------------------------------------
    % run the data analysis program
    %clicks = Clicks(1:2000,:);
    Clicks = ARLas_hpFilter(Clicks,fs,100);

    t = (0:1:length(elicitor)-1)'/fs ;
    Elicitor = randn(length(t),1).*H(:,1);
    Elicitor = Elicitor / max(abs(Elicitor));
    % based on group data (found using batchMEMR_groupData.m), the mean rms
    % SPL for elicitor is 109 (calculated in 50 ms chunks). This equals
    % 1.98 pascals. But there will be peaks of up to 20 pascals (120 dB SPL)
    Elicitor = Elicitor * 25;
    % based on group data, the mean click level was 96.26 pSPL, which is
    % 1.3 Pascals
    Clicks = Clicks(:,4);
    Clicks = Clicks / (max(abs(Clicks)));
    Clicks = Clicks * 1.3;

 figure(1) % two subplots showing clicks and noise in Pa ------------------
   subplot(2,1,1)   
    plot(t,Clicks,'k')
    xlabel('Time (s)','FontSize',14)
    ylabel('Pressure (Pa)','FontSize',14)
    set(gca,'fontsize',11)
    xlim([0,8])
    ylim([-0.5 2])
   subplot(2,1,2)
    plot(t,Elicitor,'k')
    xlabel('Time (s)','FontSize',14)
    ylabel('Pressure (Pa)','FontSize',14)
    set(gca,'fontsize',11)
    xlim([0,8])
    ylim([-20 20])
    grid on

    tbA = annotation('textbox',...
    [0.15 0.65 0.3 0.15],...
    'String',{'A'},...
    'FontSize',18,...
    'FontName','Arial',...
    'LineStyle','none',...
    'EdgeColor',[0 0 0],...
    'Position',[0.1286    0.8571    0.0580    0.0810]);

    tbB = annotation('textbox',...
    [0.15 0.65 0.3 0.15],...
    'String',{'B'},...
    'FontSize',18,...
    'FontName','Arial',...
    'LineStyle','none',...
    'EdgeColor',[0 0 0],...
    'Position',[0.1286    0.3809    0.0580    0.0810]);





keyboard

    parentDrive = 'C'; % use 'C' or 'D'
    dataPathName = [parentDrive,':\MEMR_Analysis1\']; % location of raw data
    fileName = 'MEM24_Run1_Analysis1.mat';
    dummy = load([dataPathName,fileName]); 
    MEMR_mem = dummy.MEMR_mem;
    d1 = MEMR_mem.d1;
    D1 = 20*log10(d1);

    h3 = figure;
    h3.Position = [360.0000  198.0000  560.0000  303.6667];
    plot(D1,'b.-','MarkerSize',10)
    xticks([0,20,40,60,80,100,120,140,160])
    xticklabels({'40','58','75','93','110','93','75','58','40'})
    %levels = 17.5*0.05([0,20,40,60,80,60,40,20,0]')+40
    xlabel('Elicitor Level (dB SPL)','FontSize',14)
    ylabel('\Delta Pressure (dB)','FontSize',14)
    set(gca,'fontsize',11)
    xlim([0,160])
    %ylim([-20 20])
    grid on

    tbC = annotation('textbox',...
    [0.15 0.65 0.3 0.15],...
    'String',{'C'},...
    'FontSize',18,...
    'FontName','Arial',...
    'LineStyle','none',...
    'EdgeColor',[0 0 0],...
    'Position',[0.1286    0.8571    0.0580    0.0810]);
    






%---------------------------------------------------------------------

    % show the windowed clicks
    % get the first click
    [~,indx] = min(abs(t-0.05));
    tt = t(1:indx);
    xx = Clicks(1:indx);
    Indx = header.userInfo.clickIndx; % location of the clicks


    [~,indx1] = min(abs(t-0.002));
    [~,indx2] = min(abs(t-0.005));
    tt = tt(1:indx2);
    xx = xx(1:indx2);

    winN = 14; % number of samples to reduce (before reflection comes back)
    clickN = winN*2+1; % number of analysis samples in each (windowed click)
    h = hann(clickN); % make a hann window to window the edges 
    qqq = xx(Indx(1)-winN:Indx(1)+winN,:);
    qqq = qqq .* h;    
    ttt = tt(Indx(1)-winN:Indx(1)+winN,:);

    clickN = 100 + winN; % number of analysis samples in each (windowed click)
    h = hann(winN*2); % make a hann window to window the edges 
    h = h(1:winN);
    qqqq = xx(Indx(1):Indx(1)+clickN-1,:); 
    tttt = tt(Indx(1):Indx(1)+clickN-1,:); 
    qqqq(1:winN,:) = qqqq(1:winN,:) .* h;
    qqqq = flipud(qqqq);
    qqqq(1:winN,:) = qqqq(1:winN,:) .* h;
    qqqq = flipud(qqqq);


 figure(2) % plot the windowed clicks for incident and reflected pressures
    plot(tt,xx,'Color',[.7 .7 .7],'LineWidth',2)
    hold on
    plot(ttt,qqq,'b')
    plot(tttt,qqqq,'r')

figure
plot(xx,'k')
xlim([280,320])
keyboard



    % get the stimulus levels ---------------------------------------------
    [rows,cols] = size(elicitor);
    chunkSize = round(fs*0.05);
    E = reshape(elicitor,chunkSize,rows/chunkSize);
    C = reshape(Clicks(:,1),chunkSize,rows/chunkSize);
    T = reshape(t,chunkSize,rows/chunkSize);
    [rows,cols] = size(E);
    for ii=1:cols
        RMS(ii,1) = sqrt(mean(E(:,ii).^2));
        RMSC(ii,1) = sqrt(mean(C(:,ii).^2));
        RMST(ii,1) = mean(T(:,ii));
    end
    pSPL = 20*log10(max(xx)/.00002);
    rmsSPL = 20*log10(RMSC/.00002);

 figure(3)
    plot(RMST,20*log10(RMS/.00002))
    hold on
    line([0,8],[pSPL,pSPL],'Color',[1 0 0])
    line([0,8],[rmsSPL,rmsSPL],'Color',[0 1 0])

figure(1) % make this the current figure

keyboard


 %-------------------------------------------------------------------------
    % analyze MEMR
    stabilityCheck = 0;
    analysisPart1_mem = runme(Clicks,nClicks,indx,time,fs,stabilityCheck);
    % analyze ear canal stability
    stabilityCheck = 1;
    analysisPart1_inc = runme(Clicks,nClicks,indx,time,fs,stabilityCheck);

    analysisPart1_mem.elicitor = elicitor;
    analysisPart1_inc.elicitor = elicitor;
   %keyboard
  %   plot1(analysisPart1_inc,analysisPart1_mem)


end
% INTERNAL FUNCTIONS ------------------------------------------------------
function [analysisPart1] = runme(Clicks,nClicks,indx,time,fs,stabilityCheck)
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
        % 8-14 samples
    winN = 14; % number of samples to reduce (before reflection comes back)
    if stabilityCheck == 1
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
    else
        clickN = 100 + winN; % number of analysis samples in each (windowed click)
        h = hann(winN*2); % make a hann window to window the edges 
        h = h(1:winN);
        H = repmat(h,1,nSweeps);
        C = zeros(clickN,nSweeps,nClicks); % initialize matrix of windowed clicks (C for clicks)
        Cn = zeros(clickN+winN,nSweeps,nClicks);
        for ii=1:nClicks % loop across each click position (1-160)
            q = Clicks(indx(ii):indx(ii)+clickN-1,:); % contains all the sweeps
            qn = Clicks(indx(ii)-winN:indx(ii)+clickN-1,:); % window for noise estimate
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
        FT = fft(cc,nfft); % contains all sweeps at all frequencies for a given click position
        FT = FT(indx1:indx2,:); % keep only the frequencies of interest
        FT = FT ./ scale; % scale the magnitude appropriately
        CFT(:,ii,:) = complex(FT); % Fourier transform of C, cut to frequenc6ies of interest

        % for noise 
        ccn = squeeze(Cn(:,ii,:));
        FT = fft(ccn,nfft); % contains all sweeps at all frequencies for a given click position
        FT = FT(indx1:indx2,:); % keep only the frequencies of interest
        CFTn(:,ii,:) = complex(FT); % Fourier transform of C, cut to frequenc6ies of interest
    end

    if stabilityCheck ~=5

        %q = squeeze(CFT(5,:,:));
        %q = q.';
        %plot(abs(q(:)))
        %hold on
        
        [CFT,TREND,W] = memrDetrend(CFT);
        [CFTn,TRENDn,Wn] = memrDetrend(CFTn);

        for jj=1:nFreqs
            dummy = squeeze(TREND(jj,:,:));
            dummy = dummy.';
            dummy = dummy(:);
            dummy = dummy ./ dummy(1);
            dummy = abs(dummy - 1) + 1;
            dummy = 20*log10(dummy);
            Trend(:,jj) = dummy;

            dummy = squeeze(TRENDn(jj,:,:));
            dummy = dummy.';
            dummy = dummy(:);
            dummy = dummy ./ dummy(1);
            dummy = abs(dummy - 1) + 1;
            dummy = 20*log10(dummy);
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
            sm = 0.0001; % was 0.00001;  smoothing factor (smaller numbers are more smooth)
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
             t = x * 0.05;
%             dy = gradient(imag(z_sm)./gradient(t));
%             dx = gradient(real(z_sm)./gradient(t));
%             d = sqrt(dy.^2 + dx.^2);
%             [dmax,dmaxIndx] = max(d);
%             k = ones(size(d));
%             k(dmaxIndx+1:end) = -1;
%             D = cumsum(d.*k);

            d = z_sm ./ mean([z_sm(1:3),z_sm(end-2:end)]);
            d1 = abs(d-1)+1; % the combined mag+phase change
            d2 = abs(abs(d)-1)+1; % magnitude only change


            %z_sm = z_sm ./ mean([z_sm(1),z_sm(end)]);
            Z_sm(:,jj) = z_sm;
            
            %D(:,jj) = D;
            D1(:,jj) = d1; % use this!
            D2(:,jj) = d2;



        end

        % Standard magnitude analysis only -----------------------------------------
        for ii=1:160 % loop across time
            XX = squeeze(CFT(:,:,ii));
            WW = squeeze(W(:,:,ii));
            for jj=1:nFreqs % loop across frequencies
                Xk = XX(jj,:);
                Wk = WW(jj,:);
                badIndx = find(Wk==0);
                if ~isempty(badIndx)
                    %keyboard
                    Xk(badIndx) = [];
                end
                K = length(Xk);
                Xbar = (1/K) * sum(Xk,2); % signal is the coherent mean
                Xbar2 = abs(Xbar) .^2; % signal energy
                Xbar = (1/K) * sum(Xk,2); % signal is the coherent mean
                Xbar2 = abs(Xbar) .^2; % signal energy
                XBAR = repmat(Xbar,1,K); % replicate to matrix (to avoid running a loop)
                S2 = (1/(K-1)) * sum((Xk - XBAR) .* conj(Xk - XBAR),2); % variance (this is the noise floor)
                Se2 = (1/K) * S2; % energy of the standard error
                ref = 0.00002;
                signal(ii,jj) = 10*log10(Xbar2/(ref^2));
                noiseFloor(ii,jj) = 10*log10(Se2/(ref^2));
            end
        end


        % apply smoothing -------------------------------------------------
        sm = 0.95; % smoothing factor (smaller numbers are more smooth)
        n = 160; % number of clicks
        x = linspace(0,8,n)'; % x-axis for smoothing (click number)
        xxx = linspace(0,8*3,n*3)';
        for ii=1:size(signal,2) % loop across frequencies
            
            q = signal(:,ii); % get the mean (across frequencies) growth function
            w = ones(size(x)); % weighting factor
                doPlotAR = 0;
                multiplier = 1.5;
                [rejectIndx,~] = newAR2(q',multiplier,doPlotAR); 
                w(rejectIndx) = 0;
            qqq = [q;q;q];
            www = [w;w;w];
            pp = csaps(xxx,qqq,sm,[],www); % piecewise polynomial object
            signal_sm(:,ii) = ppval(pp,x); 
            
            q = noiseFloor(:,ii); % get the mean (across frequencies) growth function
            w = ones(size(q));
                doPlotAR = 0;
                multiplier = 1.5;
                [rejectIndx,~] = newAR2(q',multiplier,doPlotAR); 
                w(rejectIndx) = 0;
            qqq = [flipud(q);q;q];
            www = [flipud(w);w;w];
            pp = csaps(xxx,qqq,sm,[],www); % piecewise polynomial object
            noiseFloor_sm(:,ii) = ppval(pp,x);             
        end

        % shift starting point to zero for all frequencies
        aveNum = 10;
        for ii=1:size(signal,2) % loop across frequencies
            q = signal_sm(:,ii);
            offset = mean([q(1:aveNum);q(end-aveNum+1:end)]);
            q = q - offset;
            signal_sm(:,ii) = q;
            q = signal(:,ii);
            q = q - offset;
            signal(:,ii) = q;

            q = noiseFloor_sm(:,ii);
            offset = mean([q(1:aveNum);q(end-aveNum+1:end)]);
            q = q - offset;
            noiseFloor_sm(:,ii) = q;
            q = noiseFloor(:,ii);
            q = q - offset;
            noiseFloor(:,ii) = q;
        end
        
        % sort smallest to largest
        [~,sortIndx] = sort(nanmean(abs(signal_sm)));
        signal_sorted = signal_sm(:,sortIndx);
        freqSig_sorted = frequency(sortIndx);
        sF = mean(freqSig_sorted(end));
%         trend = squeeze(TREND(sortIndx(end),:,:));
%         trend = trend.';
%         trend = trend(:);
%         trend = abs(trend);
%         trend = trend ./ trend(1);
%         trend = 20*log10(trend);
% 
%         trendn = squeeze(TRENDn(sortIndx(end),:,:));
%         trendn = trendn.';
%         trendn = trendn(:);
%         trendn = abs(trendn);
%         trendn = trendn ./ trendn(1);
%         trendn = 20*log10(trendn);
        
        [~,sortIndx] = sort(nanmean(abs(noiseFloor_sm)));
        noiseFloor_sorted = noiseFloor_sm(:,sortIndx);
        freqNf_sorted = frequency(sortIndx);
        nF = mean(freqNf_sorted(end));
        sigHat = abs(signal_sorted(:,end));
        %sigHat = mean(signal_sorted(:,24:26),2);
        %sigHat = mean(abs(signal_sorted),2);
        sigHatNF = mean(signal_sorted(:,1:3),2);
        %nfHat = mean(noiseFloor_sorted(:,24:26),2);
        nfHat = mean(abs(noiseFloor_sorted),2);
        nfHatNF = mean(noiseFloor_sorted(:,1:3),2);

        % [~,maxSigIndx] = max(nanmean(abs(signal_sm)));
        % [~,maxNfIndx] = max(nanmean(abs(noiseFloor_sm)));
        % sF = frequency(maxSigIndx);
        % nF = frequency(maxNfIndx);
        % [~,minSigIndx] = min(nanmean(abs(signal_sm)));
        % [~,minNfIndx] = min(nanmean(abs(noiseFloor_sm)));

        % find threshold by finding where largest response crosses the
        % noise, defined by the largest value in the first 20 samples
        cutoff = max(sigHat(1:20));
        indx = 80;
        done = 0;
        counter = 1;
        while done == 0
            q = sigHat(indx);
            if q > cutoff
                indx = indx - 1;
            else
                done = 1;
            end
            if indx < 1
                done = 1;
            else
                counter = counter + 1;
            end
        end
        thdOnIndx = indx;
        indx = 80;
        done = 0;
        counter = 1;
        while done == 0
            q = sigHat(indx);
            if q > cutoff
                indx = indx + 1;
            else
                done = 1;
            end
            if indx > 160
                done = 1;
            else
                counter = counter + 1;
            end
        end        
        thdOffIndx = indx;

        figure
        subplot(2,1,1)
        %plot(signal,'Color',[.7 .7 .7])
        hold on
        plot(signal_sm,'Color',[0 0 1])
        %plot(abs(signal_sm(:,maxSigIndx)),'g','LineWidth',2)
        %plot(abs(signal_sm(:,minSigIndx)),'k','LineWidth',1)
        plot(abs(sigHat),'g','LineWidth',2)
        plot(abs(sigHatNF),'k','LineWidth',1)
        title([num2str(sF),' Hz'])
        xlim([1 160])
        ymax = max(sigHat);
        ymax = ymax * 1.1;
        ymin = -ymax;
        ylim([ymin,ymax])
        line([thdOnIndx,thdOnIndx],[ymin,ymax],'Color',[1 0 0],'LineWidth',1)
        line([thdOffIndx,thdOffIndx],[ymin,ymax],'Color',[1 0 0],'LineWidth',1)

        subplot(2,1,2)
        
        %plot(noiseFloor,'Color',[.7 .7 .7])
        hold on
        plot(abs(noiseFloor_sm),'Color',[1 0 0])
        %plot(abs(noiseFloor_sm(:,maxNfIndx)),'g','LineWidth',2)
        %plot(abs(noiseFloor_sm(:,minNfIndx)),'k','LineWidth',1)
        plot(abs(nfHat),'g','LineWidth',2)
        plot(abs(nfHatNF),'k','LineWidth',1)        
        title([num2str(nF),' Hz'])
        xlim([1 160])

%         figure
%         Q = repmat(sigHat,1,15);
%         Q = Q(:);
%         plot(Q+trend)
%         hold on
%         plot(trend,'r')

%         keyboard

%         analysisPart1.Z = Z;
%         analysisPart1.Zn = Zn;
%         analysisPart1.TREND = TREND;
%         analysisPart1.TRENDn = TRENDn;
%         analysisPart1.nFreqs = nFreqs;
            
        analysisPart1.signal_sm = signal_sm;
        analysisPart1.sigHat = sigHat;
        analysisPart1.sigHatNF = sigHatNF;
        analysisPart1.sF = sF;
        analysisPart1.thdOnIndx = thdOnIndx;
        analysisPart1.thdOffIndx = thdOffIndx;
        analysisPart1.noiseFloor_sm = noiseFloor_sm;
        analysisPart1.nfHat = nfHat;
        analysisPart1.nfHatNF = nfHatNF;
        analysisPart1.nF = nF;
        % stuff we intend to use
        analysisPart1.Trend = Trend; % is in dB
        analysisPart1.Trendn = Trendn;
        analysisPart1.D1 = D1;
        analysisPart1.D2 = D2;
        analysisPart1.x = x;
        analysisPart1.t = t;
        analysisPart1.freq = freq;
        analysisPart1.Z = Z;
        analysisPart1.Z_sm = Z_sm;
        

        return

    end

end
function [] = plot1(analysisPart1_inc,analysisPart1_mem)

    z = analysisPart1_inc.Z;   
    TREND = analysisPart1_inc.TREND;
    nFreqs = analysisPart1_inc.nFreqs;

    kk=19;
    mr = analysisPart1_inc.modeR(:,kk);
    mi = analysisPart1_inc.modeI(:,kk);
    % ciUR = analysisPart1.hdi95UR(:,kk)-analysisPart1.modeR(:,kk) + analysisPart1.hdi95UR(1,kk);
    % ciUI = analysisPart1.hdi95UI(:,kk)-analysisPart1.modeI(:,kk) + analysisPart1.hdi95UI(1,kk);
    ciUR = analysisPart1_inc.hdi95UR(:,kk)-analysisPart1_inc.modeR(:,kk) + analysisPart1_inc.modeR(1,kk);
    ciUI = analysisPart1_inc.hdi95UI(:,kk)-analysisPart1_inc.modeI(:,kk) + analysisPart1_inc.modeI(1,kk);
    ciLR = analysisPart1_inc.hdi95LR(:,kk)-analysisPart1_inc.modeR(:,kk) + analysisPart1_inc.modeR(1,kk);
    ciLI = analysisPart1_inc.hdi95LI(:,kk)-analysisPart1_inc.modeI(:,kk) + analysisPart1_inc.modeI(1,kk);
    m = mr + 1i*mi;
    ciU = ciUR + 1i*ciUI;
    ciL = ciLR + 1i*ciLI;
    m = m ./ m(1);
    ciU = ciU ./ ciU(1);
    ciL = ciL ./ ciL(1);
    
    % nLvl = [nLvl;flipud(nLvl)];
    % nLvl = linspace(45,115,80)';
    x = (1:1:160)';
    
    plot(x,20*log10(abs(ciU)),'r')
    hold on
    plot(x,20*log10(abs(ciL)),'k')
    plot(x,20*log10(abs(m)),'b')
    xIndx = (1:18:160)';
    xticks(x(xIndx))
    xticklabels({'45', '61', '77','93', '109','106','90', '74', '58'})
    xlim([x(1),x(end)])
    [ymax,yMaxIndx] = max(20*log10(abs(m)));
    ymin = -0.5;
    line([80,80],[ymin,ymax],'LineStyle','-','Color',[0 0 0],'LineWidth',1)
    delay = yMaxIndx - 80;
    line([yMaxIndx,yMaxIndx],[ymin,ymax],'LineStyle','--','Color',[0 0 0],'LineWidth',1)
    
    % keyboard
    
    % 
    % plot(analysisPart1.modeR(:,1))
    % hold on
    % plot(analysisPart1.hdi95UR(:,1)-analysisPart1.modeR(:,1) + analysisPart1.hdi95UR(1,1))
    % plot(analysisPart1.hdi95UR(:,1),'k--')

    figure 01
    subplot(2,1,1)
    plot(abs(z))
    xlabel('Time (Click Number)')
    ylabel('Mag Change (linear)')
    hold on
    subplot(2,1,2)
    ylim([-5, 5]); 
    plot(20*log10(abs(z)))
    xlabel('Time (samples)')
    ylabel('Mag Change (dB)')

    clear trend
    for ii=1:nFreqs
        t = squeeze(TREND(ii,:,:)).';
        t = abs(t(:));
        t = t ./ t(1);
        trend(:,ii) = t;
    end
    figure 02
    plot(20*log10(trend))
    xlabel('Time (samples)')
    ylabel('Trend (dB)')


    figure
    plot(z)


    keyboard


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
        w(rejectIndx) = 0;
        
        smoothing = 0.00000000001; % smoothing factor--0 is a straight line, 1 is cubic spline
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


% OLD CODE ------------------------------
