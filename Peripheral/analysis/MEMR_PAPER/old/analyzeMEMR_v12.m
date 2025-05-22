function [MEMR] = analyzeMEMR_v12(pathName,fileNameL,fileNameR,subjectName,runNumber)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [MEMR] = analyzeMEMR_v9(pathName,fileNameL,fileNameR,subjectName,runName)
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
% Last Updated: October 3, 2023 -- ssg -- version 12. Still to do: Analyze
%               click peak to show that it is stable.
% Last Updated: October 20, 2023 -- shows that change in ZL is the same
% result as change in SPL.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
    % ---------------------------------------------------------------------
%     pathName = 'D:\MEMR_DATA\';
%     subjectName = 'MEM611\';
%     runNumber = ;
% 
%     fileNameL = 'Ch3_ER10xA_memr_0001.mat';
%     fileNameR = 'Ch4_ER10xB_memr_0001.mat';
%     
    pathName = 'D:\MEMR_DATA\';
    subjectName = 'MEM38\';
    runNumber = 2;

    fileNameL = 'Ch3_ER10xA_memr_0001.mat';
    fileNameR = 'Ch4_ER10xB_memr_0001.mat';
    

    doJitterFix = 1;
    % ---------------------------------------------------------------------

    d = dir([pathName,subjectName,'*run',num2str(runNumber)]);
    octStepSize = 1/2; % octave step size for analysis
    fmin = 125; % minimum frequency to examine
    fmax = 8000; % maximum frequency to examine
    % frequencies of interest:
    peakFreqs = round(2.^(log2(fmin):octStepSize:log2(fmax))'); % analysis center frequencies (Hz)

    % extract the raw recordings -----
    dummy = load([pathName,subjectName,d.name,'\',fileNameL]); % the clicks
    header = dummy.header;
    Clicks = dummy.data;
    clear dummy
    micCorrection = header.userInfo.micCorrectionL;
    Clicks = applyMicCorrection(Clicks,micCorrection);
    %[pl,Pl,phi,other,wf] = ARLas_convertPL(Clicks(:,1),header.userInfo.iscS1L);

    dummy = load([pathName,subjectName,d.name,'\',fileNameR]); % the noise
    headerNoise = dummy.header;
    Noise = dummy.data;
    clear dummy 
    micCorrection = header.userInfo.micCorrectionR;
    %Noise = applyMicCorrection(Noise,micCorrection);
    
    % Jitter Fix! ---------------------------------------------------------
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
    % ----> what is special about this???? 
    % Nothing, we just chose an exemplar for plotting


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

    clicks = Clicks(1:2000,:);
    cc = mean(clicks,2);

%     [pl,Pl,phi,other,wf] = ARLas_convertPL(cc,header.userInfo.iscS1L);
    % figure(100) 
    % plot(pl.pr)
    % 
    % figure(101)
    % plot(pl.spl,'k')
    % hold on
    % plot(pl.fpl,'g')
    % plot(pl.rpl,'r')
    % 
    % figure(102)
    % plot(wf.spl,'k')
    % hold on
    % plot(wf.fpl,'g')
    % plot(wf.rpl,'r')

    Clicks = ARLas_hpFilter(Clicks,fs,100);

    % % analyze ear canal stability
    % stabilityCheck = 1;
    % analysisPart1_inc = runme(Clicks,nSamples,nClicks,indx,time,fs,peakFreqs,noiseLvl,noiseRate,...
    %     subjectName,runNumber,elicitor,fileNameL,fileNameR,pathName,d.name,stabilityCheck,header);


    % analyze MEMR
    stabilityCheck = 0;
    analysisPart1_mem = runme(Clicks,nSamples,nClicks,indx,time,fs,peakFreqs,noiseLvl,noiseRate,...
        subjectName,runNumber,elicitor,fileNameL,fileNameR,pathName,d.name,stabilityCheck,header);

  keyboard
    plot1(analysisPart1_inc,analysisPart1_mem)


end
% INTERNAL FUNCTIONS ------------------------------------------------------
function [analysisPart1] = runme(Clicks,nSamples,nClicks,indx,time,fs,peakFreqs,noiseLvl,noiseRate,...
    subjectName,runName,elicitor,fileNameL,fileNameR,pathName,dname,stabilityCheck,header)
    %elicitorTime = time * 8; % time vector for each sweep
    %F = (0:1:fs-1)'*(fs/fs); % full frequency vector (Hz)
    % we're not examining the entire frequency range, so find the cut points
    %[~,fIndx1] = min(abs(min(peakFreqs)-F));
    %[~,fIndx2] = min(abs(max(peakFreqs)-F));
    %F = F(fIndx1:fIndx2);
    %T = T(fIndx1:fIndx2);

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

    % % % winN = 100;
    % % % wn = 20;
    % % % for ii=1:nClicks % loop across each click position (1-160)
    % % %     q = Clicks(indx(ii)-winN:indx(ii)+winN*2,:);
    % % %     if mod(size(q,1),2)~=0
    % % %         q = q(1:end-1,:);
    % % %     end
    % % % 
    % % %     h = hann(wn*2); % make a hann window to window the edges 
    % % %     h = h(1:wn);
    % % %     H = repmat(h,1,nSweeps);
    % % %     q(1:wn,:) = q(1:wn,:) .* H;
    % % %     q = flipud(q);
    % % %     q(1:wn,:) = q(1:wn,:) .* H;
    % % %     q = flipud(q);
    % % % 
    % % %     [pl,Pl,phi,other,wf] = ARLas_convertPL(mean(q,2),header.userInfo.iscS1L);
    % % %     Q(:,ii) = mean(q,2);
    % % %     RL(:,ii) = pl.pr; %$other.RL;
    % % % end


    winN = 14; % number of samples to reduce (before reflection comes back)
    if stabilityCheck == 1
        clickN = winN*2+1; % number of analysis samples in each (windowed click)
        h = hann(clickN); % make a hann window to window the edges 
        H = repmat(h,1,nSweeps);
        C = zeros(clickN,nSweeps,nClicks); % initialize matrix of windowed clicks (C for clicks)
        %fudgeFactor = 6;
        for ii=1:nClicks % loop across each click position (1-160)
            q = Clicks(indx(ii)-winN:indx(ii)+winN,:);
            %q = Clicks(indx(ii)-winN-fudgeFactor:indx(ii)+winN-fudgeFactor,:); % contains all the sweeps
            q = q .* H;
            C(:,:,ii) = q;
        end
    else
        clickN = 100 + winN; % number of analysis samples in each (windowed click)
        h = hann(winN*2); % make a hann window to window the edges 
        h = h(1:winN);
        H = repmat(h,1,nSweeps);
        C = zeros(clickN,nSweeps,nClicks); % initialize matrix of windowed clicks (C for clicks)
        for ii=1:nClicks % loop across each click position (1-160)
            q = Clicks(indx(ii):indx(ii)+clickN-1,:); % contains all the sweeps
            %[pl,Pl,phi,other,wf] = ARLas_convertPL(mean(q,2),header.userInfo.iscS1L);
            %zl = other.ZL;
            
            %S = fft(stimulus,nfft); % stimulus vector
            nfft = header.userInfo.iscS1R.nfft;
            R = fft(mean(q,2),nfft); % recordings vector
            %PL = R ./ S; % load pressure
            ff = (0:1:nfft-1)'*(fs/nfft);
            [~,indx1] = min(abs(ff-header.userInfo.iscS1R.fmin));
            [~,indx2] = min(abs(ff-header.userInfo.iscS1R.fmax));
            PL = R(indx1:indx2);
            %PL = PL(1:obj.fmaxIndx); % cut to highest calibrated frequency
            %PL = PL(obj.fminIndx:end);
            ZS = header.userInfo.iscS1R.ZS;
            PS = header.userInfo.iscS1R.PS;
            ZL(:,ii) = (ZS .* PL) ./ (PS - PL); % calculate load impedance
            
            
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

            C(:,:,ii) = q;
            % % %ccc(:,ii) = mean(q,2);
           
            %ZL(:,ii) = zl;

        end
        % % % pr = PR(:,1);
        % % % pr = repmat(pr,1,size(PR,2));
        % % % PR = PR ./ pr;
        % % % plot(pl.f,20*log10(PR))

    end

 fmin = 500; % minimum frequency to analyze (Hz)
 fmax = 3000; % maximym frequency to analyze (Hz)
 freq = header.userInfo.iscS1R.freq;
    [~,indx1] = min(abs(fmin-freq)); % index of minimum frequency
    [~,indx2] = min(abs(fmax-freq)); % index of maximum frequency   
    ZL = ZL(indx1:indx2,:);
    freq = freq(indx1:indx2,:);
    zl = ZL(:,1);
    zl = repmat(zl,1,size(ZL,2));
    Q = ZL ./ zl;
    figure
    plot(freq,abs(Q))
    xlim([0 5000])
    pause(0.01)

    % Convert clicks into the frequency domain ----------------------------
    %nfft = 960; % size of fft--chosen to give 100-Hz bin width
    frequency = (0:1:nfft-1)'*(fs/nfft); % frequency vector
    %fmin = 500; % minimum frequency to analyze (Hz)
    %fmax = 3000; % maximym frequency to analyze (Hz)
    [~,indx1] = min(abs(fmin-frequency)); % index of minimum frequency
    [~,indx2] = min(abs(fmax-frequency)); % index of maximum frequency
    nFreqs = indx2-indx1+1; % number of frequencies to analyze
    CFT = zeros(nFreqs,nSweeps,nClicks); % initialize matrix of Fourier transformed clicks
    TREND = zeros(nFreqs,nSweeps,nClicks); % initialize saved matrix of trend line
    W = zeros(nFreqs,nSweeps,nClicks); % initialize weighting matrix (downweight bad samples)
    for ii=1:nSweeps % loop across 15 sweeps
        cc = squeeze(C(:,ii,:));
        FT = fft(cc,nfft); % contains all sweeps at all frequencies for a given click position
        FT = FT(indx1:indx2,:); % keep only the frequencies of interest
        CFT(:,ii,:) = complex(FT); % Fourier transform of C, cut to frequenc6ies of interest
    end
    % calculate low-frequency drift across recordings
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
    % this just to double check we did it right
    % for ii=1:nFreqs
    %     for jj=1:nClicks
    %     cc = squeeze(CFT(jj,:,jj));
    %     ccend = squeeze(CFT(ii,:,jj));
    %     cc1 = squeeze(CFT(ii,:,1));
    %     plot(cc1)
    %     hold on
    %     plot(ccend)        
    % end

    % Bootstrap spline fit to data ----------------------------------------
    nIterations = 200; % number of bootstrap iterations
    %W = ones(size(CFT));
    x = (1:1:nClicks);
    mr_sm = zeros(nClicks,nIterations);
    mi_sm = zeros(nClicks,nIterations);
    for jj=1:nFreqs
        disp(['analyzing ',num2str(jj),' of ',num2str(nFreqs)])
        m = squeeze(CFT(jj,:,:)); % size m = 160 x 15
        w = squeeze(W(jj,:,:));
        Mr = real(m);
        Mi = imag(m);
        warning off
        for kk=1:nIterations 
            if mod(kk,50)== 0
                disp(['  analyzing ',num2str(kk),' of ',num2str(nIterations)])
            end
            indxkk = ceil(rand(nSweeps,1)*nSweeps);
            % take the weighted mean
            ww = w(indxkk,:);
            sumw = sum(ww,1);
            mr = sum(Mr(indxkk,:).*ww,1)./ sumw;
            mi = sum(Mi(indxkk,:).*ww,1)./ sumw;
            sm = 0.00001;
            ppr = csaps(x,mr,sm,[],sumw); % piecewise polynomial real coefficients
            ppi = csaps(x,mi,sm,[],sumw); % piecewise polynomial imaginary coefficients
            mr_sm(:,kk) = ppval(ppr,x);
            mi_sm(:,kk) = ppval(ppi,x);
        end
        for mm=1:size(mr_sm,1)
            [post] = calculatePosterior_v2(mr_sm(mm,:)');
            hdi95LR(mm,jj) = post.hdi95(1);
            hdi95UR(mm,jj) = post.hdi95(2);
            hdi68LR(mm,jj) = post.hdi68(1); % 68% hdi, corresponds to 1 sd, normally distributed
            hdi68UR(mm,jj) = post.hdi68(2);
            modeR(mm,jj) = post.mod;
        end
        for mm=1:size(mi_sm,1) % loop across 160 clicks
            [post] = calculatePosterior_v2(mi_sm(mm,:)');
            hdi95LI(mm,jj) = post.hdi95(1);
            hdi95UI(mm,jj) = post.hdi95(2);
            hdi68LI(mm,jj) = post.hdi68(1); % 68% hdi, corresponds to 1 sd, normally distributed
            hdi68UI(mm,jj) = post.hdi68(2);
            if jj==1
                 modeI(mm,jj) = 0;
            else
                modeI(mm,jj) = post.mod;
            end
        end
        warning on
        
        Z = mr_sm+1i*mi_sm;
        Z = Z / Z(1);
        z(:,jj) = mean(Z,2);

    end

figure
plot(abs(z(:,100)))
hold on
plot(abs(Q(100,:)))


    analysisPart1.z = z;
    analysisPart1.hdi95UR = hdi95UR;
    analysisPart1.hdi95LR = hdi95LR;
    analysisPart1.hdi68UR = hdi68UR;
    analysisPart1.hdi68LR = hdi68LR;
    analysisPart1.modeR = modeR;
    analysisPart1.hdi95UI = hdi95UI;
    analysisPart1.hdi95LI = hdi95LI;
    analysisPart1.hdi68UI = hdi68UI;
    analysisPart1.hdi68LI = hdi68LI;
    analysisPart1.modeI = modeI;
    analysisPart1.z = z;
    analysisPart1.TREND = TREND;
    analysisPart1.nFreqs = nFreqs;

    % figure
    % m = analysisPart1.modeR + 1i*analysisPart1.modeI;
    % ciL = analysisPart1.hdi95LI + 1i*analysisPart1.hdi95LI;
    % ciU = analysisPart1.hdi95UI + 1i*analysisPart1.hdi95UI;
    % 
    % for mm=1:size(m,2)
    %     m(:,mm) = m(:,mm) ./ m(1,mm);
    %     ciL(:,mm) = ciL(:,mm) ./ m(1,mm);
    %     ciU(:,mm) = ciU(:,mm) ./ m(1,mm);
    % end
    % m = abs(m);
    % ciL = abs(ciL);
    % ciU = abs(ciU);

end
function [] = plot1(analysisPart1_inc,analysisPart1_mem)

    z = analysisPart1_inc.z;   
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
function [] = runme2()
    % Normalize each sweep by the start of that sweep
    % use the average of baselineN clicks to do the normalization
    baselineN = 3; % how many clicks to include in baseline (from start and from finish)
    for jj=1:nFreqs
        m = M(:,:,jj);
        w = W(:,:,jj);
        baseline1 = m(1:baselineN,:);
        baseline2 = m(end-baselineN+1:end,:);

        w1 = w(1:baselineN,:);
        w2 = w(end-baselineN+1:end,:);
        N = sum(w(:));
        baseline = sum([baseline1(:);baseline2(:)],1) / N; 

        %baseline = mean([baseline;baseline2],1); % take from the end, too
        %mnorm = m ./ baseline; % normalized m 
        mnorm = m - baseline; % normalized m
        Mnorm(:,:,jj) = mnorm;

        mnorm2 = complex(abs(m)) ./ complex(mean(abs(baseline))); % use this to NOT normalize
        Mnorm2(:,:,jj) = mnorm2;
    end

    for kk=2:2:20
         m = squeeze(Mnorm(:,:,kk));
         figure(210)
         plot(mean(m,2))
         hold on
    
         figure(211)
         plot(abs(mean(m,2)))
         hold on
         pause(0.001)
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
    % estimate by taking the weighted least-squares fit, weighting by snr values
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
        m2 = Mnorm2(:,:,jj);
        type = 1; % this is for cumulative (complex) arc lenghts (CAL)
        [CAL(:,jj),CALn(:,jj),Swirls_sm(:,jj),Swirls_smn(:,jj),MAG(:,jj),MAG2(:,jj)] = crawl(timeChunk,m,type,timeShift,subjectName,runName,m2);

        %type = 2; % this is for magnitude only (MAG)
        %[MAG(:,jj),MAGn(:,jj)] = crawl(timeChunk,m,type,timeShift,subjectName,runName);
    end
    CAL = abs(CAL); % signal must be 0 or greater

    % analyze each frequency band: 
snrCut = 1; % cutoff for saying signal is present above the noise
    %   threshold up
    %   threshold down
    %   peak activation level
    %   amount of hysteresis (area under the curve, or threshold difference
    pw1 = 78; % window where the peak should be
    pw2 = 88;
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

        h = figure; % plotting thresholds and max
        h.Position = [360   198   551   420];
        sp1 = subplot(3,1,1);
        sp2 = subplot(3,1,2);
        sp3 = subplot(3,1,3);
        %sp1.Position = [0.1300    0.1468    0.6456    0.7782];
        %sp2.Position = [0.6443    0.5309    0.3444    0.3472];
        sp1.Position = [0.1300    0.1230    0.7611    0.6166];
        sp2.Position = [0.7253    0.4429    0.2361    0.2662];
        sp3.Position = [0.1306    0.8000    0.7557    0.1408];
        
        axes(sp1)
        mult = 10^(snrCut/20);
        plot(time2,CAL(:,kk),'r')
        hold on
        plot(time2,MAG(:,kk),'r:')
        plot(time2,CALn(:,kk)*mult,'k')
        ymax = max([max(CAL(:,kk)),max(CALn(:,kk)*mult)]);
        ymax = ymax + (ymax * 0.1);
        line([4,4],[1,ymax],'Color',[0 0 0],'LineStyle','-','LineWidth',0.5)
        if present(kk,1) == 1
            plot(thdUp_time(kk,1),CAL(thdUp_indx(kk,1),kk),'+b')
            plot(thdDn_time(kk,1),CAL(thdDn_indx(kk,1),kk),'+b')
            plot(peakTime(kk,1),peakAmp(kk,1),'*b')
            line([peakTime(kk,1),peakTime(kk,1)],[0,ymax],'Color',[0 0 1],'LineStyle','-','LineWidth',0.5)
            line([thdUp_time(kk,1),thdUp_time(kk,1)],[1,ymax],'Color',[0 0 1],'LineStyle','--','LineWidth',0.5)
            line([thdDn_time(kk,1),thdDn_time(kk,1)],[1,ymax],'Color',[0 0 1],'LineStyle','--','LineWidth',0.5)
        end
        xlabel('Time (s)')
        ylabel('Relative Change')
        set(sp1,'XTick',[1,2,3,4,5,6,7,8],'xTickLabel',{'1','2','3','4','5','6','7','8'},'FontSize',10)
        grid on
        legend('Arc Length','Magnitude','Noise Floor','Location','NorthWest')
        ylim([1,ymax])
        %line([4,4],[-20,20],'Color',[0 0 0],'LineStyle','-','LineWidth',0.5)

        axes(sp3)
        plot(elicitorTime,elicitor,'k')
        %plot(elicitorTime,20*log10(abs(elicitor)/.00002),'k')
        timeTargets = [0,1,2,3,4,5,6,7,8];
        %LVL = linspace(45,115,length(elicitorTime)/2)';
        %LVL = [LVL;flipud(LVL)];
        %for ii=1:length(timeTargets)
        %    [~,indx] = min(abs(timeTargets(ii)-elicitorTime));
        %    lvlTargets(ii,1) = LVL(indx);
        %end
        lvlLabels = {'45', '63', '80', '98', '115', '98', '80', '63', '45'};
        set(sp3,'XTick',timeTargets,'xTickLabel',lvlLabels,'FontSize',10)
        xlabel('Elicitor Levels (dB SPL)')
        ylim([-20 20])
        ylabel('Ampl (Pa)')
        title([subjectName(1:end-1),'     ',runName(1:end-1),'     ',num2str(peakFreqs(kk)),' kHz'])
        grid on
        line([4,4],[-20,20],'Color',[0 0 0],'LineStyle','-','LineWidth',0.5)
        if present(kk,1) == 1
            line([peakTime(kk,1),peakTime(kk,1)],[-20,20],'Color',[0 0 1],'LineStyle','-','LineWidth',0.5)
            line([thdUp_time(kk,1),thdUp_time(kk,1)],[-20,20],'Color',[0 0 1],'LineStyle','--','LineWidth',0.5)
            line([thdDn_time(kk,1),thdDn_time(kk,1)],[-20,20],'Color',[0 0 1],'LineStyle','--','LineWidth',0.5)
        end

        axes(sp2)
        qqq = Swirls_sm(:,kk); 
        qqqn = Swirls_smn(:,kk);
        %qqq = qqq + 1;
        %qqqn = qqqn + 1;
        plot(qqq,'r')
        hold on
        plot(qqqn,'k')
        xymax = max([max(abs(real(qqq))),max(abs(imag(qqq))),max(abs(real(qqqn))),max(abs(imag(qqqn)))]);
        xymax = xymax + (xymax*.05);
        if xymax < 1.1
            xymax = 1.1;
        end
        line([-1,1],[0,0],'Color',[0 0 0],'LineWidth',0.5,'LineStyle','-')
        line([0 0],[-1,1],'Color',[0 0 0],'LineWidth',0.5,'LineStyle','-')
        xlim([-xymax,xymax])
        ylim([-xymax,xymax])
        box on
        grid on
        axis square
        xlabel('Real')
        ylabel('Imag')
        theta = linspace(0,2*pi,1000); % Linearly-spaced vector 
        UC = cos(theta)+1i*sin(theta); % Create unit circle
        plot(UC,'k')

        %FN = [pathName,['ArcLength_',num2str(kk)]];
        PN = [pathName,subjectName,dname,'\',dname,'_analysis\'];
        FN = [PN,'ArcLength_',num2str(kk)];
        saveas(h,FN,'jpg')
        pause(0.01)
        close(h)
    end

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
    MEMR.time = time2; % time across each sweep (8 seconds)
    MEMR.noiseLvl = noiseLvl2; % noise level as a function of time
    MEMR.noiseLvlCorrected = noiseLvl2 - (timeShift * noiseRate); % noise level corrected for shift
    MEMR.Swirls = Swirls;
    MEMR.Swirls_sm = Swirls_sm;

    % create a measure of overall activation
    total = 0;
    dl = median(diff(MEMR.noiseLvl));
    for ii=1:length(MEMR.peakFreqs)
        if MEMR.present(ii) == 1
            new = sum(MEMR.CAL(MEMR.thdUp_indx(ii):MEMR.thdDn_indx(ii),ii)) * dl;
            total = total + new;
        end
    end
    totalActivation = total / length(MEMR.peakFreqs) / 2;
    % total = divide by number of frequencies divided by 2 (becase on and off)
    if isempty(totalActivation)
        totalActivation = 0;
    end
    MEMR.totalActivation = totalActivation;
    MEMR.thdUpEstimate = sum(MEMR.thdUp_lvl .* MEMR.snr) / sum(MEMR.snr);
    MEMR.thdDnEstimate = sum(MEMR.thdDn_lvl .* MEMR.snr) / sum(MEMR.snr);
    FN = [PN,'MEMRanalysis'];
    save(FN,'MEMR')

    h = figure; % all swirls --------------------------
    hold on
    LG = [];
    for ii=1:length(MEMR.peakFreqs)
        if MEMR.present(ii) == 1
            plot(Swirls_sm(:,ii))
            LG = [LG,{[num2str(peakFreqs(ii)),' Hz']}];
        end
    end
    xymax = max([max(abs(real(Swirls_sm(:))));max(abs(imag(Swirls_sm(:))))]);
    xymax = xymax + (xymax*.05);
    if xymax < 1.1
        xymax = 1.1;
    end
    line([-1,1],[0,0],'Color',[0 0 0],'LineWidth',1,'LineStyle','-')
    line([0 0],[-1,1],'Color',[0 0 0],'LineWidth',1,'LineStyle','-')
    xlim([-xymax,xymax])
    ylim([-xymax,xymax])
    box on
    %grid on
    axis square
    xlabel('Real')
    ylabel('Imag')
    theta = linspace(0,2*pi,1000); % Linearly-spaced vector 
    UC = cos(theta)+1i*sin(theta); % Create unit circle
    plot(UC,'k','LineWidth',1)
    plot(1,0,'k.','MarkerSize',10)
    if any(MEMR.present(ii) == 1)
        legend(LG)
    end
    title([subjectName(1:end-1),'     ',runName(1:end-1),'    TotalAct=',num2str(round(MEMR.totalActivation,1)),'    ThdUp=',num2str(round(MEMR.thdUpEstimate,1)),' dB SPL'])
    FN = [PN,'Swirls'];
    saveas(h,FN,'jpg')
    pause(0.01)
    close(h)
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
function [delta,nf,swirl,swirl_n,absLen,m2] = crawl(t,x,type,timeShift,subjectName,runName,m2)
    % calculate the amount of change by finding the length of the smoothed curve
    % t = timeChunk (1 x 160 vector)
    % x = complex vector of changes across time
    % arcLen_sm = the smoothed cumulative arc length
         medSm = 5; % number of samples (to both left and right) over which to smooth for median
        meanSm = 12; % same, for mean
    [xrs,xis,xrn,xin] = doSmoothing(x,type,medSm,meanSm);
    swirl = mean(xrs,2)+1i*mean(xis,2);
    swirl_n = mean(xrn,2)+1i*mean(xin,2);
    offset = mean(swirl(1:3));
    swirl = swirl - offset + 1;
    offset = mean(swirl_n(1:3));
    swirl_n = swirl_n - offset + 1;

    [xrs2,xis2,xrn2,xin2] = doSmoothing(m2,type,medSm,meanSm);
    m2 = mean(xrs2,2)+1;

    Q = [];
    Qn = [];

    dt = median(gradient(t));
    [~,midPoint] = min(abs(t-4));
    foldIndx = midPoint + round(timeShift/0.05);
    foldIndx2 = (80 +  80 - foldIndx) + 1;
    
    for iii=1:1 % no loop
        sr = mean(xrs,2)';
        si = mean(xis,2)';
        nr = mean(xrn,2)';
        ni = mean(xin,2)';
        nn = length(sr);

        % CALCULATE FOR SIGNAL --------------------------------------------
        % integrate from the left
        if type == 1
        baseR = mean(sr(1:3)); % real signal
        baseI = mean(si(1:3)); % imgaingary signal
        sr = sr - baseR;
        si = si - baseI;

        %absLen = abs(complex(sr+1i*si)+1);
        absLen = abs(sr+1i*si)+1;
            grad = sqrt((  (gradient(sr,t)).^2   + (gradient(si,t)).^2  )  );
            k = ones(size(grad));
            k(foldIndx:end) = -k(foldIndx:end);
            q1 = cumsum(grad.*k) * dt;
            sr = fliplr(sr);
            si = fliplr(si);
        baseR = mean(sr(1:3)); % real signal
        baseI = mean(si(1:3)); % imgaingary signal
        sr = sr - baseR;
        si = si - baseI;
            grad = sqrt((  (gradient(sr,t)).^2   + (gradient(si,t)).^2  )  );
            k = ones(size(grad));
            k(foldIndx2:end) = -k(foldIndx2:end);
            q2 = fliplr(cumsum(grad.*k) * dt);
            wf1 = linspace(1,0,nn); % weighting function
            wf2 = linspace(0,1,nn);
            arcLen = ((q1.*wf1)+(q2.*wf2)) ./ (wf1+wf2); % weighted sum 
            
        elseif type == 2
            arcLenL = abs(complex(sr+1i*si));
        elseif type == 3
            arcLenL = abs(complex((sr+1i*si)-1));
        end
        % intergrate from the right
        if type == 1
        else
            arcLen = arcLenL;
        end
        meanSm = 2; 
        arcLen = meanSmoother(arcLen',meanSm)';

        % CALCULATE FOR NOISE ---------------------------------------------
        % integrate from the left
        if type == 1
        baseR = mean(nr(1:3)); % real signal
        baseI = mean(ni(1:3)); % imgaingary signal
        nr = nr - baseR;
        ni = ni - baseI;            
            grad = sqrt((  (gradient(nr,t)).^2   + (gradient(ni,t)).^2  )  );
            k = ones(size(grad));
            k(foldIndx:end) = -k(foldIndx:end);
            q1 = cumsum(grad.*k) * dt;
            nr = fliplr(nr);
            ni = fliplr(ni);
        baseR = mean(nr(1:3)); % real signal
        baseI = mean(ni(1:3)); % imgaingary signal
        nr = nr - baseR;
        ni = ni - baseI;            
            grad = sqrt((  (gradient(nr,t)).^2   + (gradient(ni,t)).^2  )  );
            k = ones(size(grad));
            k(foldIndx2:end) = -k(foldIndx2:end);
            q2 = fliplr(cumsum(grad.*k) * dt);
            wf1 = linspace(1,0,nn); % weighting function
            wf2 = linspace(0,1,nn);
            arcLen_n = ((q1.*wf1)+(q2.*wf2)) ./ (wf1+wf2); % weighted sum 

        elseif type == 2
            arcLenLn = abs(complex(nr+1i*ni));
        elseif type == 3
            arcLenLn = abs(complex(nr+1i*ni)-1);
        end
        % flip over a the center point to make the descending I/O function
        if type == 1
        else
            arcLen_n = arcLenLn;
        end
        arcLen_n = meanSmoother(arcLen_n',meanSm)';

        Q = [Q;arcLen]; 
        Qn = [Qn;arcLen_n];
    end
    if type == 1
        %Delta = Qn';
        %Delta = (Q-Qn)';
        delta = mean(Q',2);
        nf = mean(Qn',2);

        %hat = max(delta(60:100));
        %hatn = max(nf(60:100));
        %targetLvl = hat - hatn;
        %k = targetLvl / hat;
        %delta = delta * k;
        %nf = nf * k;

        %hat = max(delta(60:100));
        %delta = delta - deltan;
        %tah = max(delta);
        %delta = (delta ./ tah) * hat;

        %nf1 = std(Q',[],2) / sqrt(size(Q',2));
        %nf2 = std(Qn',[],2) / sqrt(size(Qn',2));
        %nf = mean([nf1,nf2],2);
        %nf = (nf ./ tah) * hat;
    else
        Delta = Q';
        DeltaN = Qn';
        delta = mean(Delta,2);
        nf = mean(DeltaN,2);
    end
    delta = delta(:);
    nf = nf(:);
    absLen = absLen(:);

    meanSm = 20;
    nf = meanSmoother(nf,meanSm);
    meanSm = 8;
    absLen = meanSmoother(absLen,meanSm);

    hat = max(absLen(40:120));
    %hat = max(abs(m2(40:120)));
    tah = max(delta(40:120));
    delta = (delta ./ tah) * hat;
    nf = (nf ./ tah) * hat;

    % one more shift here
    hat = hat - 1;
    tah = max(delta(40:120));
    delta = (delta ./ tah) * hat;
    nf = (nf ./ tah) * hat;
    delta = delta + 1;
    nf = nf + 1;

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
function [delta,nf,swirl,swirl_n,absLen] = crawl2(t,x,type,timeShift,subjectName,runName)
    % calculate the amount of change by finding the length of the smoothed curve
    % t = timeChunk (1 x 160 vector)
    % x = complex vector of changes across time
    % arcLen_sm = the smoothed cumulative arc length
         medSm = 5; % number of samples (to both left and right) over which to smooth for median
        meanSm = 20; % same, for mean
    [xrs,xis,xrn,xin] = doSmoothing(x,type,medSm,meanSm);
    swirl = mean(xrs,2)+1i*mean(xis,2);
    swirl_n = mean(xrn,2)+1i*mean(xin,2);
    offset = mean(swirl(1:3));
    swirl = swirl - offset;
    offset = mean(swirl_n(1:3));
    swirl_n = swirl_n - offset;

    Q = [];
    Qn = [];

    dt = median(gradient(t));
    [~,midPoint] = min(abs(t-4));
    foldIndx = midPoint + round(timeShift/0.05);
    foldIndx2 = (80 +  80 - foldIndx) + 1;
    
    for iii=1:size(xrs,2) % loop over the 7 buffers
        sr = xrs(:,iii)';
        si = xis(:,iii)';
        nr = xrn(:,iii)';
        ni = xin(:,iii)';
        nn = length(sr);

        % CALCULATE FOR SIGNAL --------------------------------------------
        % integrate from the left
        if type == 1
        baseR = mean(sr(1:3)); % real signal
        baseI = mean(si(1:3)); % imgaingary signal
        sr = sr - baseR;
        si = si - baseI;

        absLen = abs(complex(sr+1i*si));

            grad = sqrt((  (gradient(sr,t)).^2   + (gradient(si,t)).^2  )  );
            k = ones(size(grad));
            k(foldIndx:end) = -k(foldIndx:end);
            q1 = cumsum(grad.*k) * dt;
            sr = fliplr(sr);
            si = fliplr(si);
        baseR = mean(sr(1:3)); % real signal
        baseI = mean(si(1:3)); % imgaingary signal
        sr = sr - baseR;
        si = si - baseI;
            grad = sqrt((  (gradient(sr,t)).^2   + (gradient(si,t)).^2  )  );
            k = ones(size(grad));
            k(foldIndx2:end) = -k(foldIndx2:end);
            q2 = fliplr(cumsum(grad.*k) * dt);
            wf1 = linspace(1,0,nn); % weighting function
            wf2 = linspace(0,1,nn);
            arcLen = ((q1.*wf1)+(q2.*wf2)) ./ (wf1+wf2); % weighted sum 
            
        elseif type == 2
            arcLenL = abs(complex(sr+1i*si));
        elseif type == 3
            arcLenL = abs(complex((sr+1i*si)-1));
        end
        % intergrate from the right
        if type == 1
        else
            arcLen = arcLenL;
        end
        meanSm = 2; 
        arcLen = meanSmoother(arcLen',meanSm)';

        % CALCULATE FOR NOISE ---------------------------------------------
        % integrate from the left
        if type == 1
        baseR = mean(nr(1:3)); % real signal
        baseI = mean(ni(1:3)); % imgaingary signal
        nr = nr - baseR;
        ni = ni - baseI;            
            grad = sqrt((  (gradient(nr,t)).^2   + (gradient(ni,t)).^2  )  );
            k = ones(size(grad));
            k(foldIndx:end) = -k(foldIndx:end);
            q1 = cumsum(grad.*k) * dt;
            nr = fliplr(nr);
            ni = fliplr(ni);
        baseR = mean(nr(1:3)); % real signal
        baseI = mean(ni(1:3)); % imgaingary signal
        nr = nr - baseR;
        ni = ni - baseI;            
            grad = sqrt((  (gradient(nr,t)).^2   + (gradient(ni,t)).^2  )  );
            k = ones(size(grad));
            k(foldIndx2:end) = -k(foldIndx2:end);
            q2 = fliplr(cumsum(grad.*k) * dt);
            wf1 = linspace(1,0,nn); % weighting function
            wf2 = linspace(0,1,nn);
            arcLen_n = ((q1.*wf1)+(q2.*wf2)) ./ (wf1+wf2); % weighted sum 

        elseif type == 2
            arcLenLn = abs(complex(nr+1i*ni));
        elseif type == 3
            arcLenLn = abs(complex(nr+1i*ni)-1);
        end
        % flip over a the center point to make the descending I/O function
        if type == 1
        else
            arcLen_n = arcLenLn;
        end
        arcLen_n = meanSmoother(arcLen_n',meanSm)';

        Q = [Q;arcLen]; 
        Qn = [Qn;arcLen_n];
    end
    if type == 1
        %Delta = Qn';
        %Delta = (Q-Qn)';
        delta = mean(Q',2);
        nf = mean(Qn',2);

        %hat = max(delta(60:100));
        %hatn = max(nf(60:100));
        %targetLvl = hat - hatn;
        %k = targetLvl / hat;
        %delta = delta * k;
        %nf = nf * k;

        %hat = max(delta(60:100));
        %delta = delta - deltan;
        %tah = max(delta);
        %delta = (delta ./ tah) * hat;

        %nf1 = std(Q',[],2) / sqrt(size(Q',2));
        %nf2 = std(Qn',[],2) / sqrt(size(Qn',2));
        %nf = mean([nf1,nf2],2);
        %nf = (nf ./ tah) * hat;
    else
        Delta = Q';
        DeltaN = Qn';
        delta = mean(Delta,2);
        nf = mean(DeltaN,2);
    end
    delta = delta(:);
    nf = nf(:);
    absLen = absLen(:);

    meanSm = 20;
    nf = meanSmoother(nf,meanSm);
    meanSm = 8;
    absLen = meanSmoother(absLen,meanSm);

    hat = max(absLen(40:120));
    tah = max(delta(40:120));
    delta = (delta ./ tah) * hat;
    nf = (nf ./ tah) * hat;    

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
function [post] = calculatePosterior_v2(chain)
    % calculate the posterior distribution using a kernel method
    dist = fitdist(chain,'kernel');
    minx = min(chain);
    maxx = max(chain);
    N = 1000;
    xx = linspace(minx,maxx,N);
    yy = dist.pdf(xx);
    
    alpha = 0.05;
    [hdi95,HDIdensity,inOut95] = getHDI(alpha,xx,yy);
    alpha = 0.68;
    [hdi68,HDIdensity,inOut68] = getHDI(alpha,xx,yy);
    med = dist.median;
    [~,medIndx] = min(abs(med-xx));
    medY = yy(medIndx);
    modY = max(yy);
    indxMod = find(yy==modY);
    mod = xx(indxMod);

    post.dist = dist; % distribution object
    post.hdi95 = hdi95; % hdi interval for 95%
    post.hdi68 = hdi68; % hid interval for 50%
    post.med = med; % median of distribution
    post.medY = medY; % height of median of distribution
    post.xx = xx; % x-axis of distribution
    post.yy = yy; % y-axis of distribution
    post.mod = mod; % mode of distribution
    post.modY = modY; % height of mode of distribution
    post.inOut95 = inOut95;
    post.inOut68 = inOut68;
end
function [HDI,HDIdensity,inOut] = getHDI(alpha,xx,yy)
    % Calculate highest density interval on alpha
    [p,I] = sort(yy,'ascend'); % sort densities in ascending order
    cutValue = sum(p)*alpha; % find the alpha% cut value
    cutIndx = min(find(cumsum(p)>cutValue)); % find the location of the cut value
    waterline = p(cutIndx); % this is the cutoff density
    [goodValues,goodIndx] = find(yy >= waterline); % locate all values > cut
    HDI = [xx(min(goodIndx));xx(max(goodIndx))]; % determine the interval
    HDIdensity = waterline;
    inOut = (yy >= waterline);
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
    if strcmp(subjectName,'MEM10\') & runNumber == 1
        badIndx = [];
    end
    if strcmp(subjectName,'MEM37\') & runNumber == 1
        badIndx = [];
    end
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
% OLD CODE ----------------------------------------------------------------


% Xk = Q;
% Xbar = (1/K) * sum(Xk,2); % signal is the coherent mean
% Xbar2 = abs(Xbar) .^2; % signal energy
% phase = angle(Xbar); % signal phase
%     XBAR = repmat(Xbar,1,K); % replicate to matrix (to avoid running a loop)
%     S2 = (1/(K-1)) * sum((Xk - XBAR) .* conj(Xk - XBAR),2); % variance (this is the HA noise floor)
%         % Se = sqrt((1/K) * S2); % standard error of mean (calculation not necessary)
%     Se2 = (1/K) * S2; % energy of the standard error
% 

     % figure
        % plot(Z,'Color',[.7 .7 .7])
        % hold on
        % plot(z,'b')
        % pause
        % 
        % 
        % plot(mymodeR + 1i*mymodeI,'b','LineWidth',1)
        % hold on
        % plot(hdi1R + 1i*hdi1I,'r:')
        % plot(hdi2R + 1i*hdi2I,'r:')
        % plot(hdi3R + 1i*hdi3I,'r--')
        % plot(hdi4R + 1i*hdi4I,'r--')
        
%         mnorm = m - baseline; % normalized m
%         Mnorm(:,:,jj) = mnorm;
% 
%         mnorm2 = complex(abs(m)) ./ complex(mean(abs(baseline))); % use this to NOT normalize
%         Mnorm2(:,:,jj) = mnorm2;

% baseline = Z(1,:);
% base = mean(baseline);
% 
% T = (t/160)*8;
% 
% plot(T,gradient(mean(abs(Z./base),2))./gradient(T))
% plot(T,mean(abs(Z./base),2))
% plot(T,20*log10(mean(abs(Z./base),2)))

%     % find the location of the desired analysis frequencies
%     for ii=1:length(peakFreqs)
%         [~,peakIndx(ii,1)] = min(abs(F-peakFreqs(ii)));
%     end


    %nFreqs = indx(ii)+nSamples - indx(ii)+1;
    % we will make a multidimensional matrix to store the analysis 
    % nClicks x nSweeps x nFreqs 
    % so each column is a single sweep, each row contains the click change for a single frequency
    % M is size 160 clicks x 15 sweeps x n frequencies
    % nFreqs = 194/2; %length(peakIndx); % number of frequencies to analyze
    % M = zeros(nClicks,nSweeps,nFreqs);  % unnormalized analysis
    % W = ones(nClicks,nSweeps,nFreqs);  % weights
    %Mnorm = zeros(nClicks,nSweeps,nFreqs); % normalized analysis
