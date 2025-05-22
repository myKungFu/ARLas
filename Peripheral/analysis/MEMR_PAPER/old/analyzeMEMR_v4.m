function [h10,h11] = analyzeMEMR_v4(Data,time,fs,indx)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [h10,h11] = analyzeMEMR_v4(Data,time,fs,indx)
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
%
% NOTES: Still to do: 
% choose peak frequencies by picking the larges values
%   over a range. Currently we just use fixed values, which are not always
%   the best choice (i.e., our current frequencies will not usually yield
%   the largest value). This won't match with the Mepani method.
%
% Add a data structre to save the analysis. Include template, templateSNR, etc.  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    

    time = time / time(end); % express time from 0 to 1

    % frequencies of interest in 1/2 octave steps
    peakFreqs = [125, 177, 250, 354, 500, 707, 1000, 1414, 2000, 2828]';

    % When considering which portions of each time window to analyze:
    % Within each click window, we want to include several round trip
    % travel times of the click in the canal, but no low-frequency OAEs. We
    % should consider only the first 2 ms of the click, which is 
    % round(0.002*fs) = 192 samples. 
    % We used a noise activator with a rise time of 4 seconds and a fall
    % time of 4 seconds. This gives 8 seconds / .05 = 160 clicks.
    nSamples = 193; % number of time samples in each click analysis window
    nClicks = length(indx); % number of total clicks to analyze


    % compare changes in clicks by making a "Template".
    % Get the template by averaging the first 3 clicks and last 3 clicks.
    % NOTE: The MEMR has hysteresis, so that you can't also average across
    % multiple sweeps! Sweeps can show drift. Also, you aren't guaranteed
    % to return to baseline by the end of the sweep, so you should probably
    % avoid just averaging the start and end also?
    % Each click window is 50 ms long, so the average of 3 is 150 ms.
    templateN = 3; % number of clicks (from both beginning and end) to include in the template
    for jj=1:templateN 
        chunk = Data(indx(jj):indx(jj)+nSamples,:);
        % artifact rejection:
        hat = max(chunk,[],1); % get the maximum value in each sweep
        doPlot = 0;
        [indxAR,nRejects] = AR(hat(:)','mild',doPlot);
        if ~isempty(indxAR) % if there are rejections to make
            chunk(:,indxAR) = [];
        end
        template(:,jj) = mean(chunk,2);
    end
    counter = 4;
    for jj=nClicks-2:nClicks % loop across click
        chunk = Data(indx(jj):indx(jj)+nSamples,:);
        % artifact rejection:
        hat = max(chunk,[],1); % get the maximum value in each sweep
        doPlot = 0;
        [indxAR,nRejects] = AR(hat(:)','mild',doPlot);
        if ~isempty(indxAR) % if there are rejections to make
            chunk(:,indxAR) = [];
        end
        template(:,counter) = mean(chunk,2);
        counter = counter + 1;
    end
    templateNoise = std(template,[],2)/sqrt(size(template,2));
    template = mean(template,2);
    % average SNR in the template window
    templateSNR = 20*log10(mean(abs(template)./templateNoise));

    % run the data analysis program
    [h10,h11] = runme(template,Data,nSamples,nClicks,indx,time,fs,peakFreqs);


end
% INTERNAL FUNCTIONS ------------------------------------------------------

function [h10,h11] = runme(template,Data,nSamples,nClicks,indx,time,fs,peakFreqs)

    T = fft(template,fs); % put template into the frequency domain
    % note that we are not scaling, because we will take a ratio later, so it doesn't matter
    F = (0:1:fs-1)'*(fs/fs); % full frequency vector (Hz)
    % we're not examining the entire frequency range, so find the cut points
    [~,fIndx1] = min(abs(min(peakFreqs)-F));
    [~,fIndx2] = min(abs(max(peakFreqs)-F));
    F = F(fIndx1:fIndx2);
    T = T(fIndx1:fIndx2);

    for jj=1:nClicks % loop across each click position in the sweep
        timeChunk(1,jj) = time(indx(jj)); % the temporal postion of each click position
    end
    
    % find the location of the desired analysis frequencies
    for ii=1:length(peakFreqs)
        [~,peakIndx(ii,1)] = min(abs(F-peakFreqs(ii)));
    end
    nFreqs = length(peakIndx); % number of frequencies to analyze
    nSweeps = size(Data,2); % number of sweeps
    
    % we will make a multidimensional matrix to store the analysis 
    % nClicks x nSweeps x nFreqs 
    % so each column is a single sweep, each row contains the click change for a single frequency
    % M is size 160 x 15 x 10
    M = zeros(nClicks,nSweeps,nFreqs);  % unnormalized analysis
    Mnorm = zeros(nClicks,nSweeps,nFreqs); % normalized analysis

    for ii=1:nClicks % loop across each click position
        Chunk = Data(indx(ii):indx(ii)+nSamples,:); % contains all the sweeps
        FT = fft(Chunk,fs); % contains all sweeps at all frequencies for a given click position
        FT = FT(fIndx1:fIndx2,:); % cut down to size
        %R = FT ./ T; % take the ratio of the clicks relative to the template
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
        baseline = mean(baseline,1);
        %mnorm = m ./ baseline; % normalized m 
        mnorm = m;
        Mnorm(:,:,jj) = mnorm;
    end

    % analyze on a frequency-by-frquency basis
    % get the cumulative arc lengths
    for jj=1:nFreqs
jj=9;        
        m = Mnorm(:,:,jj);
        %m = ARfix(m);
        [CAL(:,ii),CALn(:,ii)] = crawl(timeChunk,m);
        
        
        for ii=1:size(m,2)
            mm = m(:,ii);
            [CAL(:,ii),CALn(:,ii)] = crawl(timeChunk,mm);
        end
        CALd = CAL - CALn;
        cald = mean(CALd,2);
        %caln = mean(CALn,2);
        cal_nf = std(CAL,[],2) / sqrt(size(CAL,2));
        cald_nf = std(CALd,[],2) / sqrt(size(CALd,2));
        caln_nf = std(CALn,[],2) / sqrt(size(CALn,2));
        

        %plot(caln+2*caln_nf,'b')
        %hold on
        %plot(caln-2*caln_nf,'b')
        %plot(caln,'b--')
        %plot(cal,'r')
        plot(cald,'b')
        hold on
        plot(cal_nf*2,'k--')
        %plot(caln,'r')
        plot(caln_nf,'k')
        plot(caln_nf*2,'g')
        hold off

%         plot(CAL,'Color',[.7 .7 .7])
%         hold on
%         plot(cal,'r','LineWidth',2)
%         plot(cal_nf,'k')
        keyboard

    end


    for jj=1:nClicks % loop across each click position in the sweep
        Chunk = Data(indx(jj):indx(jj)+nSamples,:); % current chunk to analyze (contains the click)
        FT = fft(Chunk,fs);
        %S(:,jj) = fft(mean(chunk,2),fs);
    end
    
    Template = [S(:,1:3),S(:,nClicks-2:nClicks)]; % template in the frequency domain (template = first 3 and last 3 clicks)
    Template = mean(Template,2);
    for jj=1:length(indx) 
        Delta(:,jj) = S(:,jj)./Template;
    end
    freq = (0:1:size(chunk,1)-1)*(fs/size(chunk,1));




fff = (1:1:fs)-1'; % full frequency vector
% find the locatino of the desired center frequencies
for ii=1:length(peakFreqs)
    [~,peakIndx(ii)] = min(abs(fff-peakFreqs(ii)));
end
%peakIndx = [349,933,1619,3312, 4039,7294,14109];    

%-------- keyboard
% for jj=1:length(peakIndx)-1
%         D = complex(Delta(peakIndx(jj):peakIndx(jj+1),:));
%         f = fff(peakIndx(jj):peakIndx(jj+1));
%         [~,fhatIndx] = max(mean(abs(D-1),2));
%         peakFreqs2(jj) = f(fhatIndx);
% end





%--------------------------------------------------------------------------
end

% INTERNAL FUNCTIONS ------------------------------------------------------
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
function [m] = ARfix(m)
    % fix artifact problems
    % check for and remove single-sample artifacts
    mr = real(m);
    mi = imag(m);
    [rows,cols] = size(mr);
    mr = mr(:).';
    mi = mi(:).';
    nn = length(mr);
    doPlot = 0;
    [indxAR,nRejects] = AR(mr','severe',doPlot);
    if ~isempty(indxAR) % if there are rejections to make
        for rindx = 1:nRejects
            if indxAR(rindx) == 1
                mr(indxAR) = mean([mr(indxAR+1),mr(indxAR+2)]);
            elseif indxAR(rindx) == nn
                mr(indxAR) = mean([mr(indxAR-1),mr(indxAR-2)]);
            else
                mr(indxAR) = mean([mr(indxAR-1),mr(indxAR+1)]);
            end
        end
    end
    [indxAR,nRejects] = AR(mi','mild',doPlot);
    if ~isempty(indxAR) % if there are rejections to make
        for rindx = 1:nRejects
            if indxAR(rindx) == 1
                mi(indxAR) = mean([mi(indxAR+1),mi(indxAR+2)]);
            elseif indxAR(rindx) == nn
                mi(indxAR) = mean([mi(indxAR-1),mi(indxAR-2)]);
            else
                mi(indxAR) = mean([mi(indxAR-1),mi(indxAR+1)]);
            end
        end
    end
    m = mr + 1i*mi;
    m = reshape(m,rows,cols);
end
function [arcLen_sm,arcLen_n] = crawl(t,x)


%  [frequency,signal,noiseFloor] = ARLas_fda(x,96000,0.00002);
% plot(frequency,signal)
% hold on
% plot(frequency,noiseFloor,'k')
% plot(frequency,noiseFloor+6,'k')

    % calculate the amount of change by finding the length of the smoothed curve
    % t = timeChunk (1 x 160 vector)
    % x = complex vector of changes across time
    % arcLen_sm = the smoothed cumulative arc length
    dt = median(gradient(t));
    %x = x-1;
    xr = real(x);
    xi = imag(x);
    %xr = medianSmoother(xr,2);
    %xi = medianSmoother(xi,2);
    %smoothing = 0.999; % 0.999999;
    %ppr = csaps(t,xr,smoothing);
    %ppi = csaps(t,xi,smoothing);
    %xr_sm = ppval(ppr,t);
    %xi_sm = ppval(ppi,t);

[rows,cols] = size(xr);
xr = xr(:);
xi = xi(:);
xr = xr(:);
xi = xi(:);

    % artifact rejection
    multiplier = 4;
    xr = newAR(xr',multiplier)';
    xi = newAR(xi',multiplier)';

    xr_sm = medianSmoother(xr,5);
    xr_sm = meanSmoother(xr_sm,2);
    xr_n = xr - xr_sm;
    xi_sm = medianSmoother(xi,5);
    xi_sm = meanSmoother(xi_sm,2);
    xi_n = xi - xi_sm;
    xr_sm = xr_sm';
    xr_n = xr_n';
    xi_sm = xi_sm';
    xi_n = xi_n';
   %xr_sm = xr';
   %xi_sm = xi';
    % so we're comparing "raw" xr (which is called xr_sm)
    % to the noise with the smoothing subtracted off (which is called xr_n)
    % artifact rejection
    multiplier = 1.5;
    xr_n = newAR(xr_n,multiplier);
    xi_n = newAR(xi_n,multiplier);


% plot(xr,'b')
% hold on
% plot(xr_sm,'r')
% plot(xr_n,'k')
% plot(xi,'c')
% plot(xi_sm,'m')
% plot(xi_n,'k--')

% % xr_n = xr.' - xr_sm;
% % xi_n = xi.' - xi_sm;
% % xr_sm = xr.';
% % xi_sm = xi.';

%     % to estimate the noise, disrupt the time dependence in this time
%     % series (i.e., randomize the vector)
%     n = length(xr);
%     indxR = rand(n,1);
%     [~,indx] = sort(indxR);
%     noise = x(indx);
%     %nr = real(noise);
%     %ni = imag(noise);
%     % this doesn't really work:
%     %nr = xr.' - xr_sm;
%     %ni = xi.' - xi_sm;
%     ppr = csaps(t,nr,smoothing);
%     ppi = csaps(t,ni,smoothing);
%     nr_sm = ppval(ppr,t);
%     ni_sm = ppval(ppi,t);


    %int = cumsum(( (gradient(xr,t))   + (gradient(xi,t))  )) * dt;
    %int_smL = cumsum(( (gradient(xr_sm,t))   + (gradient(xi_sm,t))  )) * dt;
    %int_smR = fliplr(cumsum(fliplr(( (gradient(xr_sm,t))   + (gradient(xi_sm,t))  ))) * dt);

Q = [];
Qn = [];
pdr = makedist('normal','mu',0,'sigma',std(xr_n));
pdi = makedist('normal','mu',0,'sigma',std(xi_n));

xr_sm = reshape(xr_sm,rows,cols);
xi_sm = reshape(xi_sm,rows,cols);


for iii=1:100
    % to better estimate the noise, disrupt the time dependence in this time
    % series (i.e., randomize the vector)
    nnn = length(xr_n);
%     indxR = rand(n,1);
%     [~,indx] = sort(indxR);
%     xr_n = xr_n(indx);
%     indxR = rand(n,1);
%     [~,indx] = sort(indxR);
%     xi_n = xi_n(indx);
    
    qr_n = pdr.random(1,nnn);
    qi_n = pdi.random(1,nnn);

    qr_n = medianSmoother(qr_n',5);
    qr_n = meanSmoother(qr_n,2)';
    qi_n = medianSmoother(qi_n',5);
    qi_n = meanSmoother(qi_n,2)';
    qr_n = reshape(qr_n,rows,cols);
    qi_n = reshape(qi_n,rows,cols);

    
    % integrate from the left
        %% 
        arcLen_smL = cumsum(sqrt((  (gradient(xr_sm,t)).^2   + (gradient(xi_sm,t)).^2  )  )) * dt;
    % intergrate from the right    
        arcLen_smR = fliplr(cumsum(fliplr(sqrt((  (gradient(xr_sm,t)).^2   + (gradient(xi_sm,t)).^2  )  ))) * dt);
    % find where the two functions cross: this is the peak
        indxLT = find(arcLen_smL <= arcLen_smR);
        %[~,peakIndx2] = max(indxLT);
        peakIndx2 = 82;
        peak = arcLen_smL(peakIndx2);
    % flip the descending side
        arcLen_sm = [arcLen_smL(1:peakIndx2),arcLen_smR(peakIndx2+1:end)];
        %[~,peakIndx3] = max(arcLen_sm);
        peakIndx3 = 81;
        arcLen_sm(peakIndx3) = mean([arcLen_sm(peakIndx3-1),arcLen_sm(peakIndx3),arcLen_sm(peakIndx2+1)]);
    Q = [Q;arcLen_sm];


    % now same thing for the noise
    % integrate from the left
        arcLen_smL = cumsum(sqrt((  (gradient(qr_n,t)).^2   + (gradient(qi_n,t)).^2  )  )) * dt;
        %arcLen_smL = cumsum(((  (gradient(xr.'-xr_sm,t))   + (gradient(xi.'-xi_sm,t))  )  )) * dt;
    % intergrate from the right
        arcLen_smR = fliplr(cumsum(fliplr(sqrt(( (gradient(qr_n,t)).^2   + (gradient(qi_n,t)).^2  )  ))) * dt);
        %arcLen_smR = fliplr(cumsum(fliplr(((  (gradient(xr.'-xr_sm,t))   + (gradient(xi.'-xi_sm,t))  )  ))) * dt);
    % find where the two functions cross: this is the peak
        %indxLT = find(arcLen_smL <= arcLen_smR);
    % flip the descending side
        arcLen_n = [arcLen_smL(1:peakIndx2),arcLen_smR(peakIndx2+1:end)];
        %arcLen_n = arcLen_smL;
        arcLen_n(peakIndx3) = mean([arcLen_n(peakIndx3-1),arcLen_n(peakIndx3),arcLen_n(peakIndx2+1)]);
    Qn = [Qn;arcLen_n];
end


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

arcLen_sm = mean(Q,1);
arcLen_n = mean(Qn,1);

end


% OLD CODE ----------------------------------------------------------------
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
