function [] = batchMEMR_mepani3()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% batchMEMR_mepani3
% 
% Batch analysis of MEMR data using analyzeMEMR_vmepani.m.
%
% Author: Shawn Goodman & Ehsan Khalili
% Date: January 26, 2024
% Last Updated: January 26, 2024
% Last Updated: February 16, 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    parentDrive = 'D'; % use 'C' or 'D'
    dataPathName = [parentDrive,':\MEMR_MepaniClicks\']; % location of raw data
    fileNameL = 'Ch3_ER10xA_memr_0001.mat';
    fileNameR = 'Ch4_ER10xB_memr_0001.mat';
    savePath = [parentDrive,':\MEMR_MepaniAnalysis\']; % where to save analyzed data

    d = dir(dataPathName); % read in file directory information
    nSubjects = size(d,1)-2; % number of subject folders
    counter = 1;
    
    for ii=1:nSubjects
        disp(['Analyzing subject ',num2str(ii),' of ',num2str(nSubjects)])
        dummy = load([dataPathName,d(ii+2).name]); % get clicks
        Clicks2 = dummy.Clicks2; 
        Clicks1 = dummy.Clicks1;
        Clicks1 = Clicks1(1:999,:,:);        
        clear dummy

        % hat = 482; % where the peak of the click should be
        % figure
        % hold on

        doPlot = 0;
        [Clicks1] = fixJitter(Clicks1,doPlot);
        [Clicks2] = fixJitter(Clicks2,doPlot);

        if doPlot == 1
            pause
            c
        end

        [C_mem1,C_inc1] = analyzeMEMR_vmepani3(Clicks1,counter);
        [C_mem2,C_inc2] = analyzeMEMR_vmepani3(Clicks2,counter);
        %[C_memD,C_incD] = analyzeMEMR_vmepani3(Clicks2-Clicks1,counter);
        


        C_mem.d = C_mem2.d./C_mem1.d;
        C_mem.d1 = C_mem2.d1./C_mem1.d1;
        C_mem.d2 = C_mem2.d2./C_mem1.d2;
        C_mem.z_sm = C_mem2.Z_sm./C_mem1.Z_sm;
        C_mem.Z = C_mem2.Z./C_mem1.Z;
        C_mem.D = C_mem2.D./C_mem1.D;
        C_mem.D1 = C_mem2.D1./C_mem1.D1;
        C_mem.D2 = C_mem2.D2./C_mem1.D2;
        C_mem.x = C_mem1.x;
        C_mem.t = C_mem1.t;
        C_mem.freq = C_mem1.freq;
        C_mem.trend = C_mem2.trend./C_mem1.trend;
        
        C_inc.d = C_inc2.d./C_inc1.d;
        C_inc.d1 = C_inc2.d1./C_inc1.d1;
        C_inc.d2 = C_inc2.d2./C_inc1.d2;
        C_inc.z_sm = C_inc2.Z_sm./C_inc1.Z_sm;
        C_inc.Z = C_inc2.Z./C_inc1.Z;
        C_inc.D = C_inc2.D./C_inc1.D;
        C_inc.D1 = C_inc2.D1./C_inc1.D1;
        C_inc.D2 = C_inc2.D2./C_inc1.D2;
        C_inc.x = C_inc1.x;
        C_inc.t = C_inc1.t;
        C_inc.freq = C_inc1.freq;
        C_inc.trend = C_inc2.trend./C_inc1.trend;


        counter = counter + 1;


    % 
    % figure
    % phi = linspace(0,2*pi,2000)';
    % cc = exp(1i*phi);
    % plot(cc)
    % hold on
    % plot(C_mem.z_sm)
    % ylim([-.1 .1])
    % xlim([-.1 .1]+1)

%    pause
    c

        % c
        % figure
        % plot(abs(C_mem1.Z_sm))
        % hold on
        % plot(abs(C_mem1.d2),'b')
        % plot(abs(C_mem1.d1),'r')
        % 
        % figure
        % plot(abs(C_mem2.Z_sm))
        % hold on
        % plot(abs(C_mem2.d2),'b')
        % plot(abs(C_mem2.d1),'r')

        % figure
        % plot(20*log10(C_mem2.d1./C_mem1.d1))
        % pause



        %hold on
        %subplot(2,1,2)
        %plot(20*log10(C_memD.d1),'b')
        %title('Difference')
        %pause

        %----------------------------------
        Z = C_mem2.Z_sm ./ C_mem1.Z_sm;
        D1 = abs(Z-1)+1;
        d1 = mean(D1,2);
        % -----------------------
        C_mem.Z = Z;
        C_mem.D1 = D1;
        C_mem.d1 = d1;

        Z = C_inc2.Z_sm ./ C_inc1.Z_sm;
        D1 = abs(Z-1)+1;
        d1 = mean(D1,2);
        C_inc.Z = Z;
        C_inc.D1 = D1;
        C_inc.d1 = d1;

    % Hack job -- change names here    
    MEMR_mem = C_mem;
    MEMR_inc = C_inc;
    MEMR_mem.timeTrend = (0:288-1)';
    subjectName = d(ii+2).name;
    subjectName = subjectName(1:5);

    n = length(MEMR_mem.D1); % number of clicks
    MEMR_mem.t = linspace(0,8,n)'; % x-axis for smoothing (click number)
    MEMR_inc.t = MEMR_mem.t;


    [d1_sm,d1] = secondSmooth(MEMR_mem.t,MEMR_mem.d1);
    MEMR_mem.d1_sm = d1_sm;
    MEMR_mem.d1 = d1;

    % get measures
    [MEMR_mem] = extractMeasures(MEMR_mem);


% -------------------------------------
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
        xlabel('Time (s)','FontSize',11)
        ylabel('Change (dB)','FontSize',11)
        title([subjectName])
        grid on
      axes(sp2)
        ax = gca;
        ax.FontSize = 11; 
        plot(MEMR_mem.t,20*log10(MEMR_mem.d1),'Color',[.7 .7 .7])
        hold on
        plot(MEMR_mem.t,20*log10(MEMR_mem.d1_sm),'k')
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

        % plot(MEMR_mem.t,(MEMR_mem.elicitorLevel)/max(MEMR_mem.elicitorLevel)*ymax)


       
    catch ME
        keyboard
    end


% ----------------------------------


        fileName = d(ii+2).name;
        fileName = fileName(1:end-4);
        save([savePath,fileName,'_Analyzed'],'C_mem','C_inc','C_mem1','C_inc1','C_mem2','C_inc2','MEMR_mem')

        % h = figure(10);
        % plot(20*log10(D1),'r')
        % title(fileName)

        saveas(h1,[savePath,fileName,'.bmp'])
        pause(0.01)
        c

    end

end

% INTERNAL FUNCTIONS ------------------------------------------------------
function [d1_sm,d1] = secondSmooth(t1,d1)
    n = length(t1);
    xxx = (1:1:n*3)'-1;
    x = xxx(1:n);
    
    % smooth the response (effectively lowpass filter)
    % this assumes knot conditions at the ends. Here, it isn't
    % always true, so subtract out the linear trend prior to
    % fitting to make it true. Then add the line back in.
    n = length(d1);
    start = d1(1);
    finish = d1(end);
    dx = length(d1)-1;
    dy = finish - start;
    slope = dy/dx;
    lineR = ((0:1:n-1)'*slope)+start;
    d1 = d1 - lineR;
    
    w = ones(size(d1));

    mrrr = [d1;d1;d1];
    www = [w;w;w];
    sm = 0.5; %0.0001; % was 0.00001;  smoothing factor (smaller numbers are more smooth)
    ppr = csaps(xxx,mrrr,sm,[],www); % piecewise polynomial real coefficients
    d1_sm = ppval(ppr,x); % evaluate only at the original x values
    %d1 = d1 + lineR;
    %d1_sm = d1_sm + lineR;
    d1 = d1 + 1;
    d1_sm = d1_sm + 1;

    % plot(d1)
    % hold on
    % plot(d1_sm)
    % keyboard

end
function [Clicks1] = fixJitter(Clicks1,doPlot)
        clicks = [];
        [X,Y,Z] = size(Clicks1);
        for jj = 1:24 
            clicks = [clicks,squeeze(Clicks1(:,:,jj))];
        end
        template = median(clicks,2);
        q = xcorr(template,template);
        hat = size(clicks,1);
        for jj=1:size(clicks,2)
            qq = xcorr(template,clicks(:,jj));
            [~,indx(jj,1)] = max(qq);
        end
        offset = indx - hat;
        offset = -offset;
        for kk=1:size(clicks,2)
            if abs(offset(kk)) > 1
                dummy = clicks(:,kk);
                zpad = zeros(abs(offset(kk)),1);
                if offset(kk) >0
                    zpad = zeros(abs(offset(kk)),1);
                    dummy = [dummy(offset(kk):end-1);zpad];
                else % negative offset
                    dummy = [zpad;dummy(1:end+offset(kk))];
                end
                clicks(:,kk) = dummy;
            end
        end
        template2 = mean(clicks,2);
        hat2 = 482; % where the peak of the click should be
        [~,indx2] = max(template2);
        offset = indx2-hat2;
        if abs(offset) > 0
            Zpad = zeros(abs(offset),size(clicks,2));
            if offset >0
                clicks = [clicks(offset+1:end,:);Zpad];
            else % negative offset
                clicks = [Zpad;clicks(1:end+offset,:)];
            end
        end
        % double check
        template2 = mean(clicks,2);
        hat2 = 482; % where the peak of the click should be
        [~,indx2] = max(template2);
        offset = indx2-hat2;
        if abs(offset) > 0
            disp('hit!')
            keyboard
        end


        clicks = ARLas_hpFilter(clicks,96000,500);
        if doPlot == 1
            figure(10)
            plot(clicks)
            hold on
        end
           

        Clicks1 = reshape(clicks,X,Y,Z);
end
function [MEMR_mem, thdAmp] = extractMeasures(MEMR_mem)
    % EXTRACT THE NEEDED MEASURES ----------------------------------
    %peak delay -------------------
    %There is no peak delay in the Mepani method!
        d1 = MEMR_mem.d1_sm.';
        sm = 1; % smoothing factor (smaller numbers are more smooth)
        n = length(d1); % number of clicks
        x = linspace(0,8,n)'; % x-axis for smoothing (click number)
        w = ones(size(x)); % weighting factor
    
        pp = csaps(x,d1,sm,[],w); % piecewise polynomial object
        xx = linspace(0,8,1000)';
        yy = ppval(pp,xx);
    
        dfdx = fnder(pp); % take derivative and solve for zero slope
        peakXX = fnzeros(dfdx); % peak location in seconds
        peakXX = peakXX(1,:);
    
        try
            peakYY = ppval(pp,peakXX);
            [~,peakIndx] = max(peakYY); 
            peakX = peakXX(peakIndx);
            peakY = peakYY(peakIndx);
        catch
            keyboard
        end
    
        % thresholds ------------------
        %   use the Q concept: reduce 12 from peak
        % targetLevelNoise = (90:-5:35); % this will be the descending run
        % targetLevelNoise = targetLevelNoise(:)';
        % targetLevelNoise = [targetLevelNoise,(fliplr(targetLevelNoise)+2.5)]; % add an ascending run
        % targetLevelNoise = targetLevelNoise(:);
    
    
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
        catch ME
            thdOnsetTime = NaN;
            thdOffsetTime = NaN;
        end
        % convert threshold times to thresholds re: nominal elicitor level
%        noiseTarget = MEMR_mem.noiseTarget;
        % noise went from 37.5 to 92.5 in 5 dB steps. Then back down
        % from 90 to 35 in - 5 dB Steps. 
        % Here, we are pretending that these changes occurred over 4
        % seconds:
        % 92.5-37.5 = 55 dB range in 4 seconds = (13.75 dB/s)
        % y = mx + b
        % stimLevel = 13.75x + 37.5 for the ascending run
        % stimLevel = -13.75x + 90 for the ascending run
        % but also account for reflex delay
        slope = 13.75;
        delay = 0;
        b1 = 37.5;
        b2 = 90; 
        try
            ThdLvl(1) = slope*(Thd(1)-delay) + b1;
            thdOnsetLvl = ThdLvl(1); % onset threshold re: stim level
            ThdLvl(2) = -slope*(Thd(2)-4-delay) + b2;
            thdOffsetLvl = ThdLvl(2); % offset threshold re: stim level
            elicitorLevel1 = slope*(x-delay) + b1;
            elicitorLevel2 = -slope*(x-4-delay) + b2;
            elicitorLevel = [elicitorLevel1(1:12);elicitorLevel2(13:end)];
        catch ME
            ThdLvl(1) = NaN;
            thdOnsetLvl = NaN; % onset threshold re: stim level
            ThdLvl(2) = NaN;
            thdOffsetLvl = NaN; % offset threshold re: stim level
            elicitorLevel = 17.5*(x-delay) + 45;
        end
        % go back and get threshold amplitudes from non-scaled d1
        sm = 1; % smoothing factor (smaller numbers are more smooth)
        %n = 24; % number of clicks
        x = linspace(0,8,n)'; % x-axis for smoothing (click number)
        w = ones(size(x)); % weighting factor
        pp = csaps(x,d1,sm,[],w); % piecewise polynomial object
        thdAmp = ppval(pp,thdOnsetTime);
    
        % hysteresis ---------------------
        try
            pp = csaps(x,d1,sm,[],w); % piecewise polynomial object
            hh = fnint(pp); % integrate the smoothed spline
            A = ppval(hh,Thd(1));
            B = ppval(hh,peakX);
            C = ppval(hh,Thd(2));
            aucLeft = B-A; % area under the curve left
            aucRight = C-B; % area under the curve right
            hyst = aucRight / aucLeft; % hysteresis as a ratio of area under the curves
            hysteresis = hyst;
        catch
            hysteresis = NaN;
        end
        % Calculate the slopes ------------------
        try
            peak_index = round(peakX /(1/3));
            index_xa = round(thdOnsetTime / (1/3));
            index_xb = round(thdOffsetTime / (1/3));
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
    MEMR_mem.peakTime = peakX;
    MEMR_mem.peakAmp = peakY;
    MEMR_mem.delay = delay;
    MEMR_mem.thdOnsetTime = thdOnsetTime;
    MEMR_mem.thdOffsetTime = thdOffsetTime;
    MEMR_mem.thdOnsetLvl = thdOnsetLvl;
    MEMR_mem.thdOffsetLvl = thdOffsetLvl;
    MEMR_mem.hysteresis = hysteresis;
    MEMR_mem.slopeUp = slopeUp;
    MEMR_mem.slopeDn = slopeDn;
    MEMR_mem.thd = thd;
    MEMR_mem.thdAmp = thdAmp;
    MEMR_mem.elicitorLevel = elicitorLevel;
end




% OLD CODE ----------------------------------------------------------------
%    for jj=1:24
%         M1 = squeeze(Click1(:,:,jj));
%         m1(:,jj) = mean(M1,2);
% 
%         M2 = squeeze(Click2(:,:,jj));
%         m2(:,jj) = mean(M2,2);
% 
%     end
%     rms1 = sqrt(mean(m1.^2,1));
%     rms2 = sqrt(mean(m2.^2,1));
%     rmsd = sqrt(mean((m1-m2).^2,1));
% 
%     plot(20*log10(rms2./rms1(1)),'r')
%     hold on
%     plot(20*log10(rms2./rms1),'b')
%     legend('rms/1','rms/each')
% 
%     figure
%     plot(rmsd,'Color',[0 1 0],'LineWidth',1)
%     xlabel('test time')
%     ylabel('rms difference')
% 
% 
% 
% keyboard
% 
% 
% 
%         zz = MEMR_mem.Z_sm.';
%         nFreqs = size(MEMR_mem.Z_sm,1);
%         for jj=1:nFreqs
% 
%             %base = mean([zz(1,jj),zz(24,jj)]);
%             base = zz(1,jj);
%             zz(:,jj) = zz(:,jj) ./ base;
%         end
% 
% 
%    % figure
%     % plot(m2(:)-m1(:))
%     % 
%     % figure
%     % plot(20*log10((m2(:)-m1(:))./max(m1(:))+1))
% 
%     % figure
%     % plot(sqrt(abs(m2(:)-m1(:)).^2)/sqrt(abs(m1(:).^2)))
% 
%         % EXTRACT THE NEEDED METRICS ----------------------------------
%         %peak delay -------------------
%         %There is no peak delay in the Mepani method!
%             d1 = MEMR_mem.d1.';
%             sm = 1; % smoothing factor (smaller numbers are more smooth)
%             n = length(d1); % number of clicks
%             x = linspace(0,8,n)'; % x-axis for smoothing (click number)
%             w = ones(size(x)); % weighting factor
% 
%             pp = csaps(x,d1,sm,[],w); % piecewise polynomial object
%             xx = linspace(0,8,1000)';
%             yy = ppval(pp,xx);
% 
%             dfdx = fnder(pp); % take derivative and solve for zero slope
%             peakXX = fnzeros(dfdx); % peak location in seconds
%             peakXX = peakXX(1,:);
% 
%             try
%                 peakYY = ppval(pp,peakXX);
%                 [~,peakIndx] = max(peakYY); 
%                 peakX = peakXX(peakIndx);
%                 peakY = peakYY(peakIndx);
%             catch
%                 keyboard
%             end
% 
%             % thresholds ------------------
%             %   use the Q concept: reduce 12 from peak
%             % targetLevelNoise = (90:-5:35); % this will be the descending run
%             % targetLevelNoise = targetLevelNoise(:)';
%             % targetLevelNoise = [targetLevelNoise,(fliplr(targetLevelNoise)+2.5)]; % add an ascending run
%             % targetLevelNoise = targetLevelNoise(:);
% 
% 
%             d1_scaled = (d1-1)./(peakY-1);
%             thd = 10.^(-12/20); % threshold is "Q12", or 12 dB down from peak
%             d1_scaled = d1_scaled - thd;
%             sm = 1; % smoothing factor (smaller numbers are more smooth)
%             pp = csaps(x,d1_scaled,sm,[],w); % piecewise polynomial object    
%             z2 = fnzeros(pp);
%             Thd = z2(1,:); % threshold times
%             try
%                 ons = Thd(find(Thd<peakX));
%                 thdOnsetTime = max(ons);
%                 offs = Thd(find(Thd>peakX));
%                 thdOffsetTime = min(offs);
%             catch ME
%                 thdOnsetTime = NaN;
%                 thdOffsetTime = NaN;
%             end
%             % convert threshold times to thresholds re: nominal elicitor level
%             noiseTarget = MEMR_mem.noiseTarget;
%             % noise went from 37.5 to 92.5 in 5 dB steps. Then back down
%             % from 90 to 35 in - 5 dB Steps. 
%             % Here, we are pretending that these changes occurred over 4
%             % seconds:
%             % 92.5-37.5 = 55 dB range in 4 seconds = (13.75 dB/s)
%             % y = mx + b
%             % stimLevel = 13.75x + 37.5 for the ascending run
%             % stimLevel = -13.75x + 90 for the ascending run
%             % but also account for reflex delay
%             slope = 13.75;
%             delay = 0;
%             b1 = 37.5;
%             b2 = 90; 
%             try
%                 ThdLvl(1) = slope*(Thd(1)-delay) + b1;
%                 thdOnsetLvl = ThdLvl(1); % onset threshold re: stim level
%                 ThdLvl(2) = -slope*(Thd(2)-4-delay) + b2;
%                 thdOffsetLvl = ThdLvl(2); % offset threshold re: stim level
%                 elicitorLevel1 = slope*(x-delay) + b1;
%                 elicitorLevel2 = -slope*(x-4-delay) + b2;
%                 elicitorLevel = [elicitorLevel1(1:12);elicitorLevel2(13:end)];
%             catch ME
%                 ThdLvl(1) = NaN;
%                 thdOnsetLvl = NaN; % onset threshold re: stim level
%                 ThdLvl(2) = NaN;
%                 thdOffsetLvl = NaN; % offset threshold re: stim level
%                 elicitorLevel = 17.5*(x-delay) + 45;
%             end
%             % go back and get threshold amplitudes from non-scaled d1
%             sm = 1; % smoothing factor (smaller numbers are more smooth)
%             %n = 24; % number of clicks
%             x = linspace(0,8,n)'; % x-axis for smoothing (click number)
%             w = ones(size(x)); % weighting factor
%             pp = csaps(x,d1,sm,[],w); % piecewise polynomial object
%             thdAmp = ppval(pp,thdOnsetTime);
% 
%             % hysteresis ---------------------
%             try
%                 pp = csaps(x,d1,sm,[],w); % piecewise polynomial object
%                 hh = fnint(pp); % integrate the smoothed spline
%                 A = ppval(hh,Thd(1));
%                 B = ppval(hh,peakX);
%                 C = ppval(hh,Thd(2));
%                 aucLeft = B-A; % area under the curve left
%                 aucRight = C-B; % area under the curve right
%                 hyst = aucRight / aucLeft; % hysteresis as a ratio of area under the curves
%                 hysteresis = hyst;
%             catch
%                 hysteresis = NaN;
%             end
%             % Calculate the slopes ------------------
%             try
%                 peak_index = round(peakX /(1/3));
%                 index_xa = round(thdOnsetTime / (1/3));
%                 index_xb = round(thdOffsetTime / (1/3));
%                 part1_x = x(index_xa:peak_index);
%                 part1_y = d1(index_xa:peak_index);
%                 part2_x = x(peak_index+1:index_xb);
%                 part2_y = d1(peak_index+1:index_xb);
%                 slope_ascending = polyfit(part1_x,part1_y,1);
%                 slope_descending = polyfit(part2_x,part2_y,1);
%                 slopeUp = slope_ascending(1);
%                 slopeDn = slope_descending(1);
%             catch
%                 slopeUp = NaN;
%                 slopeDn = NaN;
%             end
%         MEMR_mem.peakTime = peakX;
%         MEMR_mem.peakAmp = peakY;
%         MEMR_mem.delay = delay;
%         MEMR_mem.thdOnsetTime = thdOnsetTime;
%         MEMR_mem.thdOffsetTime = thdOffsetTime;
%         MEMR_mem.thdOnsetLvl = thdOnsetLvl;
%         MEMR_mem.thdOffsetLvl = thdOffsetLvl;
%         MEMR_mem.hysteresis = hysteresis;
%         MEMR_mem.slopeUp = slopeUp;
%         MEMR_mem.slopeDn = slopeDn;
%         MEMR_mem.thd = thd;
%         MEMR_mem.thdAmp = thdAmp;
% 
% 
%         % -------------------------------------------------------------
%         saveName = [folderName,'_Mepani','_Analysis1.mat'];
%         save([savePath,saveName],'MEMR_inc','MEMR_mem')
%         %saveName = [folderName,'_Run',num2str(jj),'_Analysis1.bmp'];
%         %saveas(h1,[savePath,saveName],'bmp')
% 
%         %keyboard
% 

        % % fix the jitter issue
        % for jj = 1:24 
        %     clicks = squeeze(Clicks1(:,:,jj));
        %     for kk=1:12
        %         [~,indx] = max(clicks(:,kk));
        %         if abs(indx-hat) > 1
        %             dummy = clicks(:,kk);
        %             offset = indx - hat;
        %             zpad = zeros(abs(offset),1);
        %             if offset >0
        %                 dummy = [dummy(offset:end-1);zpad];
        %             else % negative offset
        %                 dummy = [zpad;dummy(1:end+offset)];
        %             end
        %             clicks(:,kk) = dummy;
        %         end
        %     end
        %     plot(clicks,'b')
        %     pause(0.001)
        %     Clicks1(:,:,jj) = clicks;
        % end
        % for jj = 1:24 
        %     clicks = squeeze(Clicks2(:,:,jj));
        %     for kk=1:12
        %         [~,indx] = max(clicks(:,kk));
        %         if abs(indx-hat) > 1
        %             dummy = clicks(:,kk);
        %             offset = indx - hat;
        %             zpad = zeros(abs(offset),1);
        %             if offset >0
        %                 dummy = [dummy(offset:end-1);zpad];
        %             else % negative offset
        %                 dummy = [zpad;dummy(1:end+offset)];
        %             end
        %             clicks(:,kk) = dummy;
        %         end
        %     end
        %     plot(clicks,'r')
        %     pause(0.001)
        %     Clicks2(:,:,jj) = clicks;
        % end
        % %    clicks2 = squeeze(Clicks2(:,:,jj));
        % %     Clicks(:,counter:counter+12) = clicks1;
        % %     Clicks(:,counter+12:counter+12+12) = clicks2;
        % % 
        % %     counter = counter + 24;
        % 
        % title(num2str(ii))
        % xlim([460,510])
        % pause
        % c
