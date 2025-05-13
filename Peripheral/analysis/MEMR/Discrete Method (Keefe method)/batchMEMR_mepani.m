function [] = batchMEMR_mepani()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% batchMEMR_mepani
% 
% Batch analysis of MEMR data using analyzeMEMR_vmepani.m.
%
% Author: Shawn Goodman & Ehsan Khalili
% Date: January 26, 2024
% Last Updated: January 26, 2024
% Last Updated: February 13, 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    parentDrive = 'C'; % use 'C' or 'D'
    dataPathName = [parentDrive,':\myWork\ARLas\MEMR_DATA\']; % location of raw data
    fileNameL = 'Ch3_ER10xA_memr_0001.mat';
    fileNameR = 'Ch4_ER10xB_memr_0001.mat';
    savePath = [parentDrive,':\myWork\ARLas\MEMR_Analysis_Final_M\']; % where to save analyzed data

    d = dir(dataPathName); % read in file directory information
    nSubjects = size(d,1)-2; % number of subject folders
    
    for ii=1:nSubjects
        disp(['Analyzing subject ',num2str(ii),' of ',num2str(nSubjects)])
        folderName = d(ii+2).name;
        d2 = dir([dataPathName,folderName,'\*mepani*']); % read in file directory information
        d3C = dir([dataPathName,folderName,'\',d2.name,'\*Ch3_ER10xA_memr*']);
        d3N = dir([dataPathName,folderName,'\',d2.name,'\*Ch4_ER10xB_memr*']);
        nFiles = size(d3C,1);

        for jj=1:nFiles % loop over elicitor levels -----------------------
            disp(['  Analyzing run ',num2str(jj),' of ',num2str(nFiles)])
            dummy = load([dataPathName,folderName,'\',d2.name,'\',d3C(jj).name]); % get clicks
            header = dummy.header;
            Clicks = dummy.data; % each elicitor level has 12 clicks
            clear dummy
            dummy = load([dataPathName,folderName,'\',d2.name,'\',d3N(jj).name]); % get noise
            headerN = dummy.header;
            Noise = dummy.data;
            clear dummy 
    
            [memr_inc,memr_mem] = analyzeMEMR_vmepani(header,Clicks,headerN,Noise,folderName);

            MEMR_inc.trend(:,jj) = memr_inc.trend; % : [24×1 double]
            MEMR_inc.timeTrend(:,jj) = memr_inc.timeTrend; %: [24×1 double]
            MEMR_inc.D(:,jj) = memr_inc.D(2,:); %: [2×40 double]
            MEMR_inc.D1(:,jj) = memr_inc.D1(2,:); %: [2×40 double]
            MEMR_inc.D2(:,jj) = memr_inc.D2(2,:); %: [2×40 double]
            MEMR_inc.d(1,jj) = memr_inc.d(2,1); %: [2×1 double]
            MEMR_inc.d1(1,jj) = memr_inc.d1(2,1); %: [2×1 double]
            MEMR_inc.d2(1,jj) = memr_inc.d2(2,1); %: [2×1 double]
            MEMR_inc.freq(:,jj) = memr_inc.freq; %: [40×1 double]
            MEMR_inc.Z(:,jj) = memr_inc.Z(2,:)'; %: [2×40 double]
            MEMR_inc.Z_sm(:,jj) = memr_inc.Z_sm(2,:)'; % [2×40 double]
            
            MEMR_mem.trend(:,jj) = memr_mem.trend; % : [24×1 double]
            MEMR_mem.timeTrend(:,jj) = memr_mem.timeTrend; %: [24×1 double]
            MEMR_mem.D(:,jj) = memr_mem.D(2,:); %: [2×40 double]
            MEMR_mem.D1(:,jj) = memr_mem.D1(2,:); %: [2×40 double]
            MEMR_mem.D2(:,jj) = memr_mem.D2(2,:); %: [2×40 double]
            MEMR_mem.d(1,jj) = memr_mem.d(2,1); %: [2×1 double]
            MEMR_mem.d1(1,jj) = memr_mem.d1(2,1); %: [2×1 double]
            MEMR_mem.d2(1,jj) = memr_mem.d2(2,1); %: [2×1 double]
            MEMR_mem.freq(:,jj) = memr_mem.freq; %: [40×1 double]
            MEMR_mem.Z(:,jj) = memr_mem.Z(2,:)'; %: [2×40 double]
            MEMR_mem.Z_sm(:,jj) = memr_mem.Z_sm(2,:)'; % [2×40 double]
            MEMR_mem.RMS(1,jj) = memr_mem.RMS; 
            MEMR_mem.RMSC(1,jj) = memr_mem.RMSC;
            MEMR_mem.pSPL(1,jj) = memr_mem.pSPL;
            MEMR_mem.noiseTarget(1,jj) = memr_mem.noiseTarget; % intended elicitor level (per header file)
        end
        % EXTRACT THE NEEDED METRICS ----------------------------------
        %peak delay -------------------
        %There is no peak delay in the Mepani method!
            d1 = MEMR_mem.d1.';
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
            noiseTarget = MEMR_mem.noiseTarget;
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


        % -------------------------------------------------------------
        saveName = [folderName,'_Mepani','_Analysis1.mat'];
        save([savePath,saveName],'MEMR_inc','MEMR_mem')
       
    end

end