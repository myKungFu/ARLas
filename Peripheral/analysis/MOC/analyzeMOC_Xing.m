function [] = analyzeMOC_Xing()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% analyzeMOC_Xing;
%
% This function analyzes MOC data from Xings's study.
% 
% Required non-standard functions:
%   fastfilter
%   AR
%   AR_engine
%   ARLas_fda
%
% Author: Shawn Goodman
% Date: December 11, 2023
% Last Updated: February 15, 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %--------------------------------------------------------------------------
% Xing, Be sure to set the directory to where the data are before you run
% this code!
% %--------------------------------------------------------------------------

whichProbe = '*A_moc*';
[ANALYSIS_MOC_Left,ANALYSIS_MEM_Left] = runMe(whichProbe);

whichProbe = '*B_moc*';
[ANALYSIS_MOC_Right,ANALYSIS_MEM_Right] = runMe(whichProbe);


disp('----- Left Ear -----')
disp('change (dB)  ChangeSNR  significant  baselineSNR')
ANALYSIS_MOC_Left
ANALYSIS_MEM_Left
disp('----- Right Ear -----')
ANALYSIS_MOC_Right
ANALYSIS_MEM_Right
disp(' ')

[ANALYSIS_MOC_Left';ANALYSIS_MEM_Left';ANALYSIS_MOC_Right';ANALYSIS_MEM_Right']
save('mocAnalyzedData','ANALYSIS_MOC_Left','ANALYSIS_MEM_Left','ANALYSIS_MOC_Right','ANALYSIS_MEM_Right')

disp('FINISHED WITH MOC ANALYSIS!')

end


% INTERNAL FUNCTIONS ------------------------------------------------------
function [ANALYSIS_MOC,ANALYSIS_MEM] = runMe(whichProbe)
    d4 = dir(whichProbe);
    nFiles = size(d4,1);
    clickWith = [];
    clickWithout = [];
    for jj=1:nFiles % loop across 5 saved mat files -------------------
        disp(['  Analyzing recording ',num2str(jj),' of ',num2str(nFiles)])
        dummy = load([d4(jj).folder,'\',d4(jj).name]);
        data = dummy.data;
        header = dummy.header;
        fs = header.fs;
        template = header.userInfo.clickTrainTemplate;
        template = template(:);
        data = data(:);
        b = HPF(fs);
        data = fastFilter(b,data);
        clickIndx = find(template == 1);
        % initially, keep 3 ms before the click and 15 ms after each click
        preSamples = round(0.003 * fs);
        postSamples = round(0.015 * fs);
        nClicks = length(clickIndx);
        
        % loop across intire block. In each block, there are 6 chunks
        % (3 without noise, 3 with noise, interleaved).
        start = 1;
        finish = 1536000;
        stepSize = finish;
        counter = 1;
        nChunks = 6;
        for kk=1:nChunks
            chunk = data(start:finish);
            Click = [];
            for mm=1:nClicks % loop across clicks -------------------------
                click = chunk(clickIndx(mm)-preSamples:clickIndx(mm)+postSamples);
                Click = [Click,click];
            end
            if mod(kk,2)==0
                clickWith = [clickWith,Click];
            else
                clickWithout = [clickWithout,Click];
            end
            start = start + stepSize;
            finish = finish + stepSize;
            counter = counter + 1;
        end
    end

    % separate out clicks and emissions for further processing.
    [~,hatIndx] = max(mean(clickWithout,2)); % peak location of click--should be 290
    if hatIndx ~= 290
        disp('Warning: hatIndx is not at expected location!')
        %keyboard
    end
    % for oae, start onset ramp 1 ms after peak and end onset ram 2 ms after peak
    start = round(0.001*fs);
    mocWith = clickWith(hatIndx+start:end,:);
    mocWithout = clickWithout(hatIndx+start:end,:);
    % use 1 ms on and off ramps
    winN = round(0.001*fs);
    h = hann(winN*2);
    h = h(1:winN);
    H = repmat(h,1,size(mocWith,2));
    mocWith(1:winN,:) = mocWith(1:winN,:) .* H;
    mocWith = flipud(mocWith);
    mocWith(1:winN,:) = mocWith(1:winN,:) .* H;
    mocWith = flipud(mocWith);
    mocWithout(1:winN,:) = mocWithout(1:winN,:) .* H;
    mocWithout = flipud(mocWithout);
    mocWithout(1:winN,:) = mocWithout(1:winN,:) .* H;
    mocWithout = flipud(mocWithout);

    % for mem, start onset ramp 1 ms before peak and end onset ram 1 ms after peak
    start = round(0.001*fs);
    memWith = clickWith(hatIndx-start:hatIndx+start,:);
    memWithout = clickWithout(hatIndx-start:hatIndx+start,:);
    % use 0.5 ms on and off ramps
    winN = round(0.0005*fs);
    h = hann(winN*2);
    h = h(1:winN);
    H = repmat(h,1,size(memWith,2));
    memWith(1:winN,:) = memWith(1:winN,:) .* H;
    memWith = flipud(memWith);
    memWith(1:winN,:) = memWith(1:winN,:) .* H;
    memWith = flipud(memWith);
    memWithout(1:winN,:) = memWithout(1:winN,:) .* H;
    memWithout = flipud(memWithout);
    memWithout(1:winN,:) = memWithout(1:winN,:) .* H;
    memWithout = flipud(memWithout);

    % lowpass filter
    b = LPF(fs);
    mocWith = fastFilter(b,mocWith);
    mocWithout = fastFilter(b,mocWithout);
    memWith = fastFilter(b,memWith);
    memWithout = fastFilter(b,memWithout);

    % artifact rejection
    rms = sqrt(mean(mocWith.^2,1));
    [indxW,~] = AR(rms,'moderate',0);
    rms = sqrt(mean(mocWithout.^2,1));
    [indxWo,~] = AR(rms,'moderate',0);
    indx = unique([indxW;indxWo]);
    mocNRejects = length(indx); % number of rejected buffers
    mocPCRejects = mocNRejects / length(rms); % percentage of rejected buffers
    mocWith = AR_engine(mocWith,indx);
    mocWithout = AR_engine(mocWithout,indx);

    rms = sqrt(mean(memWith.^2,1));
    [indxW,~] = AR(rms,'moderate',0);
    rms = sqrt(mean(memWithout.^2,1));
    [indxWo,~] = AR(rms,'moderate',0);
    indx = unique([indxW;indxWo]);
    memNRejects = length(indx); % number of rejected buffers
    memPCRejects = memNRejects / length(rms); % percentage of rejected buffers
    memWith = AR_engine(memWith,indx);
    memWithout = AR_engine(memWithout,indx);

    % run the quick analysis
    Analysis_moc = quickAnalyhsis(mocWith,mocWithout,fs);
    Analysis_mem = quickAnalyhsis(memWith,memWithout,fs);

    ii=1;

    ANALYSIS_MOC(ii,1) = Analysis_moc.change;
    ANALYSIS_MOC(ii,2) = Analysis_moc.changeSNR;
    ANALYSIS_MOC(ii,3) = Analysis_moc.sig;
    ANALYSIS_MOC(ii,4) = Analysis_moc.baselineSNR;

    ANALYSIS_MEM(ii,1) = Analysis_mem.change;
    ANALYSIS_MEM(ii,2) = Analysis_mem.changeSNR;
    ANALYSIS_MEM(ii,3) = Analysis_mem.sig;
    ANALYSIS_MEM(ii,4) = Analysis_mem.baselineSNR;
end
function [Analysis] = quickAnalyhsis(mocWith,mocWithout,fs)
    % Frequency domain analysis of MOC ------------------------------------
    ref = 0.00002; % 20 micropascals reference
    nfft = 2048; % number of samples in fft
    originalN = size(mocWith,1);
    [freq,sigWo,nfWo] = ARLas_fda(mocWithout,fs,ref,nfft,originalN);
    fmin = 1000;
    fmax = 2000;
    [~,indxMin] = min(abs(freq-fmin));
    [~,indxMax] = min(abs(freq-fmax));
    freq = freq(indxMin:indxMax);
    sigWo = sigWo(indxMin:indxMax);
    nfWo = nfWo(indxMin:indxMax);
    % we are looking in bin sizes of 46.9 Hz. Smooth by taking the mean
    % across 10 bins = 469 Hz and find the place with largest SNR in the
    % original (baseline) OAE signal.
    snrWo = sigWo-nfWo;
    snrWo_sm = meanSmoother(snrWo,10);
    %[~,hat] = min(abs(snrWo_sm-max(snrWo_sm)));
    % take 5 samples on either side (total of 11 samples = 515.9 Hz average
    snrWo = mean(snrWo_sm); % this is the baseline SNR
    sigWo_sm = meanSmoother(sigWo,10);
    %sigWo = mean(sigWo_sm); % this is the baseline signal (in dB)
    
    % Is there a significant difference in rms between the with and without?
    D = mocWith - mocWithout;
    [~,sigD,nfD] = ARLas_fda(D,fs,ref,nfft,originalN);
    sigD = 10.^(sigD/20)*ref;
    nfD = 10.^(nfD/20)*ref;
    sigD = sigD(indxMin:indxMax);
    nfD = nfD(indxMin:indxMax);
    snrD = sigD ./ nfD;
    snrD_sm = meanSmoother(snrD,10);
    snrD = mean(snrD_sm); % this is the mean Difference SNR (dB)
    snrD = 20*log10(snrD);
    sigD_sm = meanSmoother(sigD,10);
    sigD = mean(sigD_sm); % this is the mean Difference

    [~,sigA,nfA] = ARLas_fda(mocWithout,fs,ref,nfft,originalN);
    sigA = 10.^(sigA/20)*ref;
    sigA = sigA(indxMin:indxMax);
    sigA_sm = meanSmoother(sigA,10);
    sigA = mean(sigA_sm); 

    [~,sigB,~] = ARLas_fda(mocWith,fs,ref,nfft,originalN);
    sigB = 10.^(sigB/20)*ref;
    sigB = sigB(indxMin:indxMax);
    sigB_sm = meanSmoother(sigB,10);
    sigB = mean(sigB_sm); 

    mocEffect = 20*log10(sigB ./ sigA); % this is the MOC effect in dB
    if snrD >= 6
        sig = 1;
    else
        sig = 0;
    end
    %mocSNR = snrD;
    %mogSig = sig;

    Analysis.change = mocEffect;
    Analysis.changeSNR = snrD;
    Analysis.sig = sig;
    Analysis.baselineSNR = snrWo;
    
end

function b = HPF(Fs)
    Fstop = 500;             % Stopband Frequency
    Fpass = 1000;            % Passband Frequency
    Dstop = 0.01;            % Stopband Attenuation
    Dpass = 0.057501127785;  % Passband Ripple
    flag  = 'scale';         % Sampling Flag
    [N,Wn,BETA,TYPE] = kaiserord([Fstop Fpass]/(Fs/2), [0 1], [Dpass Dstop]);
    if mod(N,2)~=0
        N = N + 1;
    end
    b  = fir1(N, Wn, TYPE, kaiser(N+1, BETA), flag)';
end

function b = LPF(Fs)
    Fpass = 3000;            % Passband Frequency
    Fstop = 5000;            % Stopband Frequency
    Dpass = 0.057501127785;  % Passband Ripple
    Dstop = 0.001;           % Stopband Attenuation
    flag  = 'scale';         % Sampling Flag
    [N,Wn,BETA,TYPE] = kaiserord([Fpass Fstop]/(Fs/2), [1 0], [Dstop Dpass]);
    if mod(N,2)~=0
        N = N + 1;
    end
    b  = fir1(N, Wn, TYPE, kaiser(N+1, BETA), flag)';
end


