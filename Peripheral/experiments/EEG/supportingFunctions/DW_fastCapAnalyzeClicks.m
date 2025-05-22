function [mag,nf,snr,S,cap,time,stimulus_dBreFO,phi,signal,noise,frequency] = DW_fastCapAnalyzeClicks(header,Data)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% mag,nf,snr,freq,cap,time,stimulus_dBreFO,phi,signal,noise,frequency] = DW_fastCapAnalyzeClicks(header,Data);
%
% Analysis code for CAPs measured at mulitple click bandwidths and levels.
% This code should be called by DW_fastCapAudioClicks.m
%
% Authors: Shawn S Goodman, PhD & Jeffery T Lichtenhan PhD
% Auditory Research Lab, Dept. of Comm. Sciences & Disorders, The University of Iowa
% Auditory Physiology Lab, Dept. of Otolaryngology, Washington University School of Medicine, St. Louis
% Date: October 30, 2019
% Updated: October 30, 2019 - ssg
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% de-interleave the data
sortIndx = header.userInfo.sortIndx;
Data = Data(:,sortIndx);  % put data back into sorted order

% unpack needed information from the header
levels = header.userInfo.levels;
freqs = header.userInfo.freqs;
nReps = header.userInfo.nReps;
stim_dBreFO = header.userInfo.stim_dBreFO;
fs = header.fs;
cuts = header.userInfo.cuts;
%stimSL = stimSPL;


% Analyze Data for each level and frequency.
nLevels = length(levels); % number of levels tested
nClicks = length(freqs);   % number of frequencies tested
nSamples = size(Data,1);  % number of samples in each recorded buffer
counter = 1; 
start = 1;

for ii=1:nLevels
    for jj=1:nClicks
        ss = freqs(jj);              % current frequency
        nn = nReps(ii,jj);          % number of recordings at this combintaion
        finish = start + nn-1;
        chunk = Data(:,start:finish);
        if nn > 0                     % if data were actually recorded at this level and frequency combination
            nSubs = size(chunk,2)/6;  % number of sub-averages to perform. Always in sets of 6 to cancel electrical artifact
            step2 = 6;
            start2 = 1;
            finish2 = step2;
            counter2 = 1;
            chunk2 = zeros(nSamples,nSubs);
            for kk=1:nSubs           % loop to create sub-averages
                chunk2(:,counter2) = mean(chunk(:,start2:finish2),2);
                start2 = start2 + step2;
                finish2 = finish2 + step2;
                counter2 = counter2 + 1;
            end
            % this is the new way; it uses Fourier-based RMS magnitudes
            [Cap,time] = getCAP(chunk2,cuts,fs);
            %[FF,SPL,SL,MAG,NF,nF,nMaxLevels,peakLocations1,peakLocations2] = DW_analyzeSeries(header,Data,peakLocations1,peakLocations2);
            [mag(1,counter),nf(1,counter),snr(1,counter),phi(1,counter),cap(:,counter),signal(:,counter),noise(:,counter),frequency] = DW_capAnalysis(Cap,time,fs);
            S(1,counter) = ss;
            
            %stimulusSL(1,counter) = stimSL(ii,jj);   % stimulus sensation level
            stimulus_dBreFO(1,counter) = stim_dBreFO(ii,jj); % stimulus sound pressure level
            counter = counter + 1;
            start = finish + 1;
        end
    end
end
end % EOF: end of experiment file

% INTERNAL FILES ----------------------------------------------------------
function [Cap,time] = getCAP(chunk2,cuts,fs)
    D1 = chunk2(cuts(1,1):cuts(2,1)-1,:); % pip
    D2 = chunk2(cuts(2,1):end,:);         % invertd pip
    Cap = (D1 + D2)/2;                  % compound action potential
    Cap = Cap * 1000000; % cap in uV
    nSamples = size(Cap,1);           % number of samples in each response
    time = (0:1:nSamples-1)'/fs;      % time vector associated with each CAP waveform
    time = time *1000;                % time in ms                
end
            
function b = lpf(fs,Fpass) % lowpass filter coefficients
    Fstop = 1.7 * Fpass;    % Stopband Frequency
    Dpass = 0.057501127785;  % Passband Ripple
    Dstop = 0.0001;          % Stopband Attenuation
    flag  = 'scale';         % Sampling Flag
    [N,Wn,BETA,TYPE] = kaiserord([Fpass Fstop]/(fs/2), [1 0], [Dstop Dpass]);
    if mod(N,2)~= 0
        N = N + 1;
    end
    b  = fir1(N, Wn, TYPE, kaiser(N+1, BETA), flag)';
end

% OLD CODE
%snrAdj = snr;
%snrAdj(snrAdj<0) = 0; % adjusted snr does not allow snr to be less than 0 dB
%criterion = 6; % 6 dB SNR criterion must be exceeded to say that cap is present
%present = snr > criterion;


% % prepare analyzed data for plotting
% F = unique(freq); % get the unique frequency values
% nF = length(F);   % number of unique frequencies tested
% % find the larges number of stimulus levels for any frequency
% for ii=1:nF
%    indx = find(freq==F(ii));
%    nVals(ii) = length(indx);
% end
% nL = max(nVals); % max number of stimulus levels tested
% P2P = nan(nL,nF); % initialize output matrix
% LVL = nan(nL,nF); % initialize output matrix
% for ii=1:nF
%    indx = find(freq==F(ii));
%    for jj=1:length(indx)
%        P2P(jj,ii) = mag(indx(jj));
%        NF(jj,ii) = nf(indx(jj));
%        LVL(jj,ii) = stimulusSPL(indx(jj));
%    end
% end
% FF = repmat(F,nL,1);
% offset = linspace(0,500,nL)';

