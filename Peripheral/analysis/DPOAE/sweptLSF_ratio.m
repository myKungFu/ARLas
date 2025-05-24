function [signal,nf,snr,targets,zbar] = sweptLSF_ratio(x,Y,f,xmin,xmax,step,fs)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [signal,nf,snr,tarets,pdN] = sweptLSF_ratio(x,Y,f,xmin,xmax,step,fs);
%
% Perform least squares fits on matrices of data, sweeping across chunks of
% the data witha moving hann window. 
% THIS VERSION FOR LEVEL RATIO SWEEPS
% Required functions: LSF.m, bswSNR.m
%
% x = vector of x values for the sweep (level or frequency)
% Y = Data matrix on which to perform LSFs. It is assumed that each column of Y
%       is a repeated recording (sweep/buffer).
% f = a scalar frequency (in Hz) to be fit. Usually fdp.
% xmin = minimum value in the X vector to target
% xmax = maximum value in the X vector to target
% step = step size of targets
% fs = sampling rate (samples / sec).
%
% signal = weighted mean magnitude, taken as abs of the complex mean across columns of Y.
% nf = noise floor, taken as the weighted standard error of the mean.
% snr = signal-to-noise ratio in dB
% targets = list of targets in X for downsampled vector
% zbar = complex mean
%
% noiseSampes = a vector of complex coefficients, one from each column of
%       Y. The mean has been subtracted, leaving the variability around the mean.
%       This can be used in other code to esitmate the noise distribution of the
%       recordings.
%
% Author: Shawn Goodman, PhD
% Auditory Research Lab, the University of Iowa
% Last Updated: February 1-7, 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    frameLen = 0.5; % analysis frame length (s)
    % for a 9 second sweep from 0-70
    % 248 ms gives 50% overlap
    % 300 ms gives 58% overlap
    % 500 ms gives 75% overlap
    % 750 ms gives 84% overlap
    frameN = round(fs*frameLen); % number of samples in each frame
    
    targets = (xmin:step:xmax)'; % list of target FPL levels for downsampled vector
    nIterations = length(targets);
    indices = zeros(nIterations,1);
    for ii=1:nIterations % find sample location of each target
        [~,indices(ii,1)] = min(abs(targets(ii)-x));
    end
    
    [nSamples,nSweeps] = size(Y); % size of the data
    
    frameN2 = ceil(frameN/2); % half-frame size (number of samples)
    w = hann(frameN);
    W = repmat(w,1,nSweeps);
    
    % bad results happen when the frames are too empty. Remove frames less than half filled
    [~,lossIndx] = min(abs(indices - frameN2));
    indices = indices(lossIndx:end-lossIndx);
    targets = targets(lossIndx:end-lossIndx);
    nIterations = length(indices);
    pcOverlap = 1 - (median(diff(indices)) / frameN); % percent overlap in adjacent frames
    
    for ii=1:nIterations % loop across frames
        if mod(ii,10)==0
            disp(['Analyzing frame ',num2str(ii),' of ',num2str(nIterations)])
        end
        start = indices(ii)-frameN2+1; % starting sample of current frame
        finish = indices(ii)+frameN2; % ending sample of current frame
        if start < 0 % handle frame run on and run off of the data matrix
            finish = indices(ii)+frameN2;
            Chunk = Y(1:finish,:);
            Pad1 = zeros(abs(start)+1,nSweeps);
            Chunk = [Pad1;Chunk];
        elseif finish > nSamples
            extra = finish - nSamples;
            Chunk = Y(start:end,:);
            Pad1 = zeros(extra,nSweeps);
            Chunk = [Chunk;Pad1];
        else
            Chunk = Y(start:finish,:);
        end
        Chunk = Chunk .* W;
        if ii==1 % the first time, have LSF calculate solution matrix X
            [signal(ii,1),nf(ii,1),snr(ii,1),noiseSamples(ii,:),X,zbar(ii,1)] = LSF(Chunk,f);
        else % after that, reuse solution matrix X
            [signal(ii,1),nf(ii,1),snr(ii,1),noiseSamples(ii,:),~,zbar(ii,1)] = LSF(Chunk,f,X);
        end
    end
    
    % correct group delay to account for the LSF shifting
    mag = abs(zbar);
    phi = angle(zbar);
    samples = median(diff(indices));
    seconds = samples/fs;
    cycles = seconds / (1/f);
    radians = cycles * 2*pi;
    correction = (0:1:length(zbar)-1)'*radians;
    phi = unwrap(phi)-correction;
    zbar = mag.*exp(1i*phi);
    % phase delay is getting earlier (phase getting more negative) as level
    % increases. Remeber this is not group delay! But can you still think
    % of this as a time shift?
    % plot(angle(zbar)/(2*pi) / f * 1000)
    % suggests coming back faster by about 0.1 ms?


end

% OLD CODE ----------------------------------------------------------------

    %k = sqrt(mean(w.^2)); % correct for the hann window. added by ssg 3/18/2025
    %w = w / k;
    %w = w / 1.1;

    % process the noise samples to create a noise distribution
    %NoiseSamples = noiseSamples;
    %noiseSamples = noiseSamples(:); % reshape into vector
    %[indx,nRejects] = AR(abs(noiseSamples.'),'mild',0); % get rid of extreme values
    %noiseSamples(indx) = [];
    %pdN = fitdist([real(noiseSamples);imag(noiseSamples)],'kernel'); % create kernel-based PDF for noise
