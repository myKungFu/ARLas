function [signal,noise,zbar,snr,noiseSamples] = bswSNR(Z)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [signal,noise,zbar,snr,noiseSamples] = bswSNR(Z);
%
% bisquares-weighted signal-to-noise ratio
% Calculate signal and noise levels from vector of complex coefficients (Z).
% Uses bisquares re-weighting to handle outliers (see notes at end).
% Z = a vector of complex coefficients obtained from LSF or FFT
%
% signal = weighted magnitude (dB SPL)
% noise = weighted sem (dB SPL)
% zbar = complex mean (linear units)
% snr = signal to noise ratio (dB)
% noiseSamples = residuals (Z after you subtract xbar)
%
% Author: Shawn Goodman
% Date: December 20, 2024
% Last Updated: February 1, 2025 -- ssg -- updated code to be more compatct
%               and to use variable names consistent with manuscript/poster
% Last Updated: February 7, 2025 -- ssg -- now called by LSF.m.
% Last Updated: April 5, 2025 -- ssg -- returns complex mean also
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ref = 0.00002; % reference (20 micropascals)
    Z = Z(:).'; % force coefficients to a column vector
    M = length(Z); % number of coefficients
    
    maxIterations = 50; % max iterations for determinging bisquares weights
    [Weights,iterations] = BS(Z,maxIterations); % get weights for bisquares re-weighting
    % make all the weights sum to 1
    W = Weights ./ sum(Weights);

    zbar = sum(Z.*W); % coherent mean 
    %sem = ((1/N)/(1-(1/N))) * sum(W.*((Z-zbar).*conj(Z-zbar))); % standard error of the mean
    sem = sqrt((M^-1/(1-M^-1)) * sum(W.*((Z-zbar).*conj(Z-zbar)))); % standard error of the mean
    
    signal = 20*log10(abs(zbar)/ref);
    noise = 20*log10(sem/ref);
    snr = signal - noise; % signal to noise ratio in dB
    noiseSamples = Z - zbar; % noise samples are the coefficients with the mean subtracted off

end

% INTERNAL FUNCTIONS ------------------------------------------------------
function [Weights,iterations] = BS(X,maxIterations) 
    % calculate bisquares weighting
    Weights = ones(size(X));
    [Residuals,epsilon] = weightedMean(X,Weights); % initial, unweighted mean
    k = tuningConstant(Residuals);
    delta = 1; % difference between current and previous iteration errors
    minDelta = 1E-12; % discontinue iterations when improvement is less than this
    ii=1;
    while (ii <= maxIterations) && (abs(delta) > minDelta)
        ii = ii+1;
        Weights = calculateWeights(Residuals,k); % new weighting based on previously-computed residuals
        [Residuals,epsilonNew] = weightedMean(X,Weights);
        delta = max(epsilon) - max(epsilonNew);
        epsilon = epsilonNew;
    end
    iterations = ii;
end
function [Residuals,epsilon] = weightedMean(X,Weights)  
    mu = sum(X .* Weights) ./ (sum(Weights)+eps); % complexe coherent weighted mean
    mu = abs(mu); % magnitude
    Mu = repmat(mu,1,size(X,2));
    Residuals = abs(X) - Mu; % unweighted residuals
    epsilon = sum((Weights .* Residuals).^2,2); % weighted sum of squared errors
end
function [k] = tuningConstant(Residuals)
    mar = median(abs(Residuals),2); % median absolute residuals
    sigma = mar / 0.6745; % robust estimate of the standard deviation...
    %sigma = ones(size(sigma)) * median(sigma); % ...across the entire input matrix
    k = (4.685 * sigma) + eps; % tuning constant
end
function [Weights] = calculateWeights(Residuals,k)
    K = repmat(k,1,size(Residuals,2));
    Weights = (1 - (Residuals./K).^2).^2; % calculate the weights for each sample
    Mask = abs(Residuals) < K; % find the weights that are "outliers"...
    Weights = Weights .* Mask; % ...and set them to zero
end
function [signal,noise,snr] = calculateSNRvanilla(b)
    % standard (non-weighted) SNR; here for comparison/troubleshooting purposes
    Xk = b(:)';
    K = size(Xk,2);
    Xbar = (1/K) * sum(Xk,2); % signal is the coherent mean
    Xbar2 = abs(Xbar) .^2; % signal energy
    XBAR = repmat(Xbar,1,K); % replicate to matrix (to avoid running a loop)
    S2 = (1/(K-1)) * sum((Xk - XBAR) .* conj(Xk - XBAR),2); % variance (this is the noise floor)
    Se2 = (1/K) * S2; % energy of the standard error
    ref = 0.00002;
    signal = 10*log10(Xbar2/(ref^2));
    noise = 10*log10(Se2/(ref^2));
    snr = signal - noise;
end
function [signal,noise] = getSNw(X,Weights)
    signal = sum(X .* Weights,2) ./ (sum(Weights,2)+eps); % final weighted mean
    Signal = repmat(signal,1,size(X,2));
    noise = sum(Weights .* (X - Signal).^2,2) ./ (sum(Weights,2)+eps); % weighted variance (single estimate)    
end


% Bisquare Mean NOTES:

% The value k is called a tuning constant; smaller values of k produce
% more resistance to outliers, but at the expense of lower efficiency when 
% the errors are normally distributed. The tuning constant is generally picked 
% to give reasonably high efficiency in the normal case; in particular,
% k = 4.685? for the bisquare (where ? is the standard deviation of the errors)
% produce 95-percent efficiency when the errors are normal, and still offer 
% protection against outliers. In an application, we need an estimate of the 
% standard deviation of the errors to use these results. Usually
% a robust measure of spread is employed in preference to the standard deviation 
% of the residuals. For example, a common approach is to take ? = MAR/0.6745, 
% where MAR is the median absolute residual.

% M-estimators can be vulnerable to high-leverage observations. A key concept
% in assessing influence is the breakdown point of an estimator: The breakdown point is the fraction of ‘bad’
% data that the estimator can tolerate without being affected to an arbitrarily large extent. For example, in
% the context of estimating the center of a distribution, the mean has a breakdown point of 0, because even
% one bad observation can change the mean by an arbitrary amount; in contrast the median has a breakdown
% point of 50 percent.
% There are also regression estimators that have breakdown points of nearly 50 percent. One such boundedinfluence
% estimator is least-trimmed squares (LTS) regression.
% Let us order the squared residuals from smallest to largest:
% The LTS estimator chooses the regression coefficients b to minimize the sum of the smallest m of the squared
% residuals,

% OLD CODE ----------------------------------------------------------------

    % no weighting
    % W = ones(size(W))/N;
    % zbar = sum(Z.*W); % coherent mean 
    % %sem = ((1/N)/(1-(1/N))) * sum(W.*((Z-zbar).*conj(Z-zbar))); % standard error of the mean
    % sem = (N^-1/(1-N^-1)) * sum(W.*((Z-zbar).*conj(Z-zbar))); % standard error of the mean
    % signal = 10*log10(abs(zbar^2)/ref^2)
    % noise = 10*log10(sem/ref^2)    

    % ORIGINAL CODE FOLLOWS:
 
    % 
    % b = Z;
    % 
    % Xk = b(:).'; % force coefficients to a column vector
    % 
    % %maxIterations = 50; % max iterations for determinging bisquares weights
    % %[Weights,iterations] = BS(Xk,maxIterations); % get weights for bisquares re-weighting
    % % make all the weights sum to 1
    % %wk0 = Weights ./ sum(Weights); % for N
    % wk0 = W;
    % wk1 = sum(wk0) - mean(wk0); % for N-1
    % K = size(Xk,2); % number of coefficients
    % 
    % 
    % Xbar = sum(Xk.*wk0); % signal is the (weighted) coherent mean
    %     % Xbar = (1/K) * sum(Xk,2); % signal is the coherent mean
    % Xbar2 = abs(Xbar) .^2; % signal energy
    % XBAR = repmat(Xbar,1,K); % replicate to matrix (to avoid running a loop)
    % S2 = (1/wk1) * sum(((Xk - XBAR) .* conj(Xk - XBAR)) .* wk0); % variance (this is the noise floor)
    % Se2 = (1/K) * S2; % energy of the standard error
    %     % S2 = (1/(K-1)) * sum((Xk - XBAR) .* conj(Xk - XBAR),2); % variance (this is the noise floor)
    %     % Se2 = (1/wk1) * S2; % energy of the standard error
    % 
    % signal = 10*log10(Xbar2/(ref^2)); % coherent mean in dB SPL
    % noise = 10*log10(Se2/(ref^2)); % noise is standard error of the mean
    % snr = signal - noise; % signal to noise ratio in dB
    % noiseSamples = b - Xbar; % noise samples are the coefficients with the mean subtracted off
    % 
