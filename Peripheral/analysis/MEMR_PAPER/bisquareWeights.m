function [signal,noise,Weights] = bisquareWeights(X,maxIterations) 
% function to calculate bisquares weighting from a matrix of data.

Weights = ones(size(X));
if size(X,2) == 1
    [signal,noise] = getSNw(X,Weights);
    return
end
if nargin == 1
    maxIterations = 50;
end

[Residuals,epsilon] = weightedMean(X,Weights); % initial, unweighted mean
k = tuningConstant(Residuals);

delta = 1; % difference between current and previous iteration errors
minDelta = 1E-10; % discontinue iterations when improvement is less than this
ii=1;
while (ii <= maxIterations) && (abs(delta) > minDelta)
    ii = ii+1;
    Weights = calculateWeights(Residuals,k); % new weighting based on previously-computed residuals
    [Residuals,epsilonNew] = weightedMean(X,Weights);
    delta = max(epsilon) - max(epsilonNew);
    epsilon = epsilonNew;
end
[signal,noise] = getSNw(X,Weights);

% Internal functions ------------------------------------------------------

function [signal,noise] = getSNw(X,Weights)
signal = sum(X .* Weights,2) ./ (sum(Weights,2)+eps); % final weighted mean
Signal = repmat(signal,1,size(X,2));
noise = sum(Weights .* (X - Signal).^2,2) ./ (sum(Weights,2)+eps); % weighted variance (single estimate)    

function [Residuals,epsilon] = weightedMean(X,Weights)
mu = sum(X .* Weights,2) ./ (sum(Weights,2)+eps); % weighted mean
Mu = repmat(mu,1,size(X,2));
Residuals = X - Mu; % unweighted residuals
epsilon = sum((Weights .* Residuals).^2,2); % weighted sum of squared errors

function [k] = tuningConstant(Residuals)
mar = median(abs(Residuals),2); % median absolute residuals
sigma = mar / 0.6745; % robust estimate of the standard deviation...
%sigma = ones(size(sigma)) * median(sigma); % ...across the entire input matrix
k = (4.685 * sigma) + eps; % tuning constant

function [Weights] = calculateWeights(Residuals,k)
K = repmat(k,1,size(Residuals,2));
Weights = (1 - (Residuals./K).^2).^2; % calculate the weights for each sample
Mask = abs(Residuals) < K; % find the weights that are "outliers"...
Weights = Weights .* Mask; % ...and set them to zero

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

