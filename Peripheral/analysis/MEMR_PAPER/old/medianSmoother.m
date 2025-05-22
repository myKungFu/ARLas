function [y] = medianSmoother(x,n)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% y = medianSmoother(x,n);
%
% Smooths data by taking the running median.
% x = Input data. Should be a column vector
% n = Number of contiguous samples to use in computing the mean.
% w = weighting function (optional input)
%
% Author: Shawn Goodman
% Date: January 7, 2010
% Updated: June 18, 2020 -- ssg
%                           Added capability to do weighted means
%                           Replaced mean with nanmean
% Updated to median filter March 10, 2023
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    N = size(x,1);
    nn = round((n-1)/2); % number of samples each side of the current sample
    y = zeros(N,1); % initialize output
    
    for ii=1:nn % filter run-on
        y(ii) = median(x(1:ii+nn));
    end
    for ii=nn+1:N-nn % filter full on
        y(ii) = median(x(ii-nn:ii+nn));
    end
    for ii=N-nn+1:N % filter run-off
        y(ii) = median(x(ii-nn:N));
    end
end