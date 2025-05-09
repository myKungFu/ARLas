function [y] = meanSmoother(x,n,w)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% y = meanSmoother(x,n,w);
%
% Smooths data by taking the running mean.
% x = Input data. Should be a column vector
% n = Number of contiguous samples to use in computing the mean.
% w = weighting function (optional input)
%
% Author: Shawn Goodman
% Date: January 7, 2010
% Updated: June 18, 2020 -- ssg
%                           Added capability to do weighted means
%                           Replaced mean with nanmean
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N = size(x,1);
nn = round((n-1)/2); % number of samples each side of the current sample
y = zeros(N,1); % initialize output

if nargin == 2 % stanadard (non-weighted mean)
    for ii=1:nn % filter run-on
        y(ii) = nanmean(x(1:ii+nn));
    end
    for ii=nn+1:N-nn % filter full on
        y(ii) = nanmean(x(ii-nn:ii+nn));
    end
    for ii=N-nn+1:N % filter run-off
        y(ii) = nanmean(x(ii-nn:N));
    end
elseif nargin == 3
    for ii=1:nn % filter run-on
        y(ii) = nansum(x(1:ii+nn).*w(1:ii+nn)) / nansum(w(1:ii+nn));
    end
    for ii=nn+1:N-nn % filter full on
        y(ii) = nansum(x(ii-nn:ii+nn).*w(ii-nn:ii+nn)) / nansum(w(ii-nn:ii+nn));
    end
    for ii=N-nn+1:N % filter run-off
        y(ii) = nansum(x(ii-nn:N).*w(ii-nn:N)) / nansum(w(ii-nn:N));
    end
end