function [y] = maxSmoother(x,n)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% y = maxSmoother(x,n);
%
% Smooths data by taking the running maximum.
% x = Input data. Should be a column vector
% n = Number of contiguous samples to use in computing the mean.
%
% Author: Shawn Goodman
% Date: August 15, 2013
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N = size(x,1);
nn = round((n-1)/2); % number of samples each side of the current sample
y = zeros(N,1); % initialize output

for ii=1:nn % filter run-on
    y(ii) = max(x(1:ii+nn));
end
for ii=nn+1:N-nn % filter full on
	y(ii) = max(x(ii-nn:ii+nn));
end
for ii=N-nn+1:N % filter run-off
	y(ii) = max(x(ii-nn:N));
end
