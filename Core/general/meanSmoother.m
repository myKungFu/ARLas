function [y] = meanSmoother(x,n)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% y = meanSmoother(x,n);
%
% Smooths data by taking the running mean.
% x = Input data. Should be a column vector
% n = Number of contiguous samples to use in computing the mean.
%
% Author: Shawn Goodman
% Date: January 7, 2010
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N = size(x,1);
nn = round((n-1)/2); % number of samples each side of the current sample
y = zeros(N,1); % initialize output

for ii=1:nn % filter run-on
    y(ii) = mean(x(1:ii+nn));
end
for ii=nn+1:N-nn % filter full on
	y(ii) = mean(x(ii-nn:ii+nn));
end
for ii=N-nn+1:N % filter run-off
	y(ii) = mean(x(ii-nn:N));
end
