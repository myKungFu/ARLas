function [y] = nanmedianSmoother(x,n)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% y = nanmedianSmoother(x,n);
%
% Smooths data by taking the running mean.
% x = Input data. Should be a column vector
% n = Number of contiguous samples to use in computing the median.
%     Should be an odd number, so that an even split is obtained on 
%     each side of the current data point.
%
% Author: Shawn Goodman
% Date: January 7, 2010
% Updated: 5/28/2013
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% REMOVE SAFETY GUARDS TO MAKE IT RUN FASTER 1/4/2011
% x = x(:); % force input to a column vector
% N = length(x);
% if n > N,
%     disp('ERROR: n must be <= N.')
%     y = [];
%     return
% end
% if mod(n,2)== 0
%     disp('WARNING: n must be odd!')
%     disp('         Using next higher odd number...')
%     n = n+1;
% end


N = size(x,1);
nn = round((n-1)/2); % number of samples each side of the current sample
y = zeros(N,1); % initialize output

%h = hann(n);

for ii=1:nn % filter run-on
%    hh = hann(length(x(1:ii+nn)));
    y(ii) = nanmedian(x(1:ii+nn));
end
for ii=nn+1:N-nn % filter full on
	y(ii) = nanmedian(x(ii-nn:ii+nn));
end
for ii=N-nn+1:N % filter run-off
%    hh = hann(length(x(ii-nn:N)));
	y(ii) = nanmedian(x(ii-nn:N));
end
 

