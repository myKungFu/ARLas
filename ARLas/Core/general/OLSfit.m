function [coeff,anova,modelSummary] = OLSfit(X,y,alfa,w)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [coeff,anova,modelSummary] = OLSfit(X,y,alfa,w);
%
% Ordinary Least-Squares Regression for the model
% y = a + b1*X1 + b2*X2 + ... + bk*Xk.
% X = matrix of independent variables.
% y = vector of dependent variables.
% alfa = level of significance for confidence intervals.
% w = vector of weighting values (optional input)
%
% coeff is a struct array with fields:
%       b = model coefficients
%       se = standard error of the coefficients
%       t = t-value of the coefficients
%       p = significance of the coefficients
%       ci = confidence interval for the coefficients (specified by alfa)
% anova is a struct array with fields:
%      dfReg = regression degrees of freedom
%       dfRes = residual degrees of freedom
%       dfTot = total degrees of freedom
%       SSreg = sum of squares regression
%       SSres = sum of squares residual
%       SStot = sum of squares total
%       MSreg = mean squares regression
%       MSres = mean squares residual
%       F = F-ratio
%       p = significance of the F-ratio
% modelSummary is a struct array with fields:
%       R = multiple R
%       R2 = multiple R-squared
%       R2adj = adjusted R-squared
%       rmse = root mean square error (standard error of the estimate)
%
% 
% BUT YOU NEED TO TAKE CARE OF THE WEIGHTING
% BEFORE CALCULATING RMSE, ETC..
%
% Author: Shawn Goodman
% ©  May 21, 2005
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin > 4,
    disp('ERROR: too many input arguments specified.')
    return;
elseif nargin == 3, % set default values if not specified
    w = [];
elseif nargin == 2,
    alfa = 0.05;
    w = [];
elseif nargin < 2,
    disp('ERROR: input arguments x and y must be specified.')
    return
end
% check for legal input and ensure that inputs are column vectors
[X,y,w,alfa,returnCode,switchback] = checkInput(X,y,w,alfa);
if returnCode ~= 0,
    coeff = []; anova = []; modelSummary = [];
    return
end

[rows,columns] = size(X); % size of the solution matrix
k = columns-1; % k is the number of IV (the column of ones for y-intercept/dc offset is not counted)
n = rows; % number of observations

if ~isempty(w) % if a weighting vector exists...
    w = w(:); % vectorize as a column vector...
    sqrtw = sqrt(w); % ...and take the square root
    X = repmat(sqrtw,1,columns).*X; % replicate sqrtw to matrix the size of X; J is the weighted version of X
    y = sqrtw.*y; % Dy is the weighted version of y
end

% Solve the least squares problem
[Q,R] = qr(X,0); % orthogonal-triangular decomposition; R is the Cholesky factor of the X matrix
                 % the input argument 0 makes an "economy sized" decomposition, so that
                 % [nSamples,nIVs] = size(Q), and 
                 % [nIvs,nIVs] = size(R).
b = full(R\(Q'*y)); % Same as p = D*X\(D*y); b is a vector of coefficients; also same as b = (X'*X)\(X'*Y).
yhat = X*b;      % Predicted responses at each data point.
residuals = y - yhat;

% Find a confidence interval for each component of X
nu = rows-columns;                       % Residual degrees of freedom
rmse = norm(residuals)/sqrt(nu);        % Root mean square error.
       % norm() here is the L2 norm, the sqrt of the sum of the squared values
RI = R\eye(columns); % the inverse of R
se = rmse * sqrt(sum(abs(RI).^2,2)); % standard errors of the coefficients
    % se = rmse ./ sqrt(sum(abs(R).^2,2)); % this is equivalent to the above
    % for single regressor, this is equivalent to se = rmse / sqrt(sum((X(:,2)-xbar).^2)).
    % for bivariate regressors, this is equivalent to se = rmse / sqrt(sum((X(:,2)-xbar).^2)*r2x2x3).
t = b ./ se; % t-values of the coefficients
bSig = zeros(k,1); % significance values of coefficients via t-test (2-tailed)
for ii=1:length(t),
    bSig(ii) = (1-mytcdf(t(ii),n-k-1))*2; % get the p-value, 2-tailed 
    %bSig(ii) = (1-mytcdf(t(ii),n-k-1)); % get the p-value, 1-tailed 
end

criticalT = mytinv(1-alfa/2,n-2); % critical two-tailed t-value
ci = [b-(criticalT*se), b+(criticalT*se)]; % confidence intervals on coefficients
SSreg = norm(yhat-mean(y))^2;  % Regression sum of squares.
SStotal = norm(y-mean(y))^2;   % Total sum of squares.
MSreg = SSreg / k; % mean square regression
SSres = SStotal - SSreg; % sum of squares residual
MSres = SSres / (n-k-1); % mean square residual
if MSres == 0,
    F = 0;
    p = 0;
    R2 = 0;
    R2adj = 0;
else
    F = MSreg / MSres; % F ratio
    p = 1 - myfcdf(F,k,n-k-1); % significance of the F ratio
    R2 = SSreg / SStotal;          % multiple R-square statistic.
    R2adj = 1- ((1 - R2)*((n-1)/(n-k-1))); % adjusted R squared; Calculation from Wherry (1931), cited in T&F, p 147.
end

% put output into three different structures:
modelSummary.R = sqrt(R2);
modelSummary.R2 = R2;
modelSummary.R2adj = R2adj;
modelSummary.rmse = rmse;

anova.dfReg = k;
anova.dfRes = n-k-1;
anova.dfTot = n-1;
anova.SSreg = SSreg;
anova.SSres = SSres;
anova.SStot = SStotal;
anova.MSreg = MSreg;
anova.MSres = MSres;
anova.F = F;
anova.p = p;

coeff.b = b;
coeff.se = se;
coeff.t = t;
coeff.p = bSig;
coeff.ci = ci;


%keyboard
% Beta is calculated by replacing the ordinary raw b coefficients by b times s(X)/s(Y), where
% s(Y) is the standard deviation of the dependent variable, and s(X) is the standard deviation of the
% predictor.


% old matlab code...
%s2 = rmse^2;                    % Estimator of error variance.
% F = (RSS/(p-1))/s2;       % F statistic for regression
% prob = 1 - fcdf(F,p-1,nu);   % Significance probability for regression
%rmse = sqrt(sum(residuals.^2 .* .2))/sqrt(nu); % use this to use correction factors
%xdiag=sqrt(sum((RI .* RI)',1))';




% % Find the standard errors of the residuals.
% % Get the diagonal elements of the "Hat" matrix.
% % Calculate the variance estimate obtained by removing each case (i.e. sigmai)
% % see Chatterjee and Hadi p. 380 equation 14.
% T = X*RI;
% hatdiag=sum((T .* T)',1)';
% ok = ((1-hatdiag) > sqrt(eps));
% hatdiag(~ok) = 1;
% if nu < 1, 
%   ser=rmse*ones(length(y),1);
% elseif nu > 1
%   denom = (nu-1) .* (1-hatdiag);
%   sigmai = zeros(length(denom),1);
%   sigmai(ok) = sqrt((nu*s2/(nu-1)) - (r(ok) .^2 ./ denom(ok)));
%   ser = sqrt(1-hatdiag) .* sigmai;
%   ser(~ok) = Inf;
% elseif nu == 1
%   ser = sqrt(1-hatdiag) .* rmse;
%   ser(~ok) = Inf;
% end
% 
% % Create confidence intervals for residuals.
% tval chaged to criticalT in my code...
% Z=[(r-tval*ser) (r+tval*ser)]';
% rint=Z';
% 


%---------------------------------------------
function [X,y,w,alfa,returnCode,switchback] = checkInput(X,y,w,alfa)
returnCode = 0;
switchback = 0;
[rows,columns] = size(y);
if rows + columns == 2,
    disp('ERROR: y cannot be scalar.')
    returnCode = 1;
end
if rows < columns,
    y = y';
    X = X';
    w = w';
    switchback = 1;
end
[rowsY,columnsY] = size(y);
[rowsX,columnsX] = size(X);
if rowsX ~= rowsY,
    disp('ERROR: size of y and X do not match.')
    returnCode = 2;
end
if columnsY ~= 1,
    disp('ERROR: y (dependent variable) must be a vector, not matrix.')
    returnCode = 3;
end
if ~any(sum(X)) == rowsX,
    disp('WARNING: X (matrix of independent variables) should include a vector of ones for intercept.')
    disp('         This vector is being automatically generated...')
    x = ones(rows,1);
    X = [x,X];
end
if alfa > 1 | alfa < 0,
    disp('ERROR: alfa is > 1 or < 0.')
    returnCode = 4;
end





function x = mytinv(p,v);
%TINV   Inverse of Student's T cumulative distribution function (cdf).
%   X=TINV(P,V) returns the inverse of Student's T cdf with V degrees 
%   of freedom, at the values in P.
% Initialize Y to zero, or NaN for invalid d.f.
x=zeros(size(p));
x(v <= 0) = NaN;
k = find(v == 1 & ~isnan(x));
if any(k)
  x(k) = tan(pi * (p(k) - 0.5));
end
% The inverse cdf of 0 is -Inf, and the inverse cdf of 1 is Inf.
k0 = find(p == 0 & ~isnan(x));
if any(k0)
    tmp   = Inf;
    x(k0) = -tmp(ones(size(k0)));
end
k1 = find(p ==1 & ~isnan(x));
if any(k1)
    tmp   = Inf;
    x(k1) = tmp(ones(size(k1)));
end
% For small d.f., call betainv which uses Newton's method
k = find(p >= 0.5 & p < 1 & ~isnan(x) & v < 1000);
if any(k)
    z = mybetainv(2*(1-p(k)),v(k)/2,0.5);
    x(k) = sqrt(v(k) ./ z - v(k));
end
k = find(p < 0.5 & p > 0 & ~isnan(x) & v < 1000);
if any(k)
    z = mybetainv(2*(p(k)),v(k)/2,0.5);
    x(k) = -sqrt(v(k) ./ z - v(k));
end
% For large d.f., use Abramowitz & Stegun formula 26.7.5
k = find(p>0 & p<1 & ~isnan(x) & v >= 1000);
if any(k)
   xn = mynorminv(p(k));
   df = v(k);
   x(k) = xn + (xn.^3+xn)./(4*df) + ...
           (5*xn.^5+16.*xn.^3+3*xn)./(96*df.^2) + ...
           (3*xn.^7+19*xn.^5+17*xn.^3-15*xn)./(384*df.^3) +...
           (79*xn.^9+776*xn.^7+1482*xn.^5-1920*xn.^3-945*xn)./(92160*df.^4);
end


function x = mybetainv(p,a,b);
%BETAINV Inverse of the beta cumulative distribution function (cdf).
%   X = BETAINV(P,A,B) returns the inverse of the beta cdf with 
%   parameters A and B at the values in P.
%   Initialize x to zero.
x = zeros(size(p));
%   Return NaN if the arguments are outside their respective limits.
k = find(p < 0 | p > 1 | a <= 0 | b <= 0);
if any(k),
   tmp = NaN;
   x(k) = tmp(ones(size(k))); 
end
% The inverse cdf of 0 is 0, and the inverse cdf of 1 is 1.  
k0 = find(p == 0 & a > 0 & b > 0);
if any(k0), 
    x(k0) = zeros(size(k0)); 
end
k1 = find(p==1);
if any(k1), 
    x(k1) = ones(size(k1)); 
end
% Newton's Method.
% Permit no more than count_limit interations.
count_limit = 100;
count = 0;
k = find(p > 0 & p < 1 & a > 0 & b > 0);
if isempty(k)
   return;
end
pk = p(k);
%   Use the mean as a starting guess. 
xk = a(k) ./ (a(k) + b(k));
% Move starting values away from the boundaries.
if xk == 0,
    xk = sqrt(eps);
end
if xk == 1,
    xk = 1 - sqrt(eps);
end
h = ones(size(pk));
crit = sqrt(eps); 
% Break out of the iteration loop for the following:
%  1) The last update is very small (compared to x).
%  2) The last update is very small (compared to 100*eps).
%  3) There are more than 100 iterations. This should NEVER happen. 
while(any(abs(h) > crit * abs(xk)) & max(abs(h)) > crit    ...
                                 & count < count_limit), 
    count = count+1;    
    h = (mybetacdf(xk,a(k),b(k)) - pk) ./ mybetapdf(xk,a(k),b(k));
    xnew = xk - h;
% Make sure that the values stay inside the bounds.
% Initially, Newton's Method may take big steps.
    ksmall = find(xnew <= 0);
    klarge = find(xnew >= 1);
    if any(ksmall) | any(klarge)
        xnew(ksmall) = xk(ksmall) /10;
        xnew(klarge) = 1 - (1 - xk(klarge))/10;
    end
    xk = xnew;  
end
% Return the converged value(s).
x(k) = xk;
if count==count_limit, 
    fprintf('\nWarning: BETAINV did not converge.\n');
    str = 'The last step was:  ';
    outstr = sprintf([str,'%13.8f'],h);
    fprintf(outstr);
end


function x = mynorminv(p,mu,sigma);
%NORMINV Inverse of the normal cumulative distribution function (cdf).
%   X = NORMINV(P,MU,SIGMA) finds the inverse of the normal cdf with
%   mean, MU, and standard deviation, SIGMA.
if nargin < 3, 
    sigma = 1;
end
if nargin < 2;
    mu = 0;
end
[errorcode p mu sigma] = mydistchck(3,p,mu,sigma);
% Allocate space for x.
x = zeros(size(p));
% Return NaN if the arguments are outside their respective limits.
k = find(sigma <= 0 | p < 0 | p > 1 | isnan(p));
if any(k)
    tmp  = NaN;
    x(k) = tmp(ones(size(k))); 
end
% Put in the correct values when P is either 0 or 1.
k = find(p == 0);
if any(k)
    tmp  = Inf;
    x(k) = -tmp(ones(size(k)));
end
k = find(p == 1);
if any(k)
    tmp  = Inf;
    x(k) = tmp(ones(size(k))); 
end
% Compute the inverse function for the intermediate values.
k = find(p > 0  &  p < 1 & sigma > 0);
if any(k),
    x(k) = sqrt(2) * sigma(k) .* erfinv(2 * p(k) - 1) + mu(k);
end



function [errorcode,varargout] = mydistchck(nparms,varargin)
%DISTCHCK Checks the argument list for the probability functions.
errorcode = 0;
n = nargout-1;
varargout = cell(1,n);
if nparms == 1
    varargout{1} = varargin{1};
    return;
end
% Get size of each input, check for scalars, copy to output
sz = cell(1,n);
isscalar = logical(zeros(1,n));
for j=1:n
   s = size(varargin{j});
   sz{j} = s;
   isscalar(j) = (prod(s) == 1);
   varargout{j} = varargin{j};
end
% Done if all inputs are scalars.  Otherwise fetch their common size.
if (all(isscalar)), return; end
t = sz(~isscalar);
size1 = t{1};
% Scalars receive this size.  Other arrays must have the proper size.
for j=1:n
   sizej = sz{j};
   if (isscalar(j))
      t = zeros(size1);
      t(:) = varargin{j};
      varargout{j} = t;
   elseif (~isequal(sizej,size1))
      errorcode = 1;
      return;
   end
end


function p = mybetacdf(x,a,b);
%BETACDF Beta cumulative distribution function.
%   P = BETACDF(X,A,B) returns the beta cumulative distribution
%   function with parameters A and B at the values in X.
if nargin<3, 
   error('Requires three input arguments.'); 
end
[errorcode x a b] = mydistchck(3,x,a,b);
if errorcode > 0
   error('Requires non-scalar arguments to match in size.');
end
% Initialize P to 0.
p = zeros(size(x));
k1 = find(a<=0 | b<=0);
if any(k1)
   tmp = NaN;
   p(k1) = tmp(ones(size(k1))); 
end
% If is X >= 1 the cdf of X is 1. 
k2 = find(x >= 1);
if any(k2)
   p(k2) = ones(size(k2));
end
k = find(x > 0 & x < 1 & a > 0 & b > 0);
if any(k)
   p(k) = betainc(x(k),a(k),b(k));
end
% Make sure that round-off errors never make P greater than 1.
k = find(p > 1);
p(k) = ones(size(k));



function y = mybetapdf(x,a,b)
%BETAPDF Beta probability density function.
%   Y = BETAPDF(X,A,B) returns the beta probability density 
%   function with parameters A and B at the values in X.
if nargin < 3, 
   error('Requires three input arguments.');
end
[errorcode x a b] = mydistchck(3,x,a,b);
if errorcode > 0
    error('Requires non-scalar arguments to match in size.');
end
% Initialize Y to zero.
y = zeros(size(x));
% Return NaN for parameter values outside their respective limits.
k1 = find(a <= 0 | b <= 0);
if any(k1)
    tmp = NaN;
    y(k1) = tmp(ones(size(k1))); 
end
% Return Inf for x = 0 and a < 1 or x = 1 and b < 1.
% Required for non-IEEE machines.
k2 = find((x == 0 & a < 1) | (x == 1 & b < 1));
if any(k2)
    tmp = Inf;
    y(k2) = tmp(ones(size(k2))); 
end
% Return the beta density function for valid parameters.
k = find(~(a <= 0 | b <= 0 | x <= 0 | x >= 1));
if any(k)
%    y(k) = x(k) .^ (a(k) - 1) .* (1 - x(k)) .^ (b(k) - 1) ./ beta(a(k),b(k));
     tmp(k) = (a(k) - 1).*log(x(k)) + (b(k) - 1).*log((1 - x(k))) - betaln(a(k),b(k));
     y(k) = exp(tmp(k));
end

function p = mytcdf(x,v)
%TCDF   Student's T cumulative distribution function (cdf).
%   P = TCDF(X,V) computes the cdf for Student's T distribution
%   with V degrees of freedom, at the values in X.
normcutoff = 1e7;
% Initialize P to zero.
p=zeros(size(x));
% use special cases for some specific values of v
k = find(v==1);
    % See Devroye pages 29 and 450.
    % (This is also the Cauchy distribution)
if any(k)
    p(k) = .5 + atan(x(k))/pi;
end
k = find(v>=normcutoff);
if any(k)
    p(k) = normcdf(x(k));
end
% See Abramowitz and Stegun, formulas 26.5.27 and 26.7.1
k = find(x ~= 0 & v ~= 1 & v > 0 & v < normcutoff);
if any(k),                            % first compute F(-|x|)
    xx = v(k) ./ (v(k) + x(k).^2);
    p(k) = betainc(xx, v(k)/2, 0.5)/2;
end
% Adjust for x>0.  Right now p<0.5, so this is numerically safe.
k = find(x > 0 & v ~= 1 & v > 0 & v < normcutoff);
if any(k), p(k) = 1 - p(k); end
p(x == 0 & v ~= 1 & v > 0) = 0.5;
% Return NaN for invalid inputs.
p(v <= 0 | isnan(x) | isnan(v)) = NaN;

function p = myfcdf(x,v1,v2)
%FCDF   F cumulative distribution function.
%   P = FCDF(X,V1,V2) returns the F cumulative distribution function
%   with V1 and V2 degrees of freedom at the values in X.
p = zeros(size(x)); %   Initialize P to zero.
t = (v1 <= 0 | v2 <= 0 | isnan(x));
p(t) = NaN;
% Compute P when X > 0.
k = find(x > 0 & ~t & isfinite(v1) & isfinite(v2));
if any(k), 
% use A&S formula 26.6.2 to relate to incomplete beta function 
    % Also use 26.5.2 to avoid cancellation by subtracting from 1
    xx = x(k)./(x(k) + v2(k)./v1(k));
    p(k) = betainc(xx, v1(k)/2, v2(k)/2);
end
if any(~isfinite(v1(:)) | ~isfinite(v2(:)))
   k = find(x > 0 & ~t & isfinite(v1) & ~isfinite(v2) & v2>0);
   if any(k)
      p(k) = mychi2cdf(v1(k).*x(k),v1(k));
   end
   k = find(x > 0 & ~t & ~isfinite(v1) & v1>0 & isfinite(v2));
   if any(k)
      p(k) = 1 - mychi2cdf(v2(k)./x(k),v2(k));
   end
   k = find(x > 0 & ~t & ~isfinite(v1) & v1>0 & ~isfinite(v2) & v2>0);
   if any(k)
      p(k) = (x(k)>=1);
   end
end

