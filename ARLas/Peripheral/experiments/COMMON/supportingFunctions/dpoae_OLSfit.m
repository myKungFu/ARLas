function [B1,B2,Bdp,residuals] = dpoae_OLSfit(y1c,y1s,y2c,y2s,ydpc,ydps,data,w)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [coeff,anova,modelSummary] = OLSfit(X,y,alfa,w);
%
% Ordinary Least-Squares Regression for the model
% y = a + b1*X1 + b2*X2 + ... + bk*Xk.
% X = matrix of independent variables.
% y = vector of dependent variables.
% w = vector of weighting values (optional input)
%
% Author: Shawn Goodman
% Data: June 22, 2021
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% create a solution matrix.
% rather than sine/cosine components that are pure sinusoids, use the
% stimuli themselves

if isempty(ydpc)
    X = [y1c,y1s,y2c,y2s,ones(size(y1c))];
else
    X = [y1c,y1s,y2c,y2s,ydpc,ydps,ones(size(y1c))];
end
    y = data;


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

% yhat2 = X(:,1:4)*b(1:4);
% residuals2 = y - yhat2;
% 
% yhat3 = X(:,5:6)*b(5:6);
% residuals3 = y - yhat3;

B1 = b(1) + 1i*b(2);
B2 = b(3) + 1i*b(4);
try
    Bdp = b(5) + 1i*b(6);
catch
    Bdp = [];
end


% % Find a confidence interval for each component of X
% nu = rows-columns;                       % Residual degrees of freedom
% rmse = norm(residuals)/sqrt(nu);        % Root mean square error.
%        % norm() here is the L2 norm, the sqrt of the sum of the squared values
% RI = R\eye(columns); % the inverse of R
% se = rmse * sqrt(sum(abs(RI).^2,2)); % standard errors of the coefficients
%     % se = rmse ./ sqrt(sum(abs(R).^2,2)); % this is equivalent to the above
%     % for single regressor, this is equivalent to se = rmse / sqrt(sum((X(:,2)-xbar).^2)).
%     % for bivariate regressors, this is equivalent to se = rmse / sqrt(sum((X(:,2)-xbar).^2)*r2x2x3).
% t = b ./ se; % t-values of the coefficients
% bSig = zeros(k,1); % significance values of coefficients via t-test (2-tailed)
% for ii=1:length(t),
%     bSig(ii) = (1-mytcdf(t(ii),n-k-1))*2; % get the p-value, 2-tailed 
%     %bSig(ii) = (1-mytcdf(t(ii),n-k-1)); % get the p-value, 1-tailed 
% end
% 
% criticalT = mytinv(1-alfa/2,n-2); % critical two-tailed t-value
% ci = [b-(criticalT*se), b+(criticalT*se)]; % confidence intervals on coefficients
% SSreg = norm(yhat-mean(y))^2;  % Regression sum of squares.
% SStotal = norm(y-mean(y))^2;   % Total sum of squares.
% MSreg = SSreg / k; % mean square regression
% SSres = SStotal - SSreg; % sum of squares residual
% MSres = SSres / (n-k-1); % mean square residual
% if MSres == 0,
%     F = 0;
%     p = 0;
%     R2 = 0;
%     R2adj = 0;
% else
%     F = MSreg / MSres; % F ratio
%     p = 1 - myfcdf(F,k,n-k-1); % significance of the F ratio
%     R2 = SSreg / SStotal;          % multiple R-square statistic.
%     R2adj = 1- ((1 - R2)*((n-1)/(n-k-1))); % adjusted R squared; Calculation from Wherry (1931), cited in T&F, p 147.
% end
% 
% % put output into three different structures:
% modelSummary.R = sqrt(R2);
% modelSummary.R2 = R2;
% modelSummary.R2adj = R2adj;
% modelSummary.rmse = rmse;
% 
% anova.dfReg = k;
% anova.dfRes = n-k-1;
% anova.dfTot = n-1;
% anova.SSreg = SSreg;
% anova.SSres = SSres;
% anova.SStot = SStotal;
% anova.MSreg = MSreg;
% anova.MSres = MSres;
% anova.F = F;
% anova.p = p;
% 
% coeff.b = b;
% coeff.se = se;
% coeff.t = t;
% coeff.p = bSig;
% coeff.ci = ci;


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


