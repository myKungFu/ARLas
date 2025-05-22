function [B,yhat,residuals] = ARLas_wlsf(y,X,w)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% weighted least-squares fit
% WLSF Analysis
    % Ordinary Least-Squares Regression for the model
    % y = b1*X1 + b2*X2 + ... + bk*Xk.
    % X = matrix of independent variables (your "solution matrix" of signals of interest)
    % y = vector of dependent variables (your recorded data)

%    X = [xCos1000,xSin1000,xCos2000,xSin2000,dc]; % solution matrix

    [rows,columns] = size(X); % size of the solution matrix
    %k = columns-1; % k is the number of IV (the column of ones for y-intercept/dc offset is not counted)
    %n = rows; % number of observations
    sqrtw = sqrt(w); % and take the square root of the weighting vector
    X = repmat(sqrtw,1,columns).*X; % replicate sqrtw to matrix the size of X; J is the weighted version of X
    y = sqrtw.*y; % Dy is the weighted version of y
    % Solve the least squares problem
    [Q,R] = qr(X,0); % orthogonal-triangular decomposition; R is the Cholesky factor of the X matrix
                     % the input argument 0 makes an "economy sized" decomposition, so that
                     % [nSamples,nIVs] = size(Q), and 
                     % [nIvs,nIVs] = size(R).
    B = full(R\(Q'*y)); % Same as p = D*X\(D*y); b is a vector of coefficients; also same as b = (X'*X)\(X'*Y).
    % The next two lines are not used here, but may be useful for other model fitting
    yhat = X*B;      % Predicted responses at each data point.
    residuals = y - yhat;

end
