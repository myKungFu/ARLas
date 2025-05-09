function [B1,B2,Bdp,residuals] = ARLas_dpoae_OLSfit_v1(y1c,y1s,y2c,y2s,ydpc,ydps,data,w)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [coeff,anova,modelSummary] = ARLas_dpoae_OLSfit_v1(X,y,alfa,w);
%
% Ordinary Least-Squares Regression for the model
% y = a + b1*X1 + b2*X2 + ... + bk*Xk.
% X = matrix of independent variables.
% y = vector of dependent variables.
% w = vector of weighting values (optional input)
%
% Author: Shawn Goodman, PhD
% Auditory Research Lab, the University of Iowa
% Date: March 7, 2022
% Last Updated: March 7, 2022 -- ssg
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
    
    B1 = b(1) + 1i*b(2);
    B2 = b(3) + 1i*b(4);
    try
        Bdp = b(5) + 1i*b(6);
    catch
        Bdp = [];
    end

end

