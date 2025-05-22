function [B] = ARLas_dpoae_OLSfit_v2(X,data,w)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [B] = ARLas_dpoae_OLSfit_v2(X,data,w);
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
% Last Updated: July 21, 2022 -- ssg
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % create a solution matrix.
    % rather than sine/cosine components that are pure sinusoids, use the
    % stimuli themselves
    
    X = [X,ones(size(X,1),1)];
    y = data;
    [rows,columns] = size(X); % size of the solution matrix
    %k = columns-1; % k is the number of IV (the column of ones for y-intercept/dc offset is not counted)
    %n = rows; % number of observations
    
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
    %yhat = X*b;      % Predicted responses at each data point.
    %residuals = y - yhat;
    
    b = b(1:end-1); % throw away dc offset
    % loop to put cosine and sine parts into complex form
    counter = 1;
    for ii=1:2:length(b)
        B(counter,1) = b(ii) + 1i*b(ii+1);
        counter = counter + 1;
    end
  end

