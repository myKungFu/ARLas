function [signal,nf,snr,noiseSamples,X,zbar] = LSF(Y,f,X)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [signal,nf,snr,noiseSamples,X,zbar] = LSF(Y,f,X);
%
% Required functions: bswSNR.m (bisquares weighted signal to noise ratio)
%
% Perform least squares fits on matrices of data.
% Y = matrix on which to perform LSFs. It is assumed that each column of Y
%       is a repeated recording (sweep/buffer/whatever). It is assmed that
%       a hann window has already been applied to Y.
% f = a scalar frequency (in Hz) to be fit.
% X = the solution matrix. If X is given as an input, it will be used. If X
%       is not given, it wil be created and returned.
%
% signal = mean magnitude, taken as abs of the complex mean across columns of Y (dB SPL)
% nf = noise floor, taken as the weighted standard error of the mean (dB SPL)
% snr = signal-to-noise ratio (dB)
% noiseSampes = a vector of complex coefficients, one from each column of
%       Y. The mean has been subtracted, leaving the variability around the mean.
%       This can be used in other code to esitmate the noise distribution of the
%       recordings.
% X = solution matrix
% zbar = complex mean
%
% Author: Shawn Goodman, PhD
% Auditory Research Lab, the University of Iowa
% Last Updated: February 1-7, 2025
% Last Updated: April 5, 2025 -- ssg -- now returns the complex mean zbar
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % create solution matrix for LSF
    fs = 96000; % sampling rate
    [frameN,M] = size(Y); % input size
    t = (0:1:frameN-1)'/fs; % time vector
    w = hann(frameN);

    if nargin < 3 % if solution matrix is not included as an input argument
        % create solution matrix
        nFreqs = length(f); % number of frequencies in the fit
        X = [];
        for ii=1:nFreqs
            c1 = cos(2*pi*f(ii)*t) .* w;
            s1 = -sin(2*pi*f(ii)*t) .* w;
            X = [X,c1,s1];
        end
        dc = ones(size(c1)); % add ones to account for dc offset
        X = [X,dc];
    end
    W = repmat(w,1,size(Y,2));
    Y = Y .* W;
    
    Z = zeros(M,1);
    for jj=1:M
        B = OLSfit_internal(X,Y(:,jj));
        Z(jj,1) = B(1) + 1i*B(2);
    end

    % compute weighted coherent mean, noise floor (standard error) and snr
    [signal,nf,zbar,snr,noiseSamples] = bswSNR(Z); 
  
end
% INTERNAL FUNCTIONS ------------------------------------------------------

function [b,yhat,residuals] = OLSfit_internal(X,y,w)
    % Ordinary Least-Squares Regression for the model
    % y = a + b1*X1 + b2*X2 + ... + bk*Xk.
    % X = matrix of independent variables.
    % y = vector of dependent variables.
    % w = vector of weighting values (optional input)
    if nargin == 2
        w = [];
    end
    [rows,columns] = size(X); % size of the solution matrix
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
end

