function [jar,h1,h2,h3,h4,h5] = ricianLGF_fit2(Z,Sz,X)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [jar,h1,h2,h3] = ricianLGF_fit2(Z,Sz,X);
% [jar,h1,h2,h3,h4,h5] = ricianLGF_fit2;
%
% Fit DPOAE level growth functions (LGFs) using a Rician-based statistcal model.
% The growth is modeled as a quadradic polynomial. 
% Curve Fitting toolbox needs to be installed, as well as Machine Learning
% and Statistics toolbox.
%
% If function is called without input arguments, a simulated LGF will be
%   genenerated based on a generalized logisitic nonlinearity.
%   If called without input arguments, you should have the Machine Learning
%   and Statistics toolbox installed, because the simulated data are
%   created using makedist.m. If you don't, then a fixed value will be used
%   instead.
% If function is called with input arguments, the following should be used:
%   Z = measured magnitude (dB SPL)
%   Sz = standard error of the mean (dB SPL)
%   X = L2 stimulus level (dB SPL)
% IMPORTANAT NOTE: If giving input arguments, you must cut these down so
% they only go up to the peak of the LGF!
%
% Output Arguments: jar is a structure with the following fields
% mle: maximum likelihood estimates of the coefficients for the squared and linear terms (offset is fixed at zero)
% B: the second-pass grid search for the two coefficients 
% Partials: log likelihoods for the two coefficients alone
% CI: confidence intervales for mle based on the second partial derivatives
% dydx = first derivative of fit
% d2ydx2 = second derivative of fit
% yhat_mle: best model fit based on mle (in mPa)
% Yhat_mle: best model fit based on mle (in dB SPL)
% Slope_mle: slope of the best model fit (dB)
% slope_mle: slope of the best model fit (linear)
% SlopeM: assuming no input, true slope of the simulated DP (dB) 
% slopem: assuming no input, true slope of the simulated DP (linear)
% Epsilon: assuming no input, error of the model fit (dB difference)
% epsilon: assuming no input, error of the model fit (linear ratio)
% EpsilonD: assuming no input, error of the model slope (dB difference)
% epsilond: assuming no input, error of the model slope (linear ratio)
% ref: pressure reference (fixed at 0.02, or 20 mPa)
% X: x-axis, L2 target levels (dB FPL)
% x: x-axis, L2 target levels (mPa)
% Z: y-axis, DPOAE LGF values (dB SPL)
% z: y-axis, DPOAE LGF values (mPa)
% Sz: noise floor, standard error of the mean (dB SPL)
% sz: noise floor, standard error of the mean (mPa)
% Esz: expectation (mean) of noise floor
% Zm: assuming no input, simulated DPOAE levels (dB SPL)
% zm: assuming no input, simulated DPOAE levels (mPa)
%
% Author: Shawn Goodman
% Date: March 27, 2025
% Updated: April 10-12, 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    jar = []; h1 = []; h2 = []; h2 = []; h4 = []; h5 = []; % initialize outputs
    if nargin == 0 % if no inputs, simulate from a generalized logistic model
        [Z,Sz,X,Zm,Q5,Q8,Q11] = simulateDataNL;
        % cut simulated data so only fitting below the peak of the LGF
        cut1 = 1; % start at sample 1
        cut2 = 60; % just below peak
        Zm = Zm(cut1:cut2);
        Z = Z(cut1:cut2);
        Sz = Sz(cut1:cut2);
        X = X(cut1:cut2);
    end

upSample = 100;
if upSample > 1
    originalN = length(Z);
    NN = originalN * upSample;
    XX = linspace(X(1),X(end),NN)';
    Z = interp1(X,Z,XX,'pchip');
    Sz = interp1(X,Sz,XX,'pchip');
    X = XX;
end

    % Input is in dB SPL. Convert to linear units (mPa)
    ref = 0.02; % reference is 20 mPa (not micro! ricean model expects mPa)
    sz = 10.^(Sz/20)*ref; % standard error of the mean vector
    Esz = mean(sz); % expectation (mean) of standard error of the mean vector
    z = 10.^(Z/20)*ref; % DPOAE levels (mPa)
    if nargin == 0
        zm = 10.^(Zm/20)*ref; % simulated noiseless DPOAE levels (mPa)
    end
    x = 10.^(X/20)*ref; % target L2 stimulus level (mPa)
    
    % % find initial guess of coefficients using standard LSF -------------
    beta = initialize(Z,Sz,x,z);
    % Set up grid to search over:
    b2 = makeB(beta(1));
    b1 = makeB(beta(2));
    B = [b2,b1]; % matrix of coefficients to search
    % run search ----------------------------------------------------------
    jar = searchME(B,x,Esz,z); % first pass
    jar2 = searchME(jar.B,x,Esz,z); % second pass
    jar2.B = jar.B; % put back the old coefficient search
    jar = jar2;
    % end of search -------------------------------------------------------

    % calculate derivatives
    jar.dydx = polyder([jar.mle,0]); % first derivative of fit (the slope)
    jar.d2ydx2 = polyder(jar.dydx); % second derivative of fit (single value)

    % calculate model output
    yhat_mle = polyval([jar.mle,0],x); % mle model fit (mPa)
    Yhat_mle = 20*log10(yhat_mle/ref); % mle model fit (dB SPL)
    Slope_mle = (gradient(Yhat_mle)./gradient(X)); % slope of the dB fit
    slope_mle = (gradient(yhat_mle)./gradient(x)); % slope of the linear fit
    if nargin == 0
        SlopeM = (gradient(Zm)./gradient(X)); % slope of the simulation in dB
        slopem = (gradient(zm)./gradient(x)); % slope of the simulation linear
    else
        SlopeM = [];
        slopem = [];
    end
    % package values into output structure
    jar.yhat_mle = yhat_mle; 
    jar.Yhat_mle = Yhat_mle;
    jar.Slope_mle = Slope_mle;
    jar.slope_mle = slope_mle;
    jar.SlopeM = SlopeM;
    jar.slopem = slopem;

    if nargin == 0 % if simulated data were used, can calculate errors
        % calculate error in the model and the underlying (synthetic) DPOAE
        Epsilon = Yhat_mle - Zm; % model error (model dB - simulated dB)
        epsilon = yhat_mle ./ zm; % model error (model mPa ./ simulated mPa)
        % calculate error in the model slope and the underlying (synthetic) DPOAE slope
        EpsilonD = Slope_mle - SlopeM; % model slope error (dy/dx - dy/dx, both in dB)
        epsilond = slope_mle ./ slopem; % model slope error (dy/dx - dy/dx, both in mPa)
        jar.Epsilon = Epsilon; 
        jar.epsilon = epsilon; 
        jar.EpsilonD = EpsilonD; 
        jar.epsilond = epsilond; 
    else
        jar.Epsilon = []; 
        jar.epsilon = []; 
        jar.EpsilonD = []; 
        jar.epsilond = []; 
    end        

    % pack results into a strucuture
    jar.ref = ref; % pressure reference (mPa)
    jar.X = X; % L2 stimulus levels (dB SPL)
    jar.x = x; % L2 stimulus levels (mPa)
    jar.Z = Z; % measured DPOAE magnitudes (dB SPL)
    jar.z = z; % measured DPOAE magnitudes (mPa)
    jar.Sz = Sz; % estimated noise floor, standard error of the mean (dB SPL)
    jar.sz = sz; % standard error (mPa)
    jar.Esz = Esz; % expected (mean) value of standard error
    if nargin == 0
        jar.Zm = Zm; % simulated DPOAE magnitudes (dB SPL, no noise)
        jar.zm = zm; % simulated DPOAE magnitudes (mPa, no noise)
    else
        jar.Zm = []; % simulated DPOAE magnitudes (dB SPL, no noise)
        jar.zm = []; % simulated DPOAE magnitudes (mPa, no noise)
    end

% % get low level error
% [~,indx] = min(abs(X - 5));
% jar.Epsilon(indx) % positive value says predicted is > actual
% 0.7945, 0.4856, 0.7487, 0.3829, 0.7914, 0.6417, 0.6821, 0.6328, 0.5676


    % plot results --------------------------------------------------------
    h1 = plotLL(jar);
    h2 = plotFit(jar);
    h3 = plotSlope(jar);
    if nargin == 0 % can only plot errors if simulated data were used
        h4 = plotSlopeError(jar);
        h5 = plotFitError(jar);
    end

    % keyboard
    % plot(5,Q5.ZnLL,'pb')
    % hold on
    % plot(5,Q5.NfnLL,'k.','MarkerSize',12)
    % plot(8,Q8.ZnLL,'pb')
    % plot(8,Q8.NfnLL,'k.','MarkerSize',12)
    % plot(11,Q11.ZnLL,'pb')
    % plot(11,Q11.NfnLL,'k.','MarkerSize',12)

end

% INTERNAL FUNCTIONS ------------------------------------------------------
function [beta] = initialize(Z,Sz,x,z)
    % find initial "ballpark" fit with standard LSF from curve fitting toolbox
    SNR = Z - Sz; % stimate SNR
    N = length(SNR);
    % assume the best SNR is at the largest L2 (right side of the LGF).
    % To find the cutoff, "Count backwards" until you fall into the noise floor
    counter = N;
    done = 0;
    cut = 3; % look for values above this SNR cutoff
    while done == 0
        if SNR(counter) > cut
            counter = counter-1;
        else
            done = 1;
        end
        if counter <= 1
            done = 1;
        end
    end
    % fit with a quadradic polynomial, bisquares weighted, offset fixed at zero
    [fitresult,gof] = createFit2(x(counter:end),z(counter:end));
    %ci = confint(fitresult); % get confidence intervals
    beta = [fitresult.p1,fitresult.p2]; % best initial fit, with noise, with no offset term
end
function [B] = makeB(b)
    % create the search vector of coefficients
    nLooks = 50; % set the vector length
    B = sort([b*0.8,b*1.2]); % the extremes are +/-20% of the initial guess
    B = linspace(B(1),B(2),nLooks)';
end
function [jar] = searchME(B,x,Esz,z)
    % find the maximum likelihood estimates and confidence intervals
    LL = getLLfunction2(B,x,Esz,z); % get the log likelihood function
    ml = max(LL(:)); % maximum log likelihood
    [mle,B,Partials,CI,redo,d22,d21] = getPartials(LL,B,ml); % get mles and partials
    if redo == 1 % if not centered, do it again
        LL = getLLfunction2(B,x,Esz,z); % get the log likelihood function
        ml = max(LL(:)); % maximum log likelihood
        [mle,B,Partials,CI,redo,d22,d21] = getPartials(LL,B,ml); % get mles and partials
        if redo == 1
            disp('Not centered a second time!')
            % keyboard % if not centered on second pass, we will need to look into 
                     % why this is still happening and implement a fix
        end
    end
    jar.mle = mle; % current maximum likelihood estimate
    jar.B = B; % new set of coefficients to search over
    jar.Partials = Partials; % slices of LL values for each coefficient
    jar.CI = CI; % confidence intervals for mle (second derivative of LL)
    jar.d2 = [d22,d21];
end
function [mle,B,Partials,CI,redo,d22,d21] = getPartials(LL,B,ml)
    % extract the likelihood functions for each coefficient. These are
    % vectors, hence  "partials" of the log likelihood matrix
    mle = []; % initialize outputs
    CI = [];
    redo = 0; % switch for redo if not centered

    [rr,cc] = find(LL==ml); % location of maximum likelihood
    partial1 = LL(rr,:); % partial is a "slice" of the LL matrix
    partial2 = LL(:,cc);
    beta1 = B(:,2); % vector of coefficients corresponding to the LL partial slice
    beta2 = B(:,1);    
    mle = [beta2(rr),beta1(cc)]; % current maximum likelihood estimate vector
    % for better accuracy, fit the tip of the partial with a quadradic, 
    % solve for root, get new constrained values of B to be used in next pass
    nLooks = length(partial2);
    Partials = [partial2(:),partial1(:)];
    % check to see whether the search was over an appropriate range.
    % If LL is not a quadradic, then the search range was wrong. Re-center
    % and redo
    [bhat,ci,mlhat,partialx,partialy,betaNew,d22] = polyme(partial2,beta2);
    if ~isempty(betaNew)
        b2 = betaNew;
        redo = 1;
        d22 = NaN;
    else
        b2 = linspace(partialx(1),partialx(end),nLooks)'; % quadradic
        mle(1) = bhat;
        CI(:,1) = ci;
    end
    [bhat,ci,mlhat,partialx,partialy,betaNew,d21] = polyme(partial1,beta1);
    if ~isempty(betaNew)
        b1 = betaNew;
        redo = 1;
        d21 = NaN;
    else
        b1 = linspace(partialx(1),partialx(end),nLooks)'; % linear
        mle(2) = bhat;
        CI(:,2) = ci;
    end
    B = [b2,b1]; % new matrix of coefficients to search over for next pass
end
function [LL] = getLLfunction2(B,x,Esz,z)
    % calculate the log likelihood function. This is a 2-D matrix
    N = size(B,1);
    LL = zeros(N,N); % initialized the log-likelihood matrix
    for ii=1:N % loop over B2 coefficients
        b2 = B(ii,1);
        for jj=1:N % loop over B1 coefficients
            b1 = B(jj,2);
            beta = [b2,b1,0]; % add the zero for fixed y-offset
            LL(ii,jj) = calcLL(beta,x,Esz,z); % calculate the log likelhood for current coefficients
        end
    end
end
function [LL] = calcLL(beta,x,Esz,z)
    % calculate the log likelihood via the Rician model
    y = polyval(beta,x); % apply the polynomial growth model
    indx = y<=0; % find any values less than zero and make them NaN 
                 % (because you can't take the log of them and they are
                 % nonsensical values)
    y(indx) = NaN;
    SNR = 20*log10(y ./ Esz); % estimated signal to noise ratio
    k = g(SNR); % calculate the k multiplier (see paper)
    q = h(SNR); % calculate the q multiplier (see paper)
    nu = k.*y; % estimate rician nu (location) parameter
    sigma = q.*Esz; % esimate rician sigma (spread) parameter
    % sove the rician PDF value in two stages
    part1 = (z./sigma.^2).*exp(-((z.^2+nu.^2)./(2*sigma.^2)));
    part2 = besseli(0,(z.*nu)./sigma.^2);
    % part 3 is the normal approximation for high SNRs
    part3 = (1./(sigma*sqrt(2*pi))) .* exp(-((z-nu).^2)./(2*sigma.^2)); % normal approximation
    % find SNR cutoffs for using the normal approximation
    r = nu./sigma; % ratio of nu to sigma
    cut = 10^(15/20); % cutoff for normal approximation is 15 dB "SNR"
    indx3 = find(r>= cut);
    % combine the parts to get the likelihood function
    f = part1 .* part2;
    f(indx3) = part3(indx3);
    % computer the log likelihood
    f = log(f);
    % replace any NaN or inf with zero
    f(isinf(f)) = 0;
    f(isnan(f)) = 0;
    LL = sum(f); % the log likelihood value
end
function [bhat,ci,ml,partialx,partialy,betaNew,d2] = polyme(partial,beta)
    % make polynomial fits to the partial LL vectors to find the best
    % mle estimate, create confidence intervals, and create new proposed
    % vectors of coefficients to search across
    bhat = []; % initialize outputs
    ci = [];
    ml = [];
    partialx = [];
    partialy = [];
    betaNew = [];
    d2 = [];

    N = length(partial);
    [~,indx] = max(partial); % location of the max LL
    direct = 'n'; % default to no directon
    % check to see if the range was adequate
    % fit the partial with a polynomial, caculate the peak and derivative
    nnn = 8; % number of samples to use either side of the peak
    if indx + nnn >= N
        direct = 'r';
    end
    if indx - nnn <= 1
        direct = 'l';
    end
    if ~strcmp(direct,'n') % if there is a range problem, suggest the solution
        dd = abs(median(diff(beta)));
        if strcmp(direct,'r')
            center = beta(N); 
            start = beta(N/2);
            finish = center + ((N/2)*dd);
            betaNew = linspace(start,finish,N)'; % new centered beta
            return
        elseif strcmp(direct,'l')
            center = beta(1);
            finish = beta(N/2);
            start = center - ((N/2)*dd);
            betaNew = linspace(start,finish,N)'; % new centered beta
            return
        else
            warning('Unexpected variable')
            keyboard
        end
    end
    yyy = partial(indx-nnn:indx+nnn);
    xxx = beta(indx-nnn:indx+nnn); 
    pp = polyfit(xxx,yyy,2);
    ppd = polyder(pp);
    r = roots(ppd);
    if length(r) > 1 % if there is more than one root to choose from
        keyboard
    end
    ml = polyval(pp,r);

    d2 = polyder(ppd); % second derivative of fit
    ci(1,1) = r - 1.96*(1./sqrt(-d2)); % lower bound of CI was 1.96
    ci(2,1) = r + 1.96*(1./sqrt(-d2)); % upper bound of CI 
    bhat = r;
    partialx = linspace(xxx(1),xxx(end),200);
    partialy = polyval(pp,partialx);
end
function [k,K1,K2,K3,SNR1,SNR2,SNR3] = g(SNR)
    % find k multiplier for rician fit
    % if SNR <-12, use polynomial:
    % p1*x + p2
    p1 = -0.997;
    p2 = 23.871;
    % if SNR >=-12 and <=-1, use 2 gaussians:
    % a1*exp(-((x-b1)/c1)^2) + a2*exp(-((x-b2)/c2)^2)
    a1 =  44.870;
    b1 =  -24.644;
    c1 =  23.613;
    a2 =  22.403;
    b2 =   4.522;
    c2 =  11.138;
    % if SNR > -1, use a constant:
    % k = 33.979255836496655 (in dB SPL),
    % which is 1 in linear units.

    [indx1] = find(SNR<-12);
    [indx2] = find(SNR>=-12 & SNR<=-1);
    [indx3] = find(SNR>-1);
    
    K1 = p1*SNR + p2;
    k1 = 10.^(K1/20)*0.02;

    K2 = a1*exp(-((SNR-b1)./c1).^2) + a2*exp(-((SNR-b2)./c2).^2);
    k2 = 10.^(K2/20)*0.02;

    k3 = ones(size(SNR));
    K3 = 20*log10(k3/.02);

    k = zeros(size(SNR));
    k(indx1) = k1(indx1);
    k(indx2) = k2(indx2);
    k(indx3) = k3(indx3);

    K1 = K1(indx1);
    K2 = K2(indx2);
    K3 = K3(indx3);

    SNR1 = SNR(indx1);
    SNR2 = SNR(indx2);
    SNR3 = SNR(indx3);
end
function [q,q1,q2,q3,SNR1,SNR2,SNR3] = h(SNR)
    % find q muliplier for rician fit
    % if SNR <-8, use a double exponential:
    a1 = 0.66421;
    b1 = -5.82103e-05;
    c1 = 0.14825;
    d1 = 0.16975;

    % if from -8 <= SNR < 3, use a rational function of degree 2,2
    % (p1*x^2 + p2*x + p3) / (x^2 + q1*x + q2)
    p1 = 0.70632;
    p2 = 10.08339;
    p3 = 40.03945;
    q1 = 14.21303;
    q2 = 56.21896;

    % if SNR > 3, use a double exponential:
    % a*exp(b*x) + c*exp(d*x)
    a2 = 6.98752e-04;
    b2 = -0.60987;
    c2 = 0.71056;
    d2 = -1.26337e-05;
    
    % try it out:
    [indx1] = find(SNR<=-8);
    [indx2] = find(SNR>=-8 & SNR<=3);
    [indx3] = find(SNR>=3);

    k1 = a1*exp(b1*SNR) + c1*exp(d1*SNR);
    k2 = (p1.*SNR.^2 + p2.*SNR + p3) ./ (SNR.^2 + q1.*SNR + q2);
    k3 = a2.*exp(b2.*SNR) + c2.*exp(d2.*SNR);
    
    q = zeros(size(SNR));
    q(indx1) = k1(indx1);
    q(indx2) = k2(indx2);
    q(indx3) = k3(indx3);

    q1 = k1(indx1);
    q2 = k2(indx2);
    q3 = k3(indx3);
    SNR1 = SNR(indx1);
    SNR2 = SNR(indx2);
    SNR3 = SNR(indx3);
end
function [Z,Sz,X,Zm,Q5,Q8,Q11] = simulateDataNL()
    Q5 = []; Q8 = []; Q11 = [];
    % Create a simulated LGF by passing primaries through a generalized logistic function.
    fs = 96000; % sampling rate
    ref = .002; % pressure reference for hyperbolic tan, 0.002;
    dur = 1; % input signal duration
    N = round(fs*dur); % number of samples in input
    t = (0:1:N-1)'/fs; % time vector
    f2 = 1000; % f2 frequency (Hz)
    fratio = 1.22;
    f1 = round(f2/fratio); % f1 frequency (Hz)
    fdp = 2*f1-f2; % DP frequency (Hz)
    L1 = 55; % level of f1. Yes, should ideally be 65, but to get results that look
             % like our data, you have to use this lower level
    L2 = (0:1:70)'; % vector of f2 levels
    % create the primary freuqencies at unit amplitude
    p1 = cos(2*pi*f1*t); % primary 1
    p2 = cos(2*pi*f2*t); % primary 2
    nn = length(L2); % number of L2 levels to consider
    
    % B is the growth rate of the function, the "gain".
    % lower growth moves peak to the right
    % change the slope in dB, but also leaves the upper value higher if
    % reduced. 
    % Set B outside of the loop if you want gain to be fixed
    % Here, let B be chosen randomly from a reasonable distribution
    try
        pd = makedist('norma','mu',0.75,'sigma',.1); % make a normal distribution
        B = pd.random(1,1);
    catch % if the user doesn't have the machine learning and stats toolbox installed
        B = 0.75; % just use the mean value
    end
B = ones(nn,1)*.85;    
%B = linspace(1.1,.85,nn)';

    for ii=1:nn % loop over test L2 levels
        A = -1; % left horizontal assymtote
        K = -A; % right horizontal assymtote
           % for these two, it doesn't matter which changes; the magnitude that
           % comes out is the same, you can reduce these to a single parameter.
        nu = 1; % >0, where growth occurs
        Q = 1;
        
        a1 = 10^(L1/20)*ref;
        a2 = 10^((L2(ii))/20)*ref;
        input = a1*p1 + a2*p2;
        
        nSweeps = 24; % this chosen to mirror what we actually collect
        % Note: you have a choice whether to add noise before the
        % nonlinearity (in which case it represents cochlear noise) or
        % after the nonlinearity (in which case it represents ear canal
        % noise). Here, we add nosie after the nonlinearity under the
        % assumption that the majority of noise is ear canal. No proof of
        % this though!
        % Note that we are also createding a noiseless version so we can
        % have the "true" answer about the LGF
        for jj=1:nSweeps
            % q = input + randn(size(input)) * 0; % *0
            % Z(:,jj) = generalizedLogistic(q,A,K,B,nu,Q); % create a noiseless version
            % Zn(:,jj) =  Z(:,jj) + randn(size(input))*.2; % create a version with noise

            q = input + randn(size(input)) * 0; % *0
            Z(:,jj) = generalizedLogistic(q,A,K,B(ii),nu,Q); % create a noiseless version
            Zn(:,jj) =  Z(:,jj) + randn(size(input))*.2; % create a version with noise


        end
        % having passed the signal through the nonlinearity, look at the
        % cubic distortion
        [frequency,signal,noiseFloor] = ARLas_fda(Z,fs,ref); % the noiseless version
        [~,indx] = min(abs(frequency-fdp)); % location of 2f1-f2
        Out(ii,1) = signal(indx); % Ldp in dB SPL
        Nf(ii,1) = noiseFloor(indx); % noise floor (standard error of mean)
        [frequency,signal,noiseFloor] = ARLas_fda(Zn,fs,ref); % version with noise
        Outn(ii,1) = signal(indx); % Ldp in dB SPL (noisy version)
        Nfn(ii,1) = noiseFloor(indx); % noise floor (noisy version)
    end
    % rename variables for output
    X = L2; % target L2 levels (dB FPL)
    Sz = Nfn; % noise floor (dB SPL), noisy version
    Z = Outn; % Ldp (dB SPL), noisy version
    Zm = Out; % Ldp (dB SPL), noiseless "model" version


    return
    % was using 0.19 for both
    a2 = 10^(5/20)*ref;
    [~,indx2] = min(abs(X-5));
    input = a1*p1 + a2*p2;
    nn = 2048;
    ZLL = zeros(fs,nn);
    for jj=1:nn
        q = input + randn(size(input)) * 0; % *0
        ZLL(:,jj) = generalizedLogistic(q,A,K,B(indx2),nu,Q); % create a noiseless version
        ZnLL(:,jj) =  ZLL(:,jj) + randn(size(input))*.2; % create a version with noise
    end
    [frequency,signal,noiseFloor] = ARLas_fda(ZLL,fs,ref); % the noiseless version
    ZLL = signal(indx); % Ldp in dB SPL
    NfLL = noiseFloor(indx); % noise floor (standard error of mean)
    [frequency,signal,noiseFloor] = ARLas_fda(ZnLL,fs,ref); % version with noise
    ZnLL = signal(indx); % Ldp in dB SPL (noisy version)
    NfnLL = noiseFloor(indx); % noise floor (noisy version)
    Q5.ZLL = ZLL;
    Q5.NfLL = NfLL;
    Q5.ZnLL = ZnLL;
    Q5.NfnLL = NfnLL;

    a2 = 10^(8/20)*ref;
    [~,indx2] = min(abs(X-8));
    input = a1*p1 + a2*p2;
    clear ZLL ZnLL
    nn = 2048;
    ZLL = zeros(fs,nn);
    for jj=1:nn
        q = input + randn(size(input)) * 0; % *0
        ZLL(:,jj) = generalizedLogistic(q,A,K,B(indx2),nu,Q); % create a noiseless version
        ZnLL(:,jj) =  ZLL(:,jj) + randn(size(input))*.2; % create a version with noise
    end
    [frequency,signal,noiseFloor] = ARLas_fda(ZLL,fs,ref); % the noiseless version
    ZLL = signal(indx); % Ldp in dB SPL
    NfLL = noiseFloor(indx); % noise floor (standard error of mean)
    [frequency,signal,noiseFloor] = ARLas_fda(ZnLL,fs,ref); % version with noise
    ZnLL = signal(indx); % Ldp in dB SPL (noisy version)
    NfnLL = noiseFloor(indx); % noise floor (noisy version)
    Q8.ZLL = ZLL;
    Q8.NfLL = NfLL;
    Q8.ZnLL = ZnLL;
    Q8.NfnLL = NfnLL;

    a2 = 10^(11/20)*ref;
    [~,indx2] = min(abs(X-11));
    input = a1*p1 + a2*p2;
    clear ZLL ZnLL
    nn = 2048;
    ZLL = zeros(fs,nn);
    for jj=1:nn
        q = input + randn(size(input)) * 0; % *0
        ZLL(:,jj) = generalizedLogistic(q,A,K,B(indx2),nu,Q); % create a noiseless version
        ZnLL(:,jj) =  ZLL(:,jj) + randn(size(input))*.2; % create a version with noise
    end
    [frequency,signal,noiseFloor] = ARLas_fda(ZLL,fs,ref); % the noiseless version
    ZLL = signal(indx); % Ldp in dB SPL
    NfLL = noiseFloor(indx); % noise floor (standard error of mean)
    [frequency,signal,noiseFloor] = ARLas_fda(ZnLL,fs,ref); % version with noise
    ZnLL = signal(indx); % Ldp in dB SPL (noisy version)
    NfnLL = noiseFloor(indx); % noise floor (noisy version)
    Q11.ZLL = ZLL;
    Q11.NfLL = NfLL;
    Q11.ZnLL = ZnLL;
    Q11.NfnLL = NfnLL;


end
function [Y] = generalizedLogistic(x,A,K,B,nu,Q)
    Y = A + ( (K-A) ./ ((1+Q.*exp(-B*x)).^(1/nu)) );
    %Q = 1;
    %C = 1;
    %Y = A + ((K-A)./((C+Q.*exp(-B*x))).^(1/nu));
    %Y = A + ( (K-A) ./ ((1+exp(-B*x)).^(1/nu)) );
end
function [fitresult,gof] = createFit2(x,z)
    % find best fitting quadradic polynomial with a weighted LSF
    % this requires the curve fitting toolbox
    [xData, yData] = prepareCurveData( x, z );
    ft = fittype( 'poly2' );
    opts = fitoptions( 'Method', 'LinearLeastSquares' );
    opts.Lower = [-Inf 0 0];
    opts.Robust = 'Bisquare';
    opts.Upper = [Inf 1 0];
    [fitresult, gof] = fit( xData, yData, ft, opts );
end
function [h] = plotSlope(jar)
    % plot the slopes
    h = figure;
    subplot(2,1,1)
     if ~isempty(jar.SlopeM)
         plot(jar.X,jar.SlopeM,'Color',[.7 .7 .7],'LineWidth',1)
     end
     hold on
     plot(jar.X,jar.Slope_mle,'r')
     xlim([jar.X(1),jar.X(end-1)])
     xlabel('L2 (dB)')
     ylabel('Slope (dB)')
     if ~isempty(jar.SlopeM)
        legend('true','model fit','Location','southwest')
     end
    subplot(2,1,2)
     if ~isempty(jar.slopem)
         plot(jar.x,jar.slopem,'Color',[.7 .7 .7],'LineWidth',1)
     end
     hold on
     plot(jar.x,jar.slope_mle,'r')
     xlim([jar.x(1),jar.x(end-1)])
     xlabel('L2 (mPa)')
     ylabel('slope (ratio)')
     if ~isempty(jar.slopem)
         legend('true','model fit','Location','southwest')
     end
end
function [h] = plotSlopeError(jar)
    h = figure;
    subplot(2,1,1)
    plot(jar.X,jar.EpsilonD,'m')
    hold on
    line([jar.X(1),jar.X(end)],[0,0],'Color',[0 0 0],'LineWidth',0.5,'LineStyle',':')
    xlim([jar.X(1),jar.X(end-1)])
    xlabel('L2 (dB)')
    ylabel('Error in model slope (dB)')
    subplot(2,1,2)
    plot(jar.x,jar.epsilond,'m')
    hold on
    line([jar.x(1),jar.x(end)],[1,1],'Color',[0 0 0],'LineWidth',0.5,'LineStyle',':')
    xlim([jar.x(1),jar.x(end-1)])
    xlabel('L2 (mPa)')
    ylabel('Error in model slope (ratio)')
end
function [h] = plotFitError(jar)
    h = figure;
    subplot(2,1,1)
    plot(jar.X,jar.Epsilon,'m')
    hold on
    line([jar.X(1),jar.X(end)],[0,0],'Color',[0 0 0],'LineWidth',0.5,'LineStyle',':')
    xlabel('L2 (dB)')
    ylabel('Error in model (dB)')
    subplot(2,1,2)
    plot(jar.x,jar.epsilon,'m')
    hold on
    line([jar.x(1),jar.x(end)],[1,1],'Color',[0 0 0],'LineWidth',0.5,'LineStyle',':')
    xlabel('L2 (mPa)')
    ylabel('Error in model (ratio)')
end
function [h] = plotFit(jar)
    h = figure;
    subplot(2,1,1)
     hold on
     if ~isempty(jar.Zm)
         plot(jar.X,jar.Zm,'Color',[.7 .7 .7],'LineWidth',2)
     end
     plot(jar.X,jar.Z,'b')
     plot(jar.X,jar.Sz,'k')
     plot(jar.X,jar.Yhat_mle,'r-')
     xlim([jar.X(1),jar.X(end)])
     xlabel('Target L2 (dB FPL)')
     ylabel('Ldp (dB SPL)')
     grid on
     if ~isempty(jar.Zm)
         legend('true DP','measured DP','noise floor','model fit','Location','northwest')
     else
         legend('measured DP','noise floor','model fit','Location','northwest')
     end
    subplot(2,1,2)
     hold on
     if ~isempty(jar.zm)
         plot(jar.x,jar.zm,'Color',[.7 .7 .7],'LineWidth',2)
     end
     plot(jar.x,jar.z,'b')
     plot(jar.x,jar.sz,'k')
     plot(jar.x,jar.yhat_mle,'r-')
     xlim([jar.x(1),jar.x(end)])
     xlabel('Target L2 (mPa)')
     ylabel('Ldp (mPa)')
     grid on
     if ~isempty(jar.Zm)
         legend('true DP','measured DP','noise floor','model fit','Location','northwest')
     else
         legend('measured DP','noise floor','model fit','Location','northwest')
     end
end
function [h] = plotLL(jar)
    h = figure;
    subplot(1,2,1)
     ymin = min(jar.Partials(:,1));
     ymax = max(jar.Partials(:,1));
     ymax = ymax + ((ymax-ymin)*.02);
     xmin = min([min(jar.CI(:,1)),min(jar.B(:,1))]);
     xmax = max([max(jar.CI(:,1)),max(jar.B(:,1))]);
     plot(jar.B(:,1),jar.Partials(:,1),'k');
     hold on
     plot(jar.mle(1),max(jar.Partials(:,1)),'r*')
     title([num2str(jar.mle(1)),'x^2'])
     xlabel('\beta_2')
     ylabel('Log Likelihood')
     grid on
     line([jar.CI(1,1),jar.CI(1,1)],[ymin,ymax],'Color',[1 0 0],'LineStyle','--','LineWidth',0.5)
     line([jar.CI(2,1),jar.CI(2,1)],[ymin,ymax],'Color',[1 0 0],'LineStyle','--','LineWidth',0.5)
     ylim([ymin,ymax])
     xlim([xmin,xmax])
    subplot(1,2,2)
     ymin = min(jar.Partials(:,2));
     ymax = max(jar.Partials(:,2));
     ymax = ymax + ((ymax-ymin)*.02);
     xmin = min([min(jar.CI(:,2)),min(jar.B(:,1))]);
     xmax = max([max(jar.CI(:,2)),max(jar.B(:,1))]);
     plot(jar.B(:,2),jar.Partials(:,2),'k');
     hold on
     plot(jar.mle(2),max(jar.Partials(:,2)),'r*')
     title([num2str(jar.mle(2)),'x'])
     xlabel('\beta_1')
     ylabel('Log Likelihood')
     grid on
     line([jar.CI(1,2),jar.CI(1,2)],[ymin,ymax],'Color',[1 0 0],'LineStyle','--','LineWidth',0.5)
     line([jar.CI(2,2),jar.CI(2,2)],[ymin,ymax],'Color',[1 0 0],'LineStyle','--','LineWidth',0.5)
     xlim([min(jar.B(:,2)),max(jar.B(:,2))])
     ylim([ymin,ymax])
end

% OLD CODE ----------------------------------------------------------------

% Notes associated with simulations:
% Notes: In noiseless system with curve fitting
% % --- % showing the percentage of error in estimating the slope:
% %      linear gives 20% error
% %      quadradic gives 5% error
% %      cubic gives < 0.5% error
%d = gradient(z)./gradient(x);
%D = gradient(Z)./gradient(X);
% 
% pp = polyfit(x,z,1);
% pd = polyder(pp);
% %xx = linspace(x(1),x(end),1000)';
% xx = x;
% yy = polyval(pd,xx);
% figure
% plot(xx(1:end-1),yy(1:end-1)./d(1:end-1));
% hold on
% 
% pp = polyfit(x,z,2);
% pd = polyder(pp);
% yy = polyval(pd,xx);
% plot(xx(1:end-1),yy(1:end-1)./d(1:end-1));
% 
%pp = polyfit(x,z,3);
% pd = polyder(pp);
% yy = polyval(pd,xx);
% plot(xx(1:end-1),yy(1:end-1)./d(1:end-1));
% 
% %plot(X,Sz)
% figure
% plot(x,z)
% hold on
% plot(xx,yy,'r')
% 
% figure
% plot(X,Z)
% hold on
% plot(X,Zn)
% plot(X,Szn,'k')
