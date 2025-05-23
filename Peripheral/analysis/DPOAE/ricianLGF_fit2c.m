function [jar,h1,h2,h3,h4,h5] = ricianLGF_fit2c(Z,Sz,X)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [jar,h1,h2,h3,h4,h5] = ricianLGF_fit2c(Z,Sz,X);
% [jar,h1,h2,h3,h4,h5] = ricianLGF_fit2c;
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
% Updated: May 9, 2025 -- ssg -- added new initial fit for better stability
%                         when SNR is poor. Uses updated linking functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
doPlot = 1; % make plots (1) or not (0)
doBS = 1; % do bootstrap CIs (1) or not (0)
nSweeps = 24;

    

    jar = []; h1 = []; h2 = []; h3 = []; h4 = []; h5 = []; % initialize outputs
    if nargin < 3 % if no inputs, simulate from a generalized logistic model
        if nargin == 1
            BBB = Z;
        end
        %[Z,Sz,X,Zm,Q5,Q8,Q11] = simulateDataNL;
        [Z,Sz,X] = simulateLGF(1,BBB);
        [Zm,Szm,Xm] = simulateLGF(0,BBB); % model fit
        % cut simulated data so only fitting below the peak of the LGF
        
        if doPlot == 1
        figure(10)
        plot(X,Z,'b')
        hold on
        plot(X,Sz,'k')
        end

        cut1 = 1; % start at sample 1
        Zsm = meanSmoother(Z,5);
        Zsm = meanSmoother(Zsm,5);
        [DPmax,cut2] = max(Zsm); % just below peak
         cut2 = cut2 - 5;
        %Zm = Zm(cut1:cut2);
        Z = Z(cut1:cut2);
        Sz = Sz(cut1:cut2);
        X = X(cut1:cut2);

        SNR = Zsm - mean(Sz);
        indx = find(SNR>6); % index of data points with SNR > 6
        N6 = length(indx);
        M6 = max(SNR);

        if doPlot == 1
        plot(X,Z,'b','LineWidth',2)
        end


    end

    upSample = 1;
    if upSample > 1
        originalN = length(Z);
        NN = originalN * upSample;
        XX = linspace(X(1),X(end),NN)';
        Z = interp1(X,Z,XX,'pchip');
        Sz = interp1(X,Sz,XX,'pchip');
        X = XX;
    end


    % ---------------------------------------------------------------------
    % Input is in dB SPL. Convert to linear units (mPa)
    ref = 0.02; % reference is 20 mPa (not micro! ricean model expects mPa)
    sz = 10.^(Sz/20)*ref; % standard error of the mean vector
    Esz = mean(sz); % expectation (mean) of standard error of the mean vector
    z = 10.^(Z/20)*ref; % DPOAE levels (mPa)
    % if nargin == 0
    %     zm = 10.^(Zm/20)*ref; % simulated noiseless DPOAE levels (mPa)
    % end
    x = 10.^(X/20)*ref; % target L2 stimulus level (mPa)
    
    if doPlot == 1
    figure(11)
    plot(x,z,'b')
    hold on
    plot(x,sz,'k')
    end

    % package values into output structure --------------------------------
    jar.ref = ref; % pressure reference (mPa)
    jar.X = X; % L2 stimulus levels (dB SPL)
    jar.x = x; % L2 stimulus levels (mPa)
    jar.Z = Z; % measured DPOAE magnitudes (dB SPL)
    jar.z = z; % measured DPOAE magnitudes (mPa)
    jar.Sz = Sz; % estimated noise floor, standard error of the mean (dB SPL)
    jar.sz = sz; % standard error (mPa)
    jar.Esz = Esz; % expected (mean) value of standard error
    jar.nSweeps = nSweeps; % number of sweeps in the average

    % ------------------------------------------
    % ------------------------------------------
    disp('Computing Rician fit')
    jar.nLooks = 100; % set the vector length -- make sure an even number!
    jar = fittingRoutine(jar); 
    % ------------------------------------------
    % ------------------------------------------

    % ------------------------------------------
    % ------------------------------------------
    if doBS == 1
        disp('Computing bootstrap confidence intervals')
        tic % get bootstrapped confidence intervals
        nIterations = 10;
        jar = getCIs(jar,nIterations,jar.B);
        toc/60
    end
    % ------------------------------------------
    % ------------------------------------------

    % plot results --------------------------------------------------------

    % compare rician results with lsf2 and lsf3
    jar.yhat_lsf = polyval([jar.betaLSF,0],x); % mle model fit (mPa)
    jar.Yhat_lsf = 20*log10(jar.yhat_lsf/ref);
    XX = [x.^3,x.^2,x];
    [Q,R] = qr(XX,0); % orthogonal-triangular decomposition; R is the Cholesky factor of the X matrix
    b = full(R\(Q'*z)); % Same as p = D*X\(D*y); b is a vector of coefficients; also same as b = (X'*X)\(X'*Y).
    jar.yhat_lsf3 = polyval([b(:)',0],x); % mle model fit (mPa)
    jar.Yhat_lsf3 = 20*log10(jar.yhat_lsf3/ref);
    jar.b3 = b;

if doPlot == 1
figure(10)
plot(jar.X,jar.Yhat_mle,'r','lineWidth',1)
plot(jar.X,jar.Yhat_lsf,'g--')
plot(jar.X,jar.Yhat_lsf3,'c--')
plot(jar.X,Zm(cut1:cut2),'k')

figure(11)
plot(jar.x,jar.yhat_mle,'r','lineWidth',1)
plot(jar.x,jar.yhat_lsf,'g--')
plot(jar.x,jar.yhat_lsf3,'c--')
plot(x,10.^(Zm(cut1:cut2)/20)*0.02,'k')
end

% percent mean absolute error
zm = 10.^(Zm(cut1:cut2)/20)*.02;
pmae_mle = 100*(mean(abs((jar.yhat_mle ./ zm)-1)));
pmae_lsf = 100*(mean(abs((jar.yhat_lsf ./ zm)-1)));
pmae_lsf3 = 100*(mean(abs((jar.yhat_lsf3 ./ zm)-1)));

disp(['percent mean absolute error:'])
disp(['mle: ',num2str(pmae_mle)])
disp(['lsf2: ',num2str(pmae_lsf)])
disp(['lsf3: ',num2str(pmae_lsf3)])

if doPlot == 1
h1 = plotLL(jar);
end

jar.pmae_mle = pmae_mle;
jar.pmae_lsf = pmae_lsf;
jar.pmae_lsf3 = pmae_lsf3;
jar.N6 = N6;
jar.M6 = M6;

return
keyboard

gra = gradient(10.^(Zm(cut1:cut2)/20)*.02)./gradient(x); % true gradient
Gra = 20*log10(gra);
graLSF = gradient(jar.yhat_lsf)./gradient(x); % LSF gradient
GraLSF = 20*log10(graLSF);


h1 = plotLL(jar);
%h2 = plotFit(jar);
h3 = plotSlope(jar);

figure(h3)
subplot(2,1,2)
plot(x,gra,'k')
plot(x,graLSF,'g--')
subplot(2,1,1)
plot(X,Gra,'k')
plot(X,GraLSF,'g--')




errMag_rice = sqrt(mean(jar.yhat_mle - zm).^2);
ErrMag_rice = 20*log10(errMag_rice);
errMag_lsf = sqrt(mean(jar.yhat_lsf - zm).^2);
ErrMag_lsf = 20*log10(errMag_lsf);

errSlope_rice = sqrt(mean(jar.slope_mle - gra).^2);
ErrSlope_rice = 20*log10(errSlope_rice);
errSlope_lsf = sqrt(mean(graLSF - gra).^2);
ErrSlope_lsf = 20*log10(errSlope_lsf);


    %ErrSlope_lsf - ErrSlope_rice
    % positive errors mean lsf is better


keyboard

end

% INTERNAL FUNCTIONS ------------------------------------------------------
function [jar] = fittingRoutine(jar,B)
    x = jar.x;
    z = jar.z;
    Esz = jar.Esz;
    nLooks = jar.nLooks;
    ref = jar.ref;
    if nargin == 1 % if running for the initial fit
        % % find initial guess of coefficients using standard LSF -------------
        jar.betaLSF = initialize(x,z);
        % Set up grid to search over:
        b2 = makeB2(jar.betaLSF(1),nLooks);
        b1 = makeB1(jar.betaLSF(2),nLooks);
        B = [b2,b1]; % matrix of coefficients to search
        
        % run search broad ----------------------------------------------------------
        LL = getLLfunction2(B,x,Esz,z); % get the log likelihood function
        ml = max(LL(:)); % maximum log likelihood
        [rr,cc] = find(LL==ml); % location of maximum likelihood
        partial1 = LL(rr,:); % partial is a "slice" of the LL matrix
        partial2 = LL(:,cc);
        beta1 = B(:,2); % vector of coefficients corresponding to the LL partial slice
        beta2 = B(:,1);    
        mle = [beta2(rr),beta1(cc)]; % current maximum likelihood estimate vector
        grad = gradient(partial1(:))./gradient(beta1);
        grad = gradient(grad)./gradient(beta1);
        d2_1 = grad(cc);
        grad = gradient(partial2(:))./gradient(beta2);
        grad = gradient(grad)./gradient(beta2);
        d2_2 = grad(rr);
        CI(1,1) = mle(1) - 16*(1./sqrt(-d2_2)); % lower bound of CI
        CI(2,1) = mle(1) + 16*(1./sqrt(-d2_2)); % upper bound of CI 
        CI(1,2) = mle(2) - 16*(1./sqrt(-d2_1)); % lower bound of CI
        CI(2,2) = mle(2) + 16*(1./sqrt(-d2_1)); % upper bound of CI 
       
    
        % figure(10)
        % subplot(2,1,1)
        % hold on
        % plot(beta2,partial2)
        % plot(mle(1),ml,'*m')
        % plot(CI(:,1),ml,'*c')
        % xlabel('\beta_2')
        % subplot(2,1,2)
        % hold on
        % plot(beta1,partial1)
        % plot(mle(2),ml,'*m')
        % plot(CI(:,2),ml,'*c')
        % xlabel('\beta_1')
    
    
        % prep search narrow ----------------------------------------------------------
        b2 = linspace(CI(1,1),CI(2,1),nLooks)';
        if max(b2)>0 % b2 must be negative
            b2 = b2 - max(b2);
        end
        b1 = linspace(CI(1,2),CI(2,2),nLooks)';
        if min(b1)<0
            b1 = b1 + -min(b1);
        end
        B = [b2,b1]; % matrix of coefficients to search
    end

    % run search narrow ----------------------------------------------------------
    LL = getLLfunction2(B,x,Esz,z); % get the log likelihood function
    ml = max(LL(:)); % maximum log likelihood
    [rr,cc] = find(LL==ml); % location of maximum likelihood
    partial1 = LL(rr,:); % partial is a "slice" of the LL matrix
    partial2 = LL(:,cc);
    beta1 = B(:,2); % vector of coefficients corresponding to the LL partial slice
    beta2 = B(:,1);    
    mle = [beta2(rr),beta1(cc)]; % current maximum likelihood estimate vector
    grad = gradient(partial1(:))./gradient(beta1);
    grad = gradient(grad)./gradient(beta1);
    d2_1 = grad(cc);
    grad = gradient(partial2(:))./gradient(beta2);
    grad = gradient(grad)./gradient(beta2);
    d2_2 = grad(rr);
    CI(1,1) = mle(1) - 1.96*(1./sqrt(-d2_2)); % lower bound of CI was 1.96
    CI(2,1) = mle(1) + 1.96*(1./sqrt(-d2_2)); % upper bound of CI 
    CI(1,2) = mle(2) - 1.96*(1./sqrt(-d2_1)); % lower bound of CI was 1.96
    CI(2,2) = mle(2) + 1.96*(1./sqrt(-d2_1)); % upper bound of CI 

    jar.mle = mle; % current maximum likelihood estimate
    jar.ml = ml; % maximum likelihood value
    jar.B = B; % new set of coefficients to search over
    jar.Partials = [partial2(:),partial1(:)]; % slices of LL values for each coefficient
    jar.CI = CI; % confidence intervals for mle (second derivative of LL)
    jar.d2 = [d2_2,d2_1];
    jar.LL = LL;

    % figure(10)
    % subplot(2,1,1)
    % hold on
    % plot(beta2,partial2)
    % plot(mle(1),ml,'or')
    % plot(CI(:,1),ml,'*g')
    % xlabel('\beta_2')
    % subplot(2,1,2)
    % hold on
    % plot(beta1,partial1)
    % plot(mle(2),ml,'or')
    % plot(CI(:,2),ml,'*g')
    % xlabel('\beta_1')

    if nargin == 1
        % calculate model output
        jar.yhat_mle = polyval([jar.mle,0],x); % mle model fit (mPa)
        jar.Yhat_mle = 20*log10(jar.yhat_mle/ref); % mle model fit (dB SPL)
        jar.slope_mle = (gradient(jar.yhat_mle)./gradient(x)); % slope of the linear fit
        jar.Slope_mle = 20*log10(jar.slope_mle);
        jar.dydx = polyder([jar.mle,0]); % first derivative of fit (the slope)
        jar.d2ydx2 = polyder(jar.dydx); % second derivative of fit (single value)

        % confidence intervals (bootstrapped) will be filled in later
        ciFit = [];
        CIFit = [];
        ciSlope = [];
        CISlope = [];
        ciD2 = [];
        ciMLE = [];
    end
end
function [beta] = initialize(x,z)
    % find initial fit using yhat = beta2*x^2 + beta1*x + 0
    % Solve the least squares problem
    X = [x.^2,x];
    [Q,R] = qr(X,0); % orthogonal-triangular decomposition; R is the Cholesky factor of the X matrix
    b = full(R\(Q'*z)); % Same as p = D*X\(D*y); b is a vector of coefficients; also same as b = (X'*X)\(X'*Y).
    
    %yhat = X*b;      % Predicted responses at each data point.
    %plot(x,z)
    %hold on
    %plot(x,yhat,'r')

    beta = [b(1),b(2)]; % best initial fit, with noise, with no offset term
end
function [B] = makeB2(b,nLooks)
    % create the search vector of coefficients for beta_2 * x.^2
    %B = sort([b*0.8,b*1.2]); % the extremes are +/-20% of the initial guess
    %B = sort([b*0.05,b*1.95]);

    d = abs(b) * 3;
    B = sort([b-d,0]);
    %B = sort([0,b+d]);
    %B = sort([0.1, 0.001]);
    B = linspace(B(1),B(2),nLooks)';
end
function [B] = makeB1(b,nLooks)
    % create the search vector of coefficients for beta_1 * x
    %B = sort([b*0.05,b*1.95]);
    %B = sort([b*0.01,b*2]);
    %d = abs(b) * 3;
    %B = sort([b-d,b+d]);
    %B = sort([b-d,0]);
    
    %B = sort([0.1, 0.001]);
    B = linspace(-75,-20,nLooks)';
    B = 10.^(B./20);
    %B = linspace(B(1),B(2),nLooks)';
end

function [jar] = getCIs(jar,nIterations,B)
    jar2.x = jar.x;
    jar2.nLooks = jar.nLooks;
    jar2.ref = jar.ref;

    nSweeps = jar.nSweeps;
    ref = jar.ref;
    % create confidence intervals -----------------------------------------
    % nf = sqrt(a^2 + a^2) / sqrt(nSweeps); % nf = actual noise floor (pa)
    % s = sqrt(A^2 + A^2) / sqrt(nSweeps);
    % s = sqrt((A^2 + A^2)/nSweeps)
    % s^2 = (A^2 + A^2)/nSweeps
    % s^2 * nSweeps = 2*A^2
    % (s^2 * nSweeps)/2 = A^2
    % sqrt((s^2 * nSweeps)/2) = A
    A = sqrt((jar.Esz.^2 * nSweeps)/2);
    % noise distribution
    %noiseDist = makedist('normal','mu',A,'sigma',std(sz));
    %noiseDist = makedist('tlocationscale','mu',A,'sigma',std(sz),'nu',2);
    %a = abs(noiseDist.random(1));
    a = A;
    pdN = makedist('normal','mu',0,'sigma',a);
    Y = jar.yhat_mle; % create the original DP as straight line
    y = 10.^(Y/20)*ref;
    yy = jar.yhat_mle.*exp(1i*0); % dpoae in complex form
    yy = repmat(yy,1,nSweeps);
    nLevels = length(yy);
    
    MLE = zeros(2,nIterations); % initialize output
    YHAT = zeros(length(jar.yhat_mle),nIterations);
    for ii=1:nIterations
        if mod(ii,10)==0
            disp(['   find CIs: ',num2str(ii),' of ',num2str(nIterations)])
        end
        noise = pdN.random(nLevels,nSweeps) + 1i*pdN.random(nLevels,nSweeps); % generate random noise
        zz = yy + noise;
        z2 = abs(mean(zz,2));
        sz2 = std(zz,[],2)/sqrt(nSweeps);
        %Z2 = 20*log10(z2/ref);
        Sz2 = 20*log10(sz2/ref);
        ESz2 = mean(Sz2);
        Esz2 = 10.^(ESz2/20)*ref;
        
        jar2.z = z2;
        jar2.Esz = Esz2;
        jar2 = fittingRoutine(jar2,B);

        MLE(:,ii) = jar2.mle;
        YHAT(:,ii) = polyval([jar2.mle,0],jar2.x); % mle model fit (mPa)
    end

    alfa = 0.05;
    ciFit(:,1) = quantile(YHAT,alfa/2,2);
    ciFit(:,2) = quantile(YHAT,1-(alfa/2),2);
    CIFit = 20*log10(ciFit/ref);

    % calculate derivatives
    MLE = MLE';
    SLOPE = zeros(length(jar2.x),nIterations);
    DYDX = zeros(nIterations,2);
    D2YDX2 = zeros(nIterations,1);
    for ii=1:nIterations
        DYDX(ii,:) = polyder([MLE(ii,:),0]); % first derivative of fit (the slope)
        SLOPE(:,ii) = polyval([DYDX(ii,:)],jar2.x); % mle model fit (mPa)
        D2YDX2(ii,1) = polyder(DYDX(ii,:)); % second derivative of fit (single value)
    end

    ciSlope(:,1) = quantile(SLOPE,alfa/2,2);
    ciSlope(:,2) = quantile(SLOPE,1-(alfa/2),2);
    CISlope = 20*log10(ciSlope);
    
    ciD2(1,1) = quantile(D2YDX2,alfa/2);
    ciD2(2,1) = quantile(D2YDX2,1-(alfa/2));

    ciMLE(1,1) = quantile(MLE(:,1),alfa/2);
    ciMLE(1,2) = quantile(MLE(:,1),1-(alfa/2));
    ciMLE(2,1) = quantile(MLE(:,2),alfa/2);
    ciMLE(2,2) = quantile(MLE(:,2),1-(alfa/2));

    jar.ciFit = ciFit;
    jar.CIFit = CIFit;
    jar.ciSlope = ciSlope;
    jar.CISlope = CISlope;
    jar.ciD2 = ciD2;
    jar.ciMLE = ciMLE;
    %jar.dydx = polyder([jar.mle,0]); % first derivative of fit (the slope)
    %jar.d2ydx2 = polyder(jar.dydx); % second derivative of fit (single value)
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
    snr = y ./ Esz; % estimated signal to noise ratio
    ESZ = repmat(Esz,length(snr),1);
    nu = gg(snr,ESZ,y); % estimate rician nu (location) parameter
    sigma = hh(snr,ESZ); % esimate rician sigma (spread) parameter
        
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

function [nu] = gg(snr,sz,y)
    % For nu:
    % define piecewise function nu = g(y,sz)
    %       if SNR >= 0, nu = y
    %       if SNR <= -15, nu = sz*0.3069 (which is -10.26 dB)
    %       otherwise nu = sz * (d + (a-d)/(1 + (snr/c)^b))
    % We fit k = nu/sz, so nu = k*sz
    % otherwise a 4-parameter logistic
    % g(snr) = d + (a-d)/(1 + (snr/c)^b);

    a = 0.309983012361646;
    b = 3.605490836836389;
    c = 0.828415511855433;
    d = 1.331525891400627;
    
    SNR = 20*log10(snr);

    indx1 = find(SNR>= 0);
    indx2 = find(SNR<= -20);
    indx3 = find(SNR>-20 & SNR<0);
    nu1 = y;
    nu2 = sz * 0.3069;
    % note: .3069 found this way: mean(nu(indx2)) / mean(sz(indx2))
    nu3 = sz .* (d + (a-d)./(1 + (snr/c).^b));

    nu = zeros(size(SNR));
    nu(indx1) = nu1(indx1);
    nu(indx2) = nu2(indx2);
    nu(indx3) = nu3(indx3);

end
function [sigma] = hh(snr,sz)
    % For sigma:
    % define piecewise function sigma = h(sz)
    %       if SNR >= 0, sigma = sz * 0.711
    %       if SNR <= -15, sigma = sz * 0.667 
    %       otherwise sigma = sz * (p1*snr^4 + p2*snr^3 + p3*snr^2 + p4*snr + p5)
    % We fit q = sigma/sz, so sigma = q*sz
    %  a 4th-order polyomial 
    % h(snr) = (p1*snr^4 + p2*snr^3 + p3*snr^2 + p4*snr + p5);
    
    p1 = 0.676563246959900;
    p2 = -1.462574793496384;
    p3 = 0.922065507608489;
    p4 = -0.094255444028801;
    p5 = 0.670568403772949;
    
     SNR = 20*log10(snr);

    indx1 = find(SNR>= 0);
    indx2 = find(SNR<= -20);
    indx3 = find(SNR>-20 & SNR<0);
    sigma1 = sz * 0.711;
    % note: 0.711 found this way: mean(sigma(indx1)) / mean(sz(indx1))
    sigma2 = sz * 0.667;
    % note: 0.667 found this way: mean(sigma(indx2)) / mean(sz(indx2))
    sigma3 = sz .* (p1*snr.^4 + p2*snr.^3 + p3*snr.^2 + p4*snr + p5);
    
    sigma = zeros(size(SNR));
    sigma(indx1) = sigma1(indx1);
    sigma(indx2) = sigma2(indx2);
    sigma(indx3) = sigma3(indx3);
end

function [h] = plotFit(jar)
    h = figure;
    subplot(1,2,1)
     hold on
     X = [jar.X;flipud(jar.X)];
     Y = [jar.CIFit(:,1);flipud(jar.CIFit(:,2))];
     C = [.8 .8 .8];
     fill(X,Y,C,'EdgeColor','none');
     %if ~isempty(jar.Zm)
     %    plot(jar.X,jar.Zm,'Color',[.7 .7 .7],'LineWidth',2)
     %end
     plot(jar.X,jar.Z,'b')
     plot(jar.X,jar.Sz,'k')
     plot(jar.X,jar.Yhat_mle,'r-')
     xlim([jar.X(1),jar.X(end)])
     xlabel('Target L2 (dB FPL)','FontSize',14)
     ylabel('Ldp (dB SPL)','FontSize',14)
     grid on
     set(gca,'FontSize',12)
     % if ~isempty(jar.Zm)
     %     legend('true DP','measured DP','noise floor','model fit','Location','northwest')
     % else
     %     legend('measured DP','noise floor','model fit','Location','northwest')
     % end
    subplot(1,2,2)
     hold on
     %if ~isempty(jar.zm)
     %    plot(jar.x,jar.zm,'Color',[.7 .7 .7],'LineWidth',2)
     %end
     X = [jar.x;flipud(jar.x)];
     Y = [jar.ciFit(:,1);flipud(jar.ciFit(:,2))];
     C = [.8 .8 .8];
     fill(X,Y,C,'EdgeColor','none');
     plot(jar.x,jar.z,'b')
     plot(jar.x,jar.sz,'k')
     plot(jar.x,jar.yhat_mle,'r-')
     xlim([jar.x(1),jar.x(end)])
     xlabel('Target L2 (mPa)','FontSize',14)
     ylabel('Ldp (mPa)','FontSize',14)
     grid on
     set(gca,'FontSize',12)
     %if ~isempty(jar.Zm)
     %    legend('true DP','measured DP','noise floor','model fit','Location','northwest')
     %else
     %    legend('measured DP','noise floor','model fit','Location','northwest')
     %end
    xxx = [30 40 45 50];
    xticks(10.^(xxx/20)*jar.ref);
    xl = {};
    for ii=1:length(xxx)
        xl = [xl,{num2str(xxx(ii))}];
    end
    xticklabels(xl)    
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
     title(['\beta_2_M_L_E = ',num2str(jar.mle(1)),'x^2'])
     xlabel('\beta_2','FontSize',14)
     ylabel('Log Likelihood','FontSize',14)
     grid on
     set(gca,'FontSize',14)
     % second derivative ci's (too narrow!)
     line([jar.CI(1,1),jar.CI(1,1)],[ymin,ymax],'Color',[1 0 0],'LineStyle','--','LineWidth',0.5)
     line([jar.CI(2,1),jar.CI(2,1)],[ymin,ymax],'Color',[1 0 0],'LineStyle','--','LineWidth',0.5)
     % bootstrapped ci's
     line([jar.ciMLE(1,1),jar.ciMLE(1,1)],[ymin,ymax],'Color',[0 1 0],'LineStyle','--','LineWidth',0.5)
     line([jar.ciMLE(1,2),jar.ciMLE(1,2)],[ymin,ymax],'Color',[0 1 0],'LineStyle','--','LineWidth',0.5)
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
     title(['\beta_1_M_L_E = ',num2str(jar.mle(2)),'x'])
     xlabel('\beta_1','FontSize',14)
     ylabel('Log Likelihood','FontSize',14)
     grid on
     set(gca,'FontSize',14)
     % second derivative ci's (too narrow!)
     line([jar.CI(1,2),jar.CI(1,2)],[ymin,ymax],'Color',[1 0 0],'LineStyle','--','LineWidth',0.5)
     line([jar.CI(2,2),jar.CI(2,2)],[ymin,ymax],'Color',[1 0 0],'LineStyle','--','LineWidth',0.5)
     % bootstrapped ci's
     line([jar.ciMLE(2,1),jar.ciMLE(2,1)],[ymin,ymax],'Color',[0 1 0],'LineStyle','--','LineWidth',0.5)
     line([jar.ciMLE(2,2),jar.ciMLE(2,2)],[ymin,ymax],'Color',[0 1 0],'LineStyle','--','LineWidth',0.5)
     xlim([min(jar.B(:,2)),max(jar.B(:,2))])
     ylim([ymin,ymax])
end
function [h] = plotSlope(jar)
    % plot the slopes
    h = figure; % plot on log-log scale
    subplot(2,1,1)
     %if ~isempty(jar.SlopeM)
     %    plot(jar.X,jar.SlopeM,'Color',[.7 .7 .7],'LineWidth',1)
     %end
     hold on
     X = [jar.X;flipud(jar.X)];
     Y = [jar.CISlope(:,1);flipud(jar.CISlope(:,2))];
     C = [.8 .8 .8];
     fill(X,Y,C,'EdgeColor','none');
     plot(jar.X,jar.Slope_mle,'r')
     xlim([jar.X(1),jar.X(end-1)])
     xlabel('L2 (dB)','FontSize',14)
     ylabel('Slope (dB)','FontSize',14)
     set(gca,'FontSize',12)
     %if ~isempty(jar.SlopeM)
     %   legend('true','model fit','Location','southwest')
     %end
    subplot(2,1,2) % plot on linear-linear scale
     %if ~isempty(jar.slopem)
     %    plot(jar.x,jar.slopem,'Color',[.7 .7 .7],'LineWidth',1)
     %end
     hold on
     X = [jar.x;flipud(jar.x)];
     Y = [jar.ciSlope(:,1);flipud(jar.ciSlope(:,2))];
     C = [.8 .8 .8];
     fill(X,Y,C,'EdgeColor','none');
     plot(jar.x,jar.slope_mle,'r')
     xlim([jar.x(1),jar.x(end-1)])
     xlabel('L2 (mPa)','FontSize',14)
     ylabel('slope (ratio)','FontSize',14)
     set(gca,'FontSize',12)
     %if ~isempty(jar.slopem)
     %    legend('true','model fit','Location','southwest')
     %end
    xxx = [10 30 40 45 50];
    xticks(10.^(xxx/20)*jar.ref);
    xl = {};
    for ii=1:length(xxx)
        xl = [xl,{num2str(xxx(ii))}];
    end
    xticklabels(xl)    
end

function [Z,Sz,X,zbar] = simulateLGFvInternal(B)
    % Create a simulated LGF by passing primaries through a generalized logistic function.
    ref = 0.002; %.002; % pressure reference;

    % create the nonlinearity 
    % standard B value = 2, which gives an intercept of about -37
    % B = 4 --> peak 30, intercept -21
    % B = 3 --> peak 23, intercept -26
    % B = 2 --> peak 14, intercept -36
    % B = 1 --> peak 0 dB, intercept -55
    % B = 0.75 --> peak -8, troubled fit
    % B = 0.5 --> all gone!
    %
    % Summary: slope should go between a min and max of -20 and -60 dB = -0.1 to 0.001 
    B = 2;
    A = -.0125; %-1; % left horizontal assymtote
    K = -A; % right horizontal asymtote
    C = .2;
    nu = 10; %6; % >0, where growth occurs
    Q = .4; %.75;

    % A + ( (K-A) ./ ((1+Q.*exp(-B*input)).^(1/nu)) );
    % -1 + (2./(1+exp(-B*input)));

    fs = 96000; % sampling rate
    dur = 9; % input signal duration (9 seconds)
    N = round(fs*dur); % number of samples in input
    t = (0:1:N-1)'/fs; % time vector
    f2 = 2000; % f2 frequency (Hz)
    fratio = 1.22;
    f1 = round(f2/fratio); % f1 frequency (Hz)
    fdp = 2*f1-f2; % DP frequency (Hz)
    L1max = 55; % level of f1. Yes, should ideally be 65, but to get results that look
             % like our data, you have to use this lower level
    L2min = -3;
    L2max = 73;
    L2 = linspace(L2min,L2max,N)'; % vector of f2 levels (dB)
    l2 = 10.^(L2/20)*ref; % vector of f2 levels (linear)
    L1 = ones(N,1)*L1max;
    l1 = 10.^(L1/20)*ref; % vector of f1 levels (linear)

    % create the primary freuqencies at unit amplitude
    p1 = l1.*sin(2*pi*f1*t); % primary 1
    p2 = l2.*sin(2*pi*f2*t); % primary 2
    nSweeps = 24; % this chosen to mirror what we actually collect
    noiseAmp1 = 0; %.1; %.1
    noiseAmp2 = .000004; %.075 * 10.^(dBoffset/20); %0.075 is base (-21 dB SPL noise floor)
    
    for jj=1:nSweeps
        input = p1 + p2;
        input = input + (randn(size(input))*noiseAmp1);
        Y(:,jj) = A + ( (K-A) ./ ((C+Q.*exp(-B*input)).^(1/nu)) );
        %Y(:,jj) = A + ( (K-A) ./ ((C+Q.*exp(-B*Y(:,jj))).^(1/nu)) );
        Yn(:,jj) =  Y(:,jj) + (randn(size(input))*noiseAmp2); % create a version with noise
    end
    Y = Y * (.00002/ref); % put into Pa
    Yn = Yn * (.00002/ref); % put into Pa

    step = 1;
    if whichOne == 1
        [signal,nf,snr,targets,zbar] = sweptLSF(L2,Yn,fdp,L2min,L2max,step,fs);

        % plot(targets,signal)
        % hold on
        % plot(targets,nf)
        % xlim([0 70])
    
        Z = signal;
        Sz = nf;
        X = targets;

    elseif whichOne == 0
        [signal,nf,snr,targets,zbar] = sweptLSF(L2,Y,fdp,L2min,L2max,step,fs);
        Z = signal;
        Sz = nf;
        X = targets;
    
        % plot(targets,signal)
        % hold on
        % plot(targets,nf)
        
    end
    
    Z = Z + 87;
    Sz = Sz + 87;

    phi = angle(zbar);
    z = 10.^(Z/20)*.02;
    zbar = z.*exp(1i*phi);

end


% INTERNAL FUNCTIONS ------------------------------------------------------




% OLD CODE ----------------------------------------------------------------
% Plotting for 3D Log Likelihood
    % B = jar.B;
    % LL = jar.LL;
    % mle = jar.mle;
    % figure % plot 3D likelihood function
    % h = surface(B(:,1),B(:,2),LL,'FaceAlpha',0.5);
    % hold on
    % xlabel('\beta_2','FontSize',16)
    % ylabel('\beta_1','FontSize',16)
    % zlabel('Log-Likelihood','FontSize',14)
    % view([60 15])
    % grid on
    % colormap jet
    % %h.EdgeColor = 'none';
    % xlim([B(1,1),B(end,1)])
    % ylim([B(1,2),B(end,2)])
    % zmin = -2564;
    % zmax = 250;
    % zlim([zmin,zmax])
    % set(gca,'FontSize',14)
    % Lx = [mle(1), mle(1), mle(1), mle(1)];
    % Ly = [min(B(:,2)), max(B(:,2)), max(B(:,2)), min(B(:,2))];
    % Lz = [zmin,zmin,zmax,zmax];
    % L1 = fill3(Lx,Ly,Lz,'g','FaceAlpha',.45,'EdgeColor',[0 1 0]);
    % Ly = [mle(2), mle(2), mle(2), mle(2)];
    % Lx = [min(B(:,1)), max(B(:,1)), max(B(:,1)), min(B(:,1))];
    % Lz = [zmin,zmin,zmax,zmax];
    % L1 = fill3(Lx,Ly,Lz,'g','FaceAlpha',.45,'EdgeColor',[0 1 0]);
    % plot3(mle(1),mle(2),max(LL(:)),'r.','MarkerSize',20)
    % view([64 46])
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
