function [OUT] = myLogisticReg(X,Y,pdPrior)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [OUT] = myLogisticReg(X,Y)
%
% Shawn Goodman
% January 14, 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    doPlot = 0;

    N = 200;
    b_min = -10;
    b_max = 100;
    B = linspace(b_min,b_max,N)';

    m_min = .75;
    m_max = 1.4;
    Nm = 50;
    %M = linspace(m_min,m_max,Nm)';
    M = 1;

    if nargin ==0
        X = [40 50 55 55 55 60 60 60 65 70 80]';
        Y = [0   0  0  1  0  0  1  1  1  1  1]';
    else
        X = X(:);
        Y = Y(:);
    end

    if nargin < 3
        priorsPathName = 'C:\myWork\ARLas\Peripheral\experiments\ARL\Audiometer\';
        priorsFileName = 'HarpAudioFPL.mat';
        [PRIORS] = getPriors(priorsPathName,priorsFileName,1000,30);
        pdPrior = makedist('normal','mu',PRIORS.mu,'sigma',PRIORS.sigma);
    end
    prior = pdPrior.cdf(B);

    LL = iterateLL(M,B,X,Y,pdPrior);
    [thd,sigma,CDF,pd] = calculateML(B,M,LL);

    xxx = (0:.1:80)';
    ppp1 = pd.cdf(xxx)'; % normal cdf with sigma estimate
    %yyy = 1.*xxx - thd;
    %ppp2 = exp(yyy) ./ (1+exp(yyy)); % logit

    yyy = 1.*X - thd;
    Epsilon = Y - (exp(yyy) ./ (1+exp(yyy)));
    D = X - thd;
    coherence = abs(D) .* abs(Epsilon);
    badIndx = find(coherence > 5);
    if isempty(badIndx)
        nRejects = 0;
    else
        nRejects = length(badIndx);
    end
    %nRejects


    % Plotting --------
    if doPlot == 1
        figure(10)
        %plot(xxx,ppp1,'b')
        hold on
        plot(X,Y,'*','Color',[1 0 0])
        plot(thd,0.5,'o','Color',[.7 .7 .7])
        plot(B,prior,'--','Color',[.7 .7 .7])
        plot(B,CDF,'Color',[.7 .7 .7])
        line([thd-sigma*2,thd + sigma*2],[0.5 0.5],'Color',[.7 .7 .7],'LineWidth',2)
        grid on
        title(['thd = ',num2str(thd)])
    end

    % get rid of outliers
    if nRejects > 0
        plot(X(badIndx),Y(badIndx),'*','Color',[.7 .7 .7])
    
        X(badIndx) = [];
        Y(badIndx) = [];
        LL = iterateLL(M,B,X,Y,pdPrior);
        [thd,sigma,CDF,pd] = calculateML(B,M,LL);
    end

    if doPlot == 1
        plot(thd,0.5,'ro')
        plot(B,CDF,'g')
        line([thd-sigma*2,thd + sigma*2],[0.5 0.5],'Color',[0 1 0],'LineWidth',2)
    end

    OUT.thd = thd;
    OUT.sigma = sigma;
    OUT.CDF = CDF;
    OUT.X = X;
    OUT.Y = Y;
    OUT.B = B;
    OUT.badIndx = badIndx;
    OUT.nRejects = nRejects;
    OUT.prior = prior;

 
end
% -------------------------------------------------------------------------
function [thd,sigma,CDF,pd] = calculateML(B,M,LL)
    % find the Maximum Likelihood Estimate as the theshold
    [row,col] = find(LL == max(LL(:))); % the location of thetaHat_MLE 
    m = M(row);
    b = B(col);
    thd = b / m; % threshold
    % find sigma as the estimate of 95% confidence around threshold
    indx = find(LL~=-Inf);
    LL = LL(indx);
    ll = exp(LL);
    %L = L(indx);
    %l = exp(L);
    BB = B;
    BB  = BB(indx);
    ll = ll ./ sum(ll);
    BB = BB(:);
    ll = ll(:);

    [HDI,HDIdensity,inOut] = getHDI(0.05,BB,ll);
    sigma2 = HDI(2) - HDI(1);
    sigma = sigma2 / 2 / 2;
    pd = makedist('normal','mu',thd,'sigma',sigma);
    CDF = pd.cdf(B);

end
function [LL] = calculateLL(m,b,X,Y,pd)
    y = m.*X - b;
    p = 1./(1+exp(-y)); % same as: p = exp(y) ./ (1+exp(y));
    %L(ii,jj) = prod(p.^Y .* (1-p).^(1-Y));
    %LL(ii,jj) = sum((Y.*log(p)) + (1-Y).*log(1-p));
                %sum(y.*log(p)) + sum((1-y).*log(1-p));
                
    q1 = Y.*log(p);
    q2 = (1-Y).*log(1-p); 
    thd = b / m;
    q3 = log(pd.pdf(thd));
    q4 = q3; 

    q1(find(isnan(12))) = 0;
    q2(find(isnan(q2))) = 0;
    q3(find(isnan(q3))) = 0;
    q4(find(isnan(q4))) = 0;
    %Q3(ii,jj) = sum(q3);

    %L(ii,jj) = sum(q1 + q2);
    LL = sum(q1) + sum(q2) + sum(q3) + sum(q4);
    %LL(ii,jj) = sum(    (Y.*log(p)) +     (1-Y).*log(1-p)      );
end
function [LL] = iterateLL(M,B,X,Y,pd)
    II = length(M);
    JJ = length(B);
    LL = zeros(II,JJ);
    for ii=1:II
        for jj=1:JJ
            m = M(ii);
            b = B(jj);
            LL(ii,jj) = calculateLL(m,b,X,Y,pd);
        end
    end
end
function [PRIORS] = getPriors(priorsPathName,priorsFileName,freqs,clientAge)
    % Find the priors for audiogram, based on Sumit's data
    % FPL thresholds for 357 subjects
    load([priorsPathName,priorsFileName])

    % look only over the subset needed for the current audiometer
    nFreqs = length(freqs);
    for ii=1:nFreqs
        [~,freqIndx(ii,1)] = min(abs(FREQS-freqs(ii)));
    end
    % look only at the client age, +/- 5 years in either direction
    age = clientAge;
    ageMin = 15;
    ageMax = 60;
    if age < ageMin
        warning('Client age < available prior data. Using priors from ',num2str(ageMin),' years.');
    end
    if age > ageMax
        warning('Client age > available prior data. Using priors from ',num2str(ageMax),' years.');
    end
    ageIndx = find(AGE>age-5 & AGE<age+5);

    nPoints = 5;
    minOutput = -10;
    maxOutput = 100;
    cix = linspace(median(minOutput),median(maxOutput),nPoints)';
    PRIORS.cix = cix;
    PRIORS.nPoints = nPoints;
    for jj=1:nFreqs
        q = AUDIO(ageIndx,jj);
        nanindx = find(isnan(q));
        q(nanindx) = [];
        MU = median(q);
        IQR = iqr(q);
        PRIORS.mu(1,jj) = MU;
        PRIORS.sigma(1,jj) = IQR;
        
        pd = makedist('normal','mu',MU,'sigma',IQR);
        cd = pd.cdf(cix);
        PRIORS.points(:,jj) = cd(:); % prior probabilities at each point
    end
end
function [HDI,HDIdensity,inOut] = getHDI(alpha,xx,yy)
    % Calculate highest density interval on alpha
    [p,I] = sort(yy,'ascend'); % sort densities in ascending order
    cutValue = sum(p)*alpha; % find the alpha% cut value
    cutIndx = min(find(cumsum(p)>cutValue)); % find the location of the cut value
    waterline = p(cutIndx); % this is the cutoff density
    [goodIndx] = find(yy >= waterline); % locate all values > cut
    HDI = [xx(min(goodIndx));xx(max(goodIndx))]; % determine the interval
    HDIdensity = waterline;
    inOut = (yy >= waterline);
end
function [Weights,iterations] = BS(X,maxIterations) 
    % calculate bisquares weighting
    Weights = ones(size(X));
    [Residuals,epsilon] = weightedMean(X,Weights); % initial, unweighted mean
    k = tuningConstant(Residuals);
    delta = 1; % difference between current and previous iteration errors
    minDelta = 1E-12; % discontinue iterations when improvement is less than this
    ii=1;
    while (ii <= maxIterations) && (abs(delta) > minDelta)
        ii = ii+1;
        Weights = calculateWeights(Residuals,k); % new weighting based on previously-computed residuals
        [Residuals,epsilonNew] = weightedMean(X,Weights);
        delta = max(epsilon) - max(epsilonNew);
        epsilon = epsilonNew;
    end
    iterations = ii;
end
function [Residuals,epsilon] = weightedMean(X,Weights)  
    mu = sum(X .* Weights) ./ (sum(Weights)+eps); % complexe coherent weighted mean
    mu = abs(mu); % magnitude
    Mu = repmat(mu,1,size(X,2));
    Residuals = abs(X) - Mu; % unweighted residuals
    epsilon = sum((Weights .* Residuals).^2,2); % weighted sum of squared errors
end
function [k] = tuningConstant(Residuals)
    mar = median(abs(Residuals),2); % median absolute residuals
    sigma = mar / 0.6745; % robust estimate of the standard deviation...
    %sigma = ones(size(sigma)) * median(sigma); % ...across the entire input matrix
    k = (4.685 * sigma) + eps; % tuning constant
end
function [Weights] = calculateWeights(Residuals,k)
    K = repmat(k,1,size(Residuals,2));
    Weights = (1 - (Residuals./K).^2).^2; % calculate the weights for each sample
    Mask = abs(Residuals) < K; % find the weights that are "outliers"...
    Weights = Weights .* Mask; % ...and set them to zero
end

% OLD CODE ----------------------------------------------------------------
   % for ii=1:length(M)
   %      for jj=1:length(B)
   %          m = M(ii);
   %          b = B(jj);
   %          y = m.*X - b;
   %          p = 1./(1+exp(-y)); % same as: p = exp(y) ./ (1+exp(y));
   % 
   %          % then p(y=1) = p
   %          %      p(y=0) = 1-p
   %          % and p = 1./(1+exp(-m*X+b));
   %          % when b is not already negative!
   % 
   %          L(ii,jj) = prod(p.^Y .* (1-p).^(1-Y));
   %          %LL(ii,jj) = sum((Y.*log(p)) + (1-Y).*log(1-p));
   %                      %sum(y.*log(p)) + sum((1-y).*log(1-p));
   % 
   %          q1 = Y.*log(p);
   %          q2 = (1-Y).*log(1-p); 
   %          thd = b / m;
   %          q3 = log(pd.pdf(thd));
   %          q4 = q3; 
   %          q1(find(isnan(12))) = 0;
   %          q2(find(isnan(q2))) = 0;
   %          q3(find(isnan(q3))) = 0;
   %          q4(find(isnan(q4))) = 0;
   %          Q3(ii,jj) = sum(q3);
   %          L(ii,jj) = sum(q1 + q2);
   %          LL(ii,jj) = sum(q1) + sum(q2) + sum(q3) + sum(q4);
   %          %LL(ii,jj) = sum(    (Y.*log(p)) +     (1-Y).*log(1-p)      );
   % 
   %          %dummy = m.*xxx - b;
   %          %ppp = exp(dummy) ./ (1+exp(dummy));
   %          %YYY(:,jj) = ppp(:);
   %      end
   %  end
   %  LL(find(LL==0)) = -inf;
   %  LL(find(LL==Inf)) = -inf;
   %  [row,col] = find(LL == max(LL(:))); % the location of thetaHat_MLE 
   %  m = M(row);
   %  b = B(col);
   %  thd = b / m;
   % 
   %  yyy = m.*xxx - b;
   %  ppp = exp(yyy) ./ (1+exp(yyy));
   % 
   %  plot(xxx,ppp)
   %  hold on
   %  plot(X,Y,'r*')
   %  plot(thd,0.5,'ro')
   %  plot(xxx,ppyyy,'Color',[.7 .7 .7])
   %  grid on
   %  title(['thd = ',num2str(thd)])
   % 
   %  indx = find(LL~=-Inf);
   %  LL = LL(indx);
   %  ll = exp(LL);
   %  L = L(indx);
   %  l = exp(L);
   %  BB = B;
   %  BB  = BB(indx);
   %  ll = ll ./ sum(ll);
   %  BB = BB(:);
   %  ll = ll(:);
   %  prior = pd.pdf(BB);
   % 
   %  [HDI,HDIdensity,inOut] = getHDI(0.05,BB,ll);
   %  sigma2 = HDI(2) - HDI(1);
   %  sigma = sigma2 / 2 / 2
   %  mu = thd;
   %  pd2 = makedist('normal','mu',mu,'sigma',sigma);
   %  CDF = pd2.cdf(B);
   % 
   %  plot(B,CDF,'g')
   %  line([thd-sigma*2,thd + sigma*2],[0.5 0.5],'Color',[0 1 0],'LineWidth',2)


       % 
    % for ii=1:length(M)
    %     for jj=1:length(B)
    %         m = M(ii);
    %         b = B(jj);
    %         y = m.*X - b;
    %         p = 1./(1+exp(-y)); % same as: p = exp(y) ./ (1+exp(y));
    % 
    %         % then p(y=1) = p
    %         %      p(y=0) = 1-p
    %         % and p = 1./(1+exp(-m*X+b));
    %         % when b is not already negative!
    % 
    %         L(ii,jj) = prod(p.^Y .* (1-p).^(1-Y));
    %         %LL(ii,jj) = sum((Y.*log(p)) + (1-Y).*log(1-p));
    %                     %sum(y.*log(p)) + sum((1-y).*log(1-p));
    % 
    %         q1 = Y.*log(p);
    %         q2 = (1-Y).*log(1-p); 
    %         thd = b / m;
    %         q3 = log(pd.pdf(thd));
    %         q4 = q3; 
    %         q1(find(isnan(12))) = 0;
    %         q2(find(isnan(q2))) = 0;
    %         q3(find(isnan(q3))) = 0;
    %         q4(find(isnan(q4))) = 0;
    %         Q3(ii,jj) = sum(q3);
    %         L(ii,jj) = sum(q1 + q2);
    %         LL(ii,jj) = sum(q1) + sum(q2) + sum(q3) + sum(q4);
    %         %LL(ii,jj) = sum(    (Y.*log(p)) +     (1-Y).*log(1-p)      );
    % 
    %         %dummy = m.*xxx - b;
    %         %ppp = exp(dummy) ./ (1+exp(dummy));
    %         %YYY(:,jj) = ppp(:);
    %     end
    % end
    % LL(find(LL==0)) = -inf;
    % LL(find(LL==Inf)) = -inf;

    %    [Weights,iterations] = BS(D,maxIterations);
