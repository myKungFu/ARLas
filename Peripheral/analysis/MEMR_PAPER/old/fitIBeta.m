classdef fitIBeta < handle
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Class definition fitIBeta
%
% Auditory Research Lab, The University of Iowa
% Deptartment of Communication Sciences & Disorders
% The University of Iowa
% Author: Shawn Goodman
% Date: August 9, 22
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

properties (SetAccess = private)
    scriptY % the estimated ablation curve, including the noise
    epsilon   % the error in the fit (RMSE)
    cf   % center frequency of stimulus
    Y_raw
    tau
    x
    N
    xx % dense x
    NN % number of dense samples
    sigmaY_raw
    sigmaY
    Y
    Yhat
    sigmaHat
    
    Y_raw2
    Y2
    p2
    residualMag
    residualPhi
    
    EN
    YHAT
    kMax
    kMin
    ablationStart
    ablationFinish
    indxAblStart
    indxAblFinish
    b
    d
    searchIterations
    optimizationParameters
    params %parameters
    slope
    offset
    resampledParams % jackknifed parameters
    resampledPointEst
    resampledCurves
    curves
    jackedI
    jackedAC
    tauTau
    pdfScaling
    pointEst
    interval
    I
    I_flp
    scriptA
    Ltotal
    flowRate
    dd
    predCf
    predCfCf
    nIterations % number of bootstrap iterations
    p % emission phase
end
properties (SetAccess = public)
    A    % A variable for beta function
    B    % B variable for beta function   
end
methods
    function obj = fitIBeta(time,data,A,B)
        % assign variables ------------------------------------------------
        obj.tau = time(:);         % 
        obj.Y_raw = data(:);       % magnitudes in original units (uV or uPa)
        obj.N = length(obj.Y_raw); % number of data samples to fit
        obj.A = A;                 % initial estimate of A
        obj.B = B;                 % initial estimate of B
        
        % create vector x for the beta distribution -----------------------
        % This is Eq. 5 in the paper.
        obj.slope = (obj.tau(end)  - obj.tau(1));
        obj.offset = obj.tau(1);
        obj.x = (obj.tau - obj.offset) / obj.slope;
        % Note: convert back to tau thus (Eq. 6):
        % obj.tau = obj.slope * obj.x + obj.offset;
        obj.NN = 500; % create vector xx for densly-spaced distributions
        spacing = (obj.tau(end)-obj.tau(1)) / (obj.NN-1);
        obj.tauTau = ((0:1:obj.NN-1)' * spacing) + obj.tau(1);
        obj.xx = (obj.tauTau - obj.offset) / obj.slope;

        % Scale the raw data so maximum of underlying ablation curve A = 1 
%         % This is Eq. 8-11 in the paper.
%         obj.ablationStart = ablationStart;         % time where ablation starts (visually determined)
%         [~,indxAblStart] = min(abs(obj.tau - obj.ablationStart)); % corresponding indices
%         signalStub = obj.Y_raw(1:indxAblStart);    % the stub of the vectors before ablation starts
%         noiseStub = obj.sigmaY_raw(1:indxAblStart);
%         Ybar = mean(signalStub);      % Eq. 8
%         avgN = length(noiseStub);    % number of samples in the average
%         correction = (gamma(avgN)*sqrt(avgN)) / (gamma(avgN + 0.5));    % bias correction for the MLE (max likelihood estimate)
%         sigmaHat = correction * sqrt((1/(2*avgN)) * sum(noiseStub.^2)); % Eq 9
%         EN = sigmaHat * sqrt(pi/2);     % Eq. 10  The expected noise value
%         obj.kMax = sqrt(Ybar^2 - EN^2); % Eq. 11  the scaling factor

%         % apply the scaling factor ----------------------------------------
%         obj.Y = obj.Y_raw / obj.kMax;               % scaled magnitude data -- residual absent
%         obj.sigmaY = obj.sigmaY_raw / obj.kMax;     % scaled noise floor
%         obj.Y2 = obj.Y_raw2 / obj.kMax;               % scaled magnitude data -- residual present
%         
%         % The ablation is not assumed to be complete; calcuate a second scaling factor to acount for incomplete ablation.
%         % This is Eq. 12 in the paper.
%         obj.ablationFinish = ablationFinish; % time where ablation finishes (visually determined)
%         [~,indxAblFinish] = min(abs(obj.tau - obj.ablationFinish));       
%         signalStub = obj.Y(indxAblFinish:end); 
%         noiseStub = obj.sigmaY(indxAblFinish:end);
%         avgN = length(noiseStub); % number of samples in the average
%         correction = (gamma(avgN)*sqrt(avgN)) / (gamma(avgN + 0.5));
%         Ybar = mean(signalStub);
%         sigmaHat = correction * sqrt((1/(2*avgN)) * sum(noiseStub.^2));
%         EN = sigmaHat * sqrt(pi/2);
%         if Ybar > EN
%             obj.kMin = sqrt(Ybar^2 - EN^2); 
%         else
%             obj.kMin = 0;
%         end
% 
%         % calculate the expected value of the noise floor from a Rayleigh distributuion
%         % Note: when avgN = 1, ER == sigmaY.
%         % so here, the expected average is just the noise floor itself.
%         % This is for Eq. 13 in the paper.
%         avgN = 1;  % number of samples in the average
%         correction = (gamma(avgN)*sqrt(avgN)) / (gamma(avgN + 0.5));
%         obj.sigmaHat = correction * sqrt((1/(2*avgN)) * obj.sigmaY.^2);
%         obj.EN = obj.sigmaHat * sqrt(pi/2);
        
        obj.Yhat = obj.Y; % on the first calculation, Yhat is the actual data

        % initialize optimization parameters
        obj.epsilon = 1e10; % initialize error value
        obj.searchIterations = 10000;
        obj.optimizationParameters = optimset('Display','off','MaxFunEvals',obj.searchIterations,'MaxIter',obj.searchIterations);
        obj.interval = 0.95; % HDI (highest density interval)
        
        obj.params.alpha = zeros(3,1);
        obj.params.beta = zeros(3,1);
        obj.resampledParams.alpha = zeros(obj.nIterations,1);
        obj.resampledParams.beta = zeros(obj.nIterations,1);
        
        obj.pointEst.mean = zeros(3,1);
        obj.pointEst.median = zeros(3,1);
        obj.pointEst.mode = zeros(3,1);
        obj.pointEst.sd = zeros(3,1);
        obj.pointEst.variance = zeros(3,1);
        obj.pointEst.skewness = zeros(3,1);
        obj.pointEst.kurtosis = zeros(3,1);
        obj.pointEst.Q0 = zeros(3,1);
        obj.pointEst.Q1 = zeros(3,1);
        obj.pointEst.Q2 = zeros(3,1);
        obj.pointEst.Q3 = zeros(3,1);
        obj.pointEst.Q4 = zeros(3,1);
        obj.pointEst.meanX = zeros(3,1);
        obj.pointEst.medianX = zeros(3,1);
        obj.pointEst.modeX = zeros(3,1);
        obj.pointEst.varianceX = zeros(3,1);
        obj.pointEst.Q1X = zeros(3,1);
        obj.pointEst.Q2X = zeros(3,1);
        obj.pointEst.Q3X = zeros(3,1);
        obj.pointEst.IQRX = zeros(3,1);
        obj.pointEst.Q0X = zeros(3,1);
        obj.pointEst.Q4X = zeros(3,1);
        obj.pointEst.skewnessX = zeros(3,1);
       
        obj.resampledPointEst.mean = zeros(obj.nIterations,1);
        obj.resampledPointEst.median = zeros(obj.nIterations,1);
        obj.resampledPointEst.mode = zeros(obj.nIterations,1);
        obj.resampledPointEst.sd = zeros(obj.nIterations,1);
        obj.resampledPointEst.variance = zeros(obj.nIterations,1);
        obj.resampledPointEst.skewness = zeros(obj.nIterations,1);
        obj.resampledPointEst.kurtosis = zeros(obj.nIterations,1);
        obj.resampledPointEst.Q0 = zeros(obj.nIterations,1);
        obj.resampledPointEst.Q1 = zeros(obj.nIterations,1);
        obj.resampledPointEst.Q2 = zeros(obj.nIterations,1);
        obj.resampledPointEst.Q3 = zeros(obj.nIterations,1);
        obj.resampledPointEst.Q4 = zeros(obj.nIterations,1);
        obj.resampledPointEst.meanX = zeros(obj.nIterations,1);
        obj.resampledPointEst.medianX = zeros(obj.nIterations,1);
        obj.resampledPointEst.modeX = zeros(obj.nIterations,1);
        obj.resampledPointEst.varianceX = zeros(obj.nIterations,1);
        obj.resampledPointEst.sdX = zeros(obj.nIterations,1);
        obj.resampledPointEst.Q1X = zeros(obj.nIterations,1);
        obj.resampledPointEst.Q2X = zeros(obj.nIterations,1);
        obj.resampledPointEst.Q3X = zeros(obj.nIterations,1);
        obj.resampledPointEst.IQRX = zeros(obj.nIterations,1);
        obj.resampledPointEst.Q0X = zeros(obj.nIterations,1);
        obj.resampledPointEst.Q4X = zeros(obj.nIterations,1);
        obj.resampledPointEst.skewnessX = zeros(obj.nIterations,1);
    end
    function [] = run(obj)
        obj.curveFit  % get the best fitting parameters A and B, given the data
            if obj.A<= 1 || obj.B<=1
                error('Alpha and Beta cannot be <=1 in the current setup!')
            end
            obj.params.alpha(1,1) = obj.A;
            obj.params.beta(1,1) = obj.B;
            obj.params.scriptY = (obj.scriptY); % the expected magnitudes in Y
            obj.params.scriptA = (obj.scriptA); % the ablation curve
            obj.pointEst = estimateParameters(obj,obj.pointEst,1); % point estimates
            obj.curves.cdf = flipud(betacdf(obj.xx,obj.params.alpha(1),obj.params.beta(1)));
            obj.curves.pdf = flipud(betapdf(obj.xx,obj.params.alpha(1),obj.params.beta(1)));
        obj.bootstrap
        obj.estimateIntervals
    end
    function [params] = estimateParameters(obj,params,ii)
        % the following are in terms of obj.x -----
        % they need to be re-expressed in terms of obj.tau
        params.meanX(ii,1) = obj.A / (obj.A + obj.B);
        params.medianX(ii,1) = (obj.A-(1/3)) / (obj.A + obj.B - (2/3)); % for B>1
        params.modeX(ii,1) = (obj.A - 1) / (obj.A + obj.B - 2); % for B>1
        params.varianceX(ii,1) = (obj.A * obj.B) / ((obj.A + obj.B)^2*(obj.A + obj.B + 1)); % variance
        params.sdX(ii,1) = sqrt(params.varianceX(ii,1)); % standard deviation
        params.Q1X(ii,1) = betainv(0.25,obj.A,obj.B); % first quartile
        params.Q2X(ii,1) = betainv(0.50,obj.A,obj.B); % second quartile (median)
        params.Q3X(ii,1) = betainv(0.75,obj.A,obj.B); % third quartile
        params.IQRX(ii,1) = params.Q3X(ii,1) - params.Q1X(ii,1); % inter-quartile range
        params.Q0X(ii,1) = params.Q1X(ii,1) - (params.IQRX(ii,1) * 1.5); % lower whisker
        params.Q4X(ii,1) = params.Q3X(ii,1) + (params.IQRX(ii,1) * 1.5); % upper whisker
        params.skewnessX(ii,1) = (2*(obj.B - obj.A)*sqrt(obj.A + obj.B + 1)) / ((obj.A + obj.B + 2)*sqrt(obj.A * obj.B));
        % convert the above values to time
        params.Q4(ii,1) = obj.slope * (1-params.Q0X(ii,1)) + obj.offset; 
        params.Q3(ii,1) = obj.slope * (1-params.Q1X(ii,1)) + obj.offset;
        params.Q2(ii,1) = obj.slope * (1-params.Q2X(ii,1)) + obj.offset;
        params.Q1(ii,1) = obj.slope * (1-params.Q3X(ii,1)) + obj.offset;
        params.Q0(ii,1) = obj.slope * (1-params.Q4X(ii,1)) + obj.offset;
        params.IQR(ii,1) = params.Q3(ii,1) - params.Q1(ii,1);
        % 
        params.mean(ii,1) = obj.slope * (1-params.meanX(ii,1)) + obj.offset;
        params.median(ii,1) = obj.slope * (1-params.medianX(ii,1)) + obj.offset;
        params.mode(ii,1) = obj.slope * (1-params.modeX(ii,1)) + obj.offset;
        params.sd(ii,1) = sqrt(params.varianceX(ii,1)) * abs(obj.slope);
        params.variance(ii,1) = params.sd(ii,1)^2;
        params.skewness(ii,1) = -params.skewnessX(ii,1); % account for the flip
        % these parameters are okay as they are
        params.kurtosis(ii,1) = (6*((obj.A-obj.B)^2*(obj.A+obj.B+1)-obj.A*obj.B*(obj.A+obj.B+2))) / (obj.A*obj.B*(obj.A+obj.B+2)*(obj.A+obj.B+3));
    end
    function [] = curveFit(obj)
        [obj,fval,exitFlag,output] = fminsearchBETA(@findOptimalBeta,obj,obj.optimizationParameters);
        if exitFlag ~= 1
            disp('WARNING: Fit failed to converge.')
        end    
    end
    function [] = calculate(obj)
        obj.I = betacdf(obj.x,obj.A,obj.B);    % Eq. 13, bottom (I = regularized incomplete beta function)
        obj.I_flp = flipud(obj.I);             % Eq. 13, middle (I_flp = fliped version of I) 
        obj.scriptA = (1-obj.kMin)*obj.I_flp + obj.kMin; % Eq. 13, middle (scriptA = the underlying ablation curve)
        obj.scriptY = sqrt(obj.scriptA.^2 + obj.EN.^2); % Eq. 13, top
        obj.epsilon = sum((obj.Yhat - obj.scriptY).^2); % Eq. 13, top (epsilon is the error in the fit)
    end
    function [] = bootstrap(obj)
        obj.jackedI = zeros(obj.NN,obj.nIterations); % cdf
        obj.jackedAC = zeros(obj.N,obj.nIterations); % ablation curve
        obj.YHAT = zeros(obj.N,obj.nIterations); % ablation curve
        
        counter = 1;
        disp(' ')
        disp('Starting parametric bootstrap...')
        disp('   1')
        for kk=1:obj.nIterations
            if mod(kk,500) == 0
                disp(['   ',num2str(kk)])
            end
            
            % in this fitting scheme, alpha and beta cannot, be less than
            % 1. If this happens, throw out and re-do the bootstrap
            % iteration
            ok = 0;
            while ok == 0
                obj.A = obj.params.alpha; % always start A at the best estimate
                obj.B = obj.params.beta;  % always start B at the best estimate
                % Eq. 15: generate new instance of noise 
                n1 = random('normal',0,obj.sigmaY,obj.N,1);
                n2 = random('normal',0,obj.sigmaY,obj.N,1);
                obj.Yhat = sqrt((obj.params.scriptY + n1).^2 + n2.^2);

                obj.YHAT(:,kk) = obj.Yhat; % add to the matrix of collected Yhat values.

                obj.curveFit
                  if obj.A>1 || obj.B>1
                    ok = 1;
                  end
            end
            obj.resampledParams.alpha(kk,1) = obj.A;
            obj.resampledParams.beta(kk,1) = obj.B;    
            obj.resampledPointEst = estimateParameters(obj,obj.resampledPointEst,kk);
            obj.jackedI(:,kk) = flipud(betacdf(obj.xx,obj.A,obj.B)); % cdf
            obj.jackedAC(:,kk) = flipud(obj.scriptY);
            
            counter = counter + 1; % increment the counter
            if counter > obj.N % make counter circular
                counter = 1;
            end
        end
        disp('...bootstrap finished.')
    end
    function [] = estimateIntervals(obj)
        disp('Calculating CIs...')
        
        names = fieldnames(obj.pointEst);
        N = length(names);
        for ii=1:N
            y = getfield(obj.resampledPointEst,char(names(ii)));
            hdi = obj.hdiFinder(y);
            point = getfield(obj.pointEst,char(names(ii)));
            point = point(1);
            obj.pointEst = setfield(obj.pointEst,char(names(ii)),[hdi(1);point;hdi(2)]);
        end
        
        y = obj.resampledParams.alpha;
        hdi = obj.hdiFinder(y);
        point = obj.params.alpha(1,1);
        obj.params.alpha = [hdi(1);point;hdi(2)];
        
        y = obj.resampledParams.beta;
        hdi = obj.hdiFinder(y);
        point = obj.params.beta(1,1);
        obj.params.beta = [hdi(1);point;hdi(2)];  
        
        % get middle 95% of cdf
        subA1 = obj.resampledParams.alpha;
        subB1 = obj.resampledParams.beta;
        [q1,indx] = sort(subA1.*subB1);
        subA1 = subA1(indx);
        subB1 = subB1(indx);
        [indx1,~] = AR(q1,'mild',0);
        subA1(indx1) = [];
        subB1(indx1) = [];
        nn = length(subA1);
        resampledI = zeros(length(obj.xx),nn);
        for ii=1:nn
            resampledI(:,ii) = betacdf(obj.xx,subA1(ii),subB1(ii));
        end        
        obj.resampledCurves.I = flipud(resampledI);
        %obj.resampledCurves.Ispread = max(subA1.*subB1);
        obj.resampledCurves.Ispread = var(sqrt(obj.resampledParams.alpha.* obj.resampledParams.beta));
        
        disp('...CI calculation finished.')
    end
    function [hdi] = hdiFinder(obj,y)
        % find the highest density interval
        if sum(diff(y)) == 0
            hdi = [min(y);max(y)];
            return
        end
        kern = fitdist(y,'Kernel'); % fit a kernel density function to the data
        xmin = min(y);
        xmax = max(y);
        nSteps = 1000;
        stepSize = ((xmax-xmin)) / nSteps;
        xvals = (xmin:stepSize:xmax)';
        yvals = pdf(kern,xvals); 
        yy = sort(yvals,'ascend');
        cutoff = round(length(yvals) * (1-obj.interval));
        try
            lowerBound = yy(cutoff);
        catch
            keyboard
        end
        indx = find(yvals < lowerBound);
        yvals(indx) = []; % delete all values below lower bound
        xvals(indx) = [];
        [~,indx] = sort(yvals,'ascend'); % re-sort
        xx = xvals(indx);
        hdi(1,1) = min(xx);
        hdi(1,2) = max(xx);
    end
end
end
