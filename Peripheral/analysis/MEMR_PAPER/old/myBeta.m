classdef myBeta < handle
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Class definition myBeta
% For use with DWexplore
% Used in finding best beta distribution fit to data
%
% Auditory Research Lab, The University of Iowa
% Deptartment of Communication Sciences & Disorders
% The University of Iowa
% Authors: Shawn Goodman
% Date: July 26, 2018
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

properties (SetAccess = private)
    nMeasurements
    
    amplitudes
    fs
    len
    nSamples
    signal
    
    beta_pdf  % the underlying beta probability density function, given A and B
    beta_cdf  % the underlying beta-based ablation curve (similar to a cdf)
    ablationCurve % the estimated ablation curve, including the noise
    epsilon   % the error in the fit (RMSE)
    
    weights % weighting for the jackknife approach
    pdfFamily
    cdfFamily
    ablationFamily
    epsilonFamily
    AFamily
    BFamily
    
    
end
properties (SetAccess = public)
    time
    N

    x
    data % actual recorded data. Vector of amplitudes across time
    noiseAmplitude % estimate of noise amplitude (rms)
    A    % A variable for beta function
    B    % B variable for beta function
    
    cf   % center frequency of stimulus
    iterations % number of times this object has been iterated
    walk_A
    walk_B
    walk_E
    prior
    LH
    muA
    muB
    
    Curves
    curve
    ci
    Residuals
    residuals
    dy
end
   
methods
    function b = myBeta(time,data)
        % Assume 60 minutes recording time and measurements taken every 1 minute.
        % Then there are 60 measurements.
        b.data = data(:);
        b.time = time(:);
        b.N = length(b.data);
        for ii=1:b.N
            b.x(ii,1) = ii/(b.N+1); % time vector (0 to 1)
        end
        
        %b.x = linspace(0,1,b.nMeasurements)'; % time vector (0 to 1)
        % initialize A and B values, and noise level
        b.A = 10;
        b.B = 7;
        b.noiseAmplitude = eps;
        
        b.fs = 96000;
        b.len = 0.05;
        b.cf = 1000; % center frequency (Hz)
        b.nSamples = round(b.fs*b.len);
        t = (0:1:b.nSamples-1)'/b.fs;
        b.signal  = cos(2*pi*b.cf*t); % unit amplitude signal
        b.epsilon = 1e10; % initialize error value
        b.iterations = 0;
        b.weights = ones(size(data));
    end
    
    function [] = calculateDistribution(obj)
        obj.beta_pdf = betapdf(obj.x,obj.A,obj.B); % this specifies the signal levels
            %obj.beta_cdf = cumsum(obj.beta_pdf);       % intigrate to get the cdf.
        obj.beta_cdf = betacdf(obj.x,obj.A,obj.B);
        obj.beta_cdf = flipud(obj.beta_cdf);
            %obj.beta_cdf = obj.beta_cdf / max(obj.beta_cdf); % scale
        
        %----------------------------------------------
%         obj.ablationCurve = zeros(obj.nMeasurements,1);
%         for ii=1:obj.nMeasurements
%             ss = obj.signal * obj.beta_cdf(ii); % scale unit signal 
%             nn = randn(obj.nSamples,1) * obj.noiseAmplitude; % generate a new noise floor
%             sn = ss + nn; % signal plus noise
%             obj.ablationCurve(ii,1) = sqrt(mean(sn.^2));  % signal plus noise
%         end
%         %weights = obj.ablationCurve / obj.noiseAmplitude; % signal to noise ratio
%         obj.ablationCurve = obj.ablationCurve / sqrt(.5);
        %------------------------------------------------
        B = (sqrt(0.5) * obj.beta_cdf);
        N = obj.noiseAmplitude;
        obj.ablationCurve = sqrt(B.^2 + N.^2)/sqrt(0.5);
        %-----------------------------------------------
        
        obj.epsilon = sqrt(mean((obj.ablationCurve - obj.data).^2)); % error of the fit (unweighted)
        %weights = weights .* obj.weights; % apply the jackknife
        %obj.epsilon = sqrt(sum((obj.ablationCurve - obj.data).^2 .* weights) / sum(weights)); % snr-weighted error of the fit
        
        obj.beta_pdf = flipud(obj.beta_pdf);
        [bmax,imax] = max(obj.beta_pdf); % find the location of the maximum
        target = obj.beta_cdf(imax); % target height for pdf
        obj.beta_pdf = target * (obj.beta_pdf / max(obj.beta_pdf)); % rescale
        obj.iterations = obj.iterations + 1;
    end
    
    function [] = plotResults(obj)
        plot(obj.data,'k')
        hold on
        plot(obj.ablationCurve,'g')
        plot(obj.beta_cdf,'b')
        plot(obj.beta_pdf,'r')
        xlim([1,obj.nMeasurements])
        ylim([0,1.1])
        xlabel('Time (relative)','FontSize',12)
        ylabel('Ablation (relative)','FontSize',12)
        legend('data','ablation curve','Beta cdf','Beta pdf')
    end
    
    function [] = splineAblation(obj)
        % use a weighted-smoothing spline to fit the data.
        %keyboard
        %stepSize = .01; % make a densely-spaced x-axis
        %xx = (obj.x(1):stepSize:obj.x(end))';
        
        N = length(obj.data);
        x = linspace(1,N,N)';
        x1 = [x;N+1]';
        sm = 0.5;
        cs = csaps(x,obj.data,sm); % the smoothed cubic spline
        curve = fnval(cs,x);
        residuals = curve - obj.data;
        
        % generate kernel density function of residuals
        pd = fitdist(residuals,'Kernel');
        r = random(pd,N,1);
        %xxx = linspace(-0.1,0.1,500)';
        %yyy = pdf(pd,xxx);
        
        % family of splines using original, best-fitting spline with noise 
        % added from the estimated distribution
        nIterations = 1000;
        Curves = zeros(N,nIterations);
        Residuals = zeros(N,nIterations);
        Residuals2 = zeros(N,nIterations);
        counter = 1;
        for ii=1:nIterations
            w = ones(N,1);
            w(counter) = 0;
            noise = random(pd,N,1);
            dummy = curve + noise;
            csDummy = csaps(x,dummy,sm,x1,w); % the smoothed cubic spline
            %csDummy = csaps(x,dummy,sm,x1,w); % the smoothed cubic spline
            %Curve(:,ii) = fnval(csDummy,x);
            %Residual(:,ii) = Curve(:,ii) - dummy;
            Curves(:,ii) = csDummy(1:end-1);
            Residuals(:,ii) = Curves(:,ii) - dummy;
            Residuals2(:,ii) = obj.data - dummy;
            counter = counter + 1;
            if counter > N
                counter = 1;
            end
        end
        ci = zeros(N,2);
        for ii=1:N
            q = sort(Curves(ii,:));
            indx = round(nIterations * ((0.05/2)));
            ci(ii,1) = q(indx);
            indx = round(nIterations * (1-(0.05/2)));
            ci(ii,2) =  q(indx);
        end
        
%         figure
%         subplot(2,1,1)
%             plot(obj.data,'.k')
%             hold
%             plot(curve,'r','LineWidth',2)
%             plot(ci,'b')
%             xlim([1,N])
%         subplot(2,1,2)
%             plot(residuals)
%             xlim([1,N])

        obj.Curves = Curves;
        obj.curve = mean(Curves,2);
        obj.ci = ci;
        obj.Residuals = Residuals;
        obj.residuals = residuals;
        
        cs = csaps(x,curve,sm); % the smoothed cubic spline
        d = fnder(cs); % calculate the derivative
        obj.dy = fnval(d,x); % evaluate the derivative on the densely-spaced x-axis        
        
        
return        
        
        yy = csaps(x,y,sm,xx,w); % the smoothed spline, densly-spaced estimate
        cs = csaps(x,y,sm); % the smoothed cubic spline
            %d = fnder(cs); % calculate the derivative
            %dd = fnval(d,xx); % evaluate the derivative on the densely-spaced x-axis
        
        keyboard
            
        [~,indx1] = min(abs(window(1)-xx)); % cut down to window
        [~,indx2] = min(abs(window(2)-xx));
        dy = dd(indx1:indx2);
        xxx = xx(indx1:indx2);
        yyy = yy(indx1:indx2);
        [dHat,indx] = min(dy); % find location of the minimum value

        xHat = xxx(indx); % the location of the steepest part of the slope
        yHat = yyy(indx);  % the y-value of the steepest part of the slope


        yyy2 = yyy.^2; % power in the window
        maxPower = max(yyy2);
        halfPower = maxPower / 2;
        done = 0;
        counter = 2;
        N = length(yyy2);
        try
        while done == 0
           if yyy2(counter) < halfPower
                hpIndx = counter;
                done = 1;
           else
               counter = counter+1;
           end
           if counter >= N
               hpIndx = N;
               done = 1;
           end
        end
        catch
            keyboard
        end
        xHat2 = xxx(hpIndx); % the x-value of the half-power point
        yHat2 = yyy(hpIndx);  % the y-value of the half-power point
    end
    
    function [] = chertoffFit(obj)
%         General model:
%      f(x) = A + ((K-A)/((C+Q*exp(-B*x)).^(1/nu)))
% Coefficients (with 95% confidence bounds):
%        A =         0.2  (0.1623, 0.2377)
%        B =      0.2841  (0.1444, 0.4238)
%        C =           1  (-17, 19)
%        K =      0.9322  (-7.1e+06, 7.1e+06)
%        Q =       19.99  (-558.9, 598.9)
%        nu =   1.848e-06  (-5.954e-05, 6.323e-05)
% 
% Goodness of fit:
%   SSE: 1.803
%   R-square: 0.8823
%   Adjusted R-square: 0.8761
%   RMSE: 0.1378

        
        
        obj.beta_pdf = betapdf(obj.x,obj.A,obj.B); % this specifies the signal levels
            %obj.beta_cdf = cumsum(obj.beta_pdf);       % intigrate to get the cdf.
        obj.beta_cdf = betacdf(obj.x,obj.A,obj.B);
        obj.beta_cdf = flipud(obj.beta_cdf);
            %obj.beta_cdf = obj.beta_cdf / max(obj.beta_cdf); % scale
        obj.ablationCurve = zeros(obj.nMeasurements,1);
        for ii=1:obj.nMeasurements
            ss = obj.signal * obj.beta_cdf(ii); % scale unit signal 
            nn = randn(obj.nSamples,1) * obj.noiseAmplitude; % generate a new noise floor
            sn = ss + nn; % signal plus noise
            obj.ablationCurve(ii,1) = sqrt(mean(sn.^2));  % signal plus noise
        end
        
        weights = obj.ablationCurve / obj.noiseAmplitude; % signal to noise ratio
        obj.ablationCurve = obj.ablationCurve / sqrt(.5);
            %obj.epsilon = sqrt(mean((obj.ablationCurve - obj.data).^2)); % error of the fit (unweighted)
        weights = weights .* obj.weights; % apply the jackknife
        obj.epsilon = sqrt(sum((obj.ablationCurve - obj.data).^2 .* weights) / sum(weights)); % snr-weighted error of the fit
        
        obj.beta_pdf = flipud(obj.beta_pdf);
        [bmax,imax] = max(obj.beta_pdf); % find the location of the maximum
        target = obj.beta_cdf(imax); % target height for pdf
        obj.beta_pdf = target * (obj.beta_pdf / max(obj.beta_pdf)); % rescale
        obj.iterations = obj.iterations + 1;
        
        
    end
end
end

