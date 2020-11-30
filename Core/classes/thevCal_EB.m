classdef thevCal_EB < handle
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Class definition thevCal_EB (extended bandwidth) 
% For use with ARLas (Auditory Research Laboratory auditory software)
% Used in performing Thevenin source calibration.
%
% 1) CALCULATE PL (load pressure; pressure measured in the tubes)
% 2) CALCULATE ZLi (ideal cavity impedance of the tubes)
%  2a) CALCULATE Zpw (plane wave impedance of tubes), initial estimate
% 3a) INITIAL CALCULATION OF PS & ZS (source pressure and impedance) using Zpw only
%  2b) CALCULATE zew (evanescent wave impedance of tubes)    
%  2c) CALCULATE zfl (flow loss impedance of tubes)    
%  2d) CALCULATE ZLi (total ideal load impedance of tubes)
% 3b) FINAL CALCULATION OF PS & ZS (source pressure and impedance) using total impedanc (ZLi)
% 4) CALCULATE ZL (actual cavity impedance of the tubes, given PS & ZS)
% 5) CALCULATE ERROR
%    still to do-->estimate and remove crosstalk error
% 6) PLOT RESULTS
%
% Auditory Research Lab, The University of Iowa
% Deptartment of Communication Sciences & Disorders
% The University of Iowa
% Authors: Shawn Goodman & Sawyer Goetz
% Date: March 18 - June 22, 2018
% Updated: October 12-22, 2019 -- ssg
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

properties (SetAccess = private)
    fs
    cavityTemperature
    cavityDiameter
    cavityLengths_nominal
    stimulus
    recordings
    recordingsOrig
    
    nCavities
    nfft
    freq
    fminIndx
    fmaxIndx
    PL
    PL_timeDomain
    
    z0
    Zpw
    d
    c
    omega
    
    nFreqs
    ZS
    PS

    zfl
    
    ZLi
    
    ZL
    
    epsilon
    optSwitch % switch to optimize lengths ('L') or tau ('T')
    
end

properties (SetAccess = public)
    cavityLengths_optimal
    tau
    zew
    wn
    fmin
    fmax
    stimOrig
    stimFilt
    cut1
    cut2
    applyMicCorrection
    micCorrection
    lastCalculatedUsingMicCorrection
    nMin
end

methods
    function t = thevCal_EB(recordingParams,stimulus,recordings)
        t.fs = recordingParams.fs; % sampling rate (Hz)
        t.cavityTemperature = recordingParams.cavityTemperature; % degrees C
        t.cavityDiameter = recordingParams.cavityDiameter; % (cm)
        t.cavityLengths_nominal = recordingParams.cavityLengths_nominal; % intended cavity lengths (cm)
        t.fmin = recordingParams.fmin; % minimum frequency (Hz) over which to calibrate
        t.fmax = recordingParams.fmax; % maximum frequency (Hz) over which to calibrate
        t.stimulus = stimulus;
        t.recordings = recordings;
        t.recordingsOrig = recordings; % original recordings; no mic correction
        t.stimOrig = recordingParams.stimOrig; % original, electrical drive for the stimulus 
        t.cut1 = recordingParams.cut1;
        t.cut2 = recordingParams.cut2;
        t.stimFilt = recordingParams.stimFilt;
        t.micCorrection = recordingParams.micCorrection; % fir filter microphone correction; scalar value 1 if not provided
        t.applyMicCorrection = 1;
    end
    function [] = calculate(obj)
        if obj.applyMicCorrection == 1 % if apply microphone correction requested
            obj.recordings = fastFilter(obj.micCorrection,obj.recordingsOrig); % apply microphone correction
            if obj.micCorrection == 1 % in this case, even though technically "applied", 
                                      % the mic correction does nothing, so
                                      % same as if not applied
                obj.lastCalculatedUsingMicCorrection = 0;
            else
                obj.lastCalculatedUsingMicCorrection = 1;
            end
        else
            obj.recordings = obj.recordingsOrig; % no microphone correction
            obj.lastCalculatedUsingMicCorrection = 0;
        end
        %obj.cavityLengths_nominal = [2.9; 3.6; 4.15; 5.2; 6.9]' -0.5775;
        doParts = 0;

        if doParts == 0
            obj.calculatePL;                   % 1)  CALCULATE PL (load pressure; pressure measured in the tubes)
            obj.nFreqs = size(obj.PL,1); % number of frequencies in analysis range
            obj.ZS = zeros(obj.nFreqs,1); % initialize source impedance vector
            obj.PS = zeros(obj.nFreqs,1); % initialize source pressure vector
            % DON'T DO THIS ANYMORE--CAUSES MORE PROBLEMS THAN IT SOLVES!
                %optimizeCavityLengths(obj);         %     Make a first estimate cavity lengths based on autocorrelation
            obj.d = obj.cavityTemperature - 26.85;
            obj.c = 3.4723e4 * (1 + 0.00166 * obj.d); % speed of sound at the cavity temperature (cm/s)
            %if isempty(obj.cavityLengths_optimal) % only set back to nominal if no previous optimal exists
               obj.cavityLengths_optimal = obj.cavityLengths_nominal';
            %end
            obj.zew = zeros(obj.fmaxIndx - obj.fminIndx+1,1);  % Initialize values
            obj.zfl = zeros(obj.fmaxIndx - obj.fminIndx+1,1);
            obj.ZS = zeros(obj.fmaxIndx - obj.fminIndx+1,1);
            obj.ZLi = zeros(obj.fmaxIndx - obj.fminIndx+1,obj.nCavities);
            obj.tau = zeros(1,obj.nCavities);

            op1 = optimset('Display','off','MaxIter',1000);
            %op1 = optimset('Display','off','MaxIter',500,'MaxFunEvals',500,'TolX',10e-4,'TolFun',10e-4);
            obj.optSwitch = 'L'; % calculating optimal lengths
            [obj,fval,exitFlag,output] = fminsearchOBJ(@findOptimalLengths,obj,op1);
            %        calculateZpw(obj);                  % 2a) CALCULATE Z: Zpw (plane wave impedance of tubes)
            %        calculateThevSource(obj,obj.Zpw);   % 3a) CALCULATE PS & ZS (source pressure and impedance)
            %        calculatezew(obj);                  % 2b) CALCULATE Z: zew (evanescent wave impedance of tubes)
            %        calculatezfl(obj);                  % 2c) CALCULATE Z: zfl (flow loss impedance of tubes)
            %        calculateZLi(obj);                  % 2d) CALCULATE Z: ZLi (total ideal load impedance = Zpw + zew + zfl)

                                                %     then further refine cavity lenghts using a minimization routine

            op2 = optimset('Display','off','MaxIter',1000);
            %op2 = optimset('Display','off','MaxIter',200,'MaxFunEvals',200,'TolX',10e-4,'TolFun',10e-4);
            obj.optSwitch = 'T'; % calculate optimal tau
            [obj,fval,exitFlag,output] = fminsearchOBJ(@findOptimalLengths,obj,op2);
            %obj.calculateZpw
            %obj.calculateZLi
            %obj.calculateThevSource(obj.Zpw)

            %       calculateThevSource(obj,obj.ZLi);   % 3b) CALCULATE PS & ZS (source pressure and impedance)
            calculateZL(obj);                   % 4)  CALCULATE ZL (actual cavity impedance of the tubes, given PS & ZS)
            removeCrosstalk(obj)
            
            obj.ZS = smoother(obj.freq,obj.ZS); % apply smoothing across frequency
            obj.PS = smoother(obj.freq,obj.PS);
            obj.ZL = smoother(obj.freq,obj.ZL);
            obj.PL = smoother(obj.freq,obj.PL);
            obj.ZLi = smoother(obj.freq,obj.ZLi);            
            
            calculateError(obj);                % 5)  CALCULATE ERROR 
            plotResults(obj);                   % 6)  PLOT RESULTS
        else
            starts =   [ 200,2000, 4000, 8000,12000,16000];
            finishes = [4000,8000,12000,16000,20000,32000];
            NN = length(starts);
            for kk=1:NN
                disp(['testing range ',num2str(kk),' of ',num2str(NN)])
                obj.fmin = starts(kk);
                obj.fmax = finishes(kk);

                obj.calculatePL;                   % 1)  CALCULATE PL (load pressure; pressure measured in the tubes)
                obj.nFreqs = size(obj.PL,1); % number of frequencies in analysis range
                obj.ZS = zeros(obj.nFreqs,1); % initialize source impedance vector
                obj.PS = zeros(obj.nFreqs,1); % initialize source pressure vector



                if kk== 1
                optimizeCavityLengths(obj);         %     Make a first estimate cavity lengths based on autocorrelation
                else                                          %     Further refine cavity lenghts using a minimization routine
                    obj.cavityLengths_optimal = obj.cavityLengths_nominal';
                end

                obj.zew = zeros(obj.fmaxIndx - obj.fminIndx+1,1);  % Initialize values
                obj.zfl = zeros(obj.fmaxIndx - obj.fminIndx+1,1);
                obj.ZS = zeros(obj.fmaxIndx - obj.fminIndx+1,1);
                obj.ZLi = zeros(obj.fmaxIndx - obj.fminIndx+1,obj.nCavities);
                obj.tau = zeros(1,obj.nCavities);

                op=optimset('Display','off','MaxIter',9000);
                obj.optSwitch = 'L';
                if kk==1
                [obj,fval,exitFlag,output] = fminsearchOBJ(@findOptimalLengths,obj,op);
                end

        %        calculateZpw(obj);                  % 2a) CALCULATE Z: Zpw (plane wave impedance of tubes)
        %        calculateThevSource(obj,obj.Zpw);   % 3a) CALCULATE PS & ZS (source pressure and impedance)
        %        calculatezew(obj);                  % 2b) CALCULATE Z: zew (evanescent wave impedance of tubes)
        %        calculatezfl(obj);                  % 2c) CALCULATE Z: zfl (flow loss impedance of tubes)
        %        calculateZLi(obj);                  % 2d) CALCULATE Z: ZLi (total ideal load impedance = Zpw + zew + zfl)

                                                    %     then further refine cavity lenghts using a minimization routine

                op=optimset('Display','off','MaxIter',9000);
                obj.optSwitch = 'T';
                [obj,fval,exitFlag,output] = fminsearchOBJ(@findOptimalLengths,obj,op);


        %       calculateThevSource(obj,obj.ZLi);   % 3b) CALCULATE PS & ZS (source pressure and impedance)
                calculateZL(obj);                   % 4)  CALCULATE ZL (actual cavity impedance of the tubes, given PS & ZS)
                removeCrosstalk(obj)
                calculateError(obj);                % 5)  CALCULATE ERROR 
        %        plotResults(obj);                   % 6)  PLOT RESULTS

                if kk == 1
                    PL1 = obj.PL;
                    PS1 = obj.PS;
                    ZL1 = obj.ZL;
                    ZS1 = obj.ZS;
                    ZLi1 = obj.ZLi;
                    Freq1 = obj.freq;
                elseif kk == 2
                    PL2 = obj.PL;
                    PS2 = obj.PS;
                    ZL2 = obj.ZL;
                    ZS2 = obj.ZS;
                    ZLi2 = obj.ZLi;
                    Freq2 = obj.freq;
                elseif kk == 3
                    PL3 = obj.PL;
                    PS3 = obj.PS;
                    ZL3 = obj.ZL;
                    ZS3 = obj.ZS;
                    ZLi3 = obj.ZLi;
                    Freq3 = obj.freq;
                elseif kk == 4
                    PL4 = obj.PL;
                    PS4 = obj.PS;
                    ZL4 = obj.ZL;
                    ZS4 = obj.ZS;
                    ZLi4 = obj.ZLi;
                    Freq4 = obj.freq;
                elseif kk == 5
                    PL5 = obj.PL;
                    PS5 = obj.PS;
                    ZL5 = obj.ZL;
                    ZS5 = obj.ZS;
                    ZLi5 = obj.ZLi;
                    Freq5 = obj.freq;
                elseif kk == 6
                    PL6 = obj.PL;
                    PS6 = obj.PS;
                    ZL6 = obj.ZL;
                    ZS6 = obj.ZS;
                    ZLi6 = obj.ZLi;
                    Freq6 = obj.freq;
                else
                end
            end
            figure
            subplot(2,1,1)
            plot(Freq1,abs(PL1),'r')
            hold on
            plot(Freq2,abs(PL2),'g')
            plot(Freq3,abs(PL3),'b')
            plot(Freq4,abs(PL4),'m')
            plot(Freq5,abs(PL5),'c')
            plot(Freq6,abs(PL6),'k')
            ylabel('PL')
            subplot(2,1,2)
            plot(Freq1,abs(PS1),'r')
            hold on
            plot(Freq2,abs(PS2),'g')
            plot(Freq3,abs(PS3),'b')
            plot(Freq4,abs(PS4),'m')
            plot(Freq5,abs(PS5),'c')
            plot(Freq6,abs(PS6),'k')
            ylabel('PS')

            figure
            subplot(2,1,1)
            plot(Freq1,abs(ZL1),'r')
            hold on
            plot(Freq2,abs(ZL2),'g')
            plot(Freq3,abs(ZL3),'b')
            plot(Freq4,abs(ZL4),'m')
            plot(Freq5,abs(ZL5),'c')
            plot(Freq6,abs(ZL6),'k')
            ylabel('ZL')
            subplot(2,1,2)
            plot(Freq1,abs(ZS1),'r')
            hold on
            plot(Freq2,abs(ZS2),'g')
            plot(Freq3,abs(ZS3),'b')
            plot(Freq4,abs(ZS4),'m')
            plot(Freq5,abs(ZS5),'c')
            plot(Freq6,abs(ZS6),'k')
            ylabel('ZS')

            figure
            plot(Freq1,abs(ZLi1),'r')
            hold on
            plot(Freq2,abs(ZLi2),'g')
            plot(Freq3,abs(ZLi3),'b')
            plot(Freq4,abs(ZLi4),'m')
            plot(Freq5,abs(ZLi5),'c')
            plot(Freq6,abs(ZLi6),'k')
            ylabel('ZLi')

            [Freq,ZS,PS,ZL,PL,ZLi] = joiner(Freq1,ZS1,PS1,ZL1,PL1,ZLi1,Freq2,ZS2,PS2,ZL2,PL2,ZLi2);
            [Freq,ZS,PS,ZL,PL,ZLi] = joiner(Freq,ZS,PS,ZL,PL,ZLi,Freq3,ZS3,PS3,ZL3,PL3,ZLi3);
            [Freq,ZS,PS,ZL,PL,ZLi] = joiner(Freq,ZS,PS,ZL,PL,ZLi,Freq4,ZS4,PS4,ZL4,PL4,ZLi4);
            [Freq,ZS,PS,ZL,PL,ZLi] = joiner(Freq,ZS,PS,ZL,PL,ZLi,Freq5,ZS5,PS5,ZL5,PL5,ZLi5);
            [Freq,ZS,PS,ZL,PL,ZLi] = joiner(Freq,ZS,PS,ZL,PL,ZLi,Freq6,ZS6,PS6,ZL6,PL6,ZLi6);

            % add smoothing here
            ZS = smoother(Freq,ZS);
            PS = smoother(Freq,PS);
            ZL = smoother(Freq,ZL);
            PL = smoother(Freq,PL);
            ZLi = smoother(Freq,ZLi);

            obj.ZS = ZS;
            obj.PS = PS;
            obj.ZL = ZL;
            obj.PL = PL;
            obj.ZLi = ZLi;
            obj.freq = Freq;
            obj.nFreqs = length(obj.freq);
            obj.fmin = starts(1);
            obj.fmax = finishes(end);
            calculateError(obj);                % 5)  CALCULATE ERROR 
            plotResults(obj);                   % 6)  PLOT RESULTS
        end
    end
    
    function [] = calculatePL(obj)  % 1) CALCULATE PL (load pressure; pressure measured in the tubes)
        %obj.nCavities = length(obj.cavityLengths_nominal); % number of tubes
        [obj.nfft,obj.nCavities] = size(obj.recordings);
        obj.cavityLengths_nominal = obj.cavityLengths_nominal(1:obj.nCavities);
        %obj.nfft = size(obj.recordings,1); % fft size
        obj.freq = (0:1:(obj.nfft-1))'*(obj.fs/obj.nfft); % frequency in Hertz associated with the fft
        [~,obj.fminIndx] = min(abs(obj.freq - obj.fmin)); % index of minimum frequency
        [~,obj.fmaxIndx] = min(abs(obj.freq - obj.fmax)); % index of maximum frequency
        S = repmat(fft(obj.stimulus),1,obj.nCavities); % fft of stimulus vector, replicated to a matrix
        R = fft(obj.recordings); % fft of recordings matrix
        obj.PL = R ./ S; % divide recording by stimulus to account for non-uniform magnitude and phase
        obj.PL(1:obj.fminIndx-1,:) = 0 + eps; % zero out frequencies outside of fmin and fmax
        obj.PL(obj.fmaxIndx+1:end,:) = 0 + eps; % (add eps to avoid division by zero)
        nyquistIndx = obj.nfft/2; % nyquist frequency
        obj.PL = [obj.PL(1:nyquistIndx,:);flipud(conj(obj.PL(2:nyquistIndx-1,:)))]; % add aliased part back in
        obj.PL_timeDomain = real(ifft(abs(obj.PL).^2)); % convert to time domain; used later in computing optimal cavity lengths
        obj.PL = obj.PL(obj.fminIndx:obj.fmaxIndx,:); % get rid of frequencies outside of fmin and fmax; no longer needed
        obj.freq = obj.freq(obj.fminIndx:obj.fmaxIndx);
    end
    function [] = calculateZpw(obj) % 2a) CALCULATE Zpw (plane wave impedance of tubes)
        rho = 1.1769e-3 * (1 - 0.00335 * obj.d); % density of air
        r = obj.cavityDiameter / 2; % radius of the cavity (cm)
        obj.z0 = (rho * obj.c) / (pi * r^2);  % characteristic (surge) impedence (ohms)
        
        obj.omega = 2 * pi * obj.freq; % radian (anglar) frequency (rad/s)
        obj.omega(obj.omega<eps) = eps; % minumum value is eps, to avoid division by zero
        
        % These are the Keefe equations:
        eta = 1.846e-4 * (1 + 0.0025 * obj.d);
        Rv = r*sqrt((rho*obj.omega)/eta);
        x = (1.045 + (1.080 + 0.750 ./ Rv) ./ Rv) ./ Rv;
        y = 1 + 1.045 ./ Rv;
        obj.wn = (obj.omega/obj.c) .* complex(x,y); % wave number
        R = exp(-2 * obj.wn * obj.cavityLengths_optimal(:)'); % note: force cavity lengths to be a row vector
        % Zpw = z0 .* (1 + R) ./ (1 - R);
        % NOTE: It is possible to reduce the notch depths calculating ZLi the
        % following way, where reflNum and reflDenom are reflection
        % coefficients >0 and <=1.  Reducing R in the numerator reduces the notch depth; 
        % reducing R in the denominator reduces the peak height. This version
        % allows for separate control over notch and peak sharpness.
        % -----
        %reflNum = 1; % reduce notch depth
        %reflDenom = 1; % reduce peak height
        %obj.Zpw = obj.z0 .* (1 + (R.*reflNum)) ./ (1 - (R.*reflDenom));        
        
        
        
        
        nMax = 1;
        nMin = 1; %0.98; % 0.88; nMin = 0.875; % .85 gives .42 error
        rn = linspace(nMax,nMin,size(R,1))';
        Rn = repmat(rn,1,size(R,2));
        dMax = 1;
        dMin = 1;
        dn = linspace(dMax,dMin,size(R,1))';
        Dn = repmat(dn,1,size(R,2));
        obj.Zpw = obj.z0 .* (1 + (R.*Rn)) ./ (1 - (R.*Dn));     

% added 10/13/2019 ssg        
NMin = (0.8:0.01:1.0)';
NN = length(NMin);
for ii=1:NN
    nMax = 1;
    nMin = NMin(ii);
    rn = linspace(nMax,nMin,size(R,1))';
    Rn = repmat(rn,1,size(R,2));
    dMax = 1;
    dMin = 1;
    dn = linspace(dMax,dMin,size(R,1))';
    Dn = repmat(dn,1,size(R,2));
    obj.Zpw = obj.z0 .* (1 + (R.*Rn)) ./ (1 - (R.*Dn));     

    obj.calculateZLi;
    obj.calculateZL;
    deltaZ = obj.ZL - obj.ZLi;
    N = size(deltaZ(:),1);
    eTotal(ii,1) = sum(sum(abs(deltaZ)))/N;
end
[~,indx] = min(eTotal);
nMin = NMin(indx);
obj.nMin = nMin;
rn = linspace(nMax,nMin,size(R,1))';
Rn = repmat(rn,1,size(R,2));
dMax = 1;
dMin = 1;
dn = linspace(dMax,dMin,size(R,1))';
Dn = repmat(dn,1,size(R,2));
obj.Zpw = obj.z0 .* (1 + (R.*Rn)) ./ (1 - (R.*Dn));     



        
    end
    function [] = calculateThevSource(obj,Z) % 3) CALCULATE PS & ZS (source pressure and impedance)
        % 3a) Initial calculation using Zpw only
        % 3b) FINAL CALCULATION using idealized total impedanc (ZLi)
        obj.nFreqs = size(obj.PL,1); % number of frequencies in analysis range
        obj.ZS = zeros(obj.nFreqs,1); % initialize source impedance vector
        obj.PS = zeros(obj.nFreqs,1); % initialize source pressure vector
        for ii=1:obj.nFreqs % loop across frequencies
            z = Z(ii,:).';
            p = obj.PL(ii,:).';
            A = [z -p];
            B = z .* p;
            %%%%%            
            %w = sqrt(abs(p)/sum(abs(p)));
            %A = A .* [w,w];
            %B = B .* w;
            %%%%%            
            x = A \ B;
            obj.PS(ii) = x(1);
            obj.ZS(ii) = x(2);
        end
    end
    function [] = calculatezew(obj) % 2b) CALCULATE zew (evanescent wave impedance of tubes)  
        % Estimate wave-spread inertance from measured pressure
        %   pcd, paralell compliance delay, is needed to calculate tau.
        %   This requires an initial estimate of ZS and Z0
        %   parallel-compliance delay calculation:
        izs = unwrap(imag(1./obj.ZS)); % imaginary part of ZS (estimates parallel compliance)
        % calculate as a straight-line fit across a specified range (fmin and fmax)
        polyOrder = 1; % polynomial order
        coeff = polyfit(obj.omega,izs,polyOrder); % fit data with a straight line      
        pcd = obj.z0 * coeff(1); % slope of fit times surge (characteristic) impedance       
        % Tau, wavespread estimate, is needed to calcualte impedance of evenescent wave component
        %   tau = zeros(1,nCavities); % wavespread estimate (from parallel compliance delay)
        obj.tau = -pcd / 2 / 1000; % pcd = parallel compliance delay (? change from round trip to one-way delay and put in ms from s?)
        obj.tau = repmat(obj.tau,1,obj.nCavities); % row vector
        obj.tau(1) = obj.tau(1) * 2; % ? not sure why
        obj.tau = obj.tau(:); % force to column vector
       
        phi = 0; 
        for ii=1:obj.nCavities
            phi = phi + 1i * (obj.tau(ii) * obj.omega).^ii;
        end
        obj.zew = phi * obj.z0;
        
        %figure(2001); hold on; plot(imag(obj.zew)); pause(.05);
    end
    function [] = calculatezfl(obj) % 2c) CALCULATE zfl (flow loss impedance of tubes)
        rfl = 0.02; % resistive flow loss ?? friction factor for flow loss in an air-filled tube; Darcy friction factor
        fnl = 30e6; % ??
        obj.zfl = -rfl * sqrt(obj.freq) .* (1 - (obj.freq / fnl).^4);
    end
    function [] = calculateZLi(obj) % 2d) CALCULATE ZLi (total ideal load impedance of tubes)
        Zew = repmat(obj.zew,1,obj.nCavities);
        Zfl = repmat(obj.zfl,1,obj.nCavities);
        obj.ZLi = obj.Zpw + Zew + Zfl;
    end
    function [] = calculateZL(obj)  % 4) CALCULATE ZL (actual cavity impedance of the tubes, given PS & ZS)
        obj.ZL = zeros(size(obj.PL)); % initialize
        for ii=1:obj.nCavities  % calculate load impedance for each cavity
            try
            obj.ZL(:,ii) = obj.ZS .* obj.PL(:,ii) ./ (obj.PS - obj.PL(:,ii));
            catch ME
                keyboard
            end
        end
    end
    function [] = calculateError(obj) % CALCULATE ERROR (compare calculated and measured pressure)
        s1 = 0; s2 = 0; % initialize sums of squares variable
        accumulatedPressure = zeros(obj.nFreqs,1); % calculating pressure error (Neely method)
            for ii=1:obj.nCavities % loop across cavities
               pressureDiff = obj.PL(:,ii) - obj.PS .* obj.ZLi(:,ii) ./ (obj.ZS + obj.ZLi(:,ii)); % Numerator |Pc - Ps * Zc / (Zs + Zc)|
               accumulatedPressure = accumulatedPressure + pressureDiff; % sum of pressure differences across cavities; only used if crosstalk is being considered
               s1 = s1 + sum(abs(pressureDiff).^2); % sum of squares:  square the numerator and sum across frequency;
               s2 = s2 + sum(abs(obj.PL(:,ii)).^2); % sum of squares: square the measured cavity pressure and sum across frequency
            end
        s1 = s1 - sum(abs(accumulatedPressure).^2) / obj.nCavities; % line not used if crosstalk isn't being considered
        obj.epsilon = 10000 * s1 / s2; % divide the sum of squares of the pressure difference with the sum of squares of the measured pressure and multiply by 10,000
    end
    
    function [] = optimizeCavityLengths(obj)
        %obj.cavityLengths_optimal = obj.cavityLengths_nominal;
        %return
        
        % estimate cavity length from autocorrelation peak
        %function [lk,td]=cavity_length(pc,ss,sr,lmn,lmx)
        
        
        %minLength = obj.cavityDiameter; % minimum allowable cavity length (cm)
        %if min(obj.cavityLengths_nominal) < minLength
        %    error('Minimum cavity length must be >= cavity diameter.')
        %end
        minLength = min(obj.cavityLengths_nominal) * 0.8; % minimum cavity cavity length (cm)
        maxLength = max(obj.cavityLengths_nominal) * 1.2; % maximum cavity cavity length (cm)
        % constrain where to look for autocorrelation peaks based on cavity length min and max
        obj.d = obj.cavityTemperature - 26.85;
        obj.c = 3.4723e4 * (1 + 0.00166 * obj.d); % speed of sound at the cavity temperature (cm/s)
        minIndx = round(minLength * obj.fs *2 ./ obj.c)+1;  % minimum location (samples) where delay can be located
        maxIndx = round(maxLength * obj.fs *2 ./ obj.c)+1;  % maximum location (samples) where delay can be located
        obj.cavityLengths_optimal = zeros(obj.nCavities,1); % initialize output
        for ii=1:obj.nCavities
            %p = obj.ffs(abs(obj.PL_timeDomain(:,ii)).^2); % convert measured pressure to time domain
            %p = obj.PL_timeDomain(:,ii);
            p = obj.PL_timeDomain(minIndx:maxIndx,ii);
            [~,mHat] = max(p);
            spread = 1;
            if mHat+spread > length(p)
                mHat = 1;
            end
            
            if mHat > 1
                polyOrder = 2; % polynomial order 
                x = (1:1:length(p))';
                coeff = polyfit(x(mHat-spread:mHat+spread),p(mHat-spread:mHat+spread),polyOrder); % fit data with a straight line      
                dcoeff = polyder(coeff);
                mHat = roots(dcoeff);
                mHat = mHat + minIndx - 1;
            else
                mHat = minIndx;
            end
            timeDelay = (mHat-1) / obj.fs; % round trip delay (s)
            %[~,mHat] = max(p(minIndx:maxIndx)); % locate the maximum peak within the specified range
            %mHat = mHat + minIndx - 1; % adjust min indx to the full vector; account for time 0 at sample 1
            %if (mHat > minIndx) % if the maximum is not at or below the lower edge
            %   delay =(p(mHat-1)-p(mHat+1))/(p(mHat-1)-2*p(mHat)+p(mHat+1))/2; % calculate sample delay
            %else
            %   delay = 0; % otherwise, set delay to zero
            %end
            %timeDelay = (mHat+delay-1) / obj.fs; % time delay (seconds?)
            obj.cavityLengths_optimal(ii,1) = timeDelay * obj.c / 2; % length estimate; 2 is factoring round-trip delay
        end
        
    end
    function [] = plotResults(obj)
        figure(413) % PS: Calculated Source Pressure
        subplot(2,1,1)
            plot(obj.freq/1000,20*log10(abs(obj.PS)))
            hold on
            title('PS')
            xlabel('Frequency (kHz)','FontSize',12)
            ylabel('Magnitude in dB','FontSize',12)
            xlim([obj.freq(1)/1000,obj.freq(end)/1000])
        subplot(2,1,2);
            plot(obj.freq/1000,angle(obj.PS)/(2*pi))
            hold on
            xlabel('Frequency (kHz)','FontSize',12)
            ylabel('Phase (cycles)')
            xlim([obj.freq(1)/1000,obj.freq(end)/1000])

        figure(414) % ZS: Calculated Source Impedance
        subplot(2,1,1);
            plot(obj.freq/1000,20*log10(abs(obj.ZS)))
            hold on
            title('ZS')
            xlabel('Frequency (kHz)','FontSize',12)
            ylabel('Magnitude (dB)','FontSize',12)
            xlim([obj.freq(1)/1000,obj.freq(end)/1000])
        subplot(2,1,2);
            plot(obj.freq/1000,angle(obj.ZS)/(2*pi))
            hold on
            xlabel('Frequency (kHz)','FontSize',12)
            ylabel('Phase (cycles)','FontSize',12)
            xlim([obj.freq(1)/1000,obj.freq(end)/1000])

        figure(415) % Error Estimate: ZLi vs ZL (ideal versus measured load impedance)
        subplot(2,1,1);
            plot(obj.freq/1000,20*log10(abs(obj.ZLi)),':')
            hold on
            plot(obj.freq/1000,20*log10(abs(obj.ZL)),'-')
            title(['Error Estimate: ',num2str(obj.epsilon)])
            xlabel('Frequency (kHz)','FontSize',12)
            ylabel('Magnitude (dB)','FontSize',12)
            xlim([obj.freq(1)/1000,obj.freq(end)/1000])
        subplot(2,1,2);
            plot(obj.freq/1000,angle(obj.ZLi)/(2*pi),':')
            hold on
            plot(obj.freq/1000,angle(obj.ZL)/(2*pi),'-')
            xlabel('Frequency (kHz)','FontSize',12)
            ylabel('Phase (cycles)','FontSize',12)
            xlim([obj.freq(1)/1000,obj.freq(end)/1000])
    end
    function [] = removeCrosstalk(obj)
        [nf,nc] = size(obj.PL);         % number of cavities
        px = zeros(nf,1);
        zs = repmat(obj.ZS,1,nc);
        hh = obj.ZLi ./ (zs + obj.ZLi);
        oo = ones(nc,1);
        for k=1:nf
           hk=transpose(hh(k,:));
           pk=transpose(obj.PL(k,:));
           pp=[oo hk]\pk;
           px(k)=pp(1);
           ps(k)=pp(2);
        end
        ps = ps(:);
        obj.PS = ps - px;
    end
end
end

% External Functions (not formally part of class definition) --------------

function [epsilon] = findOptimalLengths(obj)
    % Minimize the error based on ideal versus calculated PRESSURE (Neely method)
    calculateZpw(obj);  % Always calculate plane-wave impedance based on current optimal lengths
    calculatezew(obj)
    calculatezfl(obj)
    calculateZLi(obj)  % compute total impedance (Zpw + zew + zfl)
%     if strcmp(obj.optSwitch,'L')
%         calculateThevSource(obj,obj.Zpw);
%     elseif strcmp(obj.optSwitch,'T')
%         calculateThevSource(obj,obj.ZLi);
%     end
    calculateThevSource(obj,obj.ZLi);
    calculateError(obj);
    epsilon = obj.epsilon;
end
function [] = calculateZpw_RW(obj) % cavity impedance (theoretical)
    L = obj.cavityLengths_optimal/100;
    radius = 0.5*obj.cavityDiameter/100;
    w = 2*pi*obj.freq;
    rho = 1.164;
    vis = 1.886*10^(-5); % http://www.lmnoeng.com/Flow/GasViscosity.htm @ 30 degrees Celcius
    C_p = 1.005; %specific heat of air at constant pressure , http://www.engineeringtoolbox.com/air-properties-d_156.html
    K_t = 0.0264*10^(-3); %thermal conductivity , http://www.engineeringtoolbox.com/air-properties-d_156.html
    varphi_v = (w.*rho./vis).^(0.5).*radius; %Benade(1968) eqn 4
    varphi_t = ((C_p.*w.*rho./K_t).^(0.5)).*radius; %Benade(1968) eqn 5
    [R,G,M,C] = calculateRGMC(radius,varphi_v,varphi_t,obj);
    
    
    rho = 1.164;
    r = obj.cavityDiameter/200; % radius in m
    z0 = (rho * obj.c) / (pi * r^2); % z0 = Ro; characteristic impedance (z0) is constant & real;     
    
    obj.omega = 2 * pi * obj.freq; % radian (anglar) frequency (rad/s)
    obj.omega(obj.omega<eps) = eps; % minumum value is eps, to avoid division by zero
    
    eta = 1.846e-4 * (1 + 0.0025 * obj.d);
    Rv = r*sqrt((rho*obj.omega)/eta);
    x = (1.045 + (1.080 + 0.750 ./ Rv) ./ Rv) ./ Rv;
    y = 1 + 1.045 ./ Rv;
    c = 331.3 + 0.606*obj.cavityTemperature;
    obj.wn = (obj.omega/c) .* complex(x,y); % wave number    
    
    
    Z0 = (R + 1i.*w.*M)./obj.wn; 
    YL1 = ((1i.*tan(-1i*obj.wn.*L(1)))./Z0);
     YL2 = ((1i.*tan(-1i*obj.wn.*L(2)))./Z0);
      YL3 = ((1i.*tan(-1i*obj.wn.*L(3)))./Z0);
       YL4 = ((1i.*tan(-1i*obj.wn.*L(4)))./Z0);
        YL5 = ((1i.*tan(-1i*obj.wn.*L(5)))./Z0);
    YLi = [YL1 YL2 YL3 YL4 YL5];   
    obj.Zpw = 1./YLi;
end
function [R,G,M,C] = calculateRGMC(radius,varphi_v,varphi_t,obj)  %Rigid Walls 
    %R,G,M,C exact solutions
    w = 2*pi*obj.freq; 
    rho = 1.164;
    VisArray = varphi_v.*(-1i).^0.5;
    BesVisNumArray = besselj(1,VisArray);
    Var_Vis_num = 2*BesVisNumArray;
    BesVisDenArray = besselj(0,VisArray);
    Var_Vis_denom = (varphi_v.*(-1i).^0.5).*BesVisDenArray;
    Var_Vis = Var_Vis_num./Var_Vis_denom;
    F_v = abs(Var_Vis);
    Phi_v = angle(Var_Vis);
    DD = ((1-F_v.*cos(Phi_v)).^2)+((F_v.*sin(Phi_v)).^2);
    ThArray = varphi_t.*(-1i).^0.5;
    BesThNumArray = besselj(1,ThArray);
    Var_Th_num = 2*BesThNumArray;
    BesThDenArray = besselj(0,ThArray);
    Var_Th_denom = (varphi_t.*(-1i).^0.5).*BesThDenArray;
    Var_Th = Var_Th_num./Var_Th_denom;
    F_t = abs(Var_Th);
    Phi_t = angle(Var_Th);
    R = - ((w.*rho)./(pi.*radius.^2)).*((F_v.*sin(Phi_v))./(DD));
    M =  ((rho)./(pi.*radius.^2)).*(1-F_v.*cos(Phi_v))./(DD);
    G = - ((w.*pi*radius*radius)./(rho*obj.c*obj.c)).*(0.4*F_t.*sin(Phi_t));
    C = ((pi*radius*radius)./(rho*obj.c*obj.c)).*(1+0.4*F_t.*cos(Phi_t));
end

function [Xsm] = smoother(Freq,X)
    w = ones(size(Freq)); % create vector of weights for the smoothing spline
    N = size(X,2);
    smoothing = 0.00001;
    for ii=1:N
        % values are complex; smooth on mag and phase individually
        mag = abs(X(:,ii));
        ang = angle(X(:,ii));
        mag = csaps(Freq,mag,smoothing,Freq,w); % the smoothed, weighted spline
        ang = csaps(Freq,unwrap(ang),smoothing,Freq,w); % the smoothed, weighted spline
        ang = angle(cos(ang) + 1i*sin(ang)); % wrap the phase back up
        Xsm(:,ii) = mag.*cos(ang) + 1i*mag.*sin(ang);
    end    
end
function [Freq,ZS,PS,ZL,PL,ZLi] = joiner(Freq1,ZS1,PS1,ZL1,PL1,ZLi1,Freq2,ZS2,PS2,ZL2,PL2,ZLi2)
    % join everything back together
    C = intersect(Freq1,Freq2);
    indx1 = find(Freq1==C(1));
    indx2 = find(Freq2==C(end));
    %h1 = hann(length(Freq1)-indx1+1);
    h =  hann(indx2*2);
    h2 = h(1:indx2);
    h1 = flipud(h2);
    ZScommon = ((ZS1(indx1:end).*h1) + (ZS2(1:indx2).*h2)) ./ (h1 + h2);
    ZSjoint = [ZS1(1:indx1-1);ZScommon;ZS2(indx2+1:end)];
    FreqJoint = [Freq1(1:indx1-1);C;Freq2(indx2+1:end)];
    PScommon = ((PS1(indx1:end).*h1) + (PS2(1:indx2).*h2)) ./ (h1 + h2);
    PSjoint = [PS1(1:indx1-1);PScommon;PS2(indx2+1:end)];
    for jj=1:size(PL1,2)
        ZLcommon(:,jj) = ((ZL1(indx1:end,jj).*h1) + (ZL2(1:indx2,jj).*h2)) ./ (h1 + h2);
        ZLjoint(:,jj) = [ZL1(1:indx1-1,jj);ZLcommon(:,jj);ZL2(indx2+1:end,jj)];
        PLcommon(:,jj) = ((PL1(indx1:end,jj).*h1) + (PL2(1:indx2,jj).*h2)) ./ (h1 + h2);
        PLjoint(:,jj) = [PL1(1:indx1-1,jj);PLcommon(:,jj);PL2(indx2+1:end,jj)];
    end
    for jj=1:size(ZLi1,2)
        ZLicommon(:,jj) = ((ZLi1(indx1:end,jj).*h1) + (ZLi2(1:indx2,jj).*h2)) ./ (h1 + h2);
        ZLijoint(:,jj) = [ZLi1(1:indx1-1,jj);ZLicommon(:,jj);ZLi2(indx2+1:end,jj)];
    end
    Freq = FreqJoint;
    ZS = ZSjoint;
    PS = PSjoint;
    ZL = ZLjoint;
    PL = PLjoint;
    ZLi = ZLijoint;
end


% OLD CODE -----------------------------------------------------
%         obj.tau = 0; % first time through, set to zero
%         q = obj.cavityLengths_optimal;
%         q = fminsearch(@(q)pc_err(obj.cavityLengths_optimal,obj.freq,obj.PL,obj.tau),obj.cavityLengths_optimal,op);
%         [q,fval,exitFlag,output] = fminsearch(@pc_err,obj,op);

%     function err = pc_err(pa,f,pc,tau)
%         %[lc,rc] = get_lc(pa);
%         if (min(lc)<0.1) err=1e9; return; end
%         if (min(tau)<0)  err=1e9; return; end
%         %
%         nc = size(pc,2);           % number of cavities
%         zc = cavimp(f,lc,rc,tau);  % calculate cavity impedances
%         [zs,ps] = thvsrc(zc,pc);   % estimate zs & ps
%         ps = repmat(ps,1,nc);
%         zs = repmat(zs,1,nc);
%         pl = ps .* zc ./ (zs + zc);
%         pd = pc - pl;
%         s1 = sum(sum(abs(pd).^2));
%         s2 = sum(sum(abs(pc).^2));
%         err = (s1 / s2) * 1e4;
%         %err = sum(sum(abs(1 - pl ./ pc).^2)) / 100;
%     end    
%     function h=ffs(obj,H)
%         m=length(H);
%         n=2*(m-1);
%         H(1,:)=real(H(1,:));
%         H(m,:)=real(H(m,:));
%         H((m+1):n,:)=conj(H((m-1):-1:2,:));
%         h=real(ifft(H));
%     end    
      

        %exitFlag >0 Indicates that the function converged to a solution x.
        %exitFlag =0 Indicates that the maximum number of function evaluations was exceeded.
        %exitFlag <0 Indicates that the function did not converge to a solution.
        %output.algorithm The algorithm used
        %output.funcCount The number of function evaluations
        %output.iterations The number of iterations taken        




