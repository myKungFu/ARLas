classdef thevCal_new < handle
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Class definition thevCal 
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
% Updated: October 12, 2019 -- ssg
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

properties (SetAccess = private)
    fs
    cavityTemperature
    cavityDiameter
    cavityLengths_nominal
    stimulus
    recordings
    recordingsOrig
    micCorrection
    
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
end

methods
    function t = thevCal_new(recordingParams,stimulus,recordings)
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
        disp('Calculating Thevenin Source, please wait...')
        
        if obj.applyMicCorrection == 1 % if apply microphone correction requested
            obj.recordings = fastFilter(obj.micCorrection,obj.recordingsOrig); % apply microphone correction
        else
            obj.recordings = obj.recordingsOrig; % no microphone correction
        end

        calculatePL(obj);                   % 1)  CALCULATE PL (load pressure; pressure measured in the tubes)
        optimizeCavityLengths(obj);         %     Make a first estimate cavity lengths based on autocorrelation
                                            %     Further refine cavity lenghts using a minimization routine
                                            
        obj.zew = zeros(obj.fmaxIndx - obj.fminIndx+1,1);  % Initialize values
        obj.zfl = zeros(obj.fmaxIndx - obj.fminIndx+1,1);
        obj.ZS = zeros(obj.fmaxIndx - obj.fminIndx+1,1);
        obj.ZLi = zeros(obj.fmaxIndx - obj.fminIndx+1,obj.nCavities);
        obj.tau = zeros(1,obj.nCavities);
                                            
        op=optimset('Display','off','MaxIter',9000);
        obj.optSwitch = 'L';
        [obj,fval,exitFlag,output] = fminsearchOBJ(@findOptimalLengths,obj,op);
        
        
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
        %removeCrosstalk(obj)
        calculateError(obj);                % 5)  CALCULATE ERROR 
        plotResults(obj);                   % 6)  PLOT RESULTS

        % check power reflectance
% SPL = obj.PL(:,1); %fft(X,nfft); % measured sound presure level
% ZLmatrix = obj.ZL(:,1);
% FPL = (SPL/2) .* (1 + (obj.z0 ./ ZLmatrix)); % complex FPL
% RPL = (SPL/2) .* (1 - (obj.z0 ./ ZLmatrix)); % complex RPL
% PRR = RPL ./ FPL;                 % complex pressure reflectance PR
% IPL = abs(FPL) + abs(RPL);        % integrated pressure level
% PR = abs(RPL).^2 ./ abs(FPL).^2;  % power reflectance
%         

        
        
        
    end
    
    function [] = calculatePL(obj)  % 1) CALCULATE PL (load pressure; pressure measured in the tubes)
        obj.nCavities = length(obj.cavityLengths_nominal); % number of tubes
        obj.nfft = size(obj.recordings,1); % fft size
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
        nMin = 0.88;%nMin = 0.875; % .85 gives .42 error
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
            obj.ZL(:,ii) = obj.ZS .* obj.PL(:,ii) ./ (obj.PS - obj.PL(:,ii));
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
        figure % PS: Calculated Source Pressure
        subplot(2,1,1)
            plot(obj.freq/1000,20*log10(abs(obj.PS)))
            title('PS')
            xlabel('Frequency (kHz)','FontSize',12)
            ylabel('Magnitude in dB','FontSize',12)
            xlim([obj.freq(1)/1000,obj.freq(end)/1000])
        subplot(2,1,2);
            plot(obj.freq/1000,angle(obj.PS)/(2*pi))
            xlabel('Frequency (kHz)','FontSize',12)
            ylabel('Phase (cycles)')
            xlim([obj.freq(1)/1000,obj.freq(end)/1000])

        figure % ZS: Calculated Source Impedance
        subplot(2,1,1);
            plot(obj.freq/1000,20*log10(abs(obj.ZS)))
            title('ZS')
            xlabel('Frequency (kHz)','FontSize',12)
            ylabel('Magnitude (dB)','FontSize',12)
            xlim([obj.freq(1)/1000,obj.freq(end)/1000])
        subplot(2,1,2);
            plot(obj.freq/1000,angle(obj.ZS)/(2*pi))
            xlabel('Frequency (kHz)','FontSize',12)
            ylabel('Phase (cycles)','FontSize',12)
            xlim([obj.freq(1)/1000,obj.freq(end)/1000])

        figure % Error Estimate: ZLi vs ZL (ideal versus measured load impedance)
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




