classdef thevCal < handle
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Class definition thevCal 
% For use with ARLas (Auditory Research Laboratory auditory software)
% Used in performing Thevenin source calibration.
% This class replaces the an older version, theveninSource.m
% Required files: findOptimalLengths.m.
%
% USAGE.  There are two ways to create the object:
% 1) Include required data when object is created
%       >> t = thevCal(recordingParams,stimulus,recordings);
%
% 2) Create the object, then manually load the required data
%       >> t = thevCcal;
%       >> t.stimulus = xxxx; % electrical stimulus is a column vector
%       >> t.recordings = xxxx; % matrix of acoustic recordings made in response to stimulus
%       >> t.fileName = xxxx; % name of file containing the recordings
%       >> t.fs = xxxx; % sampling rate in Hz
%       >> t.cavityTemperature = xxxx; % temperature in degrees C
%       >> t.cavityDiameter = xxxx; % diameter in cm
%       >> t.cavityLengths_nominal = xxxx; % intended cavity lengths (cm)
%       >> t.fmin = xxxx; % in Hz; min frequency over which to calculate (stimulus must include this)
%       >> t.fmax = xxxx; % in Hz; max frequency over which to calculate (stimulus must include this)
%       >> t.timeStamp = xxxx; % date and time when recordings were made
%       >> t.micSN = xxxx; % serial number of the recording microphone
%
% Once the data have been loaded into the object, then perform the Thevenin
% source calibration:
%       >> t.calculate
% You can manually make plots from objects of this class:
%       >> t.plotZ will plot the mag & phase of the measured vs ideal impedance
%               of the calibration cavities. The associated error values are 
%               calculated from this difference and are shown for magnitude only (eM),
%               phase only(eP), and total including magnitude and phase (eT).
%       >> t.plotP will plot the mag & phase of the measured vs ideal pressure
%               of the calibration cavities. The associated error values are 
%               calculated from this difference and are shown for the total 
%               (both mag and phase) only.
%       >> t.plotD  will plot the impedance error magnitude and phase,
%               where the error is the difference between Z measured and Z ideal
%               in the calibration cavities.
%       >> t.plotL will plot a comparison of cavity lengths: nominal (what
%              was expected by the program), estimated (initial best guess
%              based on half-wavelength resonances), and optimal (what the
%              optimization algorithm found to be the best set of lengths.
%              You can use this to look for systematic errors in your
%              initial vs estimated values and to see how the optimization
%              performed.
%       >> t.plotZS will plot the calculated source impedance (magnitude
%              and phase).
%       >> t.plotPS will plot the calculated source pressure (magnitude
%              and phase).
%
% Note that plots that with magnitude that can be expressed in dB or linear
% units are controlled by this command:
%       >> t.plot_dB = 1; (plot magnitude in dB, the default)
%       >> t.plot_dB = 0; (plot magnitude in linear units)
%
% If you are satisfied and want to keep a calibration save it using:
%       >> t.saveThev
%
% Auditory Research Lab, The University of Iowa
% Deptartment of Communication Sciences & Disorders
% The University of Iowa
% Author: Shawn Goodman
% Based originally on code generously provided by Stephen T. Neely
% Date: July 4, 2015
% Updated: September 2, 2015 - ssg
% Updated: October 31, 2016 - ssg
% Updated: June 13, 2017 - ssg for use with er10x; see line 182
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

properties (SetAccess = private)
    PS % source pressure
    ZS % source impedance
    PL % measured load pressure in cavities (expressed relative to stimulus)
    ZL % load impedance
    PLi % ideal load pressure (based on cavity dimensions and temperature)
    ZLi % ideal load impedance
    epsilon % total calculated error (used as a measure goodness of fit)
    delta % calculated error (used for plotting)
    systematicError % mean (optimal cavity length - nominal cavity length) 
    estimationError % rms ((estimated cavity length at start of optimization - nominal cavity length) 
    fminIndx % index of mininum frequency over which to compute (stimulus must contain this frequency)
    fmaxIndx % index of maximum frequency over which to compute (stimulus must contain this frequency)
    freq  % frequency vector (Hz), from 0 to fmax
    wn % wave number exponent
    cavityLengths_optimal % cavity lengths that minimize the error (based on fminsearch)
    cavityLengths_est % estimate of cavity lengths based on resonances in PL
    timeStamp_recordings % when the recordings were made
    timeStamp_analysis % when the analysis was done
    reflNum % reflectance of cavity
    reflDenom
end
properties (SetAccess = public)
    fs % sampling rate (Hz)
    cavityLengths_nominal % intended cavity lengths (cm)
    cavityDiameter % cavity diameters (cm)
    cavityTemperature % (C); 25 = room temp; 37 = body temp
    fmin % minimum frequency over which to do thevenin source calibration
    fmax % maximum frequency over which to do thevenin source calibration
    recordings; % measured waveforms in the time domain
    stimulus  % stimulus in the time domain  
    micSN % microphone serial number
    fileName % name of file where the recordings are saved. This will be used to save the calibration as well
    thevCalPathName
    cavityRecordingsPathName
    minimizationType = 2 % Minimize the error on (1) pressure or (2) impedance
    plot_dB = 1; % if ==1, plot magnitudes in dB, else plot magnitudes in linear units
    autoPlot = 1; %turn on and off automatic plotting when finished calculating
end
properties (Dependent = true, SetAccess = private) 
    d % used by Keefe equations to calculate c and z0
    c % speed of sound in air at cavity temperature
    z0 % characteristic impedance of cavity
    nCavities % number of cavities being used
    nfft % number of points in the fft (must be >= length of stimulus)   
end
methods
    function t = thevCal(recordingParams,stimulus,recordings)
        if nargin == 0
            % set default values if none were provided by user
            t.fs = 44100; % sampling rate (Hz)
            t.cavityTemperature = 25; % degrees C
            t.cavityDiameter = 0.8; % (cm)
            t.cavityLengths_nominal = [1.85, 2.56, 4, 5.4, 8.3]; % intended cavity lengths (cm)
            t.fmin = 200;
            t.fmax = 8000;
            warning(['ALERT: Default values used for fs, cavityTemperature,',...
                ' cavityDiameter, cavityLengths_nominal, fmin, and fmax.'])
        elseif nargin >= 1
            t.fs = recordingParams.fs; % sampling rate (Hz)
            t.cavityTemperature = recordingParams.cavityTemperature; % degrees C
            t.cavityDiameter = recordingParams.cavityDiameter; % (cm)
            t.cavityLengths_nominal = recordingParams.cavityLengths_nominal; % intended cavity lengths (cm)
            t.fmin = recordingParams.fmin;
            t.fmax = recordingParams.fmax;
            t.timeStamp_recordings = recordingParams.timeStamp;
            %t.micSN = recordingParams.micSN; 
            t.fileName = recordingParams.fileName;
            t.thevCalPathName = recordingParams.thevCalPathName;
            t.cavityRecordingsPathName = recordingParams.cavityRecordingsPathName;
        end
        if nargin >= 2
            t.stimulus = stimulus;
        end
        if nargin == 3
            t.recordings = recordings;
        end
        if nargin > 3
            error('ALERT: Too many input arguments.')
        end
    end
    function set.stimulus(obj,stimulus)
       [rows,cols] = size(stimulus);
        if min([rows,cols]) > 1
            error('stimulus must be a vector, not matrix.')
        end        
        obj.stimulus = stimulus(:); % stimulus waveform (ensure a column vector)
    end   
    function set.recordings(obj,recordings)
        obj.recordings = recordings; % stimulus waveform (ensure a column vector)
    end    
    function [] = calculate(obj)
        disp('Calculating Thevenin Source, please wait...')
        tic
        freqIndices(obj);
        waveNumber(obj);
        calculate_PL(obj);
        obj.reflNum = .995; % estimated cavity reflectance
        obj.reflDenom = .995;
        % ------------- Begin iterative fitting to find optimimized cavity lengths.
        %               There are two fitting options: 1) Minimize the
        %               error between ideal and calculated PRESSURE (Neely
        %               method). 2) Minimize the error between ideal and
        %               calculated IMPEDANCE (Goodman method).
        save(['C:\myWork\ARLas\temp.mat'],'obj') % Save to an external file to be used during fminsearch
        op = optimset('Display','off','MaxIter',15000); % optimize to find best tube parameters
        % fminsearch minimizes the output of the function findOptimalLengths.m, using 
        % cavityLengths_est as a starting point for the search
        
        % NOTE: the following lines changed 6/13/2017 by ssg: this is for
        % use with the new 10X system, and for that system the estimated
        % lengths from the resonance is less accurage and better
        % performance is obtained using this modified code
        if obj.minimizationType == 1 % minimize pressure error
            obj.cavityLengths_optimal = fminsearch('findOptimalLengths1',obj.cavityLengths_est,op);
        elseif obj.minimizationType == 2 % minimize impedance error
            obj.cavityLengths_optimal = fminsearch('findOptimalLengths2',obj.cavityLengths_est,op);
        end
%         if obj.minimizationType == 1 % minimize pressure error
%             obj.cavityLengths_optimal = fminsearch('findOptimalLengths1',obj.cavityLengths_nominal,op);
%         elseif obj.minimizationType == 2 % minimize impedance error
%             obj.cavityLengths_optimal = fminsearch('findOptimalLengths2',obj.cavityLengths_nominal,op);
%         end
        % ------------- End iterative fitting to find optimimized cavity lengths.        
        % optimize for cavity reflectance
        doPlot = 0;
        rMin = 0.95; % assumed minimum reflectance
        rMax = 1.00; % assumed maximum reflectance
        stepSize = 0.001;
        refl = (rMin:stepSize:rMax)';
        nn = length(refl);
        epsilon = zeros(nn,1);
        for ii=1:nn
            obj.reflNum = refl(ii);
            obj.reflDenom = refl(ii);
            calculate_ZLi(obj.cavityLengths_optimal,obj) % calculate ideal load impedance (ZLi)
            calculate_source(obj) % calculate source thevenin parameters (PS, ZS)
            calculate_ZL(obj) % calculate load impedance (ZL), given source parameters
            calculate_PLi(obj) % calculate ideal load pressure (PLi)
            calculateError(obj)  % done on truncated freq vector
            if obj.minimizationType == 1 % minimize pressure error
                epsilon(ii,1) = obj.epsilon.pT;
            elseif obj.minimizationType == 2 % minimize impedance error
                epsilon(ii,1) = obj.epsilon.zT;
            end
        end
        [y,indx] = min(epsilon);
        if doPlot == 1
            figure
            plot(refl,epsilon,'b*-')        
            hold on
            plot(refl(indx),y,'r*-')
        end
        obj.reflNum = refl(indx);
        obj.reflDenom = refl(indx);
        calculate_ZLi(obj.cavityLengths_optimal,obj) % calculate ideal load impedance (ZLi)
        calculate_source(obj) % calculate source thevenin parameters (PS, ZS)
        calculate_ZL(obj) % calculate load impedance (ZL), given source parameters
        calculate_PLi(obj) % calculate ideal load pressure (PLi)
        calculateError(obj)  % done on truncated freq vector
        obj.systematicError = mean(obj.cavityLengths_optimal' - obj.cavityLengths_nominal);
        obj.estimationError = sqrt(mean((obj.cavityLengths_est' - obj.cavityLengths_nominal).^2));
        obj.timeStamp_analysis = datestr(clock); % time stamp for when this calibration was done
        if obj.autoPlot == 1
            obj.plotZ
            obj.plotP
            obj.plotD
            obj.plotL
            obj.plotZS
            obj.plotPS
        end
        toc
    end
    function obj = freqIndices(obj)
        dummy = sum([isempty(obj.fs),isempty(obj.nfft),isempty(obj.fmax),isempty(obj.fmin)]);
        if dummy ~= 0
            error('Frequency indices cannot be dtermined without fs, nfft, fmax, and fmin.')
        end
        obj.freq = (0:1:obj.nfft-1)'*(obj.fs/obj.nfft); % full frequency vector
        [dummy,obj.fminIndx] = min(abs(obj.freq - obj.fmin)); % min index for truncated portion
        [dummy,obj.fmaxIndx] = min(abs(obj.freq - obj.fmax)); % max indx for truncted portion
        obj.freq = obj.freq(1:obj.fmaxIndx); % truncate freq vector at fmax
    end
    function [] = waveNumber(obj) % wave number, from Keefe (1984) tube equations
        w = 2 * pi * obj.freq; % radian frequency
        w(w<eps) = eps;   % avoid w = 0 to avoid later division by zero
        r = obj.cavityDiameter / 2; % radius in cm
        rho = 1.1769e-3 * (1 - 0.00335 * obj.d);
        eta = 1.846e-4 * (1 + 0.0025 * obj.d);
        Rv = r * sqrt(rho * w / eta);
        x = (1.045 + (1.080 + 0.750 ./ Rv) ./ Rv) ./ Rv;
        y = 1 + 1.045 ./ Rv;
        obj.wn = (w / obj.c) .* complex(x,y);
    end    
    function [] = calculate_PL(obj) % sound pressure measurement expresssed as a transfer function in the frequency domain
        S = repmat(fft(obj.stimulus),1,obj.nCavities); % stimulus matrix
        R = fft(obj.recordings); % recordings matrix
        obj.PL = R ./ S;    
        obj.PL = obj.PL(1:obj.fmaxIndx,:);
    end
    function [] = calculate_ZLi(cavityLengths,obj) % cavity impedance (theoretical)
  %         R = exp(-2 * obj.wn * cavityLengths'); % wave equation exponent.  Volume velocity is 0 at the closed end of the tube.
  %         obj.ZLi = obj.z0 .* (1 + R) ./ (1 - R); % Alternative notation: Zc = -iZo*cot(kL), where k = wavenumber
        R = exp(-2 * obj.wn * cavityLengths(:)'); % wave equation exponent.  Volume velocity is 0 at the closed end of the tube.
        % zc = z0 .* (1 + R) ./ (1 - R); % ideal cavity impedance; Alternative notation: Zc = -iZo*cot(kL), where k = wavenumber
        % this formulation of zc assumes perfectly reflective cavities, which
        % is not the case with real tubes. Reducing R in the numerator reduces
        % the notch depth; reducing R in the denominator reduces the peak
        % height. Both values of R can also be reduced. Here, we assume an
        % upper reflectance bound of 1 and a lower bound of .5. 
        obj.ZLi = obj.z0 .* (1 + (R*obj.reflNum)) ./ (1 - (R*obj.reflDenom)); % this reduces the notch depth
    end        
    function [] = calculate_source(obj) % calculate Thevenin source parameters
        nFreqs = size(obj.PL,1); % number of frequencies in analysis range
        obj.ZS = zeros(nFreqs,1); % initialize source impedance vector
        obj.PS = zeros(nFreqs,1); % initialize source pressure vector
        for ii=1:nFreqs
           z = obj.ZLi(ii,:).'; % cavity impedance (theoretical); looping across rows (frequency), taking all columns (cavities) simultaneously
           p = obj.PL(ii,:).'; % cavity pressure (measured); looping across rows (frequency), taking all columns (cavities) simultaneously
           A = [z -p]; % Zc - Pc = cavity impedance - cavity pressure (measured and defined as H)
           B = z .* p; % Zc * Pc = cavity impedance * cavity pressure
           x = A \ B; % matrix division
           obj.PS(ii) = x(1); % Ps = source pressure
           obj.ZS(ii) = x(2); % Zs = source impedance
        end
    end
    function [] = calculate_ZL(obj) % cavity impedance (calculated)
        obj.ZL = zeros(size(obj.PL)); % initialize
        for ii=1:obj.nCavities  % calculate load impedance for each cavity
            obj.ZL(:,ii) = obj.ZS .* obj.PL(:,ii) ./ (obj.PS - obj.PL(:,ii));
        end
    end
    function [] = calculate_PLi(obj) % calculate load pressure (ideal)
        nCavities = obj.nCavities;
        Ps = repmat(obj.PS,1,nCavities);
        Zs = repmat(obj.ZS,1,nCavities);
        obj.PLi = Ps .* obj.ZLi ./ (Zs + obj.ZLi); % cavity pressure (calculated)
    end    
    function [] = calculateError(obj)
        f1 = obj.fminIndx;
        f2 = obj.fmaxIndx;
        nFreqs = (f2 - f1) + 1;
        s1 = 0; s2 = 0; % initialize sums of squares variable
        accumulatedPressure = zeros(nFreqs,1); % calculating pressure error (Neely method)
        for ii=1:obj.nCavities % loop across cavities
           pressureDiff = obj.PL(f1:f2,ii) - obj.PS(f1:f2) .* obj.ZLi(f1:f2,ii) ./ (obj.ZS(f1:f2) + obj.ZLi(f1:f2,ii)); % Numerator |Pc - Ps * Zc / (Zs + Zc)|
           accumulatedPressure = accumulatedPressure + pressureDiff; % sum of pressure differences across cavities; only used if crosstalk is being considered
           s1 = s1 + sum(abs(pressureDiff).^2); % sum of squares:  square the numerator and sum across frequency;
           s2 = s2 + sum(abs(obj.PL(f1:f2,ii)).^2); % sum of squares: square the measured cavity pressure and sum across frequency
        end
        s1 = s1 - sum(abs(accumulatedPressure).^2) / obj.nCavities; % line not used if crosstalk isn't being considered
        obj.epsilon.pT = 10000 * s1 / s2; % divide the sum of squares of the pressure difference with the sum of squares of the measured pressure and multiply by 10,000
        % calculating impedance error (Goodman method)
        obj.delta = obj.ZL(f1:f2,:) - obj.ZLi(f1:f2,:);
        mag = abs(obj.delta);
        phi = angle(obj.delta);
        N = size(obj.delta(:),1);
        eMag = sqrt(sum(sum(mag.^2)))/N;
        w = abs(mag).^2;
        ePhi = sum(sum(abs(phi).*w))./sum(sum(w));
        eTotal = sum(sum(abs(obj.delta)))/N;
        obj.epsilon.zM = eMag;
        obj.epsilon.zP = ePhi;
        obj.epsilon.zT = eTotal;
    end
    function [] = plotZ(obj) % plot measured versus ideal cavity impedances
        if obj.plot_dB == 1
            zcMag = 20*log10(abs(obj.ZLi));
            zlMag = 20*log10(abs(obj.ZL));
            ytxt = 'Impedance Magnitude (dB ohms)';
        else 
            zcMag = abs(obj.ZLi);
            zlMag = abs(obj.ZL);
            ytxt = 'Impedance Magnitude (ohms)';
        end 
        h = figure(10);
        freq = obj.freq(1:obj.fmaxIndx,:)/1000;
        subplot(2,1,1) % Measured and Ideal Impedance: Magnitude
        plot(freq,zcMag,':') % cavity impedance (theoretical)
        hold on
        plot(freq,zlMag) % cavity impedance (calculated)
        xlabel('Frequency (kHz)','FontSize',10)
        ylabel(ytxt,'FontSize',10)
        hold off
        nDigits = 3;
        txt = ['Error (eM,eP,eT):    ',...
            num2str(obj.epsilon.zM,nDigits),'   ',num2str(obj.epsilon.zP,nDigits),'   ',...
            num2str(obj.epsilon.zT,nDigits)];
        title(txt,'FontSize',10,'FontWeight','bold')
        xlim([freq(obj.fminIndx),freq(obj.fmaxIndx)])
        subplot(2,1,2) % Measured and Ideal Impedance: Phase
        plot(freq,angle(obj.ZLi)/(2*pi),':') % cavity impedance (theoretical)
        hold on
        plot(freq,angle(obj.ZL)/(2*pi))  % cavity impedance (calculated)
        xlabel('Frequency (kHz)','FontSize',10)
        ylabel('Impedance Phase (cyc)','FontSize',10)
        hold off
        title('Ideal=dashed; Measured=solid','FontSize',10,'FontWeight','bold')
        xlim([freq(obj.fminIndx),freq(obj.fmaxIndx)])
    end
    function [] = plotP(obj) % plot ideal versus measured cavity pressures
        if obj.plot_dB == 1
            pcMag = 20*log10(abs(obj.PLi));
            plMag = 20*log10(abs(obj.PL));
            ytxt = 'Pressure Magnitude (dB ohms)';
        else 
            pcMag = abs(obj.PLi);
            plMag = abs(obj.PL);
            ytxt = 'Pressure Magnitude (ohms)';
        end 
        h = figure(11);
        freq = obj.freq(1:obj.fmaxIndx,:)/1000;
        subplot(2,1,1) % Measured and Ideal Impedance: Magnitude
        plot(freq,pcMag,':') % cavity pressure (theoretical)
        hold on
        plot(freq,plMag) % cavity pressure (measured)
        xlabel('Frequency (kHz)','FontSize',10)
        ylabel(ytxt,'FontSize',10)
        hold off
        nDigits = 3;
        txt = ['Error: ',num2str(obj.epsilon.pT,nDigits)];
        title(txt,'FontSize',10,'FontWeight','bold')
        xlim([freq(obj.fminIndx),freq(obj.fmaxIndx)])
        subplot(2,1,2) % Measured and Ideal Impedance: Phase
        plot(freq,angle(obj.PLi)/(2*pi),':') % cavity impedance (theoretical)
        hold on
        plot(freq,angle(obj.PL)/(2*pi))  % cavity impedance (calculated)
        xlabel('Frequency (kHz)','FontSize',10)
        ylabel('Pressure Phase (cyc)','FontSize',10)
        hold off
        title('Ideal=dashed; Measured=solid','FontSize',10,'FontWeight','bold')
        xlim([freq(obj.fminIndx),freq(obj.fmaxIndx)])
    end
    function [] = plotD(obj) % plot difference (error) between measured and ideal cavity impedance
        if obj.plot_dB == 1
            dMag = 20*log10(abs(obj.delta));
            ytxt = 'Error Magnitude (dB)';
        else 
            dMag = abs(obj.delta);
            ytxt = 'Error Magnitude (linear)';
        end 
        h = figure(12);
        freq = obj.freq(obj.fminIndx:obj.fmaxIndx,:)/1000;
        subplot(2,1,1) % Error magnitude
        plot(freq,dMag) 
        xlabel('Frequency (kHz)','FontSize',10)
        ylabel(ytxt,'FontSize',10)
        nDigits = 3;
        txt = ['Error (eM,eP,eT):    ',...
            num2str(obj.epsilon.zM,nDigits),'   ',num2str(obj.epsilon.zP,nDigits),'   ',...
            num2str(obj.epsilon.zT,nDigits)];
        title(txt,'FontSize',10,'FontWeight','bold')
        xlim([freq(1),freq(end)])
        subplot(2,1,2) % Measured and Ideal Impedance: Phase
        plot(freq,angle(obj.delta)/(2*pi)) % cavity impedance (theoretical)
        xlabel('Frequency (kHz)','FontSize',10)
        ylabel('Error Phase (cyc)','FontSize',10)
        title('Ideal=dashed; Measured=solid','FontSize',10,'FontWeight','bold')
        xlim([freq(1),freq(end)])
    end
    function [] = plotL(obj) % compare cavity lengths (nominal, estimated, optimal)
        figure(13)
        plot(obj.cavityLengths_nominal,'bo-')
        hold on
        plot(obj.cavityLengths_est,'g*-')
        plot(obj.cavityLengths_optimal,'r*-')
        xlabel('Cavity Number')
        ylabel('Length(cm)','FontSize',12)
        legend('nominal','peak estimate','optimal','Location','northwest')
        title(['Systematic Error = ',num2str(obj.systematicError,3),...
            '     Estimation Error = ',num2str(obj.estimationError,3)],'FontSize',12)
        xlim([0.5,obj.nCavities+0.5])
        set(gca,'XTick',(1:1:obj.nCavities))
    end
    function [] = plotZS(obj) % plot source impedance
        if obj.plot_dB == 1
            zsMag = 20*log10(abs(obj.ZS));
            ytxt = 'Impedance (dB)';
        else 
            zsMag = abs(obj.ZS);
            ytxt = 'Impedance (linear)';
        end 
        zsMag = zsMag(obj.fminIndx:obj.fmaxIndx,:);
        zsPhi = angle(obj.ZS)/(2*pi);
        zsPhi = zsPhi(obj.fminIndx:obj.fmaxIndx,:);
        freq = obj.freq(obj.fminIndx:obj.fmaxIndx,:)/1000;
        h = figure(14);
        subplot(2,1,1) % source impedance magnitude
        plot(freq,zsMag) 
        xlabel('Frequency (kHz)','FontSize',10)
        ylabel(ytxt,'FontSize',10)
        txt = 'Source Impedance (ZS)';
        title(txt,'FontSize',12,'FontWeight','bold')
        xlim([freq(1),freq(end)])
        subplot(2,1,2) 
        plot(freq,zsPhi) % source impedance Phase
        xlabel('Frequency (kHz)','FontSize',10)
        ylabel('Phase (cyc)','FontSize',10)
        xlim([freq(1),freq(end)])
    end
    function [] = plotPS(obj) % plot source pressure
        if obj.plot_dB == 1
            psMag = 20*log10(abs(obj.PS));
            ytxt = 'Pressure (dB)';
        else 
            psMag = abs(obj.ZS);
            ytxt = 'Pressure (linear)';
        end 
        psMag = psMag(obj.fminIndx:obj.fmaxIndx,:);
        psPhi = angle(obj.PS)/(2*pi);
        psPhi = psPhi(obj.fminIndx:obj.fmaxIndx,:);
        freq = obj.freq(obj.fminIndx:obj.fmaxIndx,:)/1000;
        h = figure(15);
        subplot(2,1,1) % source impedance magnitude
        plot(freq,psMag) 
        xlabel('Frequency (kHz)','FontSize',10)
        ylabel(ytxt,'FontSize',10)
        txt = 'Source Pressure (PS)';
        title(txt,'FontSize',12,'FontWeight','bold')
        xlim([freq(1),freq(end)])
        subplot(2,1,2) 
        plot(freq,psPhi) % source impedance Phase
        xlabel('Frequency (kHz)','FontSize',10)
        ylabel('Phase (cyc)','FontSize',10)
        xlim([freq(1),freq(end)])
    end    
    function [] = saveThev(obj)
        counter = 1;
        fileName = obj.fileName;
        pathName = obj.thevCalPathName;
        while exist([pathName,fileName,'_thevCal_',num2str(counter),'.mat'],'file') == 2
            counter = counter + 1;
        end
        t = obj;
        try
            save([pathName,fileName,'_thevCal_',num2str(counter),'.mat'],'t')
            OK = 1;
        catch ME
            OK = 0;
             errorTxt = {'  Issue: Error saving data.'
                 '  Action: Save aborted.'
                 '  Location: in thevCal.saveThev.'
                };
            errorMsgARLas(errorTxt);
        end
        if OK
             alertTxt = {'  Issue: Calibration.'
                 '  Action: Calibration successfully saved.'
                 '  Location: in thevCal.saveThev.'
                };
            alertMsgARLas(alertTxt);
        end
    end
    % Dependent functions ---------------------------------------------------
    function [d] = get.d(obj)
        d = obj.cavityTemperature - 26.85; % from Keefe (1984)
    end
    function [c] = get.c(obj)
        c = 3.4723e4 * (1 + 0.00166 * obj.d); % speed of sound, adjusted for temperature, from Keefe (1984)
    end
    function [z0] = get.z0(obj) % characteristic impedance, from Keefe (1984) tube equations
        rho = 1.1769e-3 * (1 - 0.00335 * obj.d);
        r = obj.cavityDiameter / 2; % radius in cm
        z0 = (rho * obj.c) / (pi * r^2); % z0 = Ro; characteristic impedance (z0) is constant & real; 
    end    
    function [nCavities] = get.nCavities(obj)
        if isempty(obj.recordings)
            nCavities = [];
        else
            nCavities = size(obj.recordings,2); 
        end
    end
    function [nfft] = get.nfft(obj) % number of points in the fft
        if ~isempty(obj.recordings) 
            nfft = size(obj.recordings,1);
        elseif ~isempty(obj.stimulus)
            nfft = size(obj.stimulus,1);
        else
            nfft = [];
        end
    end
    function [cavityLengths_est] = get.cavityLengths_est(obj) % lengths based on location of resonant peaks
        expectedResonances = obj.c ./ (2*obj.cavityLengths_nominal); % wavelengths (cm)
        percentCushion = 0.25;  % 25% empirically determined
        fminR = expectedResonances * (1-percentCushion);
        fmaxR = expectedResonances * (1+percentCushion);
        cavityLengths_est = zeros(obj.nCavities,1); % initialize variable
        for ii=1:obj.nCavities
            [dummy,rangeMin] = min(abs(obj.freq - fminR(ii)));
            [dummy,rangeMax] = min(abs(obj.freq - fmaxR(ii)));
            freq = obj.freq(rangeMin:rangeMax);
            [dummy,indxMax] = max(abs(obj.PL(rangeMin:rangeMax,ii)));
            if (indxMax == 1) || (indxMax == length(freq)) , % if the chosen value bumps up agains the edges
                cavityLengths_est(ii,1) = obj.cavityLengths_nominal(ii); % use the nominal value
            else
                resonantFreq = freq(indxMax);
                cavityLengths_est(ii,1) = (obj.c ./ resonantFreq) / 2; % resonant frequency (cm)
            end
        end
    end
end
end
