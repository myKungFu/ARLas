function [] = ARLas_couplerRecordings_DW10x(varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ARLas_couplerRecordings_DW10x(varargin);
%
% This code is used to make microphone corrections and to verify the Thevenin calibration. 
% This version (DW10x) is for the ER10X in Jeff Lichtenhan's lab. 
% It assumes that the 10x probe is attached to an ear bar is being used 
% to deliver the sound to guinea pigs. The 10x and ear bar should be attached to one 
% end of the customized brass coupler. A 1/8" condensor reference mic should be 
% attached to the opposite end. 
%
% OPTIONAL INPUT ARGUMENT:
% inputRef = struture specifying reference input
%
% OUTPUT ARGUMENT:
% isc = in-situ calibration structure.
%
% Authors: Shawn S Goodman, PhD & Jeffery T Lichtenhan PhD
% Auditory Research Lab, Dept. of Comm. Sciences & Disorders, The University of Iowa
% Auditory Physiology Lab, Dept. of Otolaryngology, Washington University School of Medicine, St. Louis
% Date: January 9, 2019
% Updated: October 10-22, 2019 -- ssg
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    obj = varargin{1};
    params = varargin{2};
    

% ---  Adjustable Parameteters --------------------------------------------
    Lvl = 75; % level of stimulus (dB type re: 1 Hz bin width)??
    nRepsChirp = 48;
    nRepsClick = 1024; % number of times to play the stimulus
    stimGain_dB = -8; % reduce so you don't overdrive the system; for click
% -------------------------------------------------------------------------
    
    % unpack the paramter structure
    applyMicCorrection = params.applyMicCorrection;
    fmin = params.fmin;
    fmax = params.fmax;
    testProbe = params.testProbe;
    folderName_mic = params.fileName_mic;
    folderName_thev1 = params.fileName_thev1;
    folderName_thev2 = params.fileName_thev2;
    fileName_mic = params.fileName_mic;
    fileName_thev1 = params.fileName_thev1;
    fileName_thev2 = params.fileName_thev2;
    createMicCorrection = params.createMicCorrection;

    fs = obj.fs; % get the system sampling rate

    [inputs,outputs] = hardwareSetup; %hardwareSetup; % read in the saved hardware setup
    if strcmp(testProbe,'A')
        input = inputs{1};        % ER10xA microphone
        output = outputs{1};      % ER10xA loudspeakers
    elseif strcmp(testProbe,'B')
        input = inputs{2};        % ER10xB microphone
        output = outputs{2};      % ER10xB loudspeakers
    end
    inputRef = inputs{4};         % GRAS 1/8" reference microphone
    
    [pathName_mic,folderName_mic,fileName_mic] = mostRecentCalibration('mic',input.label,[]);
    [pathName_thev1,folderName_thev1,fileName_thev1] = mostRecentCalibration('thev',output.label,1);
    [pathName_thev2,folderName_thev2,fileName_thev2] = mostRecentCalibration('thev',output.label,2);
    
    % -------------------------------------------------------------------------
    if applyMicCorrection == 1
        if ~isempty(fileName_mic)
            dummy = load([pathName_mic,folderName_mic,fileName_mic]);
            micCorrection = dummy.micCorrection; % microphone correction
        else
            micCorrection = 1; % convolving by this changes nothing
        end
    else
        micCorrection = 1; % convolving by this changes nothing
    end
    %--------------------------------------------------------------------------

    for jj=1:2 % loop across both output channels
    
        % LOAD AND CHECK THE THEVENIN CALIBRATION FILE
        if jj == 1
            [t,cut1,cut2,stimulus] = loadCal(pathName_thev1,folderName_thev1,fileName_thev1,stimGain_dB,fmin,fmax,obj);
        elseif jj == 2
            [t,cut1,cut2,stimulus] = loadCal(pathName_thev2,folderName_thev2,fileName_thev2,stimGain_dB,fmin,fmax,obj);
        end

        % RECORD THE CHIRP STIMULUS ---------------------------------------------
        obj.clearRecList % clear the previously used recording list
        obj.setRecList(input.ch,input.label,input.micSens,input.gain); % load the recording info for ER10X
        obj.setRecList(inputRef.ch,inputRef.label,inputRef.micSens,inputRef.gain); % recording info for reference mic
        obj.setNReps(nRepsChirp); % number of times to play stimulus
        obj.setFilter(1); % turn on default highpass filter
        obj.clearPlayList % clear the previously used playback list
        obj.setPlayList(stimulus,output.ch(jj)); % load the currently tested output channel
        obj.objPlayrec.run % run the stimulus
        if obj.killRun
           return
        end    
        [header,data] = obj.retrieveData(input.label); % get chirp recorded in the ear canal
        data = fastFilter(micCorrection,data);
        [recording,stimf] = cleanDataDW(data,fs,stimulus,header,cut1,cut2);
        %[headerRef,dataRef] = obj.retrieveData(inputRef.label); % get chirp recorded from the reference mic

        % IN-SITU CALIBRATION -------------------------------------------------
        isc = ARLas_inSituCal(t,recording,stimf); % perform an in-situ calibration

        % CREATE CLICK STIMULUS ------------------------------------------------
        len = 0.025; %0.01; % desired stimulus length (s).
        nSamples = round(len * fs); % number of samples in stimulus
        if mod(nSamples,2) ~= 0 % force to be even number
            nSamples = nSamples + 1;
        end
        %time = (0:1:nSamples-1)'/fs; % time in s
        S1 = zeros(nSamples,1); % initialize stimulus vector
        offsetSec = 0.005; % offset the click this far from the start
        offsetN = round(offsetSec * fs); % offset in samples
        S1(offsetN) = 0.999; % click stimulus

        % Part 2B of the calibration routine 
        f = [fmin,fmax]; % designate that the stimulus is broadband
        if createMicCorrection == 1 % if purpose was to make mic calibration
            [S1flt,a,errors] = ARLas_applyISC(isc,Lvl,'ipl',f,S1);
        else
            [S1flt,a,errors] = ARLas_applyISC(isc,Lvl,'fpl',f,S1);
        end
        k = (1./(max(abs(S1flt)))) * 0.2; % 0.25;
        S1flt = S1flt * k;

    % COLLECT THE FLAT FPL DATA --------------------------------------------------------
        obj.clearRecList % clear out whatever was used previously 
        obj.setRecList(input.ch,input.label,input.micSens,input.gain); % load the recording info for ER10X
        obj.setRecList(inputRef.ch,inputRef.label,inputRef.micSens,inputRef.gain); % recording info for reference mic
        obj.setNReps(nRepsClick); % number of times to play stimulus
        obj.setFilter(1);

        obj.clearPlayList % clear out whatever was used previously
        obj.setPlayList(S1flt,output.ch(jj));

        obj.objPlayrec.userInfo = []; % clear out previous, non-structure array info
        obj.objPlayrec.userInfo.Lvl = Lvl;
        obj.objPlayrec.userInfo.fs = fs;
        obj.objPlayrec.userInfo.isc = isc;
        obj.objPlayrec.userInfo.stim = S1flt;
        obj.objPlayrec.run % playback and record
        if obj.killRun
           return
        end
        [header,data] = obj.retrieveData(input.label); % get click recorded from ER10X mic
        data = data ./ k;
        [headerRef,dataRef] = obj.retrieveData(inputRef.label); % get click recorded from Reference mic
        dataRef = dataRef ./ k;
        data = fastFilter(micCorrection,data);
        data = ARLas_artifactReject(data);
        dataRef = ARLas_artifactReject(dataRef);
        % add filtering here?
        [pl,Pl,phi,other,wf] = ARLas_convertPL(data,isc);
        
        % because reference mic is low sensitivity, it has a high noise
        % floor. To overcome this without an inordinate amount of
        % averaging, use a spline smoother.
        smoothing = .0000000001; % amount of smoothing
        yy = mean(dataRef,2); % mean reference mic recording (time domain)
        nn = length(yy); % number of samples in the recording
        nfft = fs; % number of fft samples (zero padded)
        YY = fft(yy,nfft); % (frequency domain)
        xx = (0:1:nfft-1)'*(fs/nfft); % frequency vector (Hz);
        ww = ones(size(xx));
        y = abs(YY);
        mag = csaps(xx,y,smoothing,xx,ww); % the smoothed, weighted spline, densly-spaced estimate
        y = unwrap(angle(YY));
        ang = csaps(xx,y,smoothing,xx,ww); % the smoothed, weighted spline, densly-spaced estimate
        ang = angle(cos(ang) + 1i*sin(ang));
        Y = mag.*cos(ang) +1i*mag.*sin(ang);
        Signal = interp1(xx,Y,pl.f);
        scaling = nn/2;
        splRef = 20*log10(abs(Signal) / scaling / .00002); % dB SPL of ref mic
        fRef = pl.f; % reference frequency vector

        epsilon = pl.ipl-splRef; % error (dB) at each frequency
            
        wf10X = mean(data,2)*1000;
        wfRef = mean(dataRef,2)*1000;
        time = (0:1:length(S1flt)-1)'; % create time vector
        [~,xhat] = max(abs(wf10X)); % the location of the maximum peak in the waveform
        time = time - xhat; % max peak will be time zero
        time = time / fs;
        time = time * 1000; % express time in ms

        if jj == 1
            H = [];
        end
        H = doPlot(time,S1flt,wf10X,wfRef,pl,fRef,splRef,epsilon,fmin,fmax,H);
    
        % purpose of this code is either to make a mic correction or to
        % verify the total (Thev + mic) calibration.
        if createMicCorrection == 1 % if purpose was to make mic calibration
            iplWaveform = wf.ipl;
            refWaveform =  mean(dataRef,2);
            nI = size(iplWaveform,1); % make sure the waveforms are the same size
            nR = size(refWaveform,1);
            nD = nI - nR;
            if nD > 0
                refWaveform = [refWaveform;zeros(nD,1)];
            elseif nD < 0
                iplWaveform = [iplWaveform;zeros(abs(nD),1)];
            end
            if jj==1
                %iplWF1 = iplWaveform;
                %refWF1 = refWaveform;
                micCorrection1 = ARLas_makeMicCorrection_DW10x(iplWaveform,refWaveform,fmin,fmax,fs);
            elseif jj==2
                %iplWF2 = iplWaveform;
                %refWF2 = refWaveform;
                micCorrection2 = ARLas_makeMicCorrection_DW10x(iplWaveform,refWaveform,fmin,fmax,fs);
            end
        end
    end
    %mc = ARLas_makeMicCorrection2_DW10x(iplWF1,refWF1,iplWF2,refWF2,fmin,fmax,fs);
    if createMicCorrection == 1 % ip purpose was to make mic calibration
        % mic correction saved will be the average of the two estimates.
        % Average in the frequency domain
        % THIS DOESN'T WORK!
        %MC1 = fft(micCorrection1);
        %MC2 = fft(micCorrection2);
        %MC = (MC1 + MC2)/2;
        %micCorrection = real(ifft(MC));
        micCorrection = micCorrection2; % for now, can only use one. 
        
        tt = datetime('now'); % create the folder name for saving microphone calibrations
        folderName = datestr(tt,'mm_dd_yyyy');
        folderName = [folderName,'\']; 
        if exist([pathName_mic,folderName],'dir') ~= 7 % check to make sure that the save path exists
            success = mkdir([pathName_mic,folderName]); % if not, try to create it
            if success ~= 1
                warning('Specified folder for saving mic cal does not exist. Entering debug mode.')
                keyboard
            else
                addpath(genpath(pathName_mic)) % add all directories and subdirectories
            end
        end
        fileName_mic = [input.label,'.mat'];
        saveName = ARLas_saveName([pathName_mic,folderName],fileName_mic);
        save([pathName_mic,folderName,saveName],'micCorrection','fmin','fmax','fs','testProbe','folderName')
    end    

end

% Internal Functions ------------------------------------------------------
function [H] = doPlot(time,S1flt,wf10X,wfRef,pl,fRef,splRef,epsilon,fmin,fmax,H)
    if isempty(H)
        h = figure;
        h.Position = [392.2000   80.2000  503.2000  657.6000];
        sp1 = subplot(4,1,1);
        sp2 = subplot(4,1,2);
        sp3 = subplot(4,1,3);
        sp4 = subplot(4,1,4);
        sp1.Position = [0.1300    0.8382    0.7750    0.0868];
        sp2.Position = [0.1300    0.6655    0.7750    0.1667];
        sp3.Position = [0.1300    0.2713    0.7750    0.2883];
        sp4.Position = [0.1300    0.1100    0.7750    0.1564];
        H.h = h;
        H.sp1 = sp1;
        H.sp2 = sp2;
        H.sp3 = sp3;
        H.sp4 = sp4;
    else
        h = H.h;
        sp1 = H.sp1;
        sp2 = H.sp2;
        sp3 = H.sp3;
        sp4 = H.sp4;
    end

    figure(h)
    axes(sp1)
    plot(time,S1flt,'k') % plot of electrical stimulus
    grid on
    xlim([-1 2])
    set(gca,'Xticklabel',[])
    ymax = max(abs(S1flt))*1.05;
    ylim([-ymax ymax])
    legend('Electrical Stim')
    ylabel('AMP')

    axes(sp2)
    plot(time,wf10X,'b') % plot of microphone recordings in time
    hold on
    plot(time,wfRef,'r')
    grid on
    xlim([-1 2])
    ymax = max([max(abs(wf10X)),max(abs(wfRef))])*1.05;
    ylim([-ymax ymax])
    legend('ER10X mic','Ref mic')
    xlabel('TIME (ms)')
    ylabel('AMP (mPa)')

    axes(sp3)
    plt(pl.f,pl.spl,'m')
    hold on
    plt(pl.f,pl.fpl,'b')
    plt(pl.f,pl.rpl,'r')
    plt(pl.f,pl.ipl,'k')
    plt(fRef,splRef,'c','LineWidth',2)
    grid on
    xlim([fmin,fmax])
    ymax = max([max(pl.ipl),max(pl.spl)])*1.05;
    ylim([ymax-30, ymax])
    xlabel('FREQUENCY (kHz)')
    set(gca,'Xticklabel',[])
    ylabel('MAG (dB SPL)')
    legend('ER10X SPL','ER10X FPL','ER10X RPL','ER10X IPL','REF SPL','Location','SouthWest')

    axes(sp4)
    plt(fRef,epsilon)
    hold on
    legend('ER10x IPL - Ref SPL','Location','NorthWest')
    grid on
    ylabel('ERROR (dB)')
    xlim([fmin,fmax])
    xlabel('FREQUENCY (kHz)')
end
function [t,cut1,cut2,stimulus] = loadCal(pathName,folderName,fileName,stimGain_dB,fmin,fmax,obj)
    dummy = load([pathName,folderName,fileName]);
    t = dummy.t;
    if t.fmin > fmin
        error('Requested value of fmin is less than the calibrated value.')
    end
    if t.fmax < fmax
        error('Requested value of fmax is greater than the calibrated value.')
    end
    if obj.fs ~= t.fs
        error('Current sampling rate is different from the calibrated value.')
    end
    stimulus = t.stimOrig; % get the logarithmic chirp originally used to calibrate
    stimulus = stimulus * 10^(stimGain_dB/20); % scale down so don't overdrive output
    cut1 = t.cut1;         % get the original cut samples
    cut2 = t.cut2;
end

% OLD CODE ----------------------------------------------------------------
%     [frequency,signal,noiseFloor] = ARLas_fda(dataRef(1:1000,:),fs,0.00002,fs);
%     figure
%     subplot(2,1,1)
%     plot(frequency,signal,'r')
%     hold on
%     plot(frequency,noiseFloor,'k')
%     title('Reference')
%     [frequency,signal,noiseFloor] = ARLas_fda(data(1:1000,:),fs,0.00002,fs);
%     subplot(2,1,2)
%     plot(frequency,signal,'b')
%     hold on
%     plot(frequency,noiseFloor,'k')
%     title('10X')

%     % what is the group delay measured by the reference microphone?
%     % relative to the fpl group delay?
%     
%     FPL = fft(wf.fpl);
%     REF = fft(mean(dataRef,2));
%     Xfer = FPL ./ REF; % microphone transfer function
%     
%     smoothing = .000000001;
%     N = size(Xfer,1);
%     xx = (0:1:N-1)' * (fs/N);
%     x = xx;
%     w = ones(size(xx));
%     [~,indxCut1] = min(abs(x-fmin));
%     [~,indxCut2] = min(abs(x-fmax));
%     w(1:indxCut1) = 0;
%     w(indxCut2:end) = 0;
%     
%     y = abs(Xfer);
%     mag = csaps(x,y,smoothing,xx,w); % the smoothed, weighted spline, densly-spaced estimate
%     %y = unwrap(angle(Xfer));
%     %ang = csaps(x,y,smoothing,xx,w); % the smoothed, weighted spline, densly-spaced estimate
%     %ang = angle(cos(ang) + 1i*sin(ang));
%     ang = zeros(size(mag));
%     
% %     % adjsut the phase by subtracting off constant group delay that
% %     % represents the travel time in the tube
% %     q = unwrap(ang(indxCut1:indxCut2));
% %     r = x(indxCut1:indxCut2);
% %     p = polyfit(r,q,1);
% %     gd = p(1)*r + p(2);
% %     q = q - gd;
% %     ang(indxCut1:indxCut2) = q;
% %     ang(1:indxCut1-1) = 0;
% %     ang(indxCut2+1:end) = 0;
% %     mag(1:indxCut1-1) = 0;
% %     mag(indxCut2+1:end) = 0;
%     
%     % this is the mic transfer function. Invert to turn into a filter
%     mag = 1./mag;
%     mag(1:indxCut1-1) = 0;
%     mag(indxCut2+1:end) = 0;
%     ang = -ang;
%     nyquist = (N/2)+1;
%     mag(nyquist+1:end) = flipud(mag(2:nyquist-1));
%     ang(nyquist+1:end) = -flipud(ang(2:nyquist-1));
%     Y = mag.*cos(ang) +1i*mag.*sin(ang);
%     y = real(ifft(Y));
%     yy = [y(nyquist:end);y(1:nyquist-1)];
%     m = 96; % 1 ms group delay puts you within 1 dB of target values
%     m2 = m/2;
%     zz = yy(nyquist-m2:nyquist+m2);
%     h = hann(length(zz));
%     zz = zz .* h;
%     b = zz;
%     plt(x,20*log10(abs(fft(zz,2400))))
%     hold on
%     plt(x,20*log10(mag))
%     

% old cod3e
% 
%     hold on
%     plt(pl.f,pl.fpl)
%     plt(pl.f,pl.epl,'r')
%     plt(pl.f,pl.fplnf,'r')

    %BW = fmax - fmin; % bandwidth (Hz)
    %LPC = Lvl - 10*log10(BW); % level per cycle (1-Hz binwidth)

    %             
%             Q = fft(mean(data,2));
%             Q = Q ./ (size(data,1)/2);
%             q = abs(Q);
%             q = 20*log10(q/.00002);
%             fq = (0:1:size(data,1)-1)'*(fs/size(data,1));
%             
%             Q2 = fft(mean(data,2));
%             Q2 = Q2.^2;
%             Q2 = Q2 ./ (size(data,1)^2/2);
%             q2 = abs(Q2);
%             q2 = 10*log10(q2/.00002^2);
%             fq = (0:1:size(data,1)-1)'*(fs/size(data,1));
