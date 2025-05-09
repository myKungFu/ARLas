function [] = ARLas_micCal_v3(varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ARLas_micCal_v3(varargin)
%
% Create microphone calibration using Etymotic's provided mic calibrations.
% The mic calibration excel file was provided by Steve Viranyi on Sept 29,
% 2021 in an email sent to Shawn, Jeff, Vicky, and Dan.
% This version (_v3) reads in the data from the excel file. The Excel file
% contains microphone corrections (amplitude only--no phase corrections!)
% for all three ER10X systems (Jeff's and Shawn's from Wash U and FX's).
%
% As with previous versions, creates an impulse response, saves results into 
% a structure called LTC (long tube calibration), located in 
% 'C:\myWork\ARLas\Peripheral\calibrations\LTCals\' in a sub-folder under the date of the calibration.
%
% Author: Shawn Goodman, PhD
% Original MicCal Date: June 7, 2021
% Last Updated: September 17, 2021 -- ssg -- version 2
% Last Updated: October 19, 2021 -- ssg -- version 3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    obj = varargin{1}; % get the arlas object
    params = varargin{2};
    V = '19OCT2021'; % this is the current version number of this program

%------ USER MODIFIABLE PARAMETERS ----------------------------------------
%--------------------------------------------------------------------------

% There are none.

%------ END USER MODIFIABLE PARAMETERS ------------------------------------
%--------------------------------------------------------------------------

    fmin = params.fmin;
    fmax = params.fmax;
    probe = params.testProbe;
    which10X = params.which10X;
    % look beyond the required boundaries
    fmaxOrig = fmax;
    fmax = fmax * 1.2;
    fminOrig = fmin;
    fmin = fmin * 0.8;
    
    % get probe settings
    probeL = 'A'; 
    probeR = 'B'; 
    [probeInputL,probeOutputL,probeInputR,probeOutputR,refMic] = getProbeSettings(probeL,probeR,[]);
    if strcmp(probe,'A')
        probeInput = probeInputL;
        probeOutput = probeOutputL;
    elseif strcmp(probe,'B')
        probeInput = probeInputR;
        probeOutput = probeOutputR;
    end

    params.fmin = fmin;
    params.fmax = fmax;
    params.probe = probe;
    params.folderName_mic = []; % no file names or folders typically given to start
    params.folderName_longTube1 = [];
    params.folderName_longTube2 = [];
    params.fileName_mic = []; % no file names or folders typically given to start
    params.fileName_longTube1 = [];
    params.fileName_longTube2 = [];
    
    % Column Letter:D           E       
    % Mic Serial # :141         142     
    % Name         :Probe A     Probe B 
    if strcmp(probe,'A')
        col = 'D';
    elseif strcmp(probe,'B')
        col = 'E';
    end
    
    disp(' '); disp(' '); disp(' '); disp(' ')
    alertTxt = {'Reading Etymotic Microphone Corrections.'
         ['  Calibrating probe ',probe,'.']
        };
    nn = size(alertTxt,1);
    for ii=1:nn
        cprintf([0,0,.4],[alertTxt{ii},'\n']);
    end
    
    % read in the correction values ---------------------------------------
    fsep = filesep;
    FILE = [obj.map.experiments,'CALIBRATION',fsep,'MicData.xlsx']; %'C:\myWork\ARLas\Peripheral\experiments\CALIBRATION\MicData.xlsx';
    SHEET = 'data';
    RANGE = [col,'8:',col,'305']; % 303
    [correction,~,~]=xlsread(FILE,SHEET,RANGE); % read in magnitude correction
    RANGE = ['C8:C303'];
    [freq,~,~]=xlsread(FILE,SHEET,RANGE); % read in frequency vector
    RANGE = [col,'7'];
    [probeSN,~,~]=xlsread(FILE,SHEET,RANGE); % read in probe serial number
    
    fs = obj.fs; % get the system sampling rate
    micCorrection = makeMicCorrection(freq,correction,fmin,fmax,fs);

    fmin = fminOrig; % set back to original values
    fmax = fmaxOrig;
    
    % this stuff for plotting purposes ------------------------------------
    FFT = fft(micCorrection,fs); % matrix of individual fourier transforms
    FT = 20*log10(abs(FFT)); % matrix of individual fourier transforms--magnitude
    P = unwrap(angle(FFT)) / (2*pi);
    F = (0:1:fs-1)*(fs/fs); % frequency vector
    F = F / 1000; % frequency in kHz
    mc = mean(micCorrection,2);
    MC = 20*log10(abs(fft(mc,fs)));
    %MCp = unwrap(angle(fft(mc,fs))) / (2*pi);
    
    % plot the results
    figure
    subplot(2,1,1) % impulse response of mic correction
        plt(micCorrection,'b','LineWidth',1)
        xlim([1,length(micCorrection)])
        xlabel('Time (samples)')
        ylabel('Amplitude')
        title(['ER10X',probe,' SN:',num2str(probeSN),' Mic Corrections from Etymotic'])
    subplot(2,1,2) % magnitude of mic correction
        plt(F,MC,'b','LineWidth',1)
        xlim([fmin/1000,fmax/1000])
        xlabel('Frequency (kHz)')
        ylabel('Magnitude')
    pause(0.05)
    
        
        % save new mic calibration
        [pathName_mic,~,~] = mostRecentCalibration('mic',probeInput.label,[],obj.map);
        
        tt = datetime('now'); % create the folder name for saving microphone calibrations
        folderName = datestr(tt,'mm_dd_yyyy');
        folderName = [folderName,fsep]; 
        if exist([pathName_mic,folderName],'dir') ~= 7 % check to make sure that the save path exists
            success = mkdir([pathName_mic,folderName]); % if not, try to create it
            if success ~= 1
                warning('Specified folder for saving mic cal does not exist. Entering debug mode.')
                keyboard
            else
                addpath(genpath(pathName_mic)) % add all directories and subdirectories
            end
        end
        fileName_mic = [probeInput.label,'.mat'];
        saveName = ARLas_saveName([pathName_mic,folderName],fileName_mic);
        save([pathName_mic,folderName,saveName],'micCorrection','fmin','fmax','fs','probeInput','folderName','probeSN','which10X')

end

% INTERNAL FUNCTIONS ------------------------------------------------------
function [b] = makeMicCorrection(freq,correction,fmin,fmax,fs)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% b = makeMicCorrection(freq,correction,fmin,fmax,fs);
%
% Create microphone correction. 
% Assuming frequency range from 100 - 20000 Hz. 
%
% Input Arguments:
%  waveform = measured sound pressure level in the calibration assembly.
%  refWaveform = measured reference waveform from the condensor mic.
%  fmin = minimum frequency to include (100 Hz)
%  fmax = maximum frequency to include (20000 Hz)
%  fs = sampling rate (Hz; assuming 96000)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% --------------------- Adjustable Parameters -----------------------------
    filterGroupDelay = 0.0015; % desired filter group delay (s)
%--------------------------------------------------------------------------    

    M = correction; % magnitude of mic response
    offset = mean(M(1:100)); % mean offset over low frequencies
    M = M - offset; % subtract off the offset to make low frequencies zero
    M = -M; % correction is minus the response (in dB)
    N = length(freq); % original number of samples
    FF = (freq(1):1:freq(end))'; % dense frequency vector, in steps of 1 Hz
    [~,indxCut1] = min(abs(FF-fmin)); % cutoff indices
    [~,indxCut2] = min(abs(FF-fmax));
    FF = FF(indxCut1:indxCut2);
    MM = interp1(freq,M,FF,'pchip');
    
    mag = 10.^(MM/20); % put magnitude into linear units
    ang = zeros(size(mag)); % phase correction is zero (was not measured)
    
    Mag = zeros(fs,1);
    Phi = zeros(fs,1);
    NN = length(Mag);
    Mag(indxCut1:indxCut2) = mag;
    nyquist = NN/2 + 1;
    Mag(nyquist+1:end) = flipud(Mag(2:nyquist-1));
    Y = Mag.*cos(Phi) +1i*Mag.*sin(Phi); % complex rectangular form
    y = real(ifft(Y)); % time domain waveform
    yy = [y(nyquist:end);y(1:nyquist-1)]; % make filter causal
    m = round(filterGroupDelay * fs); % filter order
    if mod(m,2)~= 0 % force group delay to be even so that can remove it after filtering
        m = m + 1;
    end
    m2 = m/2; % half group delay
    
    yy = yy(nyquist-m2:nyquist+m2); % cut down filter to desired size
    h = hann(length(yy)); % window to smooth edges
    b = yy .* h; % filter coefficients in the time domain.    
end

% OLD CODE ----------------------------------------------------------------
