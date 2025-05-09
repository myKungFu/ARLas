function [X,stim] = cleanDataDW(X,fs,stim,header,cut1,cut2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [X,stim] = cleanDataDW(X,fs,stim,header,cut1,cut2);
%
% For use with ARLas (Auditory Research Laboratory auditory software)
% Used in applying Thevenin source calibration to cavity and/or ear canal recordings.
% This version (DW) is for use in the Lichtenhan lab at Wash U. It uses a
% high-frequency cutoff of 32 kHz instead of the usual 20 kHz. 96 kHz
% sampling rate only is allowed.
%
% X (input argument) = matrix of input waveforms
% fs = sampling rate (Hz)
% stim (input argument) = waveform vector; the original electrical stimulus that was used to obtain X
% cut1 = cutoff value to remove excess zero padding at the start of the signal
% cut2 = cutoff value to remove excess zero padding at the end of the signal
%
% X (output argument) = the lowpass Butterworth filtered (22 kHz cutoff), 
%       truncated waveform w/ artifact rejection and 3 ms onset/offset windowing
% stim (output argument) the stimulus with the same filtering and windowing applied.
%
% Authors: Shawn S Goodman, PhD & Jeffery T Lichtenhan PhD
% Auditory Research Lab, Dept. of Comm. Sciences & Disorders, The University of Iowa
% Auditory Physiology Lab, Dept. of Otolaryngology, Washington University School of Medicine, St. Louis
% Date: January 8, 2019
% Updated: October 22, 2019 -- ssg
%                              Added extended bandwidth on filter. Now goes out to 42 kHz.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    [rows,cols] = size(X); % orignal size of the recording
    % it is important to do the same processing to the stimulus as to the
    % recorded waveforms.
    stim = repmat(stim,1,cols+2); % replicate stimulus to matrix the same size as recording
    stim = stim(:); % reshape into a long vector
    X = X(:);
    stim = filtfilt(header.SOS,header.G,stim); % apply the playrec highpass filter to stim (same as was for recording)
    %[SOS,G] = getLPiirDW(fs); % get coefficients for a lowpass filter, up to 32 kHz
    [SOS,G] = getLPiir_EB(fs); % include out to 42 kHz
    stim = filtfilt(SOS,G,stim); % apply the lowpass filter
    X = filtfilt(SOS,G,X);
    X = reshape(X,rows,cols); % reshape into matrices
    stim = reshape(stim,rows,cols+2);
    stim = stim(:,2:end-1);
    tolerance = 'moderate'; doPlot = 0; % perform artifact rejectioin
    rms = sqrt(mean(X.^2,1));
    [indx,nRejects] = AR(rms,tolerance,doPlot);
    %disp([num2str(nRejects),' artifacts'])
    X = AR_engine(X,indx);
    stim = AR_engine(stim,indx);
    X = mean(X,2); % return mean    
    stim = mean(stim,2);
    % cut off extra zero padding on the ends
    X = X(cut1:cut2,:);
    stim = stim(cut1:cut2);
    % apply a hann window to the edges
    % fmin = 100; the period is 1/100 = 10 ms. use this for ramp on
    cut1 = round(0.01 * fs);
    % fmax = 32000 Hz; however, this period is too short for a window.
    % Instead, use 2 ms (which is actually 64 periods of 32 kHz)
    cut2 = round(0.002 * fs);
    h = hann(cut1*2);
    h = h(1:cut1);
    H = repmat(h,1,size(X,2));
    X(1:cut1,:) = X(1:cut1,:) .* H;
    stim(1:cut1) = stim(1:cut1) .* h;

    h = hann(cut2*2);
    h = h(1:cut2);
    H = repmat(h,1,size(X,2));
    X = flipud(X);
    stim = flipud(stim);
    X(1:cut2,:) = X(1:cut2,:) .* H;
    stim(1:cut2) = stim(1:cut2) .* h;
    X = flipud(X);
    stim = flipud(stim);
end
