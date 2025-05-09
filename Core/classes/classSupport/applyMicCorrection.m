function [Data] = applyMicCorrection(Data,micCorrection)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Data = applyMicCorrection(Data,micCorrection);
%
% For use with ARLas (Auditory Research Laboratory auditory software)
% Apply microphone correction to the recorded Data matrix.
%
% Data (input argument) = matrix of input waveforms.
% micCorrection = microphone correction, as a time-domain inpulse response.
%
% Data (output argument) = the recorded data, corrected for microphone magnitude and phase. 
%
% Auditory Research Lab, The University of Iowa
% Deptartment of Communication Sciences & Disorders
% The University of Iowa
% Author: Shawn Goodman
% Date: July 1, 2021
% Last Updated: March 23, 2022 -- ssg -- added line to return without
%                                   convolving if micCorrection value = 1.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % added this to speed up processing; no need to convolve if no mic
    % correction being applied
    if micCorrection == 1
        return
    end

    [Rows,Cols] = size(Data);
    Data = Data(:);
    Data = fastFilter(micCorrection,Data);
    Data = reshape(Data,Rows,Cols);
    
end