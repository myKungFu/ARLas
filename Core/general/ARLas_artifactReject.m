function [data] = ARLas_artifactReject(data)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% data = ARLas_artifactReject(data);
%
% Artifact rejection routine for use with ARLas (Auditory Research Laboratory auditory software)
%
% Auditory Research Lab, The University of Iowa
% Deptartment of Communication Sciences & Disorders
% The University of Iowa
% Author: Shawn S. Goodman, PhD
% Date: October 27, 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This method uses quartiles, suggested by Tukey.

rms = sqrt(mean(data.^2,1)); % compute the rms amplitude of each buffer
tolerance = 'moderate';
doPlot = 0; % turn off plotting
[indx,nRejects] = AR(rms,tolerance,doPlot); % find the buffers to reject
data = AR_engine(data,indx); % actually peform the rejections