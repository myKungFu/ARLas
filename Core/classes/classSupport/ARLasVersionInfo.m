function [] = ARLasVersionInfo()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ARLasVersionInfo
%
% Contains a list of the updates made to ARLas.
% Starting date is November 16, 2016
% 
% Auditory Research Lab, The University of Iowa
% Deptartment of Communication Sciences & Disorders
% The University of Iowa
% Author: Shawn S. Goodman, PhD
% Date: November 16, 2016
% Last Updated: November 29, 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 2016.11.16
%   Changed handling of initialization defaults.
%       Now uses the desired names and looks for appropriate indices.
%       This change was precipitated by problems with the operating
%       systems frequently changing device parameters. 
%   Fixed bugs in the Abort Experiment button.
%   Fixed bugs in the deleting of initialization gui
% 2016.11.29
%   Changed handling stimTrain.
%       Now stimTrain is a structure, with a field for each output channel.
%       This allows each output to be a vector or a matrix. 
%       If matrices are present, playrecARLas.queue loops over each column.
%       This allows random noise to be played in one or more channels,
%       among other things. 
%       THERE IS NO BACKWARDS COMPATIBILITY FOR THIS UPDATE: All experiment
%       files must be changed so that they load stimTrain using the channel
%       fields. Example of old version: stimTrain = randn(44100,2);
%               Example of new version: stimTrain.Ch1 = randn(44100,1);
%                                       stimTrain.Ch2 = randn(44100,1);
% 2016.12.07
%   Added userInfo as a public property of playrecARLas.
%       This will be saved as part of the header file.
%       This allows users to save important information from within
%       experiment files. 
%               Example: f2 = 1000; % frequency in Hz
%                        obj.objPlayrec.userInfo.f2Freq = 1000; 
%                        obj.objPlayrec.run % run the stimulus                         
%   Added backwards compatibility to stimTrain handling. The old way of
%       using stimTrain (simply make a matrix, with each column being an
%       output channel can still be used. A new property, stimTrainEx
%       (Ex stands for "explicit") can also be used to specify channels
%       explicitly, and to use either vectors or matrices for each channel.
%               Example of old version: stimTrain = randn(44100,2);
%               Example of new version: stimTrainEx.Ch1 = randn(44100,1);
%                                       stimTrainEx.Ch2 = randn(44100,1);
