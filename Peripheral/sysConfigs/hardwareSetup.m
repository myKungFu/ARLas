function [inputs,outputs] = hardwareSetup()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [inputs,outputs] = hardwareSetup;
%
% Get hardware inputs and outputs for use with ARLas.
% This version simply calls the new version (hardwareSetup_DW). 
% Included for sake of backwards compatibility. 
%
% Authors: Shawn S Goodman, PhD & Jeffery T Lichtenhan PhD
% Auditory Research Lab, Dept. of Comm. Sciences & Disorders, The University of Iowa
% Auditory Physiology Lab, Dept. of Otolaryngology, Washington University School of Medicine, St. Louis
% Date: October 10, 2019
% Updated: October 10, 2019
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    try
        [inputs,outputs] = hardwareSetup_DW2; % try Shawn's proxy first
    catch
        [inputs,outputs] = hardwareSetup_DW;  % if not present, assume Jeff's setup
    end

end