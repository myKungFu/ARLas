function [pathName_speakerCal,fileName_speakerCal,fileName_earBarCal,pathName_micCal,fileName_micCal] = oldCalInfo()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [pathName_speakerCal,fileName_speakerCal,fileName_earBarCal,pathName_micCal,fileName_micCal] = oldCalInfo;
%
% With the start of the new Thevenin-based calibration, there was a need to
% retain some of the old calibrations for use with old experiments;
% however, the old calibrations needed to be moved to new locations. Old
% experiments therefore will call this function to get the new locations.
%
% Eventually, this will be phased out.
%
% Authors: Shawn S Goodman, PhD & Jeffery T Lichtenhan PhD
% Auditory Research Lab, Dept. of Comm. Sciences & Disorders, The University of Iowa
% Auditory Physiology Lab, Dept. of Otolaryngology, Washington University School of Medicine, St. Louis
% Date: November 6, 2019
% Updated: November 6, 2019 -- ssg
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% CALIBRATION re. SHAWN & LAURA VISIT ON 8/9/17 --------
pathName_speakerCal = 'C:\myWork\ARLas\Peripheral\calibrations\oldCals\'; % location of loudspeaker calibrations
pathName_micCal = 'C:\myWork\ARLas\Peripheral\calibrations\oldCals\'; % location of microphone calibrations

fileName_speakerCal = 'speakerCal_10xB_1.mat';
fileName_earBarCal = 'speakerCalEarBar_10xB_1.mat';            
fileName_micCal = 'micCal_10xB_avg12-15PLASTIC.mat';

end

% OLD CODE ----------------------------------------------------------------
%pathName_speakerCal = 'C:\myWork\ARLas\Peripheral\calibrations\speakerCals\'; % location of loudspeaker calibrations
%pathName_micCal = 'C:\myWork\ARLas\Peripheral\calibrations\micCals\'; % location of microphone calibrations
