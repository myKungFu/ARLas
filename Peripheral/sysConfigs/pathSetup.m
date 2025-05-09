function [map] = pathSetup(base)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [map] = pathSetup(base);
%
% Get file path setup for use with ARLas.
% This version is for set up for FX study.
% Save this in the folder
%       ARLas\Peripheral\sysConfigs\
%
% Auditory Research Lab, The University of Iowa
% Deptartment of Communication Sciences & Disorders
% The University of Iowa
% Author: Shawn Goodman
% Date: July 31, 2021
% Last Updated: July 31, 2021 -- ssg -- updated for non-C-drive save locations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Make sure the following folders are specified and that they exist!
    % Data folder can be elsewhere; but everything else should be together.

    %                  '  ------------------------------------ '
    %                  '   ARLas (base directory)'
    %                  '     \\Core'
    %                  '         \\classes (* Note: The file arlas.m is located here.)'
    %                  '             \\classSupport'
    %                  '         \\general'
    %                  '         \\gui'
    %                  '         \\playrec'
    %                  '     \\Peripheral'
    %                  '         \\analysis'
    %                  '         \\calibrations'
    %                  '         \\experiments'
    %                  '         \\sysConfigs'
    %
    %
    %                  '     \\Data'
    %                  '  ------------------------------------ '


    fSep = filesep; 
    
    % create map structure
    map.calibrations = [base,'Peripheral',fSep,'calibrations',fSep];
    map.experiments = [base,'Peripheral',fSep,'experiments',fSep];
    map.sysConfigs = [base,'Peripheral',fSep,'sysConfigs',fSep];
            
    % The usual data location is here:
    map.data = [base,'Data',fSep];
    % alternative place to save data. 'D' here would be location of external drive (Lacie):
        %dataBase = ['D:',fSep,'myWork',fSep,'ARLas'];
        %map.data = [dataBase,fSep,'Data',fSep];        
            

end