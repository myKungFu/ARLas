function [multiplier,dBchange] = LTC_tipSizeLvlAdjust(LTC,whichColor,whichSize)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [multiplier,dBchange] = LTC_tipSizeLvlAdjust(LTC,whichColor,whichSize)
%
% Based on the tip size chosen (color and size), get the level adjustment
% re: the long tube diameter that was used to calibrate.
% 
% This is for the sanibel rubber "mushroom" tips that are provided with the
% ER10X probes. The tips cover a nominal range of 8-15 mm (excluding the
% pediatric "christmas tree" type tips). If no adjustment is made, the
% stimulus levels will range over 12 dB. Best practice is to calibrate
% using a middle tube size (around 8-9 mm) and then adjust using this
% function. In theory, all stimului would be within +/- 6 dB of desired, at
% worst (plus standing wave resonances). 
%
% LTC = long tube calibration structure
% whichColor = string indicating tip color: 'blue', 'red','green','yellow'
% whichSize = string indicating tip size: 'big', 'small'
%
%                  DIAM       DIAM
% COLOR    SIZE    NOMINAL    ACTUAL    
% blue     big     15         12.90         
% red      big     14         12.04
% green    big     13         11.18
% yellow   big     12         10.32
% blue     small   11          9.46
% red      small   10          8.60
% green    small    9          7.74
% yellow   small    8          6.88
%
%
% Author: Shawn Goodman, PhD
% Date: September 13, 2021
% Last Updated: September 13, 2021 -- ssg
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    multiplier = [];
    dBchange = [];
    
    % Determine the tip diameter, given color and size of tip.
    % NOTE: This is for the Sanibel tips that come with the ER10X probes.
    % Sanibel Supply: 855-278-4432
    % Contact@sanibelsupply.com
    % Note2: The cavityDiameter value shown is the nominal cap size of the
    % mushroom-shaped tip. The actual diamter size that will be used with
    % each tip is approximately 86% of the nominal value. This adjustment
    % is made below.
    if strcmp(whichColor,'blue') % BLUE -----
        if strcmp(whichSize,'large')
            cavityDiameter = 15;
        elseif strcmp(whichSize,'small')
            cavityDiameter = 11;
        else
            error('Unrecognized tip size.')
        end
    elseif strcmp(whichColor,'red') % RED -----
        if strcmp(whichSize,'large')
            cavityDiameter = 14;
        elseif strcmp(whichSize,'small')
            cavityDiameter = 10;
        else
            error('Unrecognized tip size.')
        end
    elseif strcmp(whichColor,'green') % GREEN ----- 
        if strcmp(whichSize,'large')
            cavityDiameter = 13;
        elseif strcmp(whichSize,'small')
            cavityDiameter = 9;
        else
            error('Unrecognized tip size.')
        end
    elseif strcmp(whichColor,'yellow') % YELLOW -----
        if strcmp(whichSize,'large')
            cavityDiameter = 12;
        elseif strcmp(whichSize,'small')
            cavityDiameter = 8;
        else
            error('Unrecognized tip size.')
        end
    else
       error('Unrecognized probe tip color.') 
    end
    cavityDiameter = cavityDiameter / 10; % convert to cm
    cavityDiameter = cavityDiameter * 0.86; % adjust the size from nominal to actual cavity size

    % Now calculate surge impedance of cavity, given the presumed ear canal diameter:
    cavityTemperature = LTC.cavityTemperature; % use same temperature as calibration cavity
    d = cavityTemperature - 26.85; % from Keefe (1984)
    C = 3.4723e4 * (1 + 0.00166 * d); % speed of sound, adjusted for temperature, from Keefe (1984)
    rho = 1.1769e-3 * (1 - 0.00335 * d);
    r = cavityDiameter ./ 2; % radius in cm
    z0 = (rho * C) ./ (pi * r.^2); % z0 = Ro; characteristic impedance (z0) is constant & real; 

    % Adjust the output levels, given the relative surge impedances of the
    % calibration cavity versus the presumed ear canal surge impedance
    % (based on canal diameter estimate from probe tip size)
    deltaZ = z0 / LTC.z0; % ratio of surge impedance:  tip re: calibration tip
    dBchange = 20*log10(deltaZ); % db difference
    multiplier = 10^(-dBchange/20); % amount to muliply signal by to achieve desired level
    
end

