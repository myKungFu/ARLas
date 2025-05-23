function [] = myAudiometer(varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
obj = varargin{1}; % get the arlas object

% User Modifiable Parameters ----------------------------------------------

% specify which ER10X probes to use    
    probeL = 'A'; % probe that is in LEFT ear; 'A', 'B', or []
    probeR = 'B'; % probe that is in RIGHT ear; 'A', 'B', or []
                  % Note: Either one or both ears can be tested at the same time.
                  % If either probe is empty set, will not be tested. 
    subjAge = 37;
    % freqs = [];
% -------------------------------------------------------------------------                  

    subjID = obj.subjectID;
    arlas_audiometerGDM(subjID,subjAge,obj,probeL,probeR);

end

% OLD CODE ----------------------------------------------------------------
