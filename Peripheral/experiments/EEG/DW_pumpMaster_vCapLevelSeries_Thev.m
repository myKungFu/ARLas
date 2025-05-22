function [] = DW_pumpMaster_vCapLevelSeries_Thev(varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DW_pumpMaster_vCapLevelSeries_Thev(varargin)
%
% Control the pump for apical injection experiment while making CAP level
% series measurements. Basically, we want to be able to do the experiments of
% Lee - Lichtenhan on our DW rig
%
% Authors: Shawn S Goodman, PhD & Jeffery T Lichtenhan PhD
% Auditory Research Lab, Dept. of Comm. Sciences & Disorders, The University of Iowa
% Auditory Physiology Lab, Dept. of Otolaryngology, Washington University School of Medicine, St. Louis
% Date: January 14, 2019, Modified from DW_pumpMaster_v2_10x
% Updated: January 22, 2019 -- ssg
%                              Added some checks to allow aborting the
%                              program. (if obj.killRun, return; end)
% Updated: September 22-23, 2019 -- ssg
%                              Updated for new rates code and new userInfo data.
% Updated: September 25, 2019 -- ssg
%                                Upgraded to save version numbers to header
%                                file.
% Updated: November 11, 2019 -- ssg
%                                Updated to use Thevenin calibration; also
%                                to present combination of pips and clicks
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
obj = varargin{1};
V = '11.12.2019'; % this is the current version number of this program

    
    % ---------------------- MODIFIABLE BY USER -------------------------------                  
    runPump = 1; % set to 1 to actually run pump, 0 to skip the pump commands
    expDuration = 60; % total duration of the experiment (MINUTES)
    preInjectionMin = 6; % get a baseline of 5 minutes before starting to pump
    %--------------------------------------------------------------------------
    
    if runPump == 1
        p = openPortAndInitializePump; %evenutally need to make this function tell us if opening and initializing was successful
    end
    pumpRateSeries = getRates(preInjectionMin); % get a list of pump rates to use
    if size(pumpRateSeries,1) < expDuration
        error('requested experiment duration time (expDuration) exceeds the length of the pump rate series variable.')
    end
    
    loopDurationDesired = 60; % duration of each loop of the experiment (SECONDS)
    expectedIterations = ceil((expDuration*60) / loopDurationDesired);
    timeStart = clock; % start time of the experiment
    timeElapsed = 0; % vector of time since the start of the experiment
    counter = 1; 
    ET = 0; % elapsed time for display
    done = 0; % controls when the loop ends

    dummy = which('DW_capLevelSeries','-ALL'); % get all instances of chosen file name
    nFiles = length(dummy); % number of files with chosen name
    if nFiles > 1
        warnTxt = {'  Issue: More than one file with the name DW_capLevelSeries exists on the search path.'
             '  Action:  Highly reccommended that you remove duplicate files.'
             '  Locations: '
            };
        warnMsgARLas(warnTxt);
        for ii=1:nFiles
            disp(dummy(ii))
        end
    else
    end

    peakLocations1 = [];
    peakLocations2 = [];
    calParams = [];
    while ~done
        disp(' ')
        disp(['PUMPMASTER: Iteration ',num2str(counter),' of ',num2str(expectedIterations)])
        fprintf(' Experiment duration thus far: %d minutes\n',round(ET));
        %----------------------------------------------------------------------
        if runPump == 1
            fprintf(p,pumpRateSeries(counter,:)) % Update the pump
            fprintf(p,'G')
        end
        disp(['   Testing CAP Level Series #: ',num2str(counter)])  

        % Call the data collection function
        [peakLocations1,peakLocations2,calParams] = DW_capLevelSeries_Thev(obj,preInjectionMin,pumpRateSeries,counter,peakLocations1,peakLocations2,V,calParams);
        %----------------------------------------------------------------------

        timeElapsed(counter,1) = (etime(clock,timeStart)/60); % time elapsed from beginning of experiment (minutes)
        %ET = timeElapsed(counter,1);
        expectedET = counter * (loopDurationDesired/60);
        ET = (etime(clock,timeStart)/60); % time elapsed from beginning of experiment (minutes)
        leftover = expectedET - ET;

        if obj.killRun % detect whether an abort was called. If so, abort
            return  
        end

        if leftover > 0 % if there is still waiting time...
            pauseLen = leftover * 60; % waiting time in seconds
            disp(['   Pausing for ',num2str(pauseLen),' seconds...'])
            pause(pauseLen) 
        else
            disp('Error: data collection is taking too long!')
            %keyboard
        end
        if obj.killRun % detect whether an abort was called. If so, abort
            return  
        end

        ET = (etime(clock,timeStart)/60); % time elapsed from beginning of experiment (minutes)
        if ET >= expDuration
            done = 1;
        end
        counter = counter + 1;

        if runPump == 1
            fprintf(p,'H') % halt the pump
        end
    end
    if runPump == 1
        closePumpPort(p) 
    end
end % end of experiment file

% OLD CODE ----------------------------------------------------------------
%expectedET = counter * (loopDurationDesired/60);
