function [pathName,folderName,fileName] = mostRecentCalibration(calType,probeLabel,probeCh,map)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [pathName,folderName,fileName] = mostRecentCalibration(calType,probeLabel,probeCh,map);
%
% Find the location of the most recent calibration file.
% 
% Input Arguments:
%   calType = A string specifying the type of calibration files desired:
%                 'thev' (for thevenin source) or 'mic' (for microphone)
%                 location is fixed based on the type--
%                      'C:\myWork\ARLas\Peripheral\calibrations\thevCals\'
%                      'C:\myWork\ARLas\Peripheral\calibrations\micCals\'
%   probeLabel = A string specifying the probe being used. This is obtained
%                from [inputs,outputs] = hardwareSetup, and is either
%                inputs.label or outputs.label.
%   probeCh = 1 or 2. This value specifies the output channel of the
%                probe. The value is unused if calType is 'mic', because
%                there is only one microphone per probe.
%
% Output Arguments:
%   pathName = A string specifyin the local directory where desired calibration
%                 file is located.
%   folderName = A string specifying the folder within calPathName
%                 containing the most recent calibration file. Calibration 
%                 folders are all of the form 'dd_mm_yyyy'.
%   fileName = A string specifying the name of the file that is the most 
%                 recent calibration.
%
% Authors: Shawn S Goodman, PhD & Jeffery T Lichtenhan PhD
% Auditory Research Lab, Dept. of Comm. Sciences & Disorders, The University of Iowa
% Auditory Physiology Lab, Dept. of Otolaryngology, Washington University School of Medicine, St. Louis
% Date: October 10, 2019
% Last Updated: October 10, 2019 -- ssg
% Last Updated: March 23, 2022 -- ssg -- added map as an input argument.
%                                   This allows users to map MAC and other
%                                   non-standard directories for ARLas.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if nargin <4 % for backwards compatibility; 3/23/2022 -- ssg
        map = [];
    end
    folderName = [];
    fileName = [];

    fsep = filesep; % get file separator
    
    % central location of the calibration files ---------------------------
    if isempty(map)
        if strcmp(calType,'thev')
            pathName = ['C:',fsep,'myWork',fsep,'ARLas',fsep,'Peripheral',fsep,'calibrations',fsep,'thevCals',fsep];
            tag = [probeLabel,'_Ch',num2str(probeCh),'_'];
        elseif strcmp(calType,'mic')
            pathName = ['C:',fsep,'myWork',fsep,'ARLas',fsep,'Peripheral',fsep,'calibrations',fsep,'micCals',fsep];
            tag = [probeLabel,'_']; % mic cals have no channel designator, since there is only one mic per probe
        else
            warning('No calibration files of the type specified exist.')
            return
        end
    else % if map contains the structure
        if strcmp(calType,'thev')
            pathName = [map.calibrations,'thevCals',fsep];
            tag = [probeLabel,'_Ch',num2str(probeCh),'_'];
        elseif strcmp(calType,'mic')
            pathName = [map.calibrations,'micCals',fsep]; 
            tag = [probeLabel,'_']; % mic cals have no channel designator, since there is only one mic per probe
        else
            warning('No calibration files of the type specified exist.')
            return
        end
    end
    % find and use the folder whose name is the most recent date
    d = dir(pathName);
    if size(d,1) < 3 % no calibration folders exist
        return
    end
    str = computer;
%     if strcmp(str,'MACI64')
%        jar = datestr(d(4).name);
%     else 
%         jar = datestr(d(3).name);
%     end 

    % find where the structure starts
    counter = 3;
    done = 0;
    while done == 0
        if counter <= size(d,1)
             try jar = datestr(d(counter).name);
                 done = 1;
             catch
                 counter = counter + 1;
             end
        else
            done = 1;
            keyboard % got stuck! fix me
        end
    end

%     try jar = datestr(d(3).name);
%     catch ME1
%         try jar = datestr(d(4).name);
%         catch ME2            
%             keyboard
%         end
%     end 
    for ii=counter:size(d,1) % added 20 June 2023-ssg- RB setup
        candidate = datestr(d(ii).name);
        if datetime(candidate) > datetime(jar) % if new folder has a more recent date
            jar = candidate; % replace curent contents with the more recent folder
        end
    end
    folderName = [datestr(jar,'mm_dd_yyyy'),fsep]; % this is the folder with the most recent date as its name

    d = dir([pathName,folderName,tag,'*.mat']); % list of files matching the probe label
    nFiles = size(d,1);
    if nFiles == 0
        return
    end
    jar = datestr(d(1).datenum);
    jarIndx = 1;
    if nFiles > 1
        for ii=2:size(d,1)
            candidate = datestr(d(ii).datenum);
            if datetime(candidate) > datetime(jar) % if new folder has a more recent date
                jar = candidate; % replace curent contents with the more recent folder
                jarIndx = ii;
            end
        end
    end
    fileName = d(jarIndx).name;
end