function [] = ARLas_calibrationBackup()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ARLas_calibrationBackup
%
% This function backs up old calibration files, moving them to the backup
% location. It is to be run prior to obtaining a new set of calibrations.
% It is assumed that the user is running ARLas from the standard ARLas directory
% structure (see below). The function backs up the all the contents of the 
% Peripheral\calibration\ folder. 
%
% Backups will be kept in myWork\ARLas_oldCalibrations\
%   Backups will be listed with a date-of-backup extension to the name.
%   NOTE: ARLas_oldVersions should NOT be located on the Matlab search path!
%   The files in ARLas must conform to the standard ARLas directory structure
%     (see below) in order to be backed up. In other words, calibration files and
%     folders that are not located inside of Peripheral\calibrations will
%     not be backed up.
%
% Standard ARLas file structure
%           \\Core
%           \\Peripheral
%               \\calibrations
%                   \\micCals
%                   \\thevCals
%           \\Data
%
%
% Authors: Shawn S Goodman, PhD & Jeffery T Lichtenhan PhD
% Auditory Research Lab, Dept. of Comm. Sciences & Disorders, The University of Iowa
% Auditory Physiology Lab, Dept. of Otolaryngology, Washington University School of Medicine, St. Louis
% Date: October 13, 2019
% Updated: October 13, 2019 -- ssg
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    OK = 1; % okay to continue
    sep = filesep; % file separator for the operating system being used
    backupDir = ['C:',sep,'myWork\ARLas_oldCalibrations',sep]; % backup location
    base = ['C:',sep,'myWork',sep,'ARLas',sep,'Peripheral',sep,'calibrations',sep]; % location of calibrations

    alertTxt = {'Running ARLas_calibrationBackup.'
         '  This program will back up your existing calibrations by moving '
         '  the contents of C:\\myWork\\ARLas\\Peripheral\\calibrations'
         '  to C:\\myWork\\ARLas_oldCalibrations.'
         '  After running this program, you will need to obtain new calibrations.'         
        };
    nn = size(alertTxt,1);
    for ii=1:nn
        cprintf([0,0,.4],[alertTxt{ii},'\n']);
    end
    disp(' ')
    
    % Make sure the backup directory exists -------------------------------
    if exist(backupDir,'dir') ~= 7 % check to make sure that the backup directory exists
        success = mkdir(backupDir); % if not, try to create it
        if success ~= 1
            error('Backup directory does not exist and could not be created. Aborting program.')
        end
        addpath(genpath(backupDir))
    end
    
    % Create the current backup folder ------------------------------------
    formatOut = 'mm_dd_yyyy';
    tag = datestr(now,formatOut);
    backupFolder = ['cal_',tag];
    if exist([backupDir,backupFolder],'dir') == 7 % if folder already exists make a new one with a counter on it
        counter = 1;
        maxCounts = 20; % maximum number of updates allowed on a given date
        done = 0;
        while done == 0
            backupFolder = [tag,'_Copy',num2str(counter)];
            if exist([backupDir,backupFolder],'dir') == 7 % if folder already exists
                counter = counter + 1;
            else
                done = 1;
            end
            if counter > maxCounts
                done = 1;
                error('Maximum number of copies exceeded.')
            end
        end
    end
    backupFolder = [backupFolder,sep];
    success = mkdir([backupDir,backupFolder]); % create backup folder
    if success ~= 1
        error('Backup folder does not exist and could not be created. Aborting program.')
    end
    addpath(genpath([backupDir,backupFolder]))
    
    % Move calibrations to the current backup folder -----------
    D = dir(base);
    if isempty(D) % if calibration folder does not exist in ARLas
        success = mkdir(base); % create it, along witht the two sub-folders (micCals and thevCals)
        if success ~= 1
            error('calibrations folder does not exist and could not be created. Aborting program.')
        end
        addpath(genpath(base))
        
        success = mkdir([base,'micCals']); % create backup folder
        if success ~= 1
            error('micCals folder does not exist and could not be created. Aborting program.')
        end
        addpath(genpath([base,'micCals']))
        
        success = mkdir([base,'thevCals']); % create backup folder
        if success ~= 1
            error('thevCals folder does not exist and could not be created. Aborting program.')
        end
        addpath(genpath([base,'thevCals']))

        D = dir(base);
    end
    
    D = D(~cellfun('isempty', {D.date}));
    N = size(D,1);
    for ii=1:N
        if strcmp(D(ii).name,'.') % skip this
        elseif strcmp(D(ii).name,'..') % skip this
        elseif strcmp(D(ii).name,'micCals') % move this
            status = copyfile([base,D(ii).name],[backupDir,backupFolder,D(ii).name]);
            if status ~= 1
                errorTxt = {'Backup copy of micCals folder unsuccessful.'
                     '  Required Action: Find the source of this error and fix it.'
                     '  ARLas_calibrationBackup aborted.'
                    };
                nn = size(alertTxt,1);
                for ii=1:nn
                    cprintf([1,0,0],[alertTxt{ii},'\n']);
                end
                disp(' ')
                return
            else
                alertTxt = {'micCals folder successfully backed up.'
                     ['  Contents of micCals copied to ','C:\\myWork\\ARLas_oldCalibrations',sep,sep,backupFolder(1:end-1),'.']
                    };
                nn = size(alertTxt,1);
                for ii=1:nn
                    cprintf([0,0,.4],[alertTxt{ii},'\n']);
                end
                disp(' ')
            end
        elseif strcmp(D(ii).name,'thevCals') % copy this
            status = copyfile([base,D(ii).name],[backupDir,backupFolder,D(ii).name]);
            if status ~= 1
                errorTxt = {'Backup copy of thevCals folder unsuccessful.'
                     '  Required Action: Find the source of this error and fix it.'
                     '  ARLas_calibrationBackup aborted.'
                    };
                nn = size(alertTxt,1);
                for ii=1:nn
                    cprintf([1,0,0],[alertTxt{ii},'\n']);
                end
                disp(' ')
                return
            else
                alertTxt = {'thevCals folder successfully backed up.'
                     ['  Contents of thevCals copied to ','C:\\myWork\\ARLas_oldCalibrations',sep,sep,backupFolder(1:end-1),'.']
                    };
                nn = size(alertTxt,1);
                for ii=1:nn
                    cprintf([0,.0,.4],[alertTxt{ii},'\n']);
                end
                disp(' ')
            end
        else
            warnTxt = {'ARLas_calibrationBackup has detected files or folders outside of the expected ARLas directory structure.'
                 '  Required Action:  Move all calibration files and sub-folders into one of the following:'
                 '           ...\\ARLas\\Peripheral\\calibrations\\micCals'
                 '           ...\\ARLas\\Peripheral\\calibrations\\thevCals'                 
                 '  ARLas_calibrationBackup has NOT successfully run.'
                };
            nn = size(alertTxt,1);
            for ii=1:nn
                cprintf([0,0,1],[alertTxt{ii},'\n']);
            end
            disp(' ')
            OK = 0;
        end
    end
    
    
    if OK == 1 % if successfully backed up, delete old contents -----------
        D = dir(base);
        D = D(~cellfun('isempty', {D.date}));
        N = size(D,1);        
        for ii=1:N
            if strcmp(D(ii).name,'.') % skip this
            elseif strcmp(D(ii).name,'..') % skip this
            elseif strcmp(D(ii).name,'micCals') % delete contents of this
                d = dir([base,D(ii).name,sep]);
                n = size(d,1);
                for jj=1:n
                    if strcmp(d(jj).name,'.') % skip this
                    elseif strcmp(d(jj).name,'..') % skip this
                    else
                        dd = dir([base,D(ii).name,sep,d(jj).name,sep,'*']);
                        nn = size(dd,1);
                        for kk=1:nn
                            if strcmp(dd(kk).name,'.') % skip this
                            elseif strcmp(dd(kk).name,'..') % skip this
                            else
                                delete([base,D(ii).name,sep,d(jj).name,sep,dd(kk).name])
                            end
                        end
                        warning off
                        rmdir([base,D(ii).name,sep,d(jj).name]);
                        warning on
                    end
                end
            elseif strcmp(D(ii).name,'thevCals') % delete contents of this
                d = dir([base,D(ii).name,sep]);
                n = size(d,1);
                for jj=1:n
                    if strcmp(d(jj).name,'.') % skip this
                    elseif strcmp(d(jj).name,'..') % skip this
                    else
                        dd = dir([base,D(ii).name,sep,d(jj).name,sep,'*']);
                        nn = size(dd,1);
                        for kk=1:nn
                            if strcmp(dd(kk).name,'.') % skip this
                            elseif strcmp(dd(kk).name,'..') % skip this
                            else
                                delete([base,D(ii).name,sep,d(jj).name,sep,dd(kk).name])
                            end
                        end
                        warning off
                        rmdir([base,D(ii).name,sep,d(jj).name]);
                        warning on
                    end
                end
            else
            end
        end
    end
end

    % 
% 
% load('C:\myWork\ARLas\Data\testV6\testV6\ER10xA_Ch1.mat')
% load('C:\myWork\ARLas\Peripheral\calibrations\micCals\micCal_10xA.mat')
% Q = fastFilter(b,recordings);
% t = thevCal_new(recordingParams,stimulus,Q);
% %t = thevCal_new(recordingParams,stimulus,recordings);
% 
% t.fmin = 200;
% t.fmax = 10000;
% t.calculate
% 
% figure(110)
% plot(t.freq,abs(t.PS),'b')
% hold on
% 
% t.fmin = 8000;
% t.fmax = 20000;
% t.calculate
% 
% figure(110)
% plot(t.freq,abs(t.PS),'r')
% 
% t.fmin = 18000;
% t.fmax = 32000;
% t.calculate
% 
% figure(110)
% plot(t.freq,abs(t.PS),'g')
% 
% 
