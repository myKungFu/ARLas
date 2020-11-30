function [] = ARLas_installUpdates()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ARLas_installUpdates
%
% This function installs and manages updates to ARLas.
% It is assumed that the user is running ARLas from the standard ARLas directory
% structure (see below). The function backs up the all the contents of the 
% Core and Peripheral folders. The Data folder is not backed up or overwritten.
% After successful backup, new content is copied to the Core and Peripheral
% folders. The function also checks for the existance of duplicate files
% and alerts the user if any exist on the Matlab search path.
%
% Backups will be kept in myWork\ARLas_oldVersions\
%   Backups will be listed with a date-of-backup extension to the name.
%   NOTE: ARLas_oldVersions should NOT be located on the Matlab search path!
%   The files in ARLas must conform to the standard ARLas directory structure
%     (see below) in order to be backed up. In other words, files and
%     folders that are not located inside of either Core or Peripheral will
%     not be backed up.
%
% Updates should be placed in myWork\ARLas_updates
%   Updates should be placed in this folder PRIOR to running this program.
%   If the ARLas_updates folder is empty, then a backup will be done, but
%     not updates will be installed.
%   The updates placed in ARLas_updates must conform to the standard ARLas
%     directory structure (see below).
%
% Standard ARLas file structure
%           \\Core
%               \\classes
%                   \\classSupport
%               \\general
%               \\gui
%               \\playrec
%           \\Peripheral
%               \\analysis
%               \\calibrations
%               \\experiments
%               \\sysConfigs
%           \\Data
%
%
% Authors: Shawn S Goodman, PhD & Jeffery T Lichtenhan PhD
% Auditory Research Lab, Dept. of Comm. Sciences & Disorders, The University of Iowa
% Auditory Physiology Lab, Dept. of Otolaryngology, Washington University School of Medicine, St. Louis
% Date: September 26, 2019
% Updated: September 26, 2019 -- ssg
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    OK = 1; % okay to continue
    sep = filesep; % file separator for the operating system being used
    
    updateDir = ['C:',sep,'myWork\ARLas_updates',sep];
    backupDir = ['C:',sep,'myWork\ARLas_oldVersions',sep];

    alertTxt = {'Starting ARLas_installUpdates.'
         '  This program will back up your existing files and replace older files with new updates.'
         '  NOTE: The DATA folder will NOT be copied, nor updated.'         
        };
    alertMsgARLas(alertTxt);

    % Make sure the backup directory exists -------------------------------
    if exist(backupDir,'dir') ~= 7 % check to make sure that the backup directory exists
        success = mkdir(backupDir); % if not, try to create it
        if success ~= 1
            error('Backup directory does not exist and could not be created. Aborting program.')
        end
    end
    
    % Create the current backup folder ------------------------------------
    formatOut = 'mm_dd_yyyy';
    tag = datestr(now,formatOut);
    backupFolder = ['ARLas_',tag];
    if exist([backupDir,backupFolder],'dir') == 7 % if folder already exists
        counter = 1;
        maxCounts = 20; % maximum number of updates allowed on a given date
        done = 0;
        while done == 0
            backupFolder = ['ARLas_',tag,'_Copy',num2str(counter)];
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
    
    % Copy ARLas, except for Data, to the current backup folder -----------
    base = 'C:\myWork\ARLas\';
    D = dir(base);
    D = D(~cellfun('isempty', {D.date}));
    N = size(D,1);
    for ii=1:N
        if strcmp(D(ii).name,'.') % skip this
        elseif strcmp(D(ii).name,'..') % skip this
        elseif strcmp(D(ii).name,'Data') % skip this
        elseif strcmp(D(ii).name,'Core') % copy this
            status = copyfile([base,D(ii).name],[backupDir,backupFolder,D(ii).name]);
            if status ~= 1
                errorTxt = {'Backup copy of Core folder unsuccessful.'
                     '  Required Action: Find the source of this error and fix it.'
                     '  ARLas_installUpdates aborted.'
                    };
                errorMsgARLas(errorTxt);
                return
            else
                alertTxt = {'Core folder successfully backed up.'
                     ['  Contents of Core copied to ','C:\\myWork\\ARLas_oldVersions',sep,sep,backupFolder(1:end-1),'.']
                    };
                alertMsgARLas(alertTxt);
            end
        elseif strcmp(D(ii).name,'Peripheral') % copy this
            status = copyfile([base,D(ii).name],[backupDir,backupFolder,D(ii).name]);
            if status ~= 1
                errorTxt = {'Backup copy of Peripheral folder unsuccessful.'
                     '  Required Action: Find the source of this error and fix it.'
                     '  ARLas_installUpdates aborted.'
                    };
                errorMsgARLas(errorTxt);
                return
            else
                alertTxt = {'Peripheral folder successfully backed up.'
                     ['  Contents of Peripheral copied to ','C:\\myWork\\ARLas_oldVersions',sep,sep,backupFolder(1:end-1),'.']
                    };
                alertMsgARLas(alertTxt);
            end
        else
            warnTxt = {'ARLas_installUpdates has detected files or folders outside of the expected ARLas directory structure.'
                 '  Required Action:  Move all files and sub-folders into one of the following:'
                 '           ...\\ARLas\\Core\\'
                 '           ...\\ARLas\\Data\\'
                 '           ...\\ARLas\\Peripheral\\'
                 '  ARLas_installUpdates will NOT continue until this issue is resolved.'
                };
            warnMsgARLas(warnTxt);
            OK = 0;
        end
    end     
    
    % Copy updates to the ARLas folder ------------------------------------
    if OK == 0 % only continue if the backup was successful
        return
    else
        % Make sure the update directory exists ---------------------------
        if exist(updateDir,'dir') ~= 7 % check to make sure that the backup directory exists
            success = mkdir(updateDir); % if not, try to create it
            if success ~= 1
                error('Update directory does not exist and could not be created. Aborting program.')
            end
        end
        D = dir(updateDir);
        D = D(~cellfun('isempty', {D.date}));
        N = size(D,1);
        if N <3 % if the update folder is empty
            warnTxt = {'The updates folder (C:\\myWork\\ARLas_updates\\) is empty.'
                 '  Existing files backed up, but NO UPDATES WERE MADE.'
                 '  Required Action:  Place desired updates into the updates folder.'
                 '  Note:  All files to be updated must be located inside one of the following locations:'
                 '           ...\\myWork\\ARLas_updates\\Core\\'
                 '           ...\\myWork\\ARLas_updates\\Peripheral\\'                         
                };
            warnMsgARLas(warnTxt);
            OK = 0;
        else
            for ii=1:N
                if strcmp(D(ii).name,'.') % skip this
                elseif strcmp(D(ii).name,'..') % skip this
                elseif strcmp(D(ii).name,'Data') % skip this
                elseif strcmp(D(ii).name,'Core') % copy this
                    status = copyfile([updateDir,D(ii).name],[base,'Core']);
                    if status ~= 1
                        errorTxt = {'Update of Core folder unsuccessful.'
                             '  Required Action: Find the source of this error and fix it.'
                             '  ARLas_installUpdates aborted.'
                            };
                        errorMsgARLas(errorTxt);
                        return
                    else
                        alertTxt = {'Core folder successfully updated.'
                            };
                        alertMsgARLas(alertTxt);
                    end
                elseif strcmp(D(ii).name,'Peripheral') % copy this
                    status = copyfile([updateDir,D(ii).name],[base,'Peripheral']);
                    if status ~= 1
                        errorTxt = {'Update of Peripheral folder unsuccessful.'
                             '  Required Action: Find the source of this error and fix it.'
                             '  ARLas_installUpdates aborted.'
                            };
                        errorMsgARLas(errorTxt);
                        return
                    else
                        alertTxt = {'Peripheral folder successfully backed up.'
                            };
                        alertMsgARLas(alertTxt);
                    end
                else
                    warnTxt = {'ARLas_installUpdates has detected files or folders outside of the expected ARLas directory structure.'
                         '  Required Action:  Move all updates into one of the following:'
                         '           ...\\myWork\\ARLas_updates\\Core\\'
                         '           ...\\myWork\\ARLas_updates\\Peripheral\\'                         
                         '  ARLas_installUpdates will NOT continue until this issue is resolved.'
                        };
                    warnMsgARLas(warnTxt);
                    OK = 0;
                end
            end
        end
    end
    
    % Check for any duplicate files ---------------------------------------
    disp('Checking for duplicate files....')
    D = dir('C:\myWork\ARLas\**\*.m');
    N = size(D,1);
    hit = 0;
    for ii=1:N
        dummy = which(D(ii).name,'-ALL'); % get all instances of chosen file name
        nFiles = length(dummy); % number of files with chosen name
        if nFiles > 1
            hit = 1;
            warnTxt = {'  Issue: More than one file with the name DW_simultaneousANOWSF exists on the search path.'
                 '  Action:  Highly reccommended that you remove duplicate files.'
                 '  Locations: '
                };
            warnMsgARLas(warnTxt);
            for jj=1:nFiles
                disp(dummy(jj))
            end
        else
        end
    end
    if hit == 0
        disp('No duplicate m-files found on the Matlab search path.')
    end
    
    if OK == 1
        alertTxt = {'Successfully finished ARLas_installUpdates!'
            '  VERY IMPORTANT: Either delete the contents of the updates folder (C:\\myWork\\ARLas_updates\\)'
            '                  or make sure that the updates folder is NOT on the Matlab search path.'
            };
        alertMsgARLas(alertTxt);
    else
        warnTxt = {'Fnished ARLas_installUpdates, but with warnings or errors.'
            };
        warnMsgARLas(warnTxt);
    end

end


