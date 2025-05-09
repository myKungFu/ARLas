function [] = findDuplicates()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% findDuplicates
%
% Function to search for and list any duplicate files in the ARLas
% directories and sub-directories. List of any duplictes will be printed 
% to the workspace.
%
% Searches the directory (and sub-directories) specified by the variable
% "basePath" in the first line of code.
%
% Author: Shawn Goodman
% Date: August 9, 2024
% Last Update: August 9, 2024 -- ssg
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    basePath = 'C:\myWork\'; % directory in which to search for duplicates
    
    origPath = pwd; % current directory
    cd(basePath)    % change directories to the search directory
    d = dir(['**',filesep,'*.m']); % get directory structure
    nFilesTotal = size(d,1); % total number of files in the search directory
    cd(origPath)    % change back to the original directory
    
    filenames = cell([],1) ; % initialize list of file names
    for ii=1:nFilesTotal % loop to create a vector of filenames
        filenames(ii) = cell(cellstr(d(ii).name));
    end
    
    [C,IA,IC] = unique(filenames,'stable'); % find a list of uniquie file names
    
    for ii=1:length(IC) % loop to find the duplicate files
        dummy = IC;
        dummy(ii) = [];
        hits(ii) = ismember(IC(ii),dummy);
    end
    indx = find(hits); % locations of all duplicates
    

    % concentrate only on duplicates --------------
    if isempty(indx)
        nUnique = 0;
        nDuplications = 0;
    else
        d = d(indx);
        
        nFiles = size(d,1); % number of files with duplicates
        filenames = cell([],1); % initialize
        for ii=1:nFiles % loop to get vector of file names
            filenames(ii) = cell(cellstr(d(ii).name));
        end
        [C,IA,IC] = unique(filenames); % find indices of unique file names
        nUnique = length(IA); % number of unique file names
        nDuplications = length(IC) - length(IA); % number if duplications
        files = cell([],1) ; % initialize
        counter = 1;
        for ii=1:nUnique
            findme = IC(IA(ii));
            dups = find(IC==findme); % vector of duplications
            for jj=1:length(dups)
                files(counter) = cell(cellstr([d(dups(jj)).folder,filesep,d(dups(jj)).name]));
                counter = counter + 1;
            end
        end
    end

    % display results on the screen
    disp(' ')
    disp('-------------------------------------------------------------------')
    disp(['Searched ',basePath,' and all sub-directories.'])
    disp(['Found a total of ',num2str(nFilesTotal),' m-files.'])
    disp(['Found ',num2str(nUnique),' m-files with duplicate copies.'])
    disp(['There are ',num2str(nDuplications),' duplicated m-files.'])
    if nUnique > 0
        disp('Examine the following list and remove the duplicates.')
        disp('IMPORTANT: "duplicate" means only that the file names are identical.')
        disp('            The program codes are not necessarily identical!')
        disp('            Be careful about what you delete and what you rename/move.')
        disp(' ')
        disp(files')
    end
    disp(' ')
    disp('-------------------------------------------------------------------')
    disp(' ')
    
end