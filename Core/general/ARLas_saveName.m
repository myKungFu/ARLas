function [fileName] = ARLas_saveName(pathName,fileName)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fileName = ARLas_saveName(pathName,fileName);
%
% Check to ensure not overwriting another file with the same name
% pathName = the complete path name, ending with a backslash (\).
% fileName = the complete desired file name, ending with the extension.
%
% Auditory Research Lab, The University of Iowa
% Deptartment of Communication Sciences & Disorders
% The University of Iowa
% Author: Shawn S. Goodman, PhD
% Date: March 27, 2017
% Updated: September 25, 2019 -- ssg
%                                Added check to make sure pathName ends in
%                                a file separator.
% Updated: October 21, 2019 -- ssg
%                                Added check to strip off counter at the
%                                end, if it exists.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

sep = filesep;
if strcmp(pathName(end),sep) == 0
   pathName = [pathName,sep]; 
end

ext = fileName(end-3:end); % the file extension
fileName = fileName(1:end-4); % the file name proper

% check to see if counting tag has already been added. If so, remove.
indx = strfind(fileName,'_');
if ~isempty(indx)
    indx = max(indx);
    if indx < length(fileName)
        tag = fileName(indx+1:end);
        q = str2double(tag);
        if ~isnan(q)
            fileName = fileName(1:indx-1);
        end
    end
end

ok = 0;
counter = 1;
while ~ok
    newFileName = [fileName,'_',num2str(counter)]; % add a numeric extension
    if exist([pathName,newFileName,ext],'file')==2 % if file already exists
        counter = counter + 1;
    else
        ok = 1;
        fileName = [newFileName,ext];
    end
end
