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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ext = fileName(end-3:end); % the file extension
fileName = fileName(1:end-4); % the file name proper
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
