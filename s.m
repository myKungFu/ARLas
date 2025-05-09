function [] = s()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% s
% Call this function first when starting MATLAB.
% It sets formatting and adds paths to all directories
% under C:\myWork\ARLas\.
%
% Author: Shawn Goodman
% Date: April 25, 2013
% Date: April 4, 2025 -- ssg -- updated to work with myWork2 instead of
% myWork.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

format compact
format short
addpath(genpath(['C:\myWork2\ARLas'])) % add all directories and subdirectories in myWork
% pathName = (['C:\myWork\ARLas\development']);
% if exist(pathName,'dir'), % remove attic and subdirectories from the path
%     p = myRemovePath(pathName);
% end
% rmpath(pathName)  
cd('C:\myWork2\ARLas')


function p = myRemovePath(d)
% Remove toolbox paths recursively
% This function is a modified form of genpath.m
methodsep = '@';  % qualifier for overloaded method directories
p = '';           % path to be returned
% Generate path based on given root directory
files = dir(d);
if isempty(files)
  return
end
% Add d to the path even if it is empty.
p = [p d pathsep];
% set logical vector for subdirectory entries in d
isdir = logical(cat(1,files.isdir));
% Recursively descend through directories which are neither
% private nor "class" directories.
dirs = files(isdir); % select only directory entries from the current listing
for i=1:length(dirs)
   dirname = dirs(i).name;
   if    ~strcmp( dirname,'.')          && ...
         ~strcmp( dirname,'..')         && ...
         ~strncmp( dirname,methodsep,1) && ...
         ~strcmp( dirname,'private')
         rmpath(fullfile(d,dirname))
      p = [p myRemovePath(fullfile(d,dirname))]; % recursive calling of this function.
   end
end