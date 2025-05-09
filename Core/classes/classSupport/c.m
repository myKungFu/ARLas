function [] = c()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% c
%
% Function for closing all open figures EXCEPT the ARLas gui.
% 
% Author: Shawn Goodman
% Date: October 4, 2018
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figHandles = get(groot, 'Children'); 
N = size(figHandles,1);
for ii=1:N
   dummy = figHandles(ii);
   name = dummy.Name;
   if ~isempty(name) % if a name has been given
       if ~strcmp(dummy.Name(1:5),'ARLas') % check to see if it is NOT the ARLas gui
           delete(dummy.Number)
       end
   else % if no name has been given
       delete(dummy.Number)
   end
end

% OLD CODE
%close all force
%while gcf~=1 delete(gcf);end;delete(gcf);