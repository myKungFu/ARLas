function [] = plt(varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plt(varargin)
%
% Function for plotting in the current figure EXCEPT the ARLas gui.
% If no figure currently exists, will generate a new figure.
% This code is a wrapper around the Matlab plot function.
% 
% Author: Shawn Goodman
% Date: October 16, 2018
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%test
figHandles = get(groot, 'Children'); 
N = size(figHandles,1);
didPlot = 0;
for ii=1:N
   if didPlot == 1
       break
   end
   dummy = figHandles(ii);
   name = dummy.Name;
   if ~isempty(name) % if a name has been given
       if ~strcmp(dummy.Name(1:5),'ARLas') % check to see if it is NOT the ARLas gui
           h = figHandles(ii);
           figure(h); % make current figure
           plot(varargin{:}); % call the usual plot function
           didPlot = 1;
       end
   else % if no name has been given go ahead and plot in this figure
       h = figHandles(ii);
       figure(h); % make current figure
       plot(varargin{:}); % call the usual plot function
       didPlot = 1;
   end
end
if didPlot == 0
    h = figure;
    plot(varargin{:});
end
end