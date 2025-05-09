function [Y] = AR_engine(Y,indx)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Y = AR_engine(Y,indx);
% This function removes buffers with artifacts from a data set.
%
% Y = data matrix. Buffers must be in columns; for example, 
%   a data set comprising 16 buffers of 1024 samples should be passed
%   to this function so that size(Y) == [1024 16].
% indx = an array of indices showing which buffers (columns) are to be deleted. 
%   This argument is returned by AR.m.
%
% Author: Shawn Goodman
% Date: August 1, 2005
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%disp('Removing artifacts...')
Y(:,indx) = [];
%disp('Artifact rejection completed.')


% Notes:
% It takes longer to do this on a row matrix than a column matrix;
% therefore, this code is set up to reject columns rather than rows. 
%
% The following code will do the same thing in a loop. It is excruciatingly slow...
% indx = flipud(indx); % start from the highest and work toward the lowest
% for ii=1:length(indx),
%    y(:,indx(ii)) = []; % take out the artifacts
% end
