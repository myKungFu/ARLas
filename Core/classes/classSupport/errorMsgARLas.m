function [] = errorMsgARLas(errorTxt)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% errorMsgARLas(errorTxt)
%
% Customized error messages for use with ARLas.
% This function calls cprintf.m, programmed and Copyright 
% by Yair M. Altman: altmany(at)gmail.com.
%
% Expected format:
%
%         errorTxt = {'  Issue: Describe problem here.'
%              '  Action: Describe what is being done about it.'
%              '  Location: Where the error took place, function and subfunction.'
%             };
%         errorMsgARLas(errorTxt);
%
% Auditory Research Lab, The University of Iowa
% Deptartment of Communication Sciences & Disorders
% The University of Iowa
% Author: Shawn S. Goodman, PhD
% Date: October 20, 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

errorHeader = 'ARLas ERROR:';
N = size(errorTxt,1);
disp(' ')
cprintf(-[1,0,1],[errorHeader,'\n']);
for ii=1:N
    cprintf([1,0,1],[errorTxt{ii},'\n']);
end
cprintf([0,0,0],'\n');
end
