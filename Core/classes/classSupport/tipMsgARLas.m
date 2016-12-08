function [] = tipMsgARLas(tipTxt)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% alertMsgARLas(tipTxt)
%
% Customized alert messages for use with ARLas.
% This function calls cprintf.m, programmed and Copyright 
% by Yair M. Altman: altmany(at)gmail.com.
%
% Expected format:
%
%         tipTxt = {'  Issue: Describe problem here.'
%              '  Suggested Action: Describe what is being done about it.'
%              '  Location: Where the problem took place, function and subfunction.'
%             };
%         tipMsgARLas(tipTxt);
%
% Auditory Research Lab, The University of Iowa
% Deptartment of Communication Sciences & Disorders
% The University of Iowa
% Author: Shawn S. Goodman, PhD
% Date: October 20, 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

errorHeader = 'ARLas TIP:';
N = size(tipTxt,1);
disp(' ')
cprintf(-[0,.8,0],[errorHeader,'\n']);
for ii=1:N
    cprintf([0,.8,0],[tipTxt{ii},'\n']);
end
cprintf([0,0,0],'\n');
end
