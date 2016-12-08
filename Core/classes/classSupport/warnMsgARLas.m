function [] = warnMsgARLas(warnTxt)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% warnMsgARLas(warnTxt)
%
% Customized warning messages for use with ARLas.
% This function calls cprintf.m, programmed and Copyright 
% by Yair M. Altman: altmany(at)gmail.com.
%
% Expected format:
%
%         errorTxt = {'  Issue: Describe problem here.'
%              '  Suggested Action: Describe what is being done about it.'
%              '  Location: Where the problem took place, function and subfunction.'
%             };
%         warnMsgARLas(warnTxt);
%
% Auditory Research Lab, The University of Iowa
% Deptartment of Communication Sciences & Disorders
% The University of Iowa
% Author: Shawn S. Goodman, PhD
% Date: October 20, 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

errorHeader = 'ARLas WARNING:';
N = size(warnTxt,1);
disp(' ')
cprintf(-[1,0.5,0],[errorHeader,'\n']);
for ii=1:N
    cprintf([1,0.5,0],[warnTxt{ii},'\n']);
end
cprintf([0,0,0],'\n');
end
