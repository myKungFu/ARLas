function [SOS,G] = getLPiir_EB(fs)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [SOS,G] = getLPiir_EB(fs);
%
% get low-pass IIR filter coefficients, extended andwidth (_EB). 
% This version (_EB) is specifically for use in the Lichtenhan lab. It uses
% a 42 kHz cutoff. It only allows for 96 kHz sampling
% rate.
%
% filter is saved as direct-form II, second-order-sections (SOS) and gain (G). 
% filters created using Matlab's Filter Designer tool.
% Fpass = 42000; Fstop = 46000; Apass = 1; Astop = 60; 
% Using a Butterworth fitler and matching passband exactly, minimum order.
% filter orders is 7, with 4 sections stored in 6 columns.
%
% Authors: Shawn S Goodman, PhD & Jeffery T Lichtenhan PhD
% Auditory Research Lab, Dept. of Comm. Sciences & Disorders, The University of Iowa
% Auditory Physiology Lab, Dept. of Otolaryngology, Washington University School of Medicine, St. Louis
% Date: October 20, 2019
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    G = []; 
    SOS = [];
    
    if fs == 96000
      G = [0.8985
    0.7950
    0.7363
    0.8470
    1.0000];
     s1 = [1
     1
     1
     1];
     s2 = [2
     2
     2
     1];
     s3 = [1
     1
     1
     0];
     s4 = [1
     1
     1
     1];
     s5 = [1.738327383783399
   1.538160300135817
   1.424637031085486
   0.694037195298380];
    s6 = [0.855562062241978
   0.641895494027604
   0.520716092957011
                   0];
       SOS = [s1,s2,s3,s4,s5,s6];
    else
        error('The only supported sampling rate is 96000.')
    end
end