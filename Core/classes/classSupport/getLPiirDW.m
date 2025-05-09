function [SOS,G] = getLPiirDW(fs)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [SOS,G] = getLPiirDW(fs);
%
% get low-pass IIR filter coefficients. 
% This version (DW) is specifically for use in the Lichtenhan lab. It uses
% a 32 kHz cutoff, instead of 20,000 Hz. It only allows for 96 kHz sampling
% rate.
%
% filter is saved as direct-form II, second-order-sections (SOS) and gain (G). 
% filters created using Matlab's Filter Designer tool.
% Fpass = 32000; Fstop = 35200; Apass = 1; Astop = 60; 
% Using a Butterworth fitler and matching stopband exactly, minimum order.
% filter orders is 30, with 15 sections stored in 6 columns.
%
% Authors: Shawn S Goodman, PhD & Jeffery T Lichtenhan PhD
% Auditory Research Lab, Dept. of Comm. Sciences & Disorders, The University of Iowa
% Auditory Physiology Lab, Dept. of Otolaryngology, Washington University School of Medicine, St. Louis
% Date: January 8, 2019
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    G = []; 
    SOS = [];
    if fs == 96000
      G = [0.7284
    0.6713
    0.6233
    0.5828
    0.5485
    0.5196
    0.4951
    0.4746
    0.4576
    0.4436
    0.4323
    0.4236
    0.4172
    0.4130
    0.4109
    1.0000];
     s1 = [1
     1
     1
     1
     1
     1
     1
     1
     1
     1
     1
     1
     1
     1
     1];
     s2 = [2
     2
     2
     2
     2
     2
     2
     2
     2
     2
     2
     2
     2
     2
     2];
     s3 = [1
     1
     1
     1
     1
     1
     1
     1
     1
     1
     1
     1
     1
     1
     1];
     s4 = [1
     1
     1
     1
     1
     1
     1
     1
     1
     1
     1
     1
     1
     1
     1];
     s5 = [0.999138296231463
   0.920861707866799
   0.854981769669545
   0.799376963751391
   0.752375237874295
   0.712652388734952
   0.679155671884157
   0.651046235200492
   0.627655618710505
   0.608452819770372
   0.593019359640163
   0.581030479808731
   0.572241107675511
   0.566475611344022
   0.563620650216120];
    s6 = [0.914528145731043
   0.764536165500476
   0.638298390015945
   0.531749610562939
   0.441685875213144
   0.365569772975099
   0.301384056701069
   0.247521335299533
   0.202700719589987
   0.165904713284988
   0.136331436075265
   0.113358592416246
   0.096516579257685
   0.085468854390788
   0.079998236904747];
       SOS = [s1,s2,s3,s4,s5,s6];
    else
        error('The only supported sampling rate is 96000.')
    end
end