function [SOS,G] = getLPiir(fs)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [SOS,G] = getLPiir(fs);
%
% get low-pass IIR filter coefficients. 
%
% filter is saved as second-order-sections (SOS) and gain (G). 
% Filters are saved for 3 sampling rates: 44100, 48000, and 96000.
% filters created using Matlab's Filter Designer tool.
% Fpass = 20,000; Fstop = 22000; Apass = 1; Astop = 60; 
% Using a Butterworth fitler and matching stopband exactly, minimum order.
% filter orders are: 57 (96 kHz), 11 (48 kHz), and 3 (44100).
%
% Auditory Research Lab, The University of Iowa
% Deptartment of Communication Sciences & Disorders
% The University of Iowa
% Author: Shawn Goodman
% Date: October 16, 2018
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    G = []; 
    SOS = [];
    if fs == 96000
      G = [0.366595645667690
       0.348498030887723
       0.332174418047852
       0.317419737340339
       0.304058580128791
       0.291940298698328
       0.280935025159668
       0.270930419417913
       0.261828999009000
       0.253545936085770
       0.246007231598790
       0.239148195722775
       0.232912178259390
       0.227249504156104
       0.222116578200356
       0.217475129962374
       0.213291575607033
       0.209536477606638
       0.206184086914749
       0.203211954998985
       0.200600605426882
       0.198333256568337
       0.196395588510764
       0.194775548549138
       0.193463190668446
       0.192450545324948
       0.191731516591057
       0.191301804385652
       0.437217165824475
       1.000000000000000];
     s1 = [1.000000000000000
       1.000000000000000
       1.000000000000000
       1.000000000000000
       1.000000000000000
       1.000000000000000
       1.000000000000000
       1.000000000000000
       1.000000000000000
       1.000000000000000
       1.000000000000000
       1.000000000000000
       1.000000000000000
       1.000000000000000
       1.000000000000000
       1.000000000000000
       1.000000000000000
       1.000000000000000
       1.000000000000000
       1.000000000000000
       1.000000000000000
       1.000000000000000
       1.000000000000000
       1.000000000000000
       1.000000000000000
       1.000000000000000
       1.000000000000000
       1.000000000000000
       1.000000000000000];
     s2 = [   2.000000000000000
       2.000000000000000
       2.000000000000000
       2.000000000000000
       2.000000000000000
       2.000000000000000
       2.000000000000000
       2.000000000000000
       2.000000000000000
       2.000000000000000
       2.000000000000000
       2.000000000000000
       2.000000000000000
       2.000000000000000
       2.000000000000000
       2.000000000000000
       2.000000000000000
       2.000000000000000
       2.000000000000000
       2.000000000000000
       2.000000000000000
       2.000000000000000
       2.000000000000000
       2.000000000000000
       2.000000000000000
       2.000000000000000
       2.000000000000000
       2.000000000000000
       1.000000000000000];
     s3 = [ 1.000000000000000
       1.000000000000000
       1.000000000000000
       1.000000000000000
       1.000000000000000
       1.000000000000000
       1.000000000000000
       1.000000000000000
       1.000000000000000
       1.000000000000000
       1.000000000000000
       1.000000000000000
       1.000000000000000
       1.000000000000000
       1.000000000000000
       1.000000000000000
       1.000000000000000
       1.000000000000000
       1.000000000000000
       1.000000000000000
       1.000000000000000
       1.000000000000000
       1.000000000000000
       1.000000000000000
       1.000000000000000
       1.000000000000000
       1.000000000000000
       1.000000000000000
                       0];
     s4 = [ 1.000000000000000
       1.000000000000000
       1.000000000000000
       1.000000000000000
       1.000000000000000
       1.000000000000000
       1.000000000000000
       1.000000000000000
       1.000000000000000
       1.000000000000000
       1.000000000000000
       1.000000000000000
       1.000000000000000
       1.000000000000000
       1.000000000000000
       1.000000000000000
       1.000000000000000
       1.000000000000000
       1.000000000000000
       1.000000000000000
       1.000000000000000
       1.000000000000000
       1.000000000000000
       1.000000000000000
       1.000000000000000
       1.000000000000000
       1.000000000000000
       1.000000000000000
       1.000000000000000];
     s5 = [ -0.481608120584464
      -0.457832720237396
      -0.436387881506046
      -0.417004197795273
      -0.399451229314789
      -0.383531065468300
      -0.369073094763771
      -0.355929732518702
      -0.343972912983102
      -0.333091195171439
      -0.323187364227121
      -0.314176435111331
      -0.305983984693643
      -0.298544753309978
      -0.291801468571436
      -0.285703853422153
      -0.280207787731658
      -0.275274598502727
      -0.270870458410870
      -0.266965876119739
      -0.263535264833237
      -0.260556578001090
      -0.258011003108008
      -0.255882706139917
      -0.254158620707033
      -0.252828276971464
      -0.251883666523229
      -0.251319140211792
      -0.125565668351050];
    s6 = [ 0.947990703255224
       0.851824843788288
       0.765085553697455
       0.686683147156629
       0.615685549829951
       0.551292260261611
       0.492813195402444
       0.439651410190355
       0.391288909019103
       0.347274939514517
       0.307216290622282
       0.270769218002433
       0.237632697731205
       0.207542769934396
       0.180267781372859
       0.155604373271651
       0.133374090159791
       0.113420508929278
       0.095606806069865
       0.079813696115680
       0.065937686540765
       0.053889604274438
       0.043593357151063
       0.034984900336470
       0.028011383380817
       0.022630458271254
       0.018809732887457
       0.016526357754398
                       0];
       SOS = [s1,s2,s3,s4,s5,s6];
    elseif fs == 48000
        G = [0.884115938346670
           0.790003589306179
           0.722566286746239
           0.677543695777765
           0.651781204401333
           0.802122266303265
           1.000000000000000];
        s1 = [1.000000000000000
           1.000000000000000
           1.000000000000000
           1.000000000000000
           1.000000000000000
           1.000000000000000]
        s2 = [2.000000000000000
           2.000000000000000
           2.000000000000000
           2.000000000000000
           2.000000000000000
           1.000000000000000];
        s3 = [1.000000000000000
           1.000000000000000
           1.000000000000000
           1.000000000000000
           1.000000000000000
                           0];
        s4 = [1.000000000000000
           1.000000000000000
           1.000000000000000
           1.000000000000000
           1.000000000000000
           1.000000000000000];
        s5 = [1.660622051331147
           1.483852200974787
           1.357185700738689
           1.272620425838064
           1.224231114638748
           0.604244532606530];
        s6 = [0.875841702055535
           0.676162156249930
           0.533079446246268
           0.437554357272998
           0.382893702966584
                           0];
        SOS = [s1,s2,s3,s4,s5,s6];    
    elseif fs == 44100
        G = [0.964424474001281
           0.965605975484490
           1.000000000000000];
        s1 = [1.000000000000000
           1.000000000000000];
        s2 = [2.000000000000000
           1.000000000000000];
        s3 = [1.000000000000000
                           0];
        s4 = [1.000000000000000
           1.000000000000000];
        s5 = [1.926401776975857
           0.931211950968979];
        s6 = [0.931296119029265
                           0];
        SOS = [s1,s2,s3,s4,s5,s6];
    else
        error('The only supported sampling rates are 44100, 48000, and 96000.')
    end
end