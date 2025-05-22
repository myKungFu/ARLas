function [tau,tau5,tau95] = sfoaeLatency(cf)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% tau = latency(cf);
%
% Calculate the expected latency for a low-level SFOAE or TEOAE. Based on
% Shera, Guinan, & Oxenham (2002) Proc Natl Acad Sci 99:3318-3323.
%
% cf = center frequency (Hz)
% tau = mean expected group delay in seconds.
% tau5 = 5th percentile expected delay (s)
% tau95 = 95th percetile expected delay (s)
%
% Authors: Shawn S. Goodman
% Dept. of Communication Sciences & Disorders
% University of Iowa
% Date: August 10, 2015
% Updated: June 25, 2021 -- ssg
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    cf = cf / 1000; % Convert cf to kHz
    
    % mean values ----------------------------
    alpha = 0.37; % alpha coefficient
    beta = 5.5; % beta coefficient
    nBM = beta * cf.^alpha; % basilar membrane group delay (dimensionless)
    nSfoae = 2 * nBM; % SFOAE group delay (dimensionless)
    tau = (nSfoae./cf)/1000; % SFOAE latency (s)

    % 5th percentile values ----------------------------
    alpha5 = 0.3; % alpha coefficient
    beta5 = 4.9; % beta coefficient
    nBM = beta5 * cf.^alpha5; % basilar membrane group delay (dimensionless)
    nSfoae = 2 * nBM; % SFOAE group delay (dimensionless)
    tau5 = (nSfoae./cf)/1000; % SFOAE latency (s)

    % 95th percentile values ----------------------------
    alpha95 = 0.44; % alpha coefficient
    beta95 = 6.1; % beta coefficient
    nBM = beta95 * cf.^alpha95; % basilar membrane group delay (dimensionless)
    nSfoae = 2 * nBM; % SFOAE group delay (dimensionless)
    tau95 = (nSfoae./cf)/1000; % SFOAE latency (s)


end