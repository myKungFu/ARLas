function [h] = inSituCheck(probeLabel,probeChannel,header,data,fs,whichProbe,h)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% h = inSituCheck(header,data,h);
%
% Check calibration chirp to make sure it looks like a good fit.
% Creates two subplots: 1) All chirps stacked on top of each other. Use
% this to get a good idea of whether the levels are correct. Peak levels
% should fall between 1 and 2 pascals in most earss. 2) All chirps
% stretched out into a long vector (the way they were recorded). Use this
% to look for probe slippage over time or other strange behavior.
%
% Author: Shawn Goodman, PhD
%         Auditory Research Lab, University of Iowa
% Date: November 17, 2023
% Last Updated: November 17, 2023
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    N = size(data,1);
    t = (0:1:N-1)'/fs*1000;

    % Plot all the waveforms stacked on top of each other -----------------
    if isempty(h)
        h = figure(333);
    end
    if whichProbe == 1
        subplot(3,2,1)
    elseif whichProbe == 2
        subplot(3,2,2)
    end
    plot(t,data)
    hold on
    line([0 t(end)],[1 1],'LineStyle','--','LineWidth',0.5,'Color',[0 0 0])
    line([0 t(end)],[-1 -1],'LineStyle','--','LineWidth',0.5,'Color',[0 0 0])
    line([0 t(end)],[.5 .5],'LineStyle',':','LineWidth',0.5,'Color',[0 0 0])
    line([0 t(end)],[-.5 -.5],'LineStyle',':','LineWidth',0.5,'Color',[0 0 0])
    xlim([0 t(end)])
    ylim([-1.25 1.25])
    xlabel('Time (ms)')
    ylabel('Amplitude (Pa)')
    title([probeLabel,'   Channel ',num2str(probeChannel)])

    % Plot the spectra ----------------------------------------------------
    if whichProbe == 1
        subplot(3,2,5)
    elseif whichProbe == 2
        subplot(3,2,6)
    end
    [frequency,signal,noiseFloor] = ARLas_fda(data,fs,0.00002);
    plot(frequency/1000,signal,'b')
    hold on
    plot(frequency/1000,noiseFloor,'k')
    hold on
    xlim([.1 20])
    %ylim([-1.5 1.5])
    xlabel('Frequency (kHz)')
    ylabel('Magnitude (dB SPL)')
    title([probeLabel,'   Channel ',num2str(probeChannel)])

    
    % Plot all the waveforms stretched out into a single vector -----------
    data = data(:);
    N = length(data);
    t = (0:1:N-1)'/fs;
    
    if whichProbe == 1
        subplot(3,2,3)
    elseif whichProbe == 2
        subplot(3,2,4)
    end
    plot(t,data)
    hold on
    line([0 t(end)],[1 1],'LineStyle','--','LineWidth',0.5,'Color',[0 0 0])
    line([0 t(end)],[-1 -1],'LineStyle','--','LineWidth',0.5,'Color',[0 0 0])
    line([0 t(end)],[.5 .5],'LineStyle',':','LineWidth',0.5,'Color',[0 0 0])
    line([0 t(end)],[-.5 -.5],'LineStyle',':','LineWidth',0.5,'Color',[0 0 0])
    xlim([0 t(end)])
    ylim([-1.25 1.25])
    xlabel('Time (s)')
    ylabel('Amplitude (Pa)')
    title([probeLabel,'   Channel ',num2str(probeChannel)])

end


