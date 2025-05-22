function [h10,h11] = analyzeMEMR_v3(Data,time,fs,indx)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [h10,h11] = analyzeMEMR_v3(DataL,time,fs,indx)
%
% Author: Shawn Goodman
% Date: August 6, 2022
% Last Updated: August 9, 2022
% Last Updated: October 10, 2022 -- ssg -- version 2; goes with
%                       goodmanMEMR_v2.m.
% Last Updated: October 11, 2022 -- ssg -- added figure handles to output.
% Last Updated: February 27, 2023 -- ssg --
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    time = time / time(end); % express time from 0 to 1

    % frequencies of interest in 1/2 octave steps
    peakFreqs = [125
                 177
                 250
                 354
                 500
                 707
                1000
                1414
                2000
                2828];



%     sigStart = 1; % signal start index
%     sigFinish = 15; % signal finish index
%     windowN = round(fs*0.01);
% 
%     startOffset = 1; %13;
%     finishOffset = 60; % 60
%     whichOne = 1; % which fft bin
%     counter = 1;

    % get the template from the first and last 3 clicks, averaged
    % Each click window is 50 ms long, so the average of 3 is 150 ms
    % We used a noise activator with a rise time of 4 seconds and a fall
    % time of 4 seconds. This gives 8 seconds / .05 = 160 clicks.

    % Within each click window, we want to include several round trip
    % travel times of the click in the canal, but no low-frequency OAEs. We
    % should consider only the first 2 ms of the click, which is 
    % round(0.002*fs) = 192 samples. 
    nSamples = 193; % number of time samples in each click analysis window
    nClicks = length(indx); % number of total clicks to analyze

    for jj=1:3 % loop across click
        chunk = Data(indx(jj):indx(jj)+nSamples,:);
% -->        % add artifact rejection here
        template(:,jj) = mean(chunk,2);

    end
    counter = 4;
    for jj=nClicks-2:nClicks % loop across click
        chunk = Data(indx(jj):indx(jj)+nSamples,:);
        template(:,counter) = mean(chunk,2);
        counter = counter + 1;
    end
    templateNoise = std(template,[],2)/sqrt(size(template,2));
    template = mean(template,2);

    if ~exist('iscS1L','var')
        iscS1L = [];
    end



%     for ii=1:size(Data,2)
%         offset = ii*.001;
%         runme(template,Data(:,ii),nSamples,nClicks,indx,time,fs,iscS1L,offset)
%     end
    [h10,h11] = runme(template,Data,nSamples,nClicks,indx,time,fs,iscS1L,peakFreqs);


end
% INTERNAL FUNCTIONS ------------------------------------------------------

function [h10,h11] = runme(template,Data,nSamples,nClicks,indx,time,fs,iscS1L,peakFreqs)
%     tic
%     for kk=1:size(Data,2)
%         disp(['Converting ',num2str(kk),' of ',num2str(size(Data,2))])
%         [pl,Pl,phi,other,wf] = ARLas_convertPL(Data(:,kk),iscS1L);
%         FPL(:,kk) = wf.fpl;
%         RPL(:,kk) = wf.rpl;
%     end
%     PRR = RPL ./ FPL;
%     Data = PRR;
%     toc / 60

    for jj=1:nClicks % loop across each click
        Chunk = Data(indx(jj):indx(jj)+nSamples,:);
        timeChunk(1,jj) = time(indx(jj));
        
        %delta(:,jj) = mean(Chunk,2) - template;
        %[frequency,signal(:,jj),noiseFloor(:,jj)] = ARLas_fda(chunk,fs,0.00002,fs);
        %S(:,jj) = fft(mean(chunk,2));
        S(:,jj) = fft(mean(chunk,2),fs);
    end
    %rms = sqrt(mean(delta.^2,1));
    
    Template = [S(:,1:3),S(:,nClicks-2:nClicks)]; % template in the frequency domain (template = first 3 and last 3 clicks)
    Template = mean(Template,2);
    for jj=1:length(indx) 
        Delta(:,jj) = S(:,jj)./Template;
    end
    freq = (0:1:size(chunk,1)-1)*(fs/size(chunk,1));




fff = (1:1:fs)-1'; % full frequency vector
% find the locatino of the desired center frequencies
for ii=1:length(peakFreqs)
    [~,peakIndx(ii)] = min(abs(fff-peakFreqs(ii)));
end
%peakIndx = [349,933,1619,3312, 4039,7294,14109];    

%-------- keyboard
% for jj=1:length(peakIndx)-1
%         D = complex(Delta(peakIndx(jj):peakIndx(jj+1),:));
%         f = fff(peakIndx(jj):peakIndx(jj+1));
%         [~,fhatIndx] = max(mean(abs(D-1),2));
%         peakFreqs2(jj) = f(fhatIndx);
% end


t = timeChunk;
dt = median(gradient(t));
for jj=1:length(peakIndx)
    x = (complex(Delta(peakIndx(jj),:)));
    
    xr = real(x);
    xi = imag(x);
    smoothing = 0.999999;
    ppr = csaps(t,xr,smoothing);
    ppi = csaps(t,xi,smoothing);
    xr_sm = ppval(ppr,t);
    xi_sm = ppval(ppi,t);
    resid_xr = xr - xr_sm; % residuals after subtracting the smoothing spline
    resid_xi = xi - xi_sm;
    ppr = csaps(t,resid_xr,smoothing);
    ppi = csaps(t,resid_xi,smoothing);
    xr_noise = ppval(ppr,t);
    xi_noise = ppval(ppi,t);




    int = cumsum(( (gradient(xr,t))   + (gradient(xi,t))  )) * dt;
    int_smL = cumsum(( (gradient(xr_sm,t))   + (gradient(xi_sm,t))  )) * dt;
    int_smR = fliplr(cumsum(fliplr(( (gradient(xr_sm,t))   + (gradient(xi_sm,t))  ))) * dt);



    %arcLen = cumsum(sqrt((  (gradient(xr,t)).^2   + (gradient(xi,t)).^2  )  )) * dt;
    % integrate from the left
        arcLen_smL = cumsum(sqrt((  (gradient(xr_sm,t)).^2   + (gradient(xi_sm,t)).^2  )  )) * dt;
        arcLen_noiseL = cumsum(sqrt((  (gradient(xr_noise,t)).^2   + (gradient(xi_noise,t)).^2  )  )) * dt;
    % intergrate from the right    
        arcLen_smR = fliplr(cumsum(fliplr(sqrt((  (gradient(xr_sm,t)).^2   + (gradient(xi_sm,t)).^2  )  ))) * dt);
    % find where the two functions cross: this is the peak
        indxLT = find(arcLen_smL <= arcLen_smR);
        [~,peakIndx2] = max(indxLT);
        peak = arcLen_smL(peakIndx2);
    % flip the descending side
        arcLen_sm = [arcLen_smL(1:peakIndx2),arcLen_smR(peakIndx2+1:end)];
        arcLen_sm(peakIndx2) = (arcLen_sm(peakIndx2-1) + arcLen_sm(peakIndx2+1))/2;
    
    plot(arcLen_smL,'b*-')
    hold on
    plot(arcLen_smR,'r*-')
    plot(peakIndx2,peak,'go')


end

%----------------------





    h10 = figure(10); hold on % plot complex data at desired frequencies
    for jj=1:length(peakIndx)
        plot(complex(Delta(peakIndx(jj),:)),'.-')
    end
    NN = 1000; % number of samples in the ellipses
    t = linspace(0,2*pi,NN); 
    eMinus0 = [cos(t);sin(t)]'; % create a unit circle
    hold on; box on;
    plot(eMinus0(:,1),eMinus0(:,2),'LineWidth',1,'Color',[0 0 0])
    % plot the 1 to 6 dB reduction lines
    eMinus1 = eMinus0 * (10^(-1/20));
    eMinus2 = eMinus0 * (10^(-2/20));
    eMinus3 = eMinus0 * (10^(-3/20));
    eMinus4 = eMinus0 * (10^(-4/20));
    eMinus5 = eMinus0 * (10^(-5/20));
    eMinus6 = eMinus0 * (10^(-6/20));
    innerColor = [0.5 0.5 0.5];
    LW = 1; % line width
    plot(eMinus1(:,1),eMinus1(:,2),'LineStyle','--','LineWidth',LW,'Color',innerColor)
    plot(eMinus2(:,1),eMinus2(:,2),'LineStyle','--','LineWidth',LW,'Color',innerColor)
    plot(eMinus3(:,1),eMinus3(:,2),'LineStyle','--','LineWidth',LW,'Color',innerColor)
    plot(eMinus4(:,1),eMinus4(:,2),'LineStyle','--','LineWidth',LW,'Color',innerColor)
    plot(eMinus5(:,1),eMinus5(:,2),'LineStyle','--','LineWidth',LW,'Color',innerColor)
    plot(eMinus6(:,1),eMinus6(:,2),'LineStyle','--','LineWidth',LW,'Color',innerColor)
    % Format the plot
    Lh = line([-1 1],[0 0]); 
    set(Lh,'Color',[0 0 0],'LineStyle','-','LineWidth',1)
    Lv = line([0 0],[-1 1]);
    set(Lv,'Color',[0 0 0],'LineStyle','-','LineWidth',1)
    axis square
    legend(['f=',num2str(round(freq(2))),' Hz'],['f=',num2str(round(freq(3))),' Hz'],['f=',num2str(round(freq(4))),' Hz'],['f=',num2str(round(freq(5))),' Hz'],['f=',num2str(round(freq(6))),' Hz'])
    xlim([0.5 1.3])
    ylim([-0.4 0.4])
    xlabel('Real')
    ylabel('Imaginary')

    %-----------------------------------
    CCC = [0 0 0
        1 0 0
        .5 .5 0
        0 .5 .5
        0 1 0
        0 0 .5
        0 0 1];
    h11 = figure(11); hold on
    for jj=1:length(peakIndx)
        delta = (complex(Delta(peakIndx(jj),:)));
        if max(delta) >=1 %max(abs(Delta(jj,:)))>=1
            multiplier = 1;
        else
            multiplier = -1;
        end
        %plot(timeChunk,multiplier*20*log10(abs(complex(Delta(peakIndx(jj),:))-1)),'.-') %,'Color',CCC(jj,:)
        %plot(timeChunk,multiplier*(abs(complex(Delta(peakIndx(jj),:))  )),'.-')
        plot(timeChunk,multiplier*(abs(complex(delta-1))),'.-')

        multiplier
    end
    %legend(['f=',num2str(round(freq(2))),' Hz'],['f=',num2str(round(freq(3))),' Hz'],['f=',num2str(round(freq(4))),' Hz'],['f=',num2str(round(freq(5))),' Hz'],['f=',num2str(round(freq(6))),' Hz'])
    legend(num2str(peakIndx))

q = 20*log10(abs(complex(Delta(peakIndx(5),:))));
A = 1;
B = 1;
%obj = fitIBeta(timeChunk,q,A,B);
obj = [];
%--------------------------------------------------------------------------








end


% OLD CODE ----------------------------------------------------------------
% 
    % A two second rise fall has four seconds, so 80 clicks, of which 6
    % (3.75 from each edge) form the template.
    
