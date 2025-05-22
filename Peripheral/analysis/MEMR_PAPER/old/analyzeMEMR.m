function [] = analyzeMEMR()
% Author: Shawn Goodman
% Date: August 6, 2022
% Last Updated: August 9, 2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %load('Julia032sec.mat') % data file to load
    load('RW14sec.mat')
    %load('Rachel_4sec.mat') % data file to load

    
    Data = DataL;
    time = time / time(end);

    sigStart = 1; 
    sigFinish = 15;
    windowN = round(fs*0.01);

    startOffset = 1; %13;
    finishOffset = 60; % 60
    whichOne = 1; % which fft bin
    counter = 1;

    % get the template from the first and last 3 clicks, averaged
    % Each click window is 50 ms long, so the average of 3 is 150 ms
    % A two second rise fall has four seconds, so 80 clicks, of which 6
    % (3.75 from each edge) form the template.
    % Within each click window, we want to include several round trip
    % travel times of the click in the canal, but no low-frequency OAEs. We
    % should consider only the first 2 ms of the click, which is 
    % round(0.002*fs) = 192 samples. 
    nSamples = 193;
    nClicks = length(indx);

    for jj=1:3 % loop across click
        chunk = Data(indx(jj):indx(jj)+nSamples,:);
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
    runme(template,Data,nSamples,nClicks,indx,time,fs,iscS1L)


end
% INTERNAL FUNCTIONS ------------------------------------------------------

function [] = runme(template,Data,nSamples,nClicks,indx,time,fs,iscS1L,offset)
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

    for jj=1:nClicks % loop across click
        chunk = Data(indx(jj):indx(jj)+nSamples,:);


        timeChunk(1,jj) = time(indx(jj));
        delta(:,jj) = mean(chunk,2) - template;
        %[frequency,signal(:,jj),noiseFloor(:,jj)] = ARLas_fda(chunk,fs,0.00002,fs);
        %S(:,jj) = fft(mean(chunk,2));
        S(:,jj) = fft(mean(chunk,2),fs);
    end
    rms = sqrt(mean(delta.^2,1));
    
    Template = [S(:,1:3),S(:,nClicks-2:nClicks)];
    Template = mean(Template,2);
    for jj=1:length(indx) 
        Delta(:,jj) = S(:,jj)./Template;
    end
    freq = (0:1:size(chunk,1)-1)*(fs/size(chunk,1));



peakFreqs = [ 125
         177
         250
         354
         500
         707
        1000
        1414
        2000
        2828];
fff = (1:1:fs)-1';
for ii=1:length(peakFreqs)
    [~,peakIndx(ii)] = min(abs(fff-peakFreqs(ii)));
end
%peakIndx = [349,933,1619,3312, 4039,7294,14109];    

    figure(10); hold on
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
    figure(11); hold on
    for jj=1:length(peakIndx)
        if max(abs(Delta(jj,:)))>=1
            multiplier = 1;
        else
            multiplier = -1;
        end
        plot(timeChunk,multiplier*20*log10(abs(complex(Delta(peakIndx(jj),:)))),'.-') %,'Color',CCC(jj,:)
    end
    %legend(['f=',num2str(round(freq(2))),' Hz'],['f=',num2str(round(freq(3))),' Hz'],['f=',num2str(round(freq(4))),' Hz'],['f=',num2str(round(freq(5))),' Hz'],['f=',num2str(round(freq(6))),' Hz'])
    legend(num2str(peakIndx))

q = 20*log10(abs(complex(Delta(peakIndx(5),:))));
A = 1;
B = 1;
obj = fitIBeta(timeChunk,q,A,B);
%--------------------------------------------------------------------------








end


% OLD CODE ----------------------------------------------------------------
% 
% 
% 
%     signalPa = 10.^(signal/20)*0.00002;
%     templatePa = [signalPa(:,1:3),signalPa(:,nClicks-2:nClicks)];
%     templatePa = mean(templatePa,2);
%     for jj=1:length(indx)
%         dIndx = signalPa(:,jj) < templatePa;
%         dummy = signalPa(dIndx,jj);
%         temp = templatePa(dIndx);
%         rmsLess(jj,1) = sum(dummy-temp) / sum(temp);
% 
%         dIndx = signalPa(:,jj) > templatePa;
%         dummy = signalPa(dIndx,jj);
%         temp = templatePa(dIndx);
%         rmsMore(jj,1) = sum(dummy-temp) / sum(temp);
%     end
% 
%     RMS = sqrt(mean(signalPa.^2,1));
%     RMSdB = 20*log10(RMS/0.00002);
% 
%     for jj=1:length(indx) % loop across click
%         if indx(jj)-windowN/2 >= 1 & indx(jj)-windowN/2 < size(DataL,1)
%             chunk = Data(indx(jj)-windowN/2:indx(jj)+windowN/2,:);
%             timeChunk(1,jj) = time(indx(jj));
%             [~,hat] = max(abs(mean(chunk,2)));
%             cc = chunk(hat+startOffset-1:hat+finishOffset+200,:); %chunk(477:540,:);
%             [frequency,signal,noiseFloor] = ARLas_fda(cc,fs,0.00002,size(cc,1));
%             signal = signal(sigStart:sigFinish);
%             if counter == 1
%                 template = mean(cc,2);
%             end
%             CC(:,jj) = mean(cc,2);
%             Delta(:,jj) = mean(cc,2) - template;
%             RMS(:,jj) = Delta(whichOne);
%             %delta = 10.^(signal(whichOne)/20)*0.00002;
%             %RMS(:,jj) = 20*log10(delta/0.00002)
%             counter = counter + 1;
%         end
%     end
% 
%     plot(Delta);
%     
%     figure
%     plot(20*log10(abs(fft(template,fs))))
%     hold on 
%     plot(20*log10(abs(fft(Delta(:,40),fs)))) % 75
% 
% 
% 
%     keyboard
% 
%     RMS = RMS(:,2:end);
%     RMS = RMS - RMS(1);
%     figure(10)
%     hold on
%     plot(timeChunk(2:end)',mean(RMS,1))
%         
% 
%     keyboard
% end        
