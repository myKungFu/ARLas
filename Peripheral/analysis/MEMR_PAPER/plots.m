
load('C:\MEMR_Analysis1\MEM09_Run1_Analysis1.mat')
subplot(2,1,1); plot(analysisPart1_mem.sigHat); hold on; plot(analysisPart1_mem.sigHatNF,'k')
subplot(2,1,2); plot(analysisPart1_mem.nfHat); hold on; plot(analysisPart1_mem.nfHatNF,'k')
load('C:\MEMR_Analysis1\MEM09_Run2_Analysis1.mat')
subplot(2,1,1); plot(analysisPart1_mem.sigHat); hold on; plot(analysisPart1_mem.sigHatNF,'k')
subplot(2,1,2); plot(analysisPart1_mem.nfHat); hold on; plot(analysisPart1_mem.nfHatNF,'k')
load('C:\MEMR_Analysis1\MEM09_Run3_Analysis1.mat')
subplot(2,1,1); plot(analysisPart1_mem.sigHat); hold on; plot(analysisPart1_mem.sigHatNF,'k')
subplot(2,1,2); plot(analysisPart1_mem.nfHat); hold on; plot(analysisPart1_mem.nfHatNF,'k')
load('C:\MEMR_Analysis1\MEM09_Run4_Analysis1.mat')
subplot(2,1,1); plot(analysisPart1_mem.sigHat); hold on; plot(analysisPart1_mem.sigHatNF,'k')
subplot(2,1,2); plot(analysisPart1_mem.nfHat); hold on; plot(analysisPart1_mem.nfHatNF,'k')

subplot(2,1,1)
load('C:\MEMR_Analysis1\MEM09_Run1_Analysis1.mat')
plot(20*log10(mean(analysisPart1_mem.D1,2)),'b')
load('C:\MEMR_Analysis1\MEM09_Run2_Analysis1.mat')
plot(20*log10(mean(analysisPart1_mem.D1,2)),'b')
load('C:\MEMR_Analysis1\MEM09_Run3_Analysis1.mat')
plot(20*log10(mean(analysisPart1_mem.D1,2)),'b')
load('C:\MEMR_Analysis1\MEM09_Run4_Analysis1.mat')
plot(20*log10(mean(analysisPart1_mem.D1,2)),'b')


figure
plot(abs(analysisPart1_mem.signal_sm),'Color',[.8 .8 .8])
hold on
plot(abs(analysisPart1_mem.signal_sm(:,1)),'r','LineWidth',1)
plot(abs(analysisPart1_mem.signal_sm(:,3)),'m','LineWidth',1)
plot(abs(analysisPart1_mem.signal_sm(:,5)),'g','LineWidth',1)
plot(abs(analysisPart1_mem.signal_sm(:,7)),'b','LineWidth',1)
plot(abs(analysisPart1_mem.signal_sm(:,9)),'c','LineWidth',1)
title('Mag only')

figure
plot(20*log10(analysisPart1_mem.D1),'Color',[.8 .8 .8])
hold on
plot(20*log10(analysisPart1_mem.D1(:,1)),'r','LineWidth',1)
plot(20*log10(analysisPart1_mem.D1(:,3)),'m','LineWidth',1)
plot(20*log10(analysisPart1_mem.D1(:,5)),'g','LineWidth',1)
plot(20*log10(analysisPart1_mem.D1(:,7)),'b','LineWidth',1)
plot(20*log10(analysisPart1_mem.D1(:,9)),'c','LineWidth',1)
title('D1')

% 
% figure
% plot(abs(analysisPart1_mem.signal_sm(:,1:5)),'b')
% hold on
% plot(abs(analysisPart1_mem.signal_sm(:,6:11)),'r')
% legend('low freq','mid freq')
% 

figure

q = analysisPart1_mem;
plot(q.Z_sm,'k')
hold on
plot(q.Z_sm(:,1),'r')
plot(q.Z_sm(:,3),'m')
plot(q.Z_sm(:,5),'g')
plot(q.Z_sm(:,7),'b')
plot(q.Z_sm(:,9),'c')





Col = [0 1 1];
q = analysisPart1_mem;
plot(q.Z_sm(80,1:8),'Color',Col)
hold on
plot(q.Z_sm(70,1:8),'Color',Col)
plot(q.Z_sm(60,1:8),'Color',Col)
plot(q.Z_sm(50,1:8),'Color',Col)
plot(q.Z_sm(40,1:8),'Color',Col)
plot(q.Z_sm(30,1:8),'Color',Col)
plot(q.Z_sm(20,1:8),'Color',Col)
plot(q.Z_sm(10,1:8),'Color',Col)


for ii=1:size(q.Z,2)
    dummy = q.Z(:,ii);
    shift(ii,1) = mean([dummy(1),dummy(end)]);
    dummy = dummy ./ shift(ii);
    Z(:,ii) = dummy;
end


