load('C:\myWork\ARLas\Data\EARO\EARO_20SEP2024\EARO_20SEP2024_FRSS002_run1\EARO_20SEP2024_FRSS002_run1_analysis\EARO_analyzedDPOAE_Right_1.mat')
f2 = DPOAE_L.f2;
Ldp = DPOAE_L.Ldp;
L2 = DPOAE_L.L2;
L1 = DPOAE_L.L1;

hh = figure;
sp1 = subplot(2,1,1);
 plot(f2,Ldp)
 hold on
 xmin = min(Ldp(:));
 xmax = max(Ldp(:));
 xlabel('Frequency (Hz)')
 ylabel('Ldp (dB SPL)')
 title('DP Fine Structure')
sp2 = subplot(2,1,2);
 xlabel('L2 (dB SPL)')
 ylabel('Ldp (dB SPL)')
 title('LGF')

[xx,yy] = ginput(1);
target = round(xx);
[~,indx] = min(abs(f2(:,1) - target));
if target < min(f2(:))
    done = 1;
else
    done = 0;
end
while done == 0
    axes(sp1)
     hold off
     plot(f2,Ldp)
     hold on
     xmin = min(Ldp(:));
     xmax = max(Ldp(:));
     line([target,target],[xmin,xmax],'Color',[1 0 0],'LineWidth',1,'LineStyle','-')
     xlabel('Frequency (Hz)')
     ylabel('Ldp (dB SPL)')
     title('DP Fine Structure')
    axes(sp2)
     hold off
     LGFx = L2(indx,:);
     LGFy = Ldp(indx,:);
     plot(L2(indx,:),Ldp(indx,:),'ro-')    
     xlabel('L2 (dB SPL)')
     ylabel('Ldp (dB SPL)')
     title('LGF')

    [xx,yy] = ginput(1);
    target2 = round(xx);
    if target2 < min(f2(:))
        done = 1;
    else
        target = target2;
        [~,indx] = min(abs(f2(:,1) - target));
    end
end
disp(['Target frequency for LGF is ',num2str(target),' Hz'])
