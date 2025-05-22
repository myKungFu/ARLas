function [] = goodmanMEMR_analysis_v2()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    fileNameL = 'Ch3_ER10xA_memr_0001.mat';
    fileNameR = 'Ch4_ER10xB_memr_0001.mat';
    
    basePath = 'C:\myWork\ARLas\Data\MEMR_goodmanTest\';
    dSubj = dir(basePath);
    nSubj = size(dSubj,1);
try    
    for ii=18:nSubj
        disp(['Analyzing subject ',num2str(ii-2),' of ',num2str(nSubj-2)])
        dRuns = dir([basePath,dSubj(ii).name]);
        nRuns = size(dRuns,1);
        for jj=3:nRuns
            disp(['  Analyzing run ',num2str(jj-2),' of ',num2str(nRuns-2)])
            subjectName = [dSubj(ii).name,'\'];
            runName = [dRuns(jj).name,'\'];
            pathName = [basePath,subjectName,runName];
            MEMR = analyzeMEMR_v7(pathName,fileNameL,fileNameR,subjectName,runName);
        end
    end
catch ME
    keyboard
end
end

% OLD CODE ----------------------------------------------------------------

%MEMR = analyzeMEMR_v5(DataL,DataR,time,fs,clickIndx);

% load('MEM20_analyzedMEMR_1.mat')
% %load('MEM01_analyzedMEMR_1.mat')
% DataL = MEMR.DataL;
% DataR = MEMR.DataR;
% fs = MEMR.fs;
% clickIndx = MEMR.clickIndx;
% time = MEMR.time;


% 
% % The arc length of a parametric curve (x(t),y(t)) over the interval (a,b)
% % can be found by integrating:
% 
% t = timeChunk;
% dt = median(gradient(t));
% x = q10;
% 
% xr = real(x);
% xi = imag(x);
% smoothing = 0.999999;
% ppr = csaps(t,xr,smoothing);
% ppi = csaps(t,xi,smoothing);
% xr_sm = ppval(ppr,t);
% xi_sm = ppval(ppi,t);
% int = cumsum(( (gradient(real(x),t))   + (gradient(imag(x),t))  )) * dt;
% int_sm = cumsum(( (gradient(xr_sm,t))   + (gradient(xi_sm,t))  )) * dt;
% 
% arclen = cumsum(sqrt((        (gradient(xr_sm,t)).^2   + (gradient(xi_sm,t)).^2        )    )) * dt;
% 
% % check for correct amount of smoothing
% % subplot(2,1,1)
% % plot(t,xr,'b')
% % hold on
% % fnplt(ppr, 'r'); 
% % subplot(2,1,2)
% % plot(t,xi,'b')
% % hold on
% % fnplt(ppi, 'r'); 
% 
% arclen = cumsum(sqrt( (gradient(real(x),t)).^2   + (gradient(imag(x),t)).^2  )) * dt;
% int = cumsum(( (gradient(real(x),t))   + (gradient(imag(x),t))  )) * dt;
% 
% plot(t,int_sm,'r')
% hold on
% 
% % compare with these:
% multiplier = -1;
% plot(t,multiplier*(abs(complex(x-1))),'b.-')
% hold on
% plot(t,(abs(complex(x)))-1,'g.-')
% 
% legend('integrate','abs-1','abs','Location','SouthEast')
% 
