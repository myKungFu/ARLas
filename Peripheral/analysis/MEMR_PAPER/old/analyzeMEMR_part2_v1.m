function [] = analyzeMEMR_part2_v1()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% analyzeMEMR_part2_v1
%
% Author:
% Date:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%    load 'MEM382'
load 'MEM062'

    % R = analysisPart1_mem.modeR;
    % I = analysisPart1_mem.modeI;
    % UR = analysisPart1_mem.hdi95UR;
    % LR = analysisPart1_mem.hdi95LR;
    % UI = analysisPart1_mem.hdi95UI;
    % LI = analysisPart1_mem.hdi95LI;
    % R = real(analysisPart1_mem.Z);
    % I = imag(analysisPart1_mem.Z);
    % Rn = real(analysisPart1_mem.Zn);
    % In = imag(analysisPart1_mem.Zn);
    Z = analysisPart1_mem.Z;
    Zorig = Z;
    Zn = analysisPart1_mem.Zn;
    nFreqs = analysisPart1_mem.nFreqs;
    nClicks = 160;

    sm = 0.01; % larger numbers (approaching 1) are less smoothing; approaching zero = more smooth
    x = (1:1:nClicks*3)';
    for ii=1:nFreqs
        z = Z(:,ii);
        zn = Zn(:,ii);

        mr = real(z);
        offsetR = mean([mr(1),mr(end)]);
        mr = mr - offsetR;
        mi = imag(z);
        offsetI = mean([mi(1),mi(end)]);
        mi = mi - offsetI;
        
        mrn = real(zn); % noise estimates
        minn = imag(zn);
        offsetRn = mean([mrn(1),mrn(end)]);
        offsetIn = mean([minn(1),minn(end)]);
        mrn = mrn - offsetRn;
        minn = minn - offsetI;
        
        qr = [flipud(mrn);mr;mrn];
        qi = [flipud(mi);mi;mi];
        qrn = [flipud(mrn);mrn;mrn];
        qinn = [flipud(minn);minn;minn];
        %qr = medianSmoother(qr,4);
        %qi = medianSmoother(qi,4);
        %mrsm = meanSmoother(mrsm,6);

        w = ones(size(qr));
        % if any(isnan(qr))
        %     keyboard
        % end
        ppr = csaps(x,qr,sm,[],w); 
        ppi = csaps(x,qi,sm,[],w); 
        pprn = csaps(x,qrn,sm,[],w); 
        ppin = csaps(x,qinn,sm,[],w); 
        
        mr_sm = ppval(ppr,x);
        mr_sm = mr_sm(161:320);
        mrn_sm = ppval(pprn,x);
        mrn_sm = mrn_sm(161:320);

        mi_sm = ppval(ppi,x);
        mi_sm = mi_sm(161:320);
        min_sm = ppval(ppin,x);
        min_sm = min_sm(161:320);

        %PD = fitdist(q','rayleigh');
        
        mr_sm = mr_sm + offsetR;
        mi_sm = mi_sm + offsetI;
        z = mr_sm + 1i*mi_sm;

        mrn_sm = mrn_sm + offsetRn;
        min_sm = min_sm + offsetIn;
        zn = mrn_sm + 1i*min_sm;
        
        Z(:,ii) = z ./ z(1);
        Zorig(:,ii) = Zorig(:,ii) ./ z(1);
        Zn(:,ii) = zn ./ z(1);
    end

    nfCutSample = 25;
    nf = abs(Zorig(1:nfCutSample,:)-1);
    nf = nf(:);
    PD = fitdist(nf,'rayleigh');
    cut = PD.icdf(0.95);

    % figure
    % plot((abs(Zorig-1)),'Color',[.7 .7 .7],'LineWidth',0.5)
    % hold on
    % plot((mean(abs(Zorig-1),2)),'Color',[0 0 1],'LineWidth',2)
    % line([1 160],([cut,cut]),'Color',[0 1 0],'LineWidth',2)   

    q = (nanmean(abs(Zorig),1));
    direction = zeros(size(q));
    indxPos = find(q>=1);
    direction(indxPos) = 1;
    indxNeg = find(q<1);
    direction(indxNeg) = -1;
    Direction = repmat(direction,size(Zorig,1),1);

    figure 
    subplot (2,1,1)
    plot((abs(Zorig-1)),'Color',[.7 .7 .7],'LineWidth',0.5)
    hold on
    plot((nanmean(abs(Zorig-1),2)),'Color',[0 0 1],'LineWidth',1)
    line([1 160],([cut,cut]),'Color',[0 1 0],'LineWidth',0.75)   

    % subplot(2,1,1)
    % plot(Direction.*(abs(Zorig-1)),'Color',[.7 .7 .7],'LineWidth',0.5)
    % hold on
    % line([1 160],([cut,cut]),'Color',[0 1 0],'LineWidth',2)   
    % line([1 160],(-[cut,cut]),'Color',[0 1 0],'LineWidth',2)   
    % 


% Apply smoothing to the mean growth function -----------------------------    
sm = 0.95; % smoothing factor (smaller numbers are more smooth)
stimLength = 8; % noise stimulus length (seconds)
n = 160;
x = linspace(0,8,n)'; % x-axis (click number)
y = nanmean(abs(Zorig-1),2); % get the mean growth function

w = ones(size(x)); % weighting factor
pp = csaps(x,y,sm,[],w); % piecewise polynomial object
ysm = ppval(pp,x); % for plotting, evaluate at x (ysm = y smoothed)

subplot (2,1,2)
plot(x,y,'Color',[.7 .7 .7],'LineWidth',0.5)
hold on
plot(x,ysm,'Color',[0 0 1],'LineWidth',1)
line([x(1),x(end)],([cut,cut]),'Color',[0 1 0],'LineWidth',0.75)  

% find the location of the peak -------------------
dfdx = fnder(pp); % take derivative and solve for zero slope
peakX = fnzeros(dfdx); % peak location in seconds
peakX = peakX(1); % function returns two outputs (both the same) take the first
peakY = ppval(pp,peakX); % max average activation
delay = peakX-4; % reflex delay in seconds
%dY = ppval(dfdx,x);

plot(peakX,peakY,'or')


% find threshold
pp2 = csaps(x,y-cut,sm,[],w); % piecewise polynomial real coefficients
z2 = fnzeros(pp2);
Thd = z2(1,:); % threshold times

plot(Thd(1),cut,'*r')
plot(Thd(2),cut,'*r')


% hysteresis ---------------------
% integration
hh = fnint(pp); % integrate the smoothed spline
A = ppval(hh,Thd(1));
B = ppval(hh,peakX);
C = ppval(hh,Thd(2));
aucLeft = B-A; % area under the curve left
aucRight = C-B; % area under the curve right
hyst = aucRight / aucLeft; % hysteresis as a ratio

title(['Hysteresis Ratio = ',num2str(hyst)])

keyboard

    for ii=2:8
        figure
        %nf = nanmean(abs(Zn(:,ii)*2+1));
        %nf2 = (-(nf-1))+1;
        plot(abs(Z(:,ii)-1))
        hold on
        plot(abs(Zn(:,ii)),'r')
        %line([1 160],[nf,nf],'Color',[1 0 0])
        %line([1 160],[nf2,nf2],'Color',[1 0 0])
        pause(0.001)
    end



    for ii=2:8
        figure
        nf = nanmean(abs(Zn(:,ii)*2+1));
        nf2 = (-(nf-1))+1;
        plot(abs(Z(:,ii)))
        hold on
        plot(abs(Zn(:,ii)+1),'r')
        line([1 160],[nf,nf],'Color',[1 0 0])
        line([1 160],[nf2,nf2],'Color',[1 0 0])
        pause(0.001)
    end


    %Z = R + 1i*I;
    %nAve = 3;
    %Z = analysisPart1_mem.z;
c    
pause(0.1)
    L = zeros(nClicks,nFreqs);
    for ii=1:nFreqs
        %z = Z(:,ii);
        %baseZ = [z(1:nAve);z(end-nAve+1:end)];
        %baseZ = mean(baseZ);
        %z = z ./ baseZ;
        %Z(:,ii) = abs(z);

        z = Z(:,ii);
        zn = Zn(:,ii);
        % z = re + 1i*im;
        % nf_zu = nf_ur + 1i*nf_ui;
        % nf_zl = nf_lr + 1i*nf_li;
        % 
        % baseline = z(1);
        % z = z / baseline;
        % nf_zu = nf_zu / baseline;
        % nf_zl = nf_zl / baseline;
        % 
        % Z(:,ii) = z;
        % NF_ZU(:,ii) = nf_zu;
        % NF_ZL(:,ii) = nf_zl;

        [L(:,ii),t,D(:,ii),inflection(ii,1)] = crawl(z);
        L1 = L(1,ii);
        L2 = L(end,ii);
        lineL = linspace(L1,L2,160)';
        Lfix(:,ii) = L(:,ii)-lineL;

        [Ln(:,ii),t,Dn(:,ii),inflectionn(ii,1)] = crawl(zn);
        L1 = L(1,ii);
        L2 = L(end,ii);
        lineL = linspace(L1,L2,160)';
        Lfix(:,ii) = L(:,ii)-lineL;


        figure
        plot(L(:,ii),'Color',[.7 .7 .7])
        hold on
        plot(Lfix(:,ii),'r')
        plot(abs(zn),'k')
        pause(0.01)


        % q1 = complex(z(1));
        % q2 = complex(z(end));
        % lineR = linspace(real(q1),real(q2),160)';
        % lineI = linspace(imag(q1),imag(q2),160)';
        % zz = real(z)-lineR + 1i*(imag(z)-lineI);
        % zz = zz + 1;
        % 
        % [L_left,t,D(:,ii),inflection(ii,1)] = crawl(z);
        % zz = flipud(z);
        % zz = zz ./ zz(1);
        % L_right = crawl(zz);
        % q1 = linspace(1,0,160)';
        % q2 = linspace(0,1,160)';
        % LL = (L_left.*q1 + L_right.*q2);
        % 
        % plot(L_left,'b')
        % hold
        % plot(L_right,'r')
        % plot(LL,'g')



        %[Ln(:,ii),t,Dn(:,ii),inflectionn(ii,1)] = crawl(zn);
        %[Lnfu(:,ii),t,D(:,ii),inflection(ii,1)] = crawl(nf_zu);
        %[Lnfl(:,ii),t,D(:,ii),inflection(ii,1)] = crawl(nf_zl);
    end


Lf = mean(Lfix,2);
dLf = gradient(Lf);
[~,thdIndx] = min(dLf(1:45));
yCut = Lf(thdIndx);

c
figure
plot(Lfix,'Color',[.7 .7 .7])
hold on
plot(mean(Lfix,2),'r','LineWidth',2)
line([1 160],[yCut,yCut],'Color',[0 0 1])


keyboard
    
    % for jj=1:10
    %     figure
    %     z = Z(:,jj);
    %     zn = Zn(:,jj);
    %     adjust = abs(z(1));
    %     plot(abs(z),'k')
    %     hold
    %     plot(abs(zn)+1,'r')
    %     pause(0.001)
    % end

keyboard
    figure (101)
    
    jj = 10;
    plot(abs(Z(:,jj)),'r')
    title('MEMR Wave','Threshold and Delay');
    hold
    plot(abs(NF_ZU(:,jj)),'k--');
    plot(abs(NF_ZL(:,jj)),'k--');
    delay = mean(inflection);
    xline(delay,'b');
    delay = ((mean(inflection)-80)*.05)*1000;


    figure
    subplot(1,2,1)
    plot(t,abs(Z))
    xlabel('Time (seconds)')
    ylabel('Magnitude')
    subplot(1,2,2)
    plot(t,L)
    xlabel('Time (seconds)')
    ylabel('Integrated Length')

    figure
    for jj=1:10
        h1 = L(1:inflection(jj),jj);
        h2 = L(inflection(jj):end,jj);
        padN = length(h1)-length(h2);
        h2 = [h2;-inf(padN,1)];
        h2 = flipud(h2);
        
        plot(h1)
        hold on
        plot(h2)
    end



    
    keyboard
end

% Internal Functions ------------------------------------------------------

function [L,t,D,inflection] = crawl(z)
% return the length of a curve
    N = 160; % number of clicks
    dt = 0.05; % change in time (time step is 50 ms)
    t = (0:1:160-1)'*dt; % time vector from 0 to 8 seconds
    x = real(z);
    y = imag(z);
    dx = gradient(x); % derivative of real part
    dy = gradient(y); % derivative of imaginary part

    d = sqrt((dx./dt).^2+(dy./dt).^2); % derivative of the curve
    % find the local minimum
    start = 60;
    finish = 100;
    [~,indxMin] = min(d(start:finish));
    % set the k vector
    k = ones(N,1);
    inflection = indxMin + start;
    k(inflection:end) = -1; % after the midpoint, k is negative
    L = cumsum(d .* k); % intigrate the derivative times k
        %L = cumsum(k.*(sqrt((dx./dt).^2+(dy./dt).^2)));
    D = gradient(L)/dt;
end
