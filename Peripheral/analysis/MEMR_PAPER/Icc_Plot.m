


load('D:\MEMR_GroupData\MEMR_groupData.mat')
D = thdOffsetLvl;
DD = D(:);
medd = median(D,2);
minn = min(DD)
maxx = max(DD)
midx = median(DD)
[~,indx] = sort(medd);
D = D(indx,:);
x = (1:1:30)';
X = repmat(x,1,4);

d = D(:);
PD = fitdist(d,'kernel');

Icdf = PD.icdf(x/30);
%xx = linspace(-.1,.5,1000)';
xx = linspace(40,110,1000)';
Pdf = PD.pdf(xx);

% h1 = figure;
% sp1 = subplot(2,1,1);
% sp2 = subplot(2,1,2);
% sp1.Position = [0.1300    0.1466    0.7750    0.7784];
% sp2.Position = [0.2076    0.5372    0.2690    0.3412];
% axes(sp1)
% plot(X,D,'.-','Color',[.7 .7 1])
% hold on
% plot(x,Icdf,'-b','LineWidth',2)
% axes(sp2)
% plot(xx,Pdf,'-b','LineWidth',2)
% hold on

% -------------
load('D:\MEMR_GroupData\MEMR_groupData.mat')
D = thdOnsetLvl;
DD = D(:);
medd = median(D,2);
minn = min(DD)
maxx = max(DD)
midx = median(DD)
[~,indx] = sort(medd);
D = D(indx,:);
x = (1:1:30)';
X = repmat(x,1,4);

d = D(:);
PD = fitdist(d,'kernel');

Icdf = PD.icdf(x/30);
%xx = linspace(-.1,.5,1000)';
Pdf = PD.pdf(xx);

% axes(sp1)
% plot(X,D,'.-','Color',[1 .7 .7])
% hold on
% plot(x,Icdf,'-r','LineWidth',2)
% axes(sp2)
% plot(xx,Pdf,'-r','LineWidth',2)




% ---------
load('D:\MEMR_GroupData\MEMR_groupData.mat')
D = delay;
DD = D(:);
medd = median(D,2);
minn = min(DD)
maxx = max(DD)
midx = median(DD)
[~,indx] = sort(medd);
D = D(indx,:);
x = (1:1:30)';
X = repmat(x,1,4);

d = D(:);
PD = fitdist(d,'kernel');

Icdf = PD.icdf(x/30);
%xx = linspace(-.1,.5,1000)';
Pdf = PD.pdf(xx);


load('D:\MEMR_GroupData\MEMR_groupData.mat')
D1 = thdOffsetLvl;
D2 = thdOnsetLvl;
DD1 = D1(:)
DD2 = D2(:)

medd = median(DD1-DD2);
minn = 20*log10(min(DD))
maxx = 20*log10(max(DD))
midx = 20*log10(median(DD))
[~,indx] = sort(medd);
D = D(indx,:);
x = (1:1:30)';
X = repmat(x,1,4);

d = D(:);
PD = fitdist(d,'kernel');

Icdf = PD.icdf(x/30);
%xx = linspace(-.1,.5,1000)';
Pdf = PD.pdf(xx);
xx = linspace(-.1,.5,1000)';
Pdf = PD.pdf(xx);
Cdf = PD.cdf(xx);
Icdf = PD.icdf(xx);


load('D:\MEMR_GroupData\MEMR_groupData.mat')
D = pAmp;
DD = D(:);
medd = median(D,2);
minn = 20*log10(min(DD))
maxx = 20*log10(max(DD))
midx = 20*log10(median(DD))
[~,indx] = sort(medd);
D = D(indx,:);
x = (1:1:30)';
X = repmat(x,1,4);

d = D(:);
PD = fitdist(d,'kernel');

Icdf = PD.icdf(x/30);
%xx = linspace(-.1,.5,1000)';
Pdf = PD.pdf(xx);