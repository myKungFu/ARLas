maxIterations = 50;
nn = floor(size(In,1)/10000);
DP = [];
chunk = 10000; % chunk for memory reasons
start = 1;
finish = chunk;
for ii=1:nn
    dp = bisquareWeights(In(start:finish,:),maxIterations);
    DP = [DP;dp];
    start = start + chunk;
    finish = finish + chunk;
end
dp = bisquareWeights(In(start:end,:),maxIterations);
DP = [DP;dp];
dp = DP; clear DP