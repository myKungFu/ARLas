function [e1,e2,E,e,a1,a2,A,a,time] = abrAnalysis_v1(header,V1,V2,A1,A2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% input should be as follows: 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    if nargin == 0
        path = 'C:\myWork\ARLas\Data\Shawnec2\Shawnec2_13OCT2023\Shawnec2_13OCT2023_003_run1\';% ear canal electrode

        %path = 'C:\myWork\ARLas\Data\AbrShawn\AbrShawn_13OCT2023\AbrShawn_13OCT2023_003_run1\'; % simple abr
         
        d = dir([path,'*optiAmp*']);
        dummy = load(d(1).name);
        name1 = d(1).name;
        V1 = dummy.data;
        header = dummy.header;
        try
            polarity = header.userInfo.polarity;
            offset = header.userInfo.offset;
        catch ME
            polarity = 'cond';
            
        end        
        V1 = unpackPolarity(polarity,V1);

        if size(d,1) == 2
            dummy = load(d(2).name);
            V2 = dummy.data;
            V2 = unpackPolarity(polarity,V2);
            name2 = d(2).name;
        else
            V2 = [];
        end        
        d = dir([path,'*10xA_eeg*']);
        dummy = load(d(1).name);
        A1 = dummy.data;
        A1 = unpackPolarity(polarity,A1);
    end
    
    location = 1;
    [e1c,e1r,E1,e1,time] = analyze(V1,offset,location); % analyze the voltage waveforms
    if ~isempty(V2)
        [e2c,e2r,E2,e2,time] = analyze(V2,offset,location); % analyze the voltage waveforms
    end
    location = 2;
    [a1,a2,A,a] = analyze(A1,offset,location); % analyze the aoustic waveforms
    
%Figure
%figure
plot(time,e1,'m')
hold on
plot(time,e2,'c')
% hold on
% plot(time,e2,'b')
xlim([time(1),time(end)])
xlabel('Time (ms)')
ylabel('Amplitude (\muV)')
grid on
%title 'R:Main G:E1 B:E2'
% ylim([1e-7, max([e,a])*1.1])
legend(name1,name2)

keyboard

% two polarities
% time window
% two buffers


end


% INTERNAL FUNCTIONS --------------------------------------------------------------------------
function [x1,x2,X,x,time] = analyze(X,offset,location)
X = X*10^6;
% bandpass filter
% fmin = 30;
fmax = 1500;
fs = 96000;
OctaveStep = 1/3;
b = getLPF(fmax,OctaveStep);

% plot(mean(X,2))
% hold on

X = fastFilter(b,X);

% plot(mean(X,2))
% 
% figure
% subplot(2,1,1)
% plot(b)
% subplot(2,1,2)
% plot(20*log10(abs(fft(b,fs))));
% xlim([1,fs/2])


%----------- EEG Delay Correction --------%


% time window
tmin = -1; % minimum time in ms
tmax = 15; % maximum time in ms
N = size(X,1);
time = (0:1:N-1)'/fs * 1000;
time = time - 5;
[~,indx1] = min(abs(time-tmin));
[~,indx2] = min(abs(time-tmax));
X = X(indx1:indx2 + offset,:);
time = time(indx1:indx2);

if location == 1 % for electric, take off the front
    X = X(offset+1:end,:);
elseif location == 2
    X = X(1:end-offset,:);
else
    error('Location must be 1 or 2')
end

% artifact rejection
rms = sqrt(mean(X.^2,1));
doPlot = 0;
tolerance = 'mild';
[indx,nRejects] = AR(rms,tolerance,doPlot);
X = AR_engine(X,indx); 

M = size(X,2);
indx1 = (1:2:M);
indx2 = (2:2:M);
X1 = X(:,indx1);
X2 = X(:,indx2);
x = mean(X,2);
x1 = mean(X1,2);
x2 = mean(X2,2);
end

function b = getLPF(f1,step)
Fs = 96000;  % Sampling Frequency

Fpass = f1;            % Passband Frequency
Fstop = round(2^(log2(f1)+step));           % Stopband Frequency
Dpass = 0.057501127785;  % Passband Ripple
Dstop = 0.0001;          % Stopband Attenuation
dens  = 16;              % Density Factor

% Calculate the order from the parameters using FIRPMORD.
[N, Fo, Ao, W] = firpmord([Fpass, Fstop]/(Fs/2), [1 0], [Dpass, Dstop]);

% Calculate the coefficients using the FIRPM function.
if mod(N,2)~=0
    N=N+1;
end
b  = firpm(N, Fo, Ao, W, {dens})';
end

function [X] = unpackPolarity(polarity,X)
    if ~strcmp(polarity,'alt')
        return
    end
    
    [rows,cols] = size(X);
    % force columns to be even
    if mod(cols,2)~=0
        X(:,end) = [];
        [rows,cols] = size(X);
    end
    X = reshape(X,rows*2,cols/2);
    x1 = X(1:rows,:);
    x2 = X(rows+1:end,:);
    X = x1 + x2;
end

