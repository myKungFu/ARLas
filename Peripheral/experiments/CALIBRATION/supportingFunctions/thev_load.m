function [Zo,Zw,Zx,Pf,Pr,Gs,Ze,ew]=thev_load(Zs,Ps,Pl,fdw)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% thev_load()
%
% Input Arguments:
%  Zs - source impedance
%  Ps - source pressure
%  Pl - measured load pressure
%  fdw - bandwidth / Nyquist_frequency
%
% Output Arguments:
% Zo - surge impedance
% Zw - waveguide impeance
% Zx - 
% Pf -
% Pr - 
% Gs - 
% Ze - evanescent impedance
% ew - 
%
% Calculate load impedance
% From Steve Neely to Shawn Goodman
% May 20, 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Estimate surge impedance Zx
    Zl=Zs.*Pl./(Ps-Pl); % measured load impedance
    [Zo,ew,Zx]=surge_imped(Zl,fdw); % surge impedance
    % Decompose Pl into Pf & Pr
    Gs=(Zs-Zo)./(Zs+Zo); % source reflectance
    Gw=((Zl-Zx)./(Zl+conj(Zx))); % waveguide reflectance
    Po=Ps.*(1-Gs)/2; % incident pressure
    Pf=Po./(1-Gs.*Gw); % forward pressure
    Pr=Pf.*Gw; % reflected pressure
    Zw=Zo*(1+Gw)./(1-Gw); % waveguide impedance
    Ze=Zx-Zo; % evanescent impedance
end


% Internal functions ------------------------------------------------------
function [Zo,Zewc,Zx]=surge_imped(Zl,fdw)
    % Zl - load impedance
    % fdw - usable bandwidth / Nyquist_frequency
    nf=length(Zl); % number of frequencies
    w=2*pi*(0:(nf-1))'/(nf-1); % normalized radian frequencies
    s=1i*w; % Laplace-transform variable
    % calculate Fourier-transform vectors
    v0=fd_window(ones(size(Zl)),1,fdw);
    v1=fd_window(sin(w/2),1,fdw);
    v5=fd_window(sin(w/10),1,fdw);
    v0=v0/sum(v0);
    v1=v1/sum(v1);
    v5=v5/sum(v5);
    % specify initial coefficients
    Zo=sum(real(Zl(2:end)).*v0(2:end));
    Z1=abs(sum(imag(Zl).*w)/sum(w.*w))/3;
    Z5=Zo*4e-4;
    % iterate to improve coefficients
    Zx=Zo+Z1*s+Z5*s.^5; % initial surge impedance
    for k=1:80
        Gw=(Zl-Zx)./(Zl+conj(Zx)); % waveguide reflectance
        Zo=Zo*(1+dot(real(Gw),v0));
        Z1=Z1*(1+dot(imag(Gw),v1));
        Z5=Z5*(1+dot(imag(Gw),v5));
        Zx=Zo+Z1*s+Z5*s.^5; % improved surge impedance
    end
    Zewc=[(Z1/Zo) (Z5/Zo)^0.2]; % evanescent-wave coefficients
end

function H=fd_window(H,srf,fdw)
    % H - transfer function
    % srf - sampling-rate factor
    % fdw - bandwidth / Nyquist_frequency
    if (fdw<=0), return; end
    [n,m]=size(H);
    p=pi*(0:(n-1))'/(n-1)/fdw;
    a=0.16; % Blackman window
    w=(1-a+cos(p)+a*cos(2*p))/2;
    w(p>pi)=0;
    Z=zeros((srf-1)*(n-1),m); % zero padding
    H=[H.*repmat(w,1,m);Z];
end