function [ Radar ] = Setup( Radar )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
Radar.T=1/Radar.PRF;
Radar.f0=10e9;
Radar.Pt=1000;
Radar.B=10e6;
Radar.c=3e8;
Radar.Nc=1000;
Radar.lambda=Radar.c/Radar.f0;
Radar.Nx=64;
Radar.Ny=12;
AntA=Radar.Nx*Radar.lambda/2*Radar.Ny*Radar.lambda/2;
Radar.Ae=AntA./(Radar.Nx*Radar.Ny);
Radar.sigmaTgt=10;

%Thermal Noise Power Computations
k = 1.3806488e-23;            % Boltzmann Constant in J/K.
To = 290;                     % Standard room Temperature in Kelvin.
F = 7;                        % Receiver Noise Figure in dB;
Te = To*(10^(F/10)-1);        % Effective Receiver Temperature in Kelvin.
Lr = 0;                    % System Losses on receive in dB.
Ts = 10^(Lr/10)*Te;           % Reception System Noise Temperature in Kelvin.
Nn = k*Ts;                    % Receiver Noise PSD in Watts/Hz.
Radar.Pn = Nn*Radar.B;                    % Receiver Noise Power in Watts

%%Clutter 
Radar.phi=linspace(0,180,Radar.Nc);        %Clutter Angles in degrees    
Radar.ha=9e3;
Radar.Rcik = 20*1852;                % (clutter) range of interest in meters.
dphi = 2*pi/Radar.Nc;               % Azimuth angle increment in rad.
Radar.dR = Radar.c/2/Radar.B;                   % Radar Range Resolution in meters.
Radar.theta0 = asind(Radar.ha/Radar.Rcik);          % Grazing angle at the clutter patch in rad (flat earth model).
Radar.gamma = 10^(-20/10);           % Terrain-dependent reflectivity factor.
%

%%Geometry
Radar.range=19*1852:Radar.dR:21*1852;   %Range of interest
Radar.theta=asind(Radar.ha./Radar.range);     %Range of theta values;
[Radar.Theta,Radar.Phi]=meshgrid(Radar.theta,Radar.phi);
Radar.v=sind(Radar.Theta);
Radar.u=cosd(Radar.Theta).*cosd(Radar.Phi);
Radar.v0=Radar.ha/Radar.Rcik;
Radar.u0=sqrt(1-(Radar.ha/Radar.Rcik)^2)*0;
Radar.Range=Radar.ha./sind(Radar.Theta);
%%Gain
Radar.G=4*pi/3159*Radar.Nx^2*Radar.Ny^2*(diric(pi*(Radar.u-Radar.u0),Radar.Nx).*diric(pi*(Radar.v-Radar.v0),Radar.Ny)).^2;
PatchArea = Radar.Range.*dphi*Radar.dR; %%1/cos(theta)omitted
sigma0 = Radar.gamma*sind(Radar.Theta);
Radar.sigma = sigma0.*PatchArea;

PowerR=240*Radar.Pt*Radar.Ae/((4*pi)^2).*Radar.G.*Radar.sigma./Radar.Range.^4;
Radar.VoltageR=sqrt(PowerR);
%
XEl=linspace(-32,31,64)*Radar.lambda/2;
ZEl=linspace(-6,5,12)*Radar.lambda/2;
[Radar.X,Radar.Z]=meshgrid(XEl,ZEl);
Radar.Pulses=linspace(0,Radar.M-1,Radar.M);
%Spatial Combiner
clc

Radar.TMat=zeros(Radar.Nx*Radar.Ny,Radar.ChannelNum);
[Hannx,Hanny]=meshgrid(hann(Radar.Nx),hann(Radar.Ny));
Hann=Hannx.*Hanny;
XF=Radar.X(:);
ZF=Radar.Z(:);
HannF=Hann(:);
for i=1:Radar.ChannelNum;
Select=Radar.Nx*Radar.Ny/Radar.ChannelNum*(i-1)+1:Radar.Nx*Radar.Ny/Radar.ChannelNum*(i-1)+Radar.Nx*Radar.Ny/Radar.ChannelNum;
Radar.TMat(Select,i)=exp(1i*2*pi/Radar.lambda*(XF(Select).*cosd(90)*cosd(Radar.theta0)-ZF(Select).*sind(Radar.theta0))).*HannF(Select);
end;
Radar.TMat=1/sqrt(Radar.Nx*Radar.Ny/Radar.ChannelNum).*Radar.TMat;
end

