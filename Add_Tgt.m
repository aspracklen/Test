function [ TestSig ] = Add_Tgt( Radar,TestSig )
%UNTITLED9 Summary of this function goes here
%   Detailed explanation goes here
GTgt=4*pi/3159*Radar.Nx^2*Radar.Ny^2*(diric(pi*(Radar.u0-Radar.u0),Radar.Nx).*diric(pi*(Radar.v0-Radar.v0),Radar.Ny)).^2;
PowerTgt=240*Radar.Pt*Radar.Ae/((4*pi)^2).*GTgt.*Radar.sigmaTgt./Radar.Rcik.^4;
VoltageTgt=sqrt(PowerTgt);

fsptgt=exp(1i*2*pi/Radar.lambda*(Radar.X.*cosd(90)*cosd(Radar.theta0)-Radar.Z.*sind(Radar.theta0)));
fdtgt=exp(1i*2*pi/Radar.lambda*2*Radar.vtgt*Radar.Pulses*Radar.T);
Vtgt=VoltageTgt.*kron(fdtgt(:),Radar.TMat'*fsptgt(:));

TestSig(:,round(length(Radar.range)/2))=TestSig(:,round(length(Radar.range)/2))+Vtgt;
TestSig(:,round(length(Radar.range)/2)+1)=TestSig(:,round(length(Radar.range)/2)+1)+Vtgt;
TestSig(:,round(length(Radar.range)/2)-1)=TestSig(:,round(length(Radar.range)/2)-1)+Vtgt;

end

