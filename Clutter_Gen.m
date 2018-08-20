function [ SigC ] = Clutter_Gen( Radar )
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

SigC=zeros(Radar.M*Radar.ChannelNum,length(Radar.range));
for i=1:length(Radar.range);
for k=1:Radar.Nc
fsp=exp(1i*2*pi/Radar.lambda*(Radar.X.*cosd(Radar.Phi(k,i))*cosd(Radar.Theta(k,i))-Radar.Z.*sind(Radar.Theta(k,i))));
fd=exp(1i*2*pi/Radar.lambda*2*Radar.vp*Radar.Pulses*Radar.T*cosd(Radar.Phi(k,i)).*cosd(Radar.Theta(k,i)));
Vc=exp(1i.*2*pi*rand(1,1))*Radar.VoltageR(k,i).*kron(fd(:),Radar.TMat'*fsp(:));
SigC(:,i)=SigC(:,i)+Vc;
end
i
end;
Noise=sqrt(Radar.Pn/2)*(randn(Radar.M*Radar.ChannelNum,length(Radar.range))+1i.*randn(Radar.M*Radar.ChannelNum,length(Radar.range)));
SigC=SigC+Noise;
end

