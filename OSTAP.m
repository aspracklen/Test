function [ Weights ] = OSTAP( Radar,InvRu )
%UNTITLED8 Summary of this function goes here
%   Detailed explanation goes here
vTest=linspace(-Radar.PRF/4*Radar.lambda,Radar.PRF/4*Radar.lambda,4*Radar.M+1);
Weights=zeros(Radar.M*Radar.ChannelNum,4*Radar.M+1);
        for i=1:length(vTest)
            at=exp(1i*2*pi/Radar.lambda*(Radar.X.*cosd(90)*cosd(Radar.theta0)-Radar.Z.*sind(Radar.theta0)));
            bt = exp(1i*2*pi/Radar.lambda*2*vTest(i)*Radar.Pulses*Radar.T); % Dummy Target Doppler Steering Vector
            vt = kron(bt(:),Radar.TMat'*at(:));
            w = InvRu*vt;
Weights(:,i)=w;

i
       end

end

