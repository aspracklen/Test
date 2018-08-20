function [ output_args ] = OSTAPApply( Radar,OSTAPWeights,TestSig,CovEst )
%UNTITLED14 Summary of this function goes here
%   Detailed explanation goes here
Data=zeros(length(Radar.range),length(OSTAPWeights));
for k=1:length(Radar.range);

Data(k,:)=abs(OSTAPWeights'*TestSig(:,k)).^2./diag(real(OSTAPWeights'*CovEst*OSTAPWeights));

k
end;
[Vel,RangeGates]=meshgrid(linspace(- Radar.PRF/4*Radar.lambda, Radar.PRF/4*Radar.lambda,4*Radar.M+1),linspace(1, length(Radar.range),length(Radar.range)));
pcolor(Vel,RangeGates,10.*log10(abs(Data)))
xlabel('Velocity(m/s)')
ylabel('Range Gates')
colormap jet
zlim([0,max(10.*log10(abs(Data(:))))+5])
caxis([0,max(10.*log10(abs(Data(:))))+5])
shading interp
end

