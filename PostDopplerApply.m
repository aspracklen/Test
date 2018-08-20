function [ output_args ] = PostDopplerApply( Radar,DopplerWeights,TestSig,CovEst )
%UNTITLED12 Summary of this function goes here
%   Detailed explanation goes here
Data=zeros(length(Radar.range),Radar.M);
for k=1:length(Radar.range);

Data(k,:)=abs(DopplerWeights'*TestSig(:,k)).^2./diag(real(DopplerWeights'*CovEst*DopplerWeights));
%Data(k,:)=abs(Weights'*TestSig(:,k)).^2;
k
end;
[Vel,RangeGates]=meshgrid(linspace(- Radar.PRF/4*Radar.lambda, Radar.PRF/4*Radar.lambda,Radar.M),linspace(1, length(Radar.range),length(Radar.range)));
pcolor(Vel,RangeGates,10.*log10(abs(Data)))
xlabel('Velocity(m/s)')
ylabel('Range Gates')
colormap jet
zlim([0,max(10.*log10(abs(Data(:))))+5])
caxis([0,max(10.*log10(abs(Data(:))))+5])
shading interp
end

