function [ output_args ] = SumPlot( Radar,TestSig )
%UNTITLED15 Summary of this function goes here
%   Detailed explanation goes here
SumData=zeros(length(Radar.range),Radar.M);
for i=1:length(Radar.range);
SumData(i,:)=sum(reshape(TestSig(:,i),Radar.ChannelNum,Radar.M));
end;
SumData=SumData.';
FFT=fftshift(fft(SumData,128),1);
[Vel,RangeGates]=meshgrid(linspace(- Radar.PRF/4*Radar.lambda, Radar.PRF/4*Radar.lambda,128),linspace(1, length(Radar.range),length(Radar.range)));
pcolor(Vel,RangeGates,10.*log10(abs(FFT.')))
xlabel('Velocity(m/s)')
ylabel('Range Gates')
colormap jet
% zlim([-50,max(10.*log10(abs(FFT(:))))+5])
% caxis([-50,max(10.*log10(abs(FFT(:))))+5])
shading interp
end

