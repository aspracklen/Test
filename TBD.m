clear all
clc
c = 3e8; % speed of light
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% USER INPUT SECTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

M = 10001; % slow-time dimension, samples; keep it odd.

% Lref is the range bin # of the target, on a 0:L-1 scale, at the center of
% the CPI (middle pulse)
% Lref = round(L/2); % puts target at middle range bin on the middle pulse

F0 = 10e9; % RF (Hz)
B = 1e6; % waveform bandwidth (Hz)

% sampling intervals and rates
Fsft = 5e6;
lambda=0.03;
PRF = 120e3;
L = floor(1/PRF*Fsft); % fast time dimension, samples
Lref = round(L/2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% END USER INPUT SECTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Some derived parameters
m_end = (M-1)/2;
ms = (-m_end:m_end); % slow time index labels

Tft = 1/Fsft; % fast time sampling interval
dr = c*Tft/2; % range bin spacing
Tst = 1/PRF; % slow-time sampling interval (PRI)
Dopv=ms*PRF/M*lambda/2;
v = Dopv(8335); % velocity in m/s towards the radar
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SNR=10;
SNRVal=10^((SNR-10*log10(10001))/10);
del_phi = -4*pi*(F0/c)*v*Tst; % pulse-to-pulse phase increment due to range change
for f=1:10;
for k=1:2
y=sqrt(SNRVal).*exp(1i*del_phi*ms.').*sinc( B*Tft*((0:L-1)-Lref) );
Noise=1/sqrt(2)*(randn(M,L)+1i.*randn(M,L));
y=y+Noise;
yFFT=fftshift(1/sqrt(M).*fft(y),1);
WinRef=50;
WinTest=1;
Win=[ones(WinRef,1);zeros(WinTest,1);ones(WinRef,1)];
Win=Win./sum(Win);

PFa=5e-5;
PFaTBD=5e-5;

Bins=100;
alpha=Bins*(PFa^(-1/Bins)-1);
alphaTBD=Bins*(PFaTBD^(-1/Bins)-1);
Threshold=zeros(size(y));
for i=1:L;
cfar=cconv(abs(yFFT(:,i)).^2,Win,M);
Threshold(:,i)=cfar;
end
PreAlarms(k,:,:)=(abs(yFFT).^2-alpha*Threshold)>0;
PreAlarmsTBD(k,:,:)=(abs(yFFT).^2-alphaTBD*Threshold)>0;
end;
Alarms(f,:,:)=squeeze(sum(PreAlarms))>=2;
AlarmsTBD(f,:,:)=squeeze(sum(PreAlarmsTBD))>=2;

f
end;

DetMap=zeros(M,10);
for i=1:10
[row,col]=find(squeeze(Alarms(i,:,:)));
DetMap(row,i)=1;
end;


TBDMap=zeros(M,10);
for i=1:10
[row,col]=find(squeeze(AlarmsTBD(i,:,:)));
TBDMap(row,i)=1;
end;
theta = 60:120;
[R,xp] = radon(TBDMap,theta);
sum(DetMap(1667,:))/10
Med=medfilt1(R.',5);

mesh(Med)

