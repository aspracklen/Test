clc;
clear;
fs=28.8e9;
ts=1/fs;
tp=20e-6;
t=-tp/2:ts:tp/2-ts;
Freq=linspace(-fs/2,fs/2-fs/(2*length(t)),length(t));
B=2e9;
alpha=B/(2*tp);
omega0=2*pi*9.35e9;
S=cos(1.*(1*omega0*t+2*pi*alpha*t.^2));
tdelay=0.1/3e8;
D=cos(1.*(1*omega0*(t-tdelay)+2*pi*alpha*(t-tdelay).^2));
SFFT=fftshift(fft(S));
DFFT=fftshift(fft(D));
Chirp0=Freq<9.35e9+1.5e9/2;
Chirp1=Freq>9.35e9-1.5e9/2;
Chirp=Chirp0.*Chirp1;
FilterOut0=Freq<11e9;
FilterOut1=Freq>7.7e9;
FilterOut=FilterOut0.*FilterOut1;
ChirpPos=find(Chirp);
[~,FreqDelta]=sort(abs(Freq-9.35e9));
FiltFreqS=exp(1i.*tp*1*0.0001*(Freq-9.35e9));
FiltFreqD=exp(1i.*tp*1*0.00025*(Freq-9.35e9));
%FiltFreqS=exp(1i.*tp*1*0.00025*(Freq-9.35e9)+1i.*(tp*0*0.0000001^2*(Freq-9.4e9).^2+tp*1*0.00000007^3*(Freq-9.2e9).^3));
SF=ifft(ifftshift(SFFT.*FilterOut));
DF=ifft(ifftshift(DFFT.*FilterOut));
% figure(1)
% subplot(121)
% plot(Freq,20*log10(abs(SFFT)),Freq,abs(FilterOut))
% xlim([7.5e9,11.2e9])
% subplot(122)
% plot(Freq(ChirpPos),angle(FiltFreqS(ChirpPos)),Freq(ChirpPos),angle(FiltFreqD(ChirpPos)))
% xlim([7.5e9,11.2e9])
%
clc
omega=2*pi*7.55e9;
Data=zeros(20,1);
Delay=linspace(-0.5,0.5,20);
for i=1:20;
tdelay=Delay(i)*tp;
SLO=cos(1.*(1*omega*(t-tdelay)+2*pi*alpha*(t-tdelay).^2));
SDe=SLO.*SF+0*randn(1,length(SF));
DDe=SLO.*DF+0*randn(1,length(SF));

SDeFFT=fftshift(fft(1/sqrt(length(t)).*SDe));
DDeFFT=fftshift(fft(1/sqrt(length(t)).*DDe));
% figure(2)
% subplot(121)
% plot(Freq,20.*log10(abs(SDeFFT)))
% subplot(122)
% plot(Freq,20.*log10(abs(DDeFFT)))
[~,maxS]=max(abs(SDeFFT));
[~,maxD]=max(abs(DDeFFT));
%[mod(angle(SDeFFT(maxS)./DDeFFT(maxD)),2*pi),mod(mean(unwrap(angle(FiltFreqS(ChirpPos))))-mean(unwrap(angle(FiltFreqD(ChirpPos)))),2*pi)]



%
bpFilt = designfilt('bandpassfir','FilterOrder',20, ...
         'CutoffFrequency1',1.8e9-375e6,'CutoffFrequency2',1.8e9+375e6, ...
         'SampleRate',fs);
lpFilt = designfilt('bandpassfir','FilterOrder',20, ...
         'CutoffFrequency1',0.001,'CutoffFrequency2',375e6, ...
         'SampleRate',fs);
SDeF=filter(bpFilt,SDe);
DDeF=filter(bpFilt,DDe);
dec=12;
tdec=t(1:dec:end);
Freqdec=linspace(-fs/(2*dec),fs/(2*dec)-fs/(length(tdec)),length(tdec));

omegaD=2*pi*1.8e9;
DLO=exp(-1i.*(1*omegaD*tdec));
SB=SDeF(1:dec:end).*DLO;
DB=DDeF(1:dec:end).*DLO;
%
SBFFT=fftshift(fft(1/sqrt(length(tdec)).*SB));
DBFFT=fftshift(fft(1/sqrt(length(tdec)).*DB));
% subplot(121)
% plot(Freqdec,20.*log10(abs(SBFFT)))
% subplot(122)
% plot(Freqdec,20.*log10(abs(DBFFT)))
[~,maxS]=max(abs(SBFFT));
[~,maxD]=max(abs(DBFFT));
Data(i)=angle(SBFFT(maxS)/DBFFT(maxD));
i
end;
hold on
plot(Delay,Data,'*')