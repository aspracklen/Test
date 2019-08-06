clear;
clc;
N=40;
lambda=0.03;
x=(-N/2:N/2-1)*lambda/2;
Win=taylorwin(N,5,-35)./max(taylorwin(N,5,-35));
M=40;
T=zeros(N,M);
for i=1:M;
T((i-1)*N/M+1:(i-1)*N/M+N/M,i)=Win((i-1)*N/M+1:(i-1)*N/M+N/M);
end;
Norm=sqrt(diag(T'*T));
T=T./Norm';

d0=x./Win.';
d0=([d0(1:N/2),-fliplr(d0(1:N/2))].');
d=sqrt(1.1)/sqrt(d0'*d0).*d0;
w=ones(N,1);
% diffconv=[Norm(1:32,1);-Norm(1:32,1)];
% diffconv=[Norm(1:2,1);-Norm(1:2,1)];
wdiff=inv(T'*T)*T'*(eye(N)-1/(w'*T*inv(T'*T)*T'*w)*(w*w'*T*inv(T'*T)*T'))*d;
save('wdiff.mat','wdiff')
%%
wSA=ones(M,1);
T1=Norm.'.*T;
Noise=1/sqrt(2)*(randn(N,10000)+1i.*randn(N,10000));
SA=T1'*Noise;

Sum=wSA'*SA./sqrt(sum(Norm.^2));


Diff=wdiff'*SA./sqrt(sum(abs(wdiff.*Norm).^2));
PowerSum=mean(abs(Sum).^2)
PowerDiff=mean(abs(Diff).^2)
%%
theta=linspace(-11.5,11.5,361);
Signal=exp(1i.*2*pi/lambda*x.'.*sind(theta));

SA=T1'*Signal;

Sum=wSA'*SA./sqrt(sum(Norm.^2)*7.1);
Diff=wdiff'*SA;
diffconv=[ones(2,1);-ones(2,1)];
DiffConv=diffconv'*SA./sqrt(sum(diffconv.^2));
plot(sind(theta),10*log10(abs(Sum).^2)-14.56,sind(theta),10*log10(abs(Diff).^2)-14.56,sind(theta),10*log10(abs(DiffConv).^2)-14.56)
ylim([-50,20])
%%
Gd0=[Diff;DiffConv];
Gd=max(abs(Gd0).^2);
GdSA=abs(SA(1,:)).^2;
plot(theta,10*log10(abs(Sum).^2),theta,10*log10(Gd),theta,10*log10(GdSA))
ylim([-50,20])