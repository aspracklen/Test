clear;
clc
Img=imread('Phase Taper.png');
Img=double(Img(:,:,1))<255;
% imagesc(Img)
% colormap('gray')
ImageData=Img.';
ImageData=ImageData(28:432,13:280);
[row,col]=find(ImageData);
[c,ia,ic]=unique(row);
Phase=linspace(-90,90,size(ImageData,2));
PhaseTaper=Phase(col(ia));
PhaseSym=[-fliplr(PhaseTaper(203:end)),PhaseTaper(202:end)];
plot(PhaseSym)
%%
lambda=0.03;
Xelem=length(PhaseSym);
x=(-(Xelem-1)/2:(Xelem-1)/2)*0.03/2;
T=exp(1i.*PhaseSym*pi/180*0.5).';
T0=exp(1i.*PhaseSym*pi/180*0).';
theta=linspace(-90,90,1001);
Sig=exp(1i.*2*pi/lambda*x.*sind(theta).').';
SigOut=abs(T'*Sig).^2;
SigRef=abs(T0'*Sig).^2;
plot(theta,10*log10(SigOut./max(SigRef)),theta,10*log10(SigRef./max(SigRef)))
max(10*log10(SigOut./max(SigRef)))
%%
z=linspace(-pi,pi,41)


PhaseSymSmall=PhaseSym(1:10:end);
lambda=0.03;
Xelem=length(PhaseSymSmall);
x=(-(Xelem-1)/2:(Xelem-1)/2)*0.03/2;
T=(exp(-1i.*2*pi/lambda*x.*sind(-30+1.3)).*exp(1i.*PhaseSymSmall*pi/180*0.7)).';
T0=exp(1i.*2*pi/lambda*x'.*sind(30).').*exp(1i.*PhaseSymSmall*pi/180*0).';
TTaylor=taylorwin(length(PhaseSymSmall));
TTaylor=(TTaylor./max(TTaylor)).*(exp(1i.*2*pi/lambda*x'.*sind(30)));
theta=linspace(0,90,1001);
Sig=exp(1i.*2*pi/lambda*x.*sind(theta).').';
SigOut=abs(T'*Sig).^2;
SigRef=abs(T0'*Sig).^2;
SigTaylor=abs(TTaylor'*Sig).^2;
SigComb=SigOut.*SigTaylor;
SigCombRef=SigRef.*SigTaylor;
% SigUltRef=SigTaylor.*SigTaylor;
 plot(theta,10*log10(SigRef./max(SigRef)),theta,10*log10(SigOut./max(SigOut)))
% ylim([-90,0])
% max(10*log10(SigOut./max(SigRef)))
% max(10*log10(SigOut(525:end)./max(SigOut)))
%%
plot(theta,10*log10(SigComb./max(SigComb)),theta,10*log10(SigCombRef./max(SigCombRef)),theta,10*log10(SigUltRef./max(SigUltRef)))

%%
%2D
PhaseSymSmall=PhaseSym(1:10:end);
lambda=0.03;
Xelem=length(PhaseSymSmall);
Yelem=length(PhaseSymSmall);
x=(-(Xelem-1)/2:(Xelem-1)/2)*0.03/2;
y=(-(Yelem-1)/2:(Yelem-1)/2)*0.03/2;
[X,Y]=meshgrid(x,y);
[U,V]=meshgrid(linspace(-1,1,801),linspace(-1,1,801));
% U=U(:);
% V=V(:);
T0=ones(41*41,1);
for i=1:801;
for j=1:801;
Sig=exp(1i.*2*pi/lambda.*X(:).*U(i,j)+1i.*2*pi/lambda.*Y(:).*V(i,j));
Data0(i,j)=T0'*Sig;
[i,j]
end;
end
%%


%%
PhaseSymSmall=PhaseSym(1:10:end);
lambda=0.03;
Xelem=length(PhaseSymSmall);
Yelem=length(PhaseSymSmall);
x=(-(Xelem-1)/2:(Xelem-1)/2)*0.03/2;
y=(-(Yelem-1)/2:(Yelem-1)/2)*0.03/2;
[X,Y]=meshgrid(x,y);
[U,V]=meshgrid(linspace(-1,1,801),linspace(-1,1,801));
T2D=zeros(41^2,1);
for i=1:41;
T2D((i-1)*41+1:(i-1)*41+41)=exp(1i.*PhaseSymSmall*pi/180);
end;

for i=1:801;
for j=1:801;
Sig=exp(1i.*2*pi/lambda.*X(:).*U(i,j)+1i.*2*pi/lambda.*Y(:).*V(i,j));
Data(i,j)=T2D'*Sig;
[i,j]
end;
end

%%
DataPower=abs(Data).^2;
mesh(10*log10(DataPower./max(DataPower(:))))
%%
PhaseSymSmall=PhaseSym(1:10:end);
lambda=0.03;
Xelem=length(PhaseSymSmall);
Yelem=length(PhaseSymSmall);
x=(-(Xelem-1)/2:(Xelem-1)/2)*0.03/2;
y=(-(Yelem-1)/2:(Yelem-1)/2)*0.03/2;
[X,Y]=meshgrid(x,y);
[U,V]=meshgrid(linspace(-1,1,801),linspace(-1,1,801));
T2D1=zeros(41^2,1);
for i=1:41;
T2D1((i-1)*41+1:(i-1)*41+41)=exp(1i.*PhaseSymSmall*pi/180);
end;
Test=reshape(T2D1,41,41);
Test1=Test.*Test';
Test1=Test1(:);


%%
for i=1:801;
for j=1:801;
Sig=exp(1i.*2*pi/lambda.*X(:).*U(i,j)+1i.*2*pi/lambda.*Y(:).*V(i,j));
Data1(i,j)=Test1'*Sig;
[i,j]
end;
end
%%
Data0Power=abs(Data0).^2;
DataPower=abs(Data1).^2;
mesh(10*log10(DataPower./max(Data0Power(:))))
max(10*log10(DataPower(:)./max(Data0Power(:))))
%%
phi=linspace(0,2*pi,1001);
theta=linspace(0,pi/2,1001);
[Phi,Theta]=meshgrid(phi,theta);
dphi=phi(2)-phi(1);
dtheta=theta(2)-theta(1);
sum(sin(Theta(:)).*dphi.*dtheta)
%%
for i=1:length(phi);
for j=1:length(theta);
Sig=exp(1i.*2*pi/lambda.*X(:).*sin(Theta(i,j))*cos(Phi(i,j))+1i.*2*pi/lambda.*Y(:).*sin(Theta(i,j))*sin(Phi(i,j)));
Datathetaphi(i,j)=Test1'*Sig;
[i,j]
end;
end
%%
PowerNorm=sum(sum(abs(Datathetaphi).^2.*sin(Theta)*dphi*dtheta))
%%
for i=1:length(phi);
for j=1:length(theta);
Sig=exp(1i.*2*pi/lambda.*X(:).*sin(Theta(i,j))*cos(Phi(i,j))+1i.*2*pi/lambda.*Y(:).*sin(Theta(i,j))*sin(Phi(i,j)));
Data0thetaphi(i,j)=T0'*Sig;
[i,j]
end;
end
%%
Power0Norm=sum(sum(abs(Data0thetaphi).^2.*sin(Theta)*dphi*dtheta))
%%
max(10*log10(1/Power0Norm*abs(Data0thetaphi(:)).^2))
max(10*log10(1/PowerNorm*abs(Datathetaphi(:)).^2))
%%
1/Power0Norm*sum(sum(abs(Data0thetaphi).^2.*sin(Theta)*dphi*dtheta))
