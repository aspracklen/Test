clear;
clc;
N=32^2;
M=16;
lambda=0.03;
x=(-16:15)*lambda/2;
z=(-16:15)*lambda/2;
[X,Z]=meshgrid(x,z);
XF=X(:);
ZF=Z(:);


SAPos=zeros(4,4,N);
for i=1:4;
for  j=1:4;
SACorn=[x(1)-0.03/4,z(1)-0.03/4;
x(1)-0.03/4,z(8)+0.03/4;
x(8)+0.03/4,z(8)+0.03/4;
x(8)+0.03/4,z(1)-0.03/4]+repmat([(i-1)*8*0.03/2,(j-1)*8*0.03/2],4,1);
% pos=find(inpolygon(X(:),Z(:),SACorn(:,1),SACorn(:,2)));
SAPos(i,j,:)=inpolygon(X(:),Z(:),SACorn(:,1),SACorn(:,2));
end;
end;
SAPos=reshape(SAPos,16,N);

% scatter(X(:),Z(:))
% for i=1:16
% pos=find(SAPos(i,:));
% hold on
% scatter(XF(pos),ZF(pos))
% end;
% hold off




Win=taylorwin(32)./max(taylorwin(3));
[WinX,WinY]=meshgrid(Win,Win);
Win=WinX.*WinY;
Win=Win(:);

T=zeros(N,M);
for i=1:M;
pos=find(SAPos(i,:));
T(pos,i)=Win(pos);
end;
Norm=sqrt(diag(T'*T));
T=T./Norm.';
T1=T.*Norm.';

d0=(XF./Win).';
d0=([d0(1:512),-fliplr(d0(1:512))].');
d=d0;

w=ones(N,1);
wdiff=inv(T'*T)*T'*(eye(N)-1/(w'*T*inv(T'*T)*T'*w)*(w*w'*T*inv(T'*T)*T'))*d;
wdiff=reshape(wdiff,4,4).';
wdiff=wdiff(:);
%%
% wdiff=[Norm(1:8);-fliplr(Norm(1:8))];


Noise=1/sqrt(2)*(randn(N,1000)+1i.*randn(N,1000));
SA=T1'*Noise;

wSA=ones(M,1);
Sum=wSA'*SA./sqrt(sum(Norm.^2));
Diff=wdiff'*SA./sqrt(sum((wdiff.*Norm).^2));
% mesh(reshape(d,32,32))
% xlabel('x')
% %
% mesh(reshape(wdiff,4,4))
% xlabel('x')
% hold on
% mesh(reshape(wdiff0,4,4)./10)
% hold off
%%
clearvars SumData DiffData SAData U V DiffConvData
Ang=linspace(-90,90,361);
diffconv=[ones(8,1);-ones(8,1)];


for i=1:length(Ang);
Signal=exp(1i.*2*pi/lambda*XF.*sind(Ang).*cosd(Ang(i))-1i.*2*pi/lambda*ZF(:)*sind(Ang(i)));
SA=T1'*Signal;
Sum=wSA'*SA./sqrt(sum(Norm.^2));
Diff=wdiff'*SA./sqrt(sum((wdiff.*Norm).^2));
DiffConv=diffconv'*SA./sqrt(sum(abs(diffconv.*Norm).^2));
SumData(i,:)=abs(Sum).^2;
DiffData(i,:)=abs(Diff).^2;
DiffConvData(i,:)=abs(DiffConv).^2;

SAData(i,:)=abs(SA(1,:)).^2;
U(i,:)=sind(Ang).*cosd(Ang(i));
V(i,:)=repmat(sind(Ang(i)),1,361);
i
end;

%%
mesh(U,V,(SumData-DiffData)>0)
%%
plot(V(:,182),10*log10(SumData(:,182)),V(:,182),10*log10(DiffConvData(:,182)),V(:,182),10*log10(SAData(:,182)))
ylim([-50,20])
%%
[Row,Col]=find(abs(U-V)<0.005);
for i=1:length(Row);
SumDiag(i)=SumData(Row(i),Col(i));
DiffConvDiag(i)=DiffConvData(Row(i),Col(i));
DiffDiag(i)=DiffData(Row(i),Col(i));

SADiag(i)=SAData(Row(i),Col(i));
end
plot(10*log10(SumDiag))
hold on
plot(10*log10(1/2*(DiffConvDiag+DiffDiag)))
% plot(DiffDiag)
hold off
ylim([-100,30])
