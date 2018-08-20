clear;
clc;
Radar.PRF=2300;
Radar.vp=150;
Radar.M=64;
Radar.ChannelNum=4;
Radar=Setup(Radar);

%%
Clutter=Clutter_Gen(Radar);
%%

%%
CovEst=Cov_Est(Radar,Clutter);
InvRu=inv(CovEst);
%%

OSTAPWeights=OSTAP(Radar,InvRu);
%%
TestClutter=Clutter_Gen(Radar);
%%
%
%%Tgt Gen
Radar.vtgt=3;
TestSig=Add_Tgt(Radar,TestClutter);

%%
OSTAPApply(Radar,OSTAPWeights,TestSig,CovEst)

%%
PostDopplerWeights=PostDoppler(Radar,CovEst);
%%
PostDopplerApply(Radar,PostDopplerWeights,TestSig,CovEst)

%%
SumPlot(Radar,TestSig)