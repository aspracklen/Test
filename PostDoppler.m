function [ DopplerWeights ] = PostDoppler( Radar,CovEst )
%UNTITLED10 Summary of this function goes here
%   Detailed explanation goes here
%% Doppler Filter Bank Creation:
dopplerfilterbank = linspace(-Radar.PRF/2,Radar.PRF/2,Radar.M+1);
omegadopplerbank = dopplerfilterbank/Radar.PRF;


% Doppler Filter Matrix Construction for Adjacent Bin Post-Doppler method
clc
U2 = zeros(Radar.M,Radar.M);
K=3;
P = floor(K/2);
    for m=1:Radar.M
        U2(:,m) =  1/sqrt(Radar.M)*exp(1i*2*pi*omegadopplerbank(m)*(0:Radar.M-1));  % Doppler Filter Steering Vector
    end



td30ab = chebwin(Radar.M,30);                                               % 60-dB Chebyshev Doppler Taper.


% Create Doppler Filter Bank in Fab matrix for Adjacent Bin post Doppler method:
Fab30 = diag(td30ab)*U2;
Fmab30 = zeros(Radar.M,K,Radar.M);
for m=1:Radar.M
        if (m-P>0) && (m+P<=Radar.M)
            Fmab30(:,:,m) = Fab30(:,m-P:m+P);
        elseif (m-P<=0) && (m+P<=Radar.M)
            Fmab30(:,:,m) = [Fab30(:,Radar.M+(m-P):Radar.M) Fab30(:,1:m+P)];
        elseif m+P>Radar.M
            Fmab30(:,:,m) = [Fab30(:,m-P:Radar.M) Fab30(:,1:m+P-Radar.M)];
        end   
end
%
DopplerWeights=zeros(256,Radar.M);
for m=1:Radar.M
    bdfb = exp(1i*2*pi*omegadopplerbank(m)*(0:Radar.M-1)).';
    at=exp(1i*2*pi/Radar.lambda*(Radar.X.*cosd(90)*cosd(Radar.theta0)-Radar.Z.*sind(Radar.theta0)));
    gt = kron(bdfb,Radar.TMat'*at(:));                   % Desired Vector common for both methods.
    %% Adjacent-Bin SINR Computations
    f30abm = Fmab30(:,:,m);
    R30abum = kron(f30abm,eye(Radar.ChannelNum))'*CovEst*kron(f30abm,eye(Radar.ChannelNum));
    gt30abm = kron(f30abm,eye(Radar.ChannelNum))'*gt;
    % Calculate K*N X 1 Adaptive Weight for m-th Doppler Bin.
    w60abm = R30abum\gt30abm;
    wab60 = kron(f30abm,eye(Radar.ChannelNum))*w60abm;  
DopplerWeights(:,m)=wab60;

end

end

