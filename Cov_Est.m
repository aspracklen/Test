function [ CovEst ] = Cov_Est( Radar,Signal )
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here
CovEst=zeros(Radar.M*Radar.ChannelNum,Radar.M*Radar.ChannelNum);
for i=1:length(Radar.range);
CovEst=CovEst+Signal(:,i)*Signal(:,i)';
end;
CovEst=(1/length(Radar.range).*CovEst)+Radar.Pn*eye(Radar.M*Radar.ChannelNum);

end

