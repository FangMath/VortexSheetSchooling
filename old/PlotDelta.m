function [] = PlotDelta(str_input)
% - plot the smoothing parameter on free sheet and bounded sheet
% \delta vs s
close all;

if nargin~=0
    str=str_input;
else
    str='./output/';
end

load( [str,'parameter.mat'],'Para');
% 'M','Delta1','Delta0','eta','p','R1','g','dt','s');
load( [str,'data.mat'],'tem_at','Wing');
% 'ZetaBody','ZetaF','GammaF','tem_at','Lenf');

tk=tem_at;
Delta0=Para.Delta0
Delta1=Para.Delta1

dt=Para.dt; eta=Para.eta; p=Para.p; s=Para.s; M=Para.M;

for ib=1:Para.Nw
zetabody=Wing(ib).ZetaBody(tk,:);
zetaf=Wing(ib).ZetaF{tk};
gammaf=Wing(ib).GammaF{tk};

%% calculate the arclength s on free sheet, and \delta(s) on both bounded and free sheet
K=Wing(ib).Lenf(tk);
distance(K)=0;
for k=K-1:-1:1
    space(k)=abs(zetaf(k)-zetaf(k+1));
    distance(k)=abs(zetaf(k)-zetaf(k+1))+distance(k+1);
    delta_f(k)=Delta1+(Delta0-Delta1)*((distance(k)/eta)^p/(1+(distance(k)/eta)^p));
end
delta_f(K)=Delta1; space(K)=0;

for m=1:M+1
    delta_b(m)=Delta1*exp(-(abs(s(m)-1)/Delta1)^2);
end

%% - plot the smoothing parameter on free sheet and bounded sheet
fig1=figure(ib);
plot(distance+1,delta_f,'r.-',s,delta_b,'b.-',distance+1,space,'g.-');grid on;
title(['Smoothing parameters \delta(s), dt=',num2str(dt)]);
% str2=[str,'Delta'];
% saveas(fig1,str2,'fig');
% saveas(fig1,str2,'png');
end
end
