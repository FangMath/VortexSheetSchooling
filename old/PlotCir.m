function [] = PlotCir(str_input,n)
% - plot the the vortex strength (circulation difference) 
%   shed at trailing edge at time t, vs time (figure 1)
% - plot the total circulation of the free sheet, vs time (figure 2)

close all;

if nargin~=0
    str=str_input;
else
    str='./output/';
end

load( [str,'parameter.mat'],'Para');

str=[str,'Ddata/'];
% load( [str,'tem.mat'], 'tem_at');
load( [str,'T',num2str(n),'.mat'], 'wing','sys');

% load( [str,'data.mat'],'Wing','tem_at');

% n=tem_at;
dt=Para.dt;
t=0:dt:(n-1)*dt;

h1=figure(1);
for ib=1:Para.Nw
gammaf=wing(ib).gammaf;

%% - plot the the vortex strength (circulation difference) shed at trailing edge at time t, vs time (figure 1)

Cirl_vortex=zeros(1,wing(ib).Lenf); %vortex strength shed at trailing edge
%NOTE: Gamma we defined is the minus of true circulation
for i=2:wing(ib).Lenf
    Cirl_vortex(i)=gammaf(i-1)-gammaf(i); 
end


plot(0:1:(wing(ib).Lenf-1),Cirl_vortex,'r.-',t,zeros(n,1),'k--'); grid on;
title(['Circulation of the vortex shed at time t, dt=',num2str(dt)]);
xlabel('number');ylabel('circulation');
str1=[str,'Cirl_vortex'];
% saveas(h1,str1,'fig');
% saveas(h1,str1,'png');
% 
% str3=[str,'Cirl_vortex.mat'];
% save (str3, 'Cirl_vortex');
% hold on;
% %% - plot the total circulation of the free sheet, vs time (figure 2)
% h2=figure(2);
% plot(t,-Wing(ib).GammaMax,'b.-',t,zeros(n,1),'k--'); grid on;
% title(['Total circulation on free vortex sheet, dt=',num2str(dt)]);
% xlabel('t/Period');ylabel('circulation');
% % str2=[str,'Total_cirl'];
% % saveas(h2,str2,'fig');
% % saveas(h2,str2,'png');
% hold on;
end
figure(2)
plot(wing(1).gammaf)
hold off;
end