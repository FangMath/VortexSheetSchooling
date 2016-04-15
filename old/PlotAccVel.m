function [] = PlotAccVel(str_input,apstr)
close all;

if nargin~=0
    str=str_input
else
    str='./output/';
end

load( [str,'parameter.mat'],'Para')
load( [str,'data.mat'],'Sys','tem_at');
% apstr=str(end-1);
str=[str,'AccVel/'];
if exist(str,'dir')
    rmdir(str,'s');
end
mkdir(str);

n=tem_at;
dt=Para.dt;
t=0:dt:(n-1)*dt;

AccX=(Sys.dCenterX(2:n)-Sys.dCenterX(1:n-1))/dt;
AccY=(Sys.dCenterY(2:n)-Sys.dCenterY(1:n-1))/dt;
AccO=(Sys.dOmega(2:n)-Sys.dOmega(1:n-1))/dt;


%% plot vertical position
h1=figure(1);
% plot(t,Sys.CenterX-Sys.CenterX(1),'r.-',t,Sys.CenterY,'b.-',t,Sys.Omega,'g.-'); grid on;
plot(t,-(Sys.CenterX-Sys.CenterX(1)),'r.-'); grid on;
% plot(t,-Sys.CenterY,'r.-',t,-(Sys.CenterX-Sys.CenterX(1)),'b.-',t,Sys.Omega,'g.-'); grid on;
title(['c.m. in y direction, \theta_{amp}=0.',apstr]);
% legend('X','Y','Omega');
xlabel('t/period');ylabel('y');
str1=[str,'PositionY_0',apstr];
saveas(h1,str1,'fig');
saveas(h1,str1,'png');

%% plot horizontal position and tilt angle
h2=figure(2);
plot(t,Sys.CenterY,'b.-',t,Sys.Omega,'g.-'); grid on;
title(['c.m. in X direction and the tilt angle, \theta_{amp}=0.',apstr]);
legend('X','Angle','Location','best');
xlabel('t/period');ylabel('X, Angle');
str2=[str,'PositionXO_0',apstr];
saveas(h2,str2,'fig');
saveas(h2,str2,'png');

%% plot vertical velocity
h3=figure(3);
plot(t,-Sys.dCenterX,'r.-',[t(1),t(end)],[0,0],'b-'); grid on;
title(['c.m. velocity in Y direction, \theta_{amp}=0.',apstr]);
xlabel('t/period');ylabel('Y velocity');
str3=[str,'VelocityY_0',apstr];
saveas(h3,str3,'fig');
saveas(h3,str3,'png');

%% plot horizontal velocity and tilt angular velocity
h4=figure(4);
plot(t,Sys.dCenterY,'b.-',t,Sys.dOmega,'g.-'); grid on;
title(['c.m. velocity in X direction and tilt angular velocity, \theta_{amp}=0.',apstr]);
legend('X velocity','Angular velocity','Location','best');
xlabel('t/period');
str4=[str,'VelocityXO_0',apstr];
saveas(h4,str4,'fig');
saveas(h4,str4,'png');

%% plot vertical acceleration
h5=figure(5);
plot(t(1:end-1),-AccX,'r.-',[t(1),t(end)],[0,0],'b-'); grid on;
title(['c.m. acceleration in Y direction, \theta_{amp}=0.',apstr]);
xlabel('t/period');ylabel('Y acceleration');
% axis([0 10 -4 4]);
str5=[str,'AccY_0',apstr];
saveas(h5,str5,'fig');
saveas(h5,str5,'png');

%% plot horizontal acceleration and tilt angular acceleration
h6=figure(6);
plot(t(1:end-1),AccY,'b.-',t(1:end-1),AccO,'g.-'); grid on;
title(['c.m. acceleration in X direction and tilt angular acceleration, \theta_{amp}=0.',apstr]);
legend('X acceleration','Angular acceleration','Location','best');
xlabel('t/period');
str6=[str,'AccXO_0',apstr];
saveas(h6,str6,'fig');
saveas(h6,str6,'png');


