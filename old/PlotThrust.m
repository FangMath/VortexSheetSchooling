% (applies only for heaving flapping)
% - plot the Lift (IntegrateP, integral of pressure jump over chord)
% - plot the Thrust (Leading edge suction)

function PlotThrust(str_input,apstr)
close all;
if nargin~=0
    str=str_input;
else
    str='./output/';
end

load( [str,'parameter.mat'],'Para');
load( [str,'data.mat'],'Sys','tem_at');
% apstr=str(end-1);
str=[str,'Forces/'];
if exist(str,'dir')
    rmdir(str,'s');
end
mkdir(str);

n=tem_at;
t=0:Para.dt:(n-1)*Para.dt;
h1=figure(1);
if floor(n*Para.dt)>0
for i=1:floor(n*Para.dt)
    meanThrust(i)=mean(Sys.Thrust((i-1)/Para.dt+1:i/Para.dt));
    meanLift(i)=mean(Sys.Lift((i-1)/Para.dt+1:i/Para.dt));
    meanTorque(i)=mean(Sys.Torque((i-1)/Para.dt+1:i/Para.dt));
end
plot(1:floor(n*Para.dt),-meanThrust,'r.-',[1,floor(n*Para.dt)],Para.g*[1,1],'b-');grid on;
title(['Average thrust for each stroke, and gravity, \theta_{amp}=0.',apstr]);
xlabel('t/period');
legend('Average thrust','minus gravity','Location','best');
save([str,'mean.mat'],'meanThrust','meanLift','meanTorque');
str1=[str,'AThrust_01',apstr];
saveas(h1,str1,'fig');
saveas(h1,str1,'png');
end

%% - plot the Thrust
h2=figure(2);
size(t)
size(Sys.Thrust)
plot(t,-Sys.Thrust,'r-',[t(1),t(end)],Para.g*[1,1],'b-'); grid on;
title(['Thrust and gravity,\theta_{amp}=0.',apstr]);
xlabel('t/period'); 
% axis([0 10 -15 20]);
legend('Thrust','minus gravity','Location','best');
str2=[str,'Thrust_01',apstr];
saveas(h2,str2,'fig');
saveas(h2,str2,'png');

% %% - plot the Lift (force in x direction
% h3=figure(3);
% plot(t,Sys.Lift,'r-'); grid on;
% xlabel('t/period');ylabel('force in x direction');
% xlim([0 10]);
% title(['Force in X direction, \theta_{amp}=0.',apstr]);
% str3=[str,'Lift_01',apstr];
% saveas(h3,str3,'fig');
% saveas(h3,str3,'png');
% 
% 
% % meanLift=mean(Sys.Lift((over-num)/Para.dt+1:over/Para.dt))
% % % meanThrust=mean(Sys.Thrust((over-num)/Para.dt+1:over/Para.dt))
% % 
% %% - plot the Torque 
% h4=figure(4);
% plot(t,Sys.Torque,'b-'); grid on;
% xlabel('t/period');ylabel('torque');
% title(['Torque, \theta_{amp}=0.',apstr]);
% xlim([0 10]);
% str4=[str,'Torque_01',apstr];
% saveas(h4,str4,'fig');
% saveas(h4,str4,'png');
% % % meanTorque=mean(Sys.Torque((over-num)/Para.dt+1:over/Para.dt))
end
