% calculate input output power, efficiency, and the averaged values.

function Power(str_input)
close all;
global wing Para
if nargin~=0
    str=str_input
else
    str='./output/';
end

load( [str,'parameter.mat'],'Para');
dt=Para.dt;Nw=Para.Nw; M = Para.M;


%% calculate input power = -(p+ps)*dtheta
str=[str,'Ddata/'];
load( [str,'tem.mat'], 'tem_at');

%% get p and ps, and save, or load from existing
st = tem_at - 1005
ed = tem_at
load( [str_input,'data.mat'],'Wing','Sys');
for ib = 1:Nw
    % calculate integrateP
    p(ib,:) = Wing(ib).integrateP;
    
    % calculate integratePS
    for tk=st:ed%1:ed-st
%         if mod(tk,100) == 0
%             tk
%         end
        load( [str,'T',num2str(tk),'.mat'], 'wing','sys');
        cx = sys.CenterX; cy = sys.CenterY;
        c=cx+cy*1i;
        PART1(ib,tk)=(wing(ib).L+wing(ib).tangential-c)*Wing(ib).integrateP(tk);
        ps(ib,tk) = Wing(ib).torque(tk)-real(conj(wing(ib).tangential)*(PART1(ib,tk)))-Wing(ib).les(tk)*real((wing(ib).zetabody(M+1)-c)*conj(wing(ib).normal));        
    
        dtheta(ib,tk) = wing(ib).dtheta; 
    end
    P(ib,:) = p(ib,st:ed);
    Ps(ib,:) = ps(ib,st:ed);
    Dtheta(ib,:) = dtheta(ib,st:ed);
end
% save([str_input,'power.mat'],'P','Ps','Dtheta','st','ed');

% load([str_input,'power.mat'],'P','Ps','Dtheta','st','ed');
PInput = zeros(1,ed-st+1);
for ib = 1:Nw
    pInput(ib,:) = -Dtheta(ib,:).*(P(ib,:)+Ps(ib,:));
    PInput = PInput+pInput(ib,:);
end
% figure(5)
% plot(PInput);grid on;
save([str_input,'power.mat'],'P','Ps','Dtheta','PInput','st','ed');

% load( [str_input,'data.mat'],'Sys');
Vel=-Sys.dCenterX(st:ed);
Thrust = Sys.Thrust(st:ed);

% figure(1)
% plot(Vel);grid on;
% figure(2)
% plot(Thrust);grid on;
% figure(3)
POutput = -Vel.*Thrust;
% plot(POutput);grid on;
% figure(4)
for k=1:10
    ApOut(k)=mean(POutput((k-1)*100+1:(k-1)*100+100));
    ApIn(k)=mean(PInput((k-1)*100+1:(k-1)*100+100));
end
eta = mean(ApOut./ApIn)
ApOut = mean(ApOut)
ApIn = mean(ApIn)
save([str_input,'power.mat'],'P','Ps','Dtheta','PInput','POutput','ApOut','ApIn','eta','st','ed');
% plot(ApOut);grid on; hold on;
% plot(ApIn);grid on;
% plot(ApOut./ApIn);grid on; hold on;


% figure(1)
% plot(P(1,:)+P(2,:));grid on;
% figure(4)
% plot(pInput(1,:)-pInput(2,:));grid on;
% figure(3)
% plot(Dtheta(1,:)+Dtheta(2,:));grid on;
% figure(2)
% plot(Ps(1,:)+Ps(2,:));grid on;

figure(1)
plot(P(1,:));hold on;plot(P(2,:));hold on;grid on;
% figure(3)
% plot(Dtheta(1,:));hold on;plot(Dtheta(2,:));hold on;grid on;
% figure(2)
% plot(Ps(1,:));hold on;plot(Ps(2,:));hold on;grid on;
end