function plotcxwh(str_input)
close all;
clearvars -except str_input
load([str_input,'/data.mat'], 'Wing','tem_at')
load([str_input,'/parameter.mat'], 'Para')
tem_at
length(Wing(2).Thrust)
Wing(2).Thrust(1:tem_at) = Wing(2).Thrust(1:tem_at) - Para.exforce(1:tem_at);

%N = 10; %number of periods to present
N = floor(tem_at/100)-1; %number of periods to present
st = 1
ed = tem_at
idx = [st : ed];
T = idx/100;
idxN = [floor(st/100) : 1 : floor(ed/100)-1] + 1;
TN = idxN - 0.5;

d = (Wing(2).Cx-Wing(1).Cx);
g = (Wing(2).Cx-Wing(1).Cx-2);
school = -g./Wing(1).dotCx;

for k = 1:floor(length(Wing(1).dotCx)/100)
Ap1(k) = mean(Wing(1).Cx(100*(k-1)+1:100*k));
Ap2(k) = mean(Wing(2).Cx(100*(k-1)+1:100*k));

Av1(k) = mean(Wing(1).dotCx(100*(k-1)+1:100*k));
Av2(k) = mean(Wing(2).dotCx(100*(k-1)+1:100*k));
end

for k = 1:floor(tem_at/100)
AFor1(k) = mean(Wing(1).Thrust(100*(k-1)+1:100*k));
AFor2(k) = mean(Wing(2).Thrust(100*(k-1)+1:100*k));

ASki1(k) = mean(Wing(1).Skin(100*(k-1)+1:100*k));
ASki2(k) = mean(Wing(2).Skin(100*(k-1)+1:100*k));

Ales1(k) = mean(Wing(1).Les(100*(k-1)+1:100*k));
Ales2(k) = mean(Wing(2).Les(100*(k-1)+1:100*k));
ASch(k) = mean(school(100*(k-1)+1:100*k));
end

AThr1 = AFor1 - ASki1;
AThr2 = AFor2 - ASki2;


Pin1 = -Para.HAmp(1)*2*pi*sin(2*pi*T).*Wing(1).IntegrateP(idx);
if length(Para.HAmp) == 1
Para.HAmp(2) = Para.HAmp(1); %% HACK
end
Pin2 = -Para.HAmp(2)*2*pi*sin(2*pi*T).*Wing(2).IntegrateP(idx);
for k = idxN
APin1(1,k-idxN(1)+1) = mean(Pin1(100*(k-idxN(1))+1:100*(k-idxN(1)+1)));
APin2(1,k-idxN(1)+1) = mean(Pin2(100*(k-idxN(1))+1:100*(k-idxN(1)+1)));


end

ASchV = [0,diff(ASch)];
ASchA = [0,diff(ASchV)];

sch1idx = [25:50];
sch1idx = [floor(ed/100)-1:floor(ed/100)];
A = [ASchA(sch1idx)',ASchV(sch1idx)', ASch(sch1idx)'];
klb = A\ones(length(sch1idx),1)
school1 = 1/klb(3)
%
%sch2idx = [100:130];
%A = [ASchA(sch2idx)',ASchV(sch2idx)', ASch(sch2idx)'];
%klb = A\ones(length(sch2idx),1)
%school2 = 1/klb(3)
%
%sch3idx = [170:200];
%A = [ASchA(sch3idx)',ASchV(sch3idx)', ASch(sch3idx)'];
%klb = A\ones(length(sch3idx),1)
%school3 = 1/klb(3)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ifplot = 1;
if ifplot
%%%%% wing velocity %%%%%%
f1 = figure(1);
plotset(3.5,3);
plot(T, -Wing(1).dotCx(idx),'b.-');hold on
plot(T, -Wing(2).dotCx(idx),'r.-');hold on
%plot([T(1),T(end)], 13.24*[1,1], 'k--');
plot(TN, -Av1(idxN), 'b*--');hold on; 
plot(TN, -Av2(idxN), 'r*--');hold on; 
legend({'W1','W2','single'}, 'Location', 'best', 'FontSize', 8);
ylabel('wing velocity');
xlabel('t');
title(['Wing Speed, \beta_0=',num2str(2*Para.HAmp(1)),',C_f=',num2str(Para.skin)]);
xlim([T(1),T(end)]);
grid on;


%%%%% wing Force %%%%%%
f2 = figure(2);
plotset(3.5,3);
plot(T, -Wing(1).Thrust(idx), 'b.-');hold on;
plot(T, -Wing(2).Thrust(idx), 'r.-');hold on;
plot(TN, -AFor1(idxN), 'b*--');hold on;
plot(TN, -AFor2(idxN), 'r*--');hold on;
legend({'W1','W2'}, 'Location', 'best', 'FontSize', 8);
grid on;
title(['Force, \beta_0=',num2str(2*Para.HAmp(1)),',C_f=',num2str(Para.skin), ',schl=',num2str(Para.iniSch)]);
ylabel('Force');
xlabel('t');

%%%%% skin friction %%%%%%
f3 = figure(3);
plotset(3.5,3);
plot(T, Wing(1).Skin(idx), 'b'); hold on;
plot(T, Wing(2).Skin(idx), 'r'); hold on;
plot(TN, ASki1(idxN), 'b.-'); hold on;
plot(TN, ASki2(idxN), 'r.-'); hold on;
title('Skin friction drag');
ylabel('Drag');
xlabel('t');
xlim([T(1),T(end)]);
grid on;

%%%%% thrust %%%%%%
f4 = figure(4);
plotset(3.5,3);
plot(T, -Wing(1).Thrust(idx)+Wing(1).Skin(idx), 'b'); hold on;
plot(T, -Wing(2).Thrust(idx)+Wing(2).Skin(idx), 'r'); hold on;
plot(TN, -AFor1(idxN)+ASki1(idxN), 'b.-'); hold on;
plot(TN, -AFor2(idxN)+ASki2(idxN), 'r.-'); hold on;
title('thrust');
ylabel('thrust');
xlabel('t');
xlim([T(1),T(end)]);
grid on;

%%%%% Power in %%%%%%
%f4 = figure(4);
%plotset(3.5,3);
%plot(T, -Pin1, 'b'); hold on;
%plot(T, -Pin2, 'r'); hold on;
%title('Power input');
%plot(TN, -APin1, 'b*--');hold on; 
%plot(TN, -APin2, 'r*--');hold on; 
%plot(TN, -(-AT1(idxN)+ASki1).*Av1(idxN), 'b.-'); hold on;
%plot(TN, -(-AT2(idxN)+ASki2).*Av2(idxN), 'r.-'); hold on;
%ylabel('P_{in}');
%xlabel('t');
%%xlim([TN(1),TN(end)]);

%%%%% Power out %%%%%%
%f5 = figure(5);
%plotset(3.5,3);
%plot(TN, (-AT1(idxN)+ASki1).*Av1(idxN)./APin1, 'b*--');hold on; 
%plot(TN, (-AT2(idxN)+ASki2).*Av2(idxN)./APin2, 'r*--');hold on; 
%title('efficiency');
%xlabel('t');

%%%%%%% fit Thrust vs S %%%%%%
%f7 = figure(7);
%plotset(3.5,3);
%Nst = 3; % start fits after Nst periods
%yfit = AThr2;
%fx = ASch(Nst:end);
%fy = yfit(Nst:end);
%grid on;
%
%%modelfun = @(b,x)(b(1)+b(2)*sin(b(3)*x+b(4)));
%%beta0 = [-2.5;4;2*pi;pi/2];
%%opts = statset('nlinfit');
%%beta = nlinfit(fx,fy,modelfun,beta0,opts)
%
%modelfun = @(b,x)(b(1)+b(2)*cos(2*pi*(x+b(3))));
%beta0 = [-.5;.4;-0.1];
%opts = statset('nlinfit');
%beta = nlinfit(fx,fy,modelfun,beta0,opts)
%ratio = -beta(1)/beta(2); rtio = ratio-sqrt(ratio^2-1)
%predictD = (-beta(1)+beta(2))/(1+rtio)^2
%predictSch = acos(rtio/2)/(2*pi)-beta(3)
%lambda = Av1(end)
%
%%Drg1 = AThr1(min(55,floor(ed/100)))
%%fstSch = acos((Drg1-beta(1))/beta(2))/(2*pi) - beta(3)+1
%%Drg2 = AThr1(end)
%%sndSch = acos((Drg2-beta(1))/beta(2))/(2*pi) - beta(3)+2
%Tht1 = beta(1) + beta(2)*cos(2*pi*(school1+beta(3)));
%Tht2 = beta(1) + beta(2)*cos(2*pi*(school2+beta(3)));
%Tht3 = beta(1) + beta(2)*cos(2*pi*(school3+beta(3)));
%hold on;
%fitx = 0:.01:4;
%plot(fitx,-(modelfun(beta,fitx)),'k-'); hold on;
%plot([fitx(1),fitx(end)],-Tht1*ones(1,2),'r--'); hold on;
%plot([fitx(1),fitx(end)],-Tht2*ones(1,2),'b--'); hold on;
%plot([fitx(1),fitx(end)],-Tht3*ones(1,2),'g--'); hold on;
%plot(ASch(Nst:end), -yfit(Nst:end), 'b.-'); hold on;
%legend({'sin fit','drag'},'Location','best','FontSize',9);
%xlabel('S'); ylabel('Thrust');
%title(['Sine fit of T(S), C_f=',num2str(Para.skin), ', \beta_0=',num2str(2*Para.HAmp)]);

%%%%% plot separation %%%%%%
f8=figure(8);
plotset(3.5,3);
plot((1:length(g))/100, g); hold on;
grid on;
xlabel('t');
ylabel('g');
xlim([T(1),T(end)]);

%%%%% plot schooling number %%%%%%
f9=figure(9);
plotset(3.5,3);
plot((1:length(school))/100, school); hold on;
%plot([0,length(school)]/100, finalschool*[1,1], 'k--');
plot([0,length(school)]/100, school1*[1,1], 'k--');
%plot([0,length(school)]/100, school2*[1,1], 'k--');
%plot([0,length(school)]/100, school3*[1,1], 'k--');
grid on;
xlabel('t');
ylabel('S');
xlim([T(1),T(end)]);


%%%%%%% schooling number phase diagram %%%%%%
%f5 = figure(5);
%plotset(3.5,3);
%plot(ASchA, ASch, 'b.-'); hold on;
%plot(ASch(2), ASchV(2), 'ko'); hold on;
%%plot(finalschool, 0, 'm*'); hold on; 
%plot(school1, 0, 'm*'); hold on; 
%plot(school2, 0, 'm*'); hold on; 
%plot(school3, 0, 'm*'); hold on; 
%title('Sching number phase diagram');
%xlabel('S'); ylabel('dotS'); grid on;
%
%f11 = figure(11);
%plotset(3.5,3);
%plot(T,Para.exforce(idx));
%grid on;
%title('External force on follower');
%xlabel('T');
%ylabel('External force');
%xlim([T(1),T(end)]);
end


ifsave = 0;
if ifsave
% change save path
str_input = ['../result_d_newdrag/vary_force/h',num2str(Para.HAmp(1)),',skin',num2str(Para.skin)]
saveas(f1, [str_input, 'school.fig']) 
saveas(f1, [str_input, 'school.png']) 

saveas(f2, [str_input, 'vel.fig']) 
saveas(f2, [str_input, 'vel.png']) 
%saveas(f3, [str_input, 'thrust.fig']) 
%saveas(f3, [str_input, 'thrust.png']) 

saveas(f4, [str_input, 'S_dotS.fig']) 
saveas(f4, [str_input, 'S_dotS.png']) 

saveas(f7, [str_input, 'Thrust_S.fig']) 
saveas(f7, [str_input, 'Thrust_S.png']) 

saveas(f9, [str_input, 'drag.fig']) 
saveas(f9, [str_input, 'drag.png']) 

saveas(f10, [str_input, 'thrust.fig']) 
saveas(f10, [str_input, 'thrust.png']) 

saveas(f11, [str_input, 'externalF.fig']) 
saveas(f11, [str_input, 'externalF.png']) 
end
