function Gmax(str_input,m,ifSavegammaMax)
close all;
global wing Tpos Tneg a an
if nargin~=0
    str=str_input;
else
    str='./output/';
end

if ifSavegammaMax==1
str=[str,'Ddata/'];
for tk=1:m
    if mod(tk,50)==0
        tk
    end
    load( [str,'T',num2str(tk),'.mat'], 'wing');
    gammaMaxRight(tk)=-wing(1).gammaMax;
    gammaMaxLeft(tk)=-wing(2).gammaMax;
end
save([str_input,'GMAX.mat'],'gammaMaxLeft','gammaMaxRight');
else
load([str_input,'GMAX.mat'],'gammaMaxLeft','gammaMaxRight');
end

h1=figure(1);
m
size(gammaMaxRight)
plot(1:m,gammaMaxRight(1:m),'r-',1:m,gammaMaxLeft(1:m),'b-');
legend('right sheet','left sheet','location','best');
grid on;
h2=figure(2);
CirPerTime_R=gammaMaxRight(2:end)-gammaMaxRight(1:end-1);
CirPerTime_L=gammaMaxLeft(2:end)-gammaMaxLeft(1:end-1);
plot(1:m-1,CirPerTime_R(1:m-1),'r.-',1:m-1,CirPerTime_L(1:m-1),'b--',1:m-1,zeros(m-1,1),'k-');
legend('right sheet','left sheet','location','best');
grid on;
end