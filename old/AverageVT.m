%calculate the average velocity, and converging rate
%calculate the average thrust, and converging rate
function AverageVT(str,pon,dT)
global wing cirFit cirDiff
load([str,'data.mat'],'Sys','tem_at');
% PlotThrust(str,'un');
load([str,'Forces/mean.mat'],'meanThrust');
a=5;
meanThrust=meanThrust(1+a:end);
inx=find(-meanThrust-2.5<0,1)
if ~isempty(inx)
meanThrust=meanThrust(1:inx-2);
end

L=length(meanThrust)
Trate=log(-meanThrust-2.5);
v=-Sys.dCenterX;

AVelocity=zeros(1,L);
for k=1+a:L+a
AVelocity(k-a)=mean(v((k-1)*100+1:(k-1)*100+100));
end

diffAVelocity=AVelocity(2:end)-AVelocity(1:end-1);
Vrate=log(abs(diffAVelocity));
VTerminal=AVelocity(end)
size(Vrate)
RVline=polyfit(1:L-1,Vrate,1)
RTline=polyfit(1:L,Trate,1)

% sheetCir
% m = length(v)
% SheetCir(str,m);
% 
% 
% if pon==1
% figure(1)
%     plot(1:L-1,Vrate,'r',1:L-1,polyval(RVline,1:L-1),'r',1:L,Trate,'b',1:L,polyval(RTline,1:L),'b',...
%         1:length(cirDiff),log(cirDiff),'.-',1:length(cirDiff),polyval(cirFit,1:length(cirDiff)));
% end

dH = 5;
% dT = 20;
RVline=polyfit(dH+1:L-1-dT,Vrate(dH+1:L-1-dT),1)
RTline=polyfit(dH+1:L-dT,Trate(dH+1:L-dT),1)

% sheetCir
% m = length(v)
% SheetCir(str,m);


if pon==1
figure(1)
    plot(dH+1:L-1-dT,Vrate(dH+1:L-1-dT),'r',dH+1:L-1-dT,polyval(RVline,dH+1:L-1-dT),'r.-',dH+1:L-dT,Trate(dH+1:L-dT),'b',dH+1:L-dT,polyval(RTline,dH+1:L-dT),'b',...
        1:length(cirDiff),log(cirDiff),'.-',1:length(cirDiff),polyval(cirFit,1:length(cirDiff)));
end

save([str,'Adata.mat'],'AVelocity','Vrate','Trate','L','VTerminal','RVline','RTline','cirFit');

end