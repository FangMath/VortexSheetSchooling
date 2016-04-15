function SheetCir(str_input,m)
% close all;
global wing Tpos Tneg piece_n piece_p Cirl_vortex cirFit cirDiff
if nargin~=0
    str=str_input;
else
    str='./output/';
end
load( [str,'parameter.mat'],'Para');
Nw=Para.Nw;

str=[str,'Ddata/'];
load( [str,'T',num2str(m),'.mat'], 'wing','sys');

ploton=0;
% figure(m)
% for ib=1:2
%     %     wing(ib).PointVX=[];
%     %     wing(ib).PointVCirl=[];
%     index_piece(wing,ib,ploton);
%     length(piece_p{ib})
%     for k=1:length(piece_p{ib})
%         Tpos{ib}(k)=sum(Cirl_vortex{ib}(piece_p{ib}(k,1):piece_p{ib}(k,2)));
%     end
%     
%     for k=1:length(piece_n{ib})
%         Tneg{ib}(k)=sum(Cirl_vortex{ib}(piece_n{ib}(k,1):piece_n{ib}(k,2)));
%     end
% end



length(wing(1).PointVCirl)
    if wing(1).PointVCirl(1)>0
%     aPos1 = [wing(1).PointVCirl(1:2:end),Tpos{1}];
%     aNeg1 = [wing(1).PointVCirl(2:2:end),Tneg{1}];
    aPos1 = [wing(1).PointVCirl(1:2:end)];
    aNeg1 = [wing(1).PointVCirl(2:2:end)];
    else
    aPos1 = [wing(1).PointVCirl(2:2:end)];
    aNeg1 = [wing(1).PointVCirl(1:2:end)];
%     aPos1 = [wing(1).PointVCirl(2:2:end),Tpos{1}];
%     aNeg1 = [wing(1).PointVCirl(1:2:end),Tneg{1}];
    end
%     length(aPos1)
%     length(aNeg1)
    figure(5);
%     plot(1:length(aPos1),aPos1(1:end),'.:',1:length(aNeg1),-aNeg1(1:end),'o-');
    plot(1:length(aPos1),aPos1(1:end)+aNeg1(1:end),'o-');
    legend('right positive sheet','right negative sheet','location','best');
    figure(2);
    head = 6;
    tail = 5;
    cirDiff = -aNeg1(head+1:end)-aPos1(head+1:end);
    tail = find(cirDiff<0,1)-10
if ~isempty(tail)
    cirDiff = cirDiff(1:tail);    
end
x = 1:length(cirDiff);
    cirFit = polyfit(x,log(abs(cirDiff)),1)
    plot(x,log(abs(cirDiff)),'o-',x,polyval(cirFit,x),'.-');grid on;
    
    figure(3)
    plot(x,cirDiff,'.-');grid on;

%     plot(x,cirDiff,'.-'); grid on;
    
    % TurnToPointVortices(sys.CenterX,5);
    
end