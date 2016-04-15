function MOVIEsnapshotV(str_input,apstr)
close all;
global Cirl_vortex Wing Sys tem_at Nw
if nargin~=0
    str=str_input;
else
    str='./output/';
end

load( [str,'parameter.mat'],'Para');
dt=Para.dt;M=Para.M;Nw=Para.Nw;
% load( [str,'data.mat'], 'tem_at','Wing','Sys');


% Initialize data at i=tem_at(1):
start=1;

% set axis frame
left=-10;right=10; down=-9; up=2;

%set x, y for the central line
x=[0,0]; y=[down,up];

%  plot body and free sheet at each time step
str=[str,'Ddata/'];
% if exist(str,'dir')
%     rmdir(str,'s');
% end
% mkdir(str);

%calculate circulation for each vortex
load( [str,'tem.mat'], 'tem_at');
tem_at

% vidObj = VideoWriter([str,'movie_01',apstr,'.avi']);
% vidObj.FrameRate=50;
% open(vidObj);

h=figure(1);

% tem_at=tem_at-2;
for i=1:tem_at-1
    i
    load( [str,'T',num2str(i+1),'.mat'], 'wing','sys');
    for ib=1:Nw
        Wing(ib).Nu(i,:)=wing(ib).nu;
        Wing(ib).ZetaBody(i,:)=wing(ib).zetabody;
        %         for j=1:wing(ib).K
        Wing(ib).Lenf(i)=length(wing(ib).zetaf);
        a=Wing(ib).Lenf(i)
            Wing(ib).ZetaF{i}=wing(ib).zetaf;
            Wing(ib).GammaF{i}=wing(ib).gammaf;
%         end
        Wing(ib).GammaMax(i)=wing(ib).gammaMax;
        
    end
    
    Sys.CenterX(i)=sys.CenterX;
    Sys.CenterY(i)=sys.CenterY;
    Sys.Omega(i)=sys.Omega;

    Sys.dCenterX(i)=sys.dCenterX;
    Sys.dCenterY(i)=sys.dCenterY;
    Sys.dOmega(i)=sys.dOmega;
    
    n=[];gammaf=[];Cirl_vortex=[];
    for ib=1:Nw
        n(ib)=Wing(ib).Lenf(i);
        aa=Wing(ib).GammaF{i};
        size(aa)
        n(ib)
        gammaf(ib,1:n(ib))=Wing(ib).GammaF{i};
        Cirl_vortex(ib,1)=0;
        %note that Gamma we defined is the minus of true circulation
        Cirl_vortex(ib,2:n(ib))=gammaf(ib,1:n(ib)-1)-gammaf(ib,2:n(ib));
        
    end
    plot(Sys.CenterY(i),-Sys.CenterX(i),'ko','MarkerFaceColor','k',...
        'MarkerSize',6);hold on;
    plot([imag(Wing(1).ZetaBody(i,end)),imag(Wing(2).ZetaBody(i,end))],[-real(Wing(1).ZetaBody(i,end)),-real(Wing(2).ZetaBody(i,end))],'-m');
    x1=imag(Wing(1).ZetaBody(i,end)+Wing(2).ZetaBody(i,end))/2;     y1=-real(Wing(1).ZetaBody(i,end)+Wing(2).ZetaBody(i,end))/2;
    x2=Sys.CenterY(i);  y2=-Sys.CenterX(i);
    slope=(y2-y1)/(x2-x1);
    %     x3=x1+2*sign(x2-x1)/(sqrt(slope^2+1));     y3=y1+slope*(x3-x1);
    plot([x1,x2],[y1,y2],'m-');

    plot(Sys.CenterY(1),-Sys.CenterX(1),'co','MarkerFaceColor','c',...
        'MarkerSize',6);hold on;
    plot([imag(Wing(1).ZetaBody(1,end)),imag(Wing(2).ZetaBody(1,end))],[-real(Wing(1).ZetaBody(1,end)),-real(Wing(2).ZetaBody(1,end))],'-c');
    plot([Sys.CenterY(1),imag(Wing(1).ZetaBody(1,end)+Wing(2).ZetaBody(1,end))/2],[-Sys.CenterX(1),-real(Wing(1).ZetaBody(1,end)+Wing(2).ZetaBody(1,end))/2],'-c');
    plot(imag(Wing(1).ZetaBody(1,:)),-real(Wing(1).ZetaBody(1,:)),'c-','LineWidth',3);hold on;
    plot(imag(Wing(2).ZetaBody(1,:)),-real(Wing(2).ZetaBody(1,:)),'c-','LineWidth',3);hold on;
    plot(x,y,'g--');hold on;
    
    %     figure(i)
    PreMOVIE(.2,i);
    for ib=1:Nw
%         index_piece(i,Wing(ib).Lenf(i),ib); %%%%%%%%%% plot wakes
        

        plot(imag(Wing(ib).ZetaBody(i,:)),-real(Wing(ib).ZetaBody(i,:)),'k-','LineWidth',3);hold on;
%          plot(real(Wing(ib).ZetaBody(i,:)),imag(Wing(ib).ZetaBody(i,:)),'m.-');hold on;
   end
    [pot,target]=VField(i);
    quiver(target(2,:),-target(1,:),imag(pot)/6,-real(pot)/6,0);hold on;
    hold off;
    axis([left right down up]);
    set(gcf,'double','on');
    axis equal; axis manual; grid on;
    title(['Symmetric flappying, \theta_{amp}=0.1',apstr,', t=',num2str((i-1)*dt)]);
%     title(['Symmetric flappying, \theta_{amp}=0.150, t=',num2str((i-1)*dt)]);
    xlabel('x'); ylabel('y');
    whitebg([1 1 1]);
        pause(.1)
    strs=[str,'T',num2str(i)]
%         strs=[str,'S_0150T',num2str(i)]
%     saveas(h,strs,'fig');
%     saveas(h,strs,'png');
%     currFrame = getframe(h);
%     writeVideo(vidObj,currFrame);
end
% close(vidObj);
end

% function index_piece(i,n,ib)
% global Cirl_vortex Wing
% pos=0;neg=0; %index of pos or neg piece
% piece_p=[];piece_n=[];
% for k=1:n
%     if k==1
%         if Cirl_vortex(ib,k)>=0 %determine the first piece
%             pos=1; piece_p(pos,1)=k; %start of pos
%         else
%             neg=1; piece_n(neg,1)=k; %start of neg
%         end
%     end
%     if k<n&&Cirl_vortex(ib,k)>=0&&Cirl_vortex(ib,k+1)<0
%         neg=neg+1; %pos piece ends, neg piece begins
%         piece_p(pos,2)=k; piece_n(neg,1)=k;  %end of pos, start of neg
%     else if k<n&&Cirl_vortex(ib,k)<0&&Cirl_vortex(ib,k+1)>=0
%             pos=pos+1; %neg piece ends, pos piece begins
%             piece_n(neg,2)=k; piece_p(pos,1)=k;  %end of pos, start of neg
%         else if k==n&&Cirl_vortex(ib,k)>=0 %final end is pos
%                 piece_p(pos,2)=k;
%             else if k==n&&Cirl_vortex(ib,k)<0 %final end is neg
%                     piece_n(neg,2)=k;
%                 end
%             end
%         end
%     end
% end
% msize=.01;
% for k=1:pos
%     p1=plot(imag(Wing(ib).ZetaF{i}(piece_p(k,1):piece_p(k,2))),-real(Wing(ib).ZetaF{i}(piece_p(k,1):piece_p(k,2))),...
%         'r.','MarkerSize',msize);hold on;
% %     set(p1,'color',[.5 .5 .5]);
% end
% for k=1:neg
%     p2=plot(imag(Wing(ib).ZetaF{i}(piece_n(k,1):piece_n(k,2))),-real(Wing(ib).ZetaF{i}(piece_n(k,1):piece_n(k,2))),...
%         'b.','MarkerSize',msize);hold on;
% %     set(p2,'color',[.5 .5 .5]);
% end
% end
