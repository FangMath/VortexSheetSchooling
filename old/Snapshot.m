function Snapshot(str_input,apstr)
close all;
global Cirl_vortex Wing Sys
if nargin~=0
    str=str_input;
else
    str='./output/';
end

load( [str,'parameter.mat'],'Para');
dt=Para.dt;M=Para.M;Nw=Para.Nw;
load( [str,'data.mat'], 'tem_at','Wing','Sys');
% apstr=str(end-1);
str=[str,'Snapshot/'];
if exist(str,'dir')
    rmdir(str,'s');
end
mkdir(str);

% Initialize data at tk=tem_at(1):
start=1;

% set axis frame
left=-10;right=10; down=-10; up=5;

%set x, y for the central line
x=[0,0]; y=[down,up];

%calculate circulation for each vortex

tem_at
% for ib=1:Nw
%     n(ib)=Wing(ib).Lenf(tem_at);
%
%     gammaf(ib,:)=Wing(ib).GammaF{tem_at};
%     Cirl_vortex(ib,1)=0;
%     %note that Gamma we defined is the minus of true circulation
%     Cirl_vortex(ib,2:n)=gammaf(ib,1:n-1)-gammaf(ib,2:n);
%
% end

% vidObj = VideoWriter([str,'movie_cir.avi']);
% vidObj.FrameRate=10;
% open(vidObj);

h=figure(1);
% tem_at=tem_at-2;
for i=1:tem_at-1
    i
    n=[];gammaf=[];Cirl_vortex=[];
    for ib=1:Nw
        n(ib)=Wing(ib).Lenf(i);
        aa=Wing(ib).GammaF{i};
        size(aa)
        gammaf(ib,1:n(ib))=Wing(ib).GammaF{i};
        Cirl_vortex(ib,1)=0;
        %note that Gamma we defined is the minus of true circulation
        Cirl_vortex(ib,2:n(ib))=gammaf(ib,1:n(ib)-1)-gammaf(ib,2:n(ib));
        
    end
    %     figure(i)
    plot(-Sys.CenterY(i),-Sys.CenterX(i),'m^','MarkerFaceColor','m',...
        'MarkerSize',5);hold on;
    plot([-imag(Wing(1).ZetaBody(i,end)),-imag(Wing(2).ZetaBody(i,end))],[-real(Wing(1).ZetaBody(i,end)),-real(Wing(2).ZetaBody(i,end))],'-m');
    plot([-Sys.CenterY(i),-imag(Wing(1).ZetaBody(i,1)+Wing(2).ZetaBody(i,1))/2],[-Sys.CenterX(i),-real(Wing(1).ZetaBody(i,1)+Wing(2).ZetaBody(i,1))/2],'-m');
    
    plot(-Sys.CenterY(1),-Sys.CenterX(1),'k^','MarkerFaceColor','k',...
        'MarkerSize',5);hold on;
    plot([-imag(Wing(1).ZetaBody(1,end)),-imag(Wing(2).ZetaBody(1,end))],[-real(Wing(1).ZetaBody(1,end)),-real(Wing(2).ZetaBody(1,end))],'-k');
    plot([-Sys.CenterY(1),-imag(Wing(1).ZetaBody(1,1)+Wing(2).ZetaBody(1,1))/2],[-Sys.CenterX(1),-real(Wing(1).ZetaBody(1,1)+Wing(2).ZetaBody(1,1))/2],'-k');
    plot(-imag(Wing(1).ZetaBody(1,:)),-real(Wing(1).ZetaBody(1,:)),'k-');hold on;
    plot(-imag(Wing(2).ZetaBody(1,:)),-real(Wing(2).ZetaBody(1,:)),'k-');hold on;
    plot(x,y,'g-');hold on;


    for ib=1:Nw
        index_piece(i,Wing(ib).Lenf(i),ib);
        plot(-imag(Wing(ib).ZetaBody(i,:)),-real(Wing(ib).ZetaBody(i,:)),'m.-');hold on;
    end
    hold off;
    axis([left right down up]);
    set(gcf,'double','on')
    axis equal; axis manual; grid on;
    title(['Symmetric flappying, \theta_{amp}=0.1',apstr,', t=',num2str((i-1)*dt)]);
%     title(['Symmetric flappying, \theta_{amp}=0.150, t=',num2str((i-1)*dt)]);
    xlabel('x'); ylabel('y');
    %     pause(.1)
    
    if mod(i,100)==0
        pause(1);
        strs=[str,'S_01',apstr,'T',num2str(i)]
%         strs=[str,'S_0150T',num2str(i)]
        saveas(h,strs,'fig');
        saveas(h,strs,'png');
    end
    
    %     currFrame = getframe(h);
    %     writeVideo(vidObj,currFrame);
end
% close(vidObj);
end

function index_piece(tk,n,ib)
global Cirl_vortex Wing
pos=0;neg=0; %index of pos or neg piece
piece_p=[];piece_n=[];
for k=1:n
    if k==1
        if Cirl_vortex(ib,k)>=0 %determine the first piece
            pos=1; piece_p(pos,1)=k; %start of pos
        else
            neg=1; piece_n(neg,1)=k; %start of neg
        end
    end
    if k<n&&Cirl_vortex(ib,k)>=0&&Cirl_vortex(ib,k+1)<0
        neg=neg+1; %pos piece ends, neg piece begins
        piece_p(pos,2)=k; piece_n(neg,1)=k;  %end of pos, start of neg
    else if k<n&&Cirl_vortex(ib,k)<0&&Cirl_vortex(ib,k+1)>=0
            pos=pos+1; %neg piece ends, pos piece begins
            piece_n(neg,2)=k; piece_p(pos,1)=k;  %end of pos, start of neg
        else if k==n&&Cirl_vortex(ib,k)>=0 %final end is pos
                piece_p(pos,2)=k;
            else if k==n&&Cirl_vortex(ib,k)<0 %final end is neg
                    piece_n(neg,2)=k;
                end
            end
        end
    end
end
msize=3;
for k=1:pos
    plot(-imag(Wing(ib).ZetaF{tk}(piece_p(k,1):piece_p(k,2))),-real(Wing(ib).ZetaF{tk}(piece_p(k,1):piece_p(k,2))),...
        'r-','MarkerSize',msize);hold on;
end
for k=1:neg
    plot(-imag(Wing(ib).ZetaF{tk}(piece_n(k,1):piece_n(k,2))),-real(Wing(ib).ZetaF{tk}(piece_n(k,1):piece_n(k,2))),...
        'b-','MarkerSize',msize);hold on;
end
end
