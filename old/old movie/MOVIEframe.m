function MOVIEframe(str_input,Von,apstr)
close all;
global Wing Sys tem_at Nw left right down up
if nargin~=0
    str=str_input;
else
    str='./output/';
end

load( [str,'parameter.mat'],'Para');
dt=Para.dt;Nw=Para.Nw;
load( [str,'data.mat'], 'tem_at','Wing','Sys');


% Initialize data at tk=tem_at(1):
start=1;

% set axis frame
left=-10;right=10; down=-9; up=2;

%set x, y for the central line
x=[0,0]; y=[down,up];


%calculate circulation for each vortex

tem_at


h=figure('Visible','off');
F(1:tem_at-1) = struct('cdata', [],'colormap', []);
%%
1
plot(Sys.CenterY(1),-Sys.CenterX(1),'ko','MarkerFaceColor','k',...
    'MarkerSize',6);hold on;
plot([imag(Wing(1).ZetaBody(1,end)),imag(Wing(2).ZetaBody(1,end))],[-real(Wing(1).ZetaBody(1,end)),-real(Wing(2).ZetaBody(1,end))],'-m');
x1=imag(Wing(1).ZetaBody(1,end)+Wing(2).ZetaBody(1,end))/2;     y1=-real(Wing(1).ZetaBody(1,end)+Wing(2).ZetaBody(1,end))/2;
x2=Sys.CenterY(1);  y2=-Sys.CenterX(1);
plot([x1,x2],[y1,y2],'m-');

plot(Sys.CenterY(1),-Sys.CenterX(1),'co','MarkerFaceColor','c',...
    'MarkerSize',6);hold on;
plot([imag(Wing(1).ZetaBody(1,end)),imag(Wing(2).ZetaBody(1,end))],[-real(Wing(1).ZetaBody(1,end)),-real(Wing(2).ZetaBody(1,end))],'-c');
plot([Sys.CenterY(1),imag(Wing(1).ZetaBody(1,end)+Wing(2).ZetaBody(1,end))/2],[-Sys.CenterX(1),-real(Wing(1).ZetaBody(1,end)+Wing(2).ZetaBody(1,end))/2],'-c');
plot(imag(Wing(1).ZetaBody(1,:)),-real(Wing(1).ZetaBody(1,:)),'c-','LineWidth',3);hold on;
plot(imag(Wing(2).ZetaBody(1,:)),-real(Wing(2).ZetaBody(1,:)),'c-','LineWidth',3);hold on;
plot(x,y,'g--');hold on;

PreMOVIE(.2,1);
for ib=1:Nw
    plot(imag(Wing(ib).ZetaBody(1,:)),-real(Wing(ib).ZetaBody(1,:)),'k-','LineWidth',3);hold on;
end
if Von==1
    [pot,target]=VField(1);
    quiver(target(2,:),-target(1,:),imag(pot)/3,-real(pot)/3,0);hold on;
end
hold off;
axis([left right down up]);
set(gcf,'double','on');
axis equal; axis manual; grid on;
title(['Symmetric flappying, \theta_{amp}=0.1',apstr,', t=',num2str((tk-1)*dt)]);
xlabel('x'); ylabel('y');


%%
set(gca,'NextPlot','replacechildren');
for tk=1:tem_at-1
    tk
    plot(Sys.CenterY(tk),-Sys.CenterX(tk),'ko','MarkerFaceColor','k',...
        'MarkerSize',6);hold on;
    plot([imag(Wing(1).ZetaBody(tk,end)),imag(Wing(2).ZetaBody(tk,end))],[-real(Wing(1).ZetaBody(tk,end)),-real(Wing(2).ZetaBody(tk,end))],'-m');
    x1=imag(Wing(1).ZetaBody(tk,end)+Wing(2).ZetaBody(tk,end))/2;     y1=-real(Wing(1).ZetaBody(tk,end)+Wing(2).ZetaBody(tk,end))/2;
    x2=Sys.CenterY(tk);  y2=-Sys.CenterX(tk);
    plot([x1,x2],[y1,y2],'m-');
    
    plot(Sys.CenterY(1),-Sys.CenterX(1),'co','MarkerFaceColor','c',...
        'MarkerSize',6);hold on;
    plot([imag(Wing(1).ZetaBody(1,end)),imag(Wing(2).ZetaBody(1,end))],[-real(Wing(1).ZetaBody(1,end)),-real(Wing(2).ZetaBody(1,end))],'-c');
    plot([Sys.CenterY(1),imag(Wing(1).ZetaBody(1,end)+Wing(2).ZetaBody(1,end))/2],[-Sys.CenterX(1),-real(Wing(1).ZetaBody(1,end)+Wing(2).ZetaBody(1,end))/2],'-c');
    plot(imag(Wing(1).ZetaBody(1,:)),-real(Wing(1).ZetaBody(1,:)),'c-','LineWidth',3);hold on;
    plot(imag(Wing(2).ZetaBody(1,:)),-real(Wing(2).ZetaBody(1,:)),'c-','LineWidth',3);hold on;
    plot(x,y,'g--');hold on;
    
    PreMOVIE(.2,tk);
    for ib=1:Nw
        plot(imag(Wing(ib).ZetaBody(tk,:)),-real(Wing(ib).ZetaBody(tk,:)),'k-','LineWidth',3);hold on;
    end
    if Von==1
        [pot,target]=VField(tk);
        quiver(target(2,:),-target(1,:),imag(pot)/3,-real(pot)/3,0);hold on;
    end
    hold off;
    axis([left right down up]);
    set(gcf,'double','on');
    axis equal; axis manual; grid on;
    title(['Symmetric flappying, \theta_{amp}=0.1',apstr,', t=',num2str((tk-1)*dt)]);
    xlabel('x'); ylabel('y');
    F(tk)=getframe(h);
    
    if Von==1
        save([str_input,'movieV_01',apstr,'.mat'],'F');
    else
        save([str_input,'movie_01',apstr,'.mat'],'F');
    end
end

successfulFrame=1
% if Von==1
%     movie2avi(F,[str_input,'movieV_01',apstr,'.avi'],'fps',50);
% else
%     movie2avi(F,[str_input,'movie_01',apstr,'.avi'],'fps',50);
% end
end