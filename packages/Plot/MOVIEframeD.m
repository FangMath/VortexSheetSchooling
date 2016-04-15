function MOVIEframeD(str_input, st, ed)
close all;
global Wing Sys tem_at Nw left right down up meshref wing

str=str_input;
Von = 1;

load( [str,'parameter.mat'],'Para');
dt=Para.dt;Nw=Para.Nw;

% set axis frame
plotlength = 15;
up = 2; down = -up; left = -503; right = left + plotlength; scl = 5;

lg = 1;
if lg 
plotlength = 25;
up = 2; down = -up; left = -1; right = left + plotlength;
scl = 10;
end

sm = 0;
if sm 
plotlength = 15;
up = 2; down = -up; left = -503; right = left + plotlength; scl = 5;
end


%  plot body and free sheet at each time step
str=[str,'Ddata/'];

%calculate circulation for each vortex
load( [str,'tem.mat'], 'tem_at');

h=figure(1);
    N = 1; %number of periods to present
    %ed = 8000%tem_at-2
    ed = 5225%tem_at-2
    %ed = st + N*100;
    F(1:ed-st+1) = struct('cdata', [],'colormap', []);

for tk=st:ed
plotset(10,3);
if lg
plotset(15,3);
end
    tk
    load( [str,'T',num2str(tk+1),'.mat'], 'wing');
    for ib=1:Nw
        Wing(ib).Nu(tk,:)=wing(ib).nu;
        Wing(ib).ZetaBody(tk,:)=wing(ib).zetabody;
        Wing(ib).Lenf(tk)=length(wing(ib).zetaf);
        Wing(ib).ZetaF{tk}=wing(ib).zetaf;
        Wing(ib).GammaF{tk}=wing(ib).gammaf;
        Wing(ib).GammaMax(tk)=wing(ib).gammaMax;
    end
    


%        plot(real(Wing(1).ZetaBody(st-200:tk,41)),imag(Wing(1).ZetaBody(st-200:tk,41)),'m-');hold on;
%fk = floor((tk-st)/25); st + 25*[0:fk]; plot(real(Wing(1).ZetaBody(st+25*[0:fk],41)),imag(Wing(1).ZetaBody(st+25*[0:fk],41)),'b.');hold on;



%        plot(real(Wing(2).ZetaBody(st-200:tk,41)),imag(Wing(2).ZetaBody(st-200:tk,41)),'g-');hold on;
%fk = floor((tk-st)/25); st + 25*[0:fk]; plot(real(Wing(2).ZetaBody(st+25*[0:fk],41)),imag(Wing(2).ZetaBody(st+25*[0:fk],41)),'r.');hold on;

    left=wing(1).cx-1; right=left+plotlength;
if lg
    left=wing(1).cx-1; right=left+plotlength;
end
if sm 
    left=wing(2).cx-5; right=left+plotlength;
end
Von = 0
    if Von==1
        xx = [floor(left)-1 : .3 : ceil(right)+1];
        yy = [down-.5 : .3 : up+.5];
        %xx = [ceil(Wing(1).ZetaBody(tk,1))+.5 : .05 : floor(Wing(2).ZetaBody(tk,end))-.5];
        %yy = [-.05: .01 : .05];
        target = bsxfun(@plus, xx', 1i*yy);
        size(target)
        target = reshape(target,1,[]);
        size(target)
        [pot]=VField2(1, tk, target);
        qver = 2;
        %%%%%%%% for plot %%%%%%%
        [X,Y] = meshgrid(xx,yy);
        reshapePot = reshape(pot, length(xx), []).';
        %colormap(jet);
        %srf = surf(X, Y, imag(reshapePot)); shading interp; colorbar; caxis([-2,2]); uistack(srf,'bottom');
        [C,hf]=contourf(X, Y, imag(reshapePot),100); 
        clbar = colorbar('position', [0.92 0.15 .02 .7]); 
        ylabel(clbar, 'vertical flow velocity', 'fontsize', 16)
        caxis([-1.5,1.5]); 
        set(hf,'LineColor', 'none');
        centerT = mean(reshape(target,length(xx),length(yy)),2);
        centerP = mean(reshape(pot,length(xx),length(yy)),2);
        %quiver(real(centerT),imag(centerT),real(centerP)/qver,imag(centerP)/qver,0);hold on;
        [maxP,index] = max(abs(imag(centerP)))
        %plot(real(centerT(index)), imag(centerT(index)), '*'); hold on;
        WaveTgt1(tk) = centerT(index);
        WaveMax1(tk) = maxP;
        save([str_input, 'wave.mat'], 'WaveTgt1', 'WaveMax1');

    end

    grey=[0.4,0.4,0.4];
    for ib=1:Nw
        plot(real(Wing(ib).ZetaBody(tk,:)),imag(Wing(ib).ZetaBody(tk,:)),'k-','LineWidth',5);hold on;
    end
        plot(real(Wing(2).ZetaBody(st:tk,81)),imag(Wing(1).ZetaBody(st:tk,81)),'--','Color',grey);hold on;

        plot(real(Wing(1).ZetaBody(st:tk,1)),imag(Wing(2).ZetaBody(st:tk,1)),'w--');hold on;

    ploton=1;
    index_piece(wing,1,ploton);
    index_piece(wing,2,ploton);


if lg
plotset(15,3);
end

ylim([down,up]);
xlim([left,right]);
set(gca,'position',[.03, .05, .87, .9]);
grid on;
    %title(['t=', num2str(tk/100)]);
    %ttl = title(['Simulation of freely moving flapping wings'], 'fontsize', 18);
    hold off;
    ps = 0.01;
    pause(ps)
    F(tk-st+1)=getframe(h);
end

str_inputF = 'movie';
save([str_inputF, '.mat'],'F');
MOVIEshow(str_inputF);

%if lg == 0
%str_inputF = [str_input, str_input(end-18:end-1), '_movie'];
%save([str_inputF, '.mat'],'F');
%MOVIEshow(str_inputF);
%end
%if lg
%save([str_input, str_input(end-18:end-1), '_movie.mat'],'F');
%MOVIEshow(str_input);
%end

end
