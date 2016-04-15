function [pot,target]=VField(tk)
% clearvars -except input_str;
% close all;
global tem_at Wing Sys Vp Para left right down up meshref
% str=input_str

%import parameter datas
% load( [str,'parameter.mat'],'Para');
% load( [str,'data.mat'], 'tem_at','Wing','Sys','Vp');


% set axis frame
dt=Para.dt;
% left=-2;right=9; down=-10; up=10;
% meshref=.5;
[Xgrid,Ygrid]=ndgrid(-up:meshref:-down,left:meshref:right);
sizeXY=size(Xgrid);
target = zeros(2,sizeXY(1)*sizeXY(2));
target(1,:)=reshape(Xgrid,1,sizeXY(1)*sizeXY(2));
target(2,:)=reshape(Ygrid,1,sizeXY(1)*sizeXY(2));
ntarget = length(target);
% num=ntarget
% plot(target(1,1:num),target(2,1:num),'.-')

% Z=Xgrid+Ygrid;
% surf(Xgrid,Ygrid,Z)
% (targetp, zetaf, gammaf, zetabodyw, nuw)

% for tk=1:25%tem_at-1
    tk
    zetaf1=Wing(1).ZetaF{tk};
    gammaf1=Wing(1).GammaF{tk};
    nuw1=Wing(1).Nu(tk,:);
    zetabodyw1=Wing(1).ZetaBody(tk,:);
    nsource1 = length(zetaf1)+length(zetabodyw1);
    zetaf2=Wing(2).ZetaF{tk};
    gammaf2=Wing(2).GammaF{tk};
    nuw2=Wing(2).Nu(tk,:);
    zetabodyw2=Wing(2).ZetaBody(tk,:);
    nsource2 = length(zetaf2)+length(zetabodyw2);
    
    source1=0;source2=0;
    
    source1(1,1:nsource1)=real([zetaf1,zetabodyw1]);
    source1(2,1:nsource1)=imag([zetaf1,zetabodyw1]);
    
    source2(1,1:nsource2)=real([zetaf2,zetabodyw2]);
    source2(2,1:nsource2)=imag([zetaf2,zetabodyw2]);
    
    
    M1=nsource1;
    X1=zeros(M1,1);
    N1=length(zetaf1);
    
    if N1==1
        X1=0;
    else
        X1(1)=gammaf1(2)-gammaf1(1);
        if N1>2
            X1(2:N1-1)=gammaf1(3:N1)-gammaf1(1:N1-2);
        end
        X1(N1)=gammaf1(N1)-gammaf1(N1-1);
    end
    X1=-X1/(4*1i*pi);
    X1(N1+1:M1)=nuw1/(-2*(M1-N1-1)*1i);
    X1([N1+1,M1])=0.5*X1([N1+1,M1]);
    dipstr1=X1;
    
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%
    M2=nsource2;
    X2=zeros(M2,1);
    N2=length(zetaf2);
    
    if N2==1
        X2=0;
    else
        X2(1)=gammaf2(2)-gammaf2(1);
        if N2>2
            X2(2:N2-1)=gammaf2(3:N2)-gammaf2(1:N2-2);
        end
        X2(N2)=gammaf2(N2)-gammaf2(N2-1);
    end
    X2=-X2/(4*1i*pi);
    X2(N2+1:M2)=nuw2/(-2*(M2-N2-1)*1i);
    X2([N2+1,M2])=0.5*X2([N2+1,M2]);
    dipstr2=X2;
    
    %% %%%%%%%%%%%%
    ifpot = 0;
    ifgrad = 0;
    ifhess = 0;
    ifpottarg = 1;
    ifgradtarg = 0;
    ifhesstarg = 0;
    
    % % % tic
    iprec=4;
    [U1]=zfmm2dpart(iprec,nsource1,source1,dipstr1,ifpot,ifgrad,ifhess,ntarget,target,ifpottarg,ifgradtarg,ifhesstarg);
    [U2]=zfmm2dpart(iprec,nsource2,source2,dipstr2,ifpot,ifgrad,ifhess,ntarget,target,ifpottarg,ifgradtarg,ifhesstarg);
    % % % total_time=toc
    
    pot=U1.pottarg;
    pot=pot+U2.pottarg;
    pot=pot+Para.FreeVelocity(tk);
    pot=conj(pot);
%     h=figure(1);
%     quiver(target(2,:),-target(1,:),imag(pot),-real(pot));
%     title([' t=',num2str((tk-1)*dt)]);
%     set(gcf,'double','on');
%     axis equal; axis manual;% grid on;
%     axis([left right down up]);
%     pause(1);
% end
end

