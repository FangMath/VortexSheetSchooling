% - plot the vortex sheet strength \gamma(s) on chord and free sheet, vs s
%   and time (showed in movie)
% - plot the speed of vortex sheding d\Gamma(t)/dt and CXVelocity*gamma(2)
%   vs t, tried to verify the Kutta condition

function PlotStrength(str_input)
close all;

if nargin~=0
    str=str_input;
else
    str='./output/';
end

load( [str,'parameter.mat'],'Para');
load( [str,'data.mat'],'tem_at','Wing');

dt=Para.dt; M=Para.M; s=Para.s; FreeVelocity=Para.FreeVelocity;
n=tem_at;

for tk=2:tem_at
    for ib=1:Para.Nw
    zetaf=Wing(ib).ZetaF{tk}; 
    gammaf=Wing(ib).GammaF{tk};
    nu=Wing(ib).Nu(tk,:);    gammaMax=Wing(ib).GammaMax(tk);
    K=length(zetaf);
    
    %% - plot the vortex sheet strength \gamma(s) on chord and free sheet
    
    %calculate vortex strength on bounded sheet
    for j=2:M
        s(j)=cos((j-1)*pi/M);
        gamma(j)=nu(j)/sqrt(1-s(j)^2); %vortex strength on bounded sheet
    end
    gamma(M+1)=gamma(M);gamma(1)=gamma(2);%gamma(M+1) and gamma(1) are NaN, take value of neighbours
    
    
    distance=zeros(1,K); %arclength on free sheet
    for k=K-1:-1:1
        distance(k)=abs(zetaf(k)-zetaf(k+1))+distance(k+1);
    end
    % calculate vortex strength on free sheet, i.e. d\Gamma(s)/ds, by center differencing
    gammafree=[];
    for m=2:K-1
        gammafree(m)=(gammaf(m+1)-gammaf(m-1))/(abs(zetaf(m+1)-zetaf(m-1)));
    end
    gammafree(1)=(gammaf(2)-gammaf(1))/(abs(zetaf(2)-zetaf(1)));
    gammafree(K)=(gammaf(K)-gammaf(K-1))/(abs(zetaf(K)-zetaf(K-1)));
    
    fig1=figure(1);
    if ib==1
    plot(s,gamma,'r-',distance+1,-gammafree,'.-');
    else
            plot(s,-gamma,'b-',distance+1,gammafree,'.-');
    end
    axis([-2 20 -20 20]);hold on;
    title(['t=',num2str(tk*dt)]);
    set(gcf,'double','on');
    axis manual;grid on;
    pause(0.3);
    
    %% - plot the speed of vortex sheding d\Gamma(t)/dt and CXVelocity*gamma(2)
    %DOTGAMMA(tk)=d\Gamma(tk)/dt
    if tk==2
        DOTGAMMA(tk)=DERIVATIVE(gammaMax,Wing(ib).GammaMax(tk-1),0,0,dt,111);
    else
        DOTGAMMA(tk)=DERIVATIVE(gammaMax,Wing(ib).GammaMax(tk-1),Wing(ib).GammaMax(tk-2),0,dt,10);
    end
    DOTgamma(tk)=-gammafree(K);
%     DOTgamma(tk)=-gamma(2);
    end
    hold off;
end
fig2=figure(2);
plot(0:dt:(tem_at-1)*dt,DOTGAMMA,'r.-',0:dt:(tem_at-1)*dt,-FreeVelocity*DOTgamma,'b.-'); grid on;
