% plot the pressure jump distribution on chord
% # Function B.m is called 

function PlotPressure(str_input,ib)
global Para Sys wing Wing
global cx cy omega tk dcx dcy domega
close all;


dt=Para.dt; M=Para.M; FreeVelocity=Para.FreeVelocity; s=Para.s;

if nargin~=0
    str=str_input
else
    str='./output/';
end

load( [str,'parameter.mat'],'Para');
load( [str,'data.mat'],'tem_at','Sys','Wing'); 

n=tem_at;
Pressure_dist=zeros(n,M+1);
n
for tk=2:n
    % calculate tau(j), mu(j), j=1,...M+1 and d\nu/dt, d\Gamma/dt
    K=Wing(ib).Lenf(tk);
    wing(ib).K=K;
    wing(ib).zetaf=Wing(ib).ZetaF{tk}(1:K);
    wing(ib).gammaf=Wing(ib).GammaF{tk}(1:K);
    omega=Sys.Omega(tk);
    cx=Sys.CenterX(tk);
    cy=Sys.CenterY(tk);
    domega=Sys.dOmega(tk);
    dcx=Sys.dCenterX(tk);
    dcy=Sys.dCenterY(tk);
    
    FLATPOSITION;
  
    tau=zeros(M+1,1);   mu=zeros(M+1,1);    dotnu=zeros(M+1,1);
    for j=1:M+1
  
        tau(j)=real(wing(ib).tangential*conj(wing(ib).dzetabody(j)));
        mu(j)=real(wing(ib).tangential*(B(ib,j)+Para.FreeVelocity(tk-1)));
        
        nu=Wing(ib).Nu(tk,:);    gammaMax=Wing(ib).GammaMax(tk);
        
        if tk==2
            dotnu(j)=DERIVATIVE(nu(j),Wing(ib).Nu(tk-1,j),0,0,dt,111);
            dotgammaMax=DERIVATIVE(gammaMax,Wing(ib).GammaMax(tk-1),0,0,dt,111);
        else
            dotnu(j)=DERIVATIVE(nu(j),Wing(ib).Nu(tk-1,j),Wing(ib).Nu(tk-2,j),0,dt,10);
            dotgammaMax=DERIVATIVE(gammaMax,Wing(ib).GammaMax(tk-1),Wing(ib).GammaMax(tk-2),0,dt,10);
        end
    end
    
    % calculate pressure jump distribution, by trapezoidal rule
    pressure_dist=zeros(M+1,1);
    for j=2:M
        gamma(j)=nu(j)/sqrt(1-s(j)^2);
        dotgamma(j)=dotnu(j)/sqrt(1-s(j)^2);
        for k=2:j-1
            pressure_dist(j)=pressure_dist(j)+.5*(dotgamma(k)+dotgamma(k+1))*(s(k+1)-s(k));
        end
        pressure_dist(j)=pressure_dist(j)+(mu(j)-tau(j))*gamma(j)+dotgammaMax;
        
    end
    pressure_dist(1)=0; pressure_dist(M+1)=0;
    
    Pressure_dist(tk,:)=pressure_dist;
    
    %% make plots comparing linear theory and numerics on pressure jump distribution
    figure(1);
%     x=-1:.01:1;
%     y=sqrt((1-x)./(1+x));
%     z=sqrt(1-x.^2);
%     period=2*pi;
%     plot(s,Pressure_dist(tk,:),'r.-',x,-100*pi*sin(period*tk*dt)*y-0.4*cos(period*tk*dt)*z,'b.-');
    plot(s,Pressure_dist(tk,:),'r.-');
    title(['t=',num2str(dt*tk)]);
    axis([-1 1 -50 50]);grid on;
    set(gcf,'double','on');
    axis manual;
    pause(.1);
end

save([str, num2str(ib),'Pressure_dist.mat'], 'Pressure_dist');
