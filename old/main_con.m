% This is the main function
function []=main_con(input_str)
clearvars -except input_str;
close all;
% Nf = maxNumCompThreads(2);
Ttotal=tic;
%% -define global variables
global tk cx cy omega dcx dcy domega
global Wing wing Sys Para
global Delta1 Delta0 eta p Nw s M Vp PointVVp StartNew data

% prompt = 'Use divided data? 1(yes) 0(no) ';
% UseDividedData = input(prompt)
%     
% strd=[str,'/Ddata/'];
% if UseDividedData==1
%     if ~exist(strd,'dir')
%         UseDividedData=0;
%         fprintf('Can not use divided data!');
%     end
% end

str=input_str

% -start from any time using existing data
if nargin~=0
    
    %import parameter datas
    load( [str,'parameter.mat'],'Para','ifSaveDividedData');
    load( [str,'data.mat'], 'tem_at','Wing','Sys','Vp','PointVVp','StartNew');
    strd=[str,'/Ddata/'];
    
    % Initialize data at tk=tem_at(1):
    tk=tem_at;
%     cx=Sys.CenterX(tk);dcx=Sys.dCenterX(tk);
%     cy=Sys.CenterY(tk);dcy=Sys.dCenterY(tk);
%     omega=Sys.Omega(tk);domega=Sys.dOmega(tk);
    
    for ib=1:Para.Nw
        %read in data
        wing(ib).zetabody=Wing(ib).ZetaBody(tk,:);
        wing(ib).nu=Wing(ib).Nu(tk,:);
        wing(ib).zetaf=Wing(ib).ZetaF{tk};
        wing(ib).gammaf=Wing(ib).GammaF{tk};
        wing(ib).gammaMax=Wing(ib).GammaMax(tk);
        wing(ib).K=length(wing(ib).zetaf);
        wing(ib).PointVX=Wing(ib).PointVX{tk};
        wing(ib).PointVCirl=Wing(ib).PointVCirl{tk};
        wing(ib).cx = Wing(ib).Cx(tk);
        wing(ib).dotcx = Wing(ib).dotCx(tk);
        wing(ib).cy = Wing(ib).Cy(tk);
        wing(ib).dotcy = Wing(ib).dotCy(tk);
        wing(ib).omega = Wing(ib).Omega(tk);
        wing(ib).domega = Wing(ib).dOmega(tk);
        
    end
%     Sys.CenterX(tk)=cx;
%     Sys.CenterY(tk)=cy;
%     Sys.Omega(tk)=omega;
%     
%     Sys.dCenterX(tk)=dcx;
%     Sys.dCenterY(tk)=dcy;
%     Sys.dOmega(tk)=domega;
    
    Delta1=Para.Delta1; Delta0=Para.Delta0; eta=Para.eta; p=Para.p;
    Nw=Para.Nw; s=Para.s; M=Para.M;
    
else
    error('you should use main.m!')
end

%% -start the main solver
MaxT=10005; %number of time steps
while tk<=MaxT
    step=tic;
    tex=tic;
    % --set background flow initially
%     Para.FreeVelocity(tk)=max(Para.iniV-Para.iniV/Para.numV*tk*Para.dt,0);
    Para.FreeVelocity(tk) = Para.iniV;
    tk %result would be at (tk+1)^th time step
    %% --explicit solver for free sheet
    EXPLICIT(Para);
    Texplicit=toc(tex)
%     figure(5)
%     for ib=1:Nw
% index_piece(wing,ib,1);
%     end    
    tic
    %% --implicit nonlinear solver for bounded sheet
    tk=tk+1; %update to next time step
    for ib=1:Para.Nw
        wing(ib).K=wing(ib).K+1; %add one node on free sheet
    end
    
%% --use explicit as initial guess
    if tk<4
%         cx=Sys.CenterX(1);
%         cy=Sys.CenterY(1);
%         omega=Sys.Omega(1);
%         
%         dcx=Sys.dCenterX(1);
%         dcy=Sys.dCenterY(1);
%         domega=Sys.dOmega(1);
for ib = 1:Nw
    wing(ib).cx = Wing(ib).Cx(1);
    wing(ib).cy = Wing(ib).Cy(1);
    wing(ib).omega = Wing(ib).Omega(1);
     wing(ib).dotcx = Wing(ib).dotCx(1);
    wing(ib).dotcy = Wing(ib).dotCy(1);
    wing(ib).domega = Wing(ib).dOmega(1);
end
    else
        for ib = 1:Nw
%     wing(ib).cx = Wing(ib).Cx(1);
    wing(ib).cx = Wing(ib).Cx(tk-1)+1.5*Para.dt*Wing(ib).dotCx(tk-1)-0.5*Para.dt*Wing(ib).dotCx(tk-2);
    wing(ib).cy = Wing(ib).Cy(1);
    wing(ib).omega = Wing(ib).Omega(1);
%      wing(ib).dotcx = Wing(ib).dotCx(1);

%%% wing1 and wing2 are seperate bodies
%     wing(ib).dotcx = Wing(ib).dotCx(tk-1)+1.5*Para.dt*(Wing(1).Thrust(tk-1)+Wing(2).Thrust(tk-1))/Para.Mass-0.5*Para.dt*(Wing(1).Thrust(tk-2)+Wing(2).Thrust(tk-2))/Para.Mass;

%% wing1 and wing2 is a whole body 
     wing(ib).dotcx = Wing(ib).dotCx(tk-1)+1.5*Para.dt*(Wing(1).Thrust(tk-1)+Wing(2).Thrust(tk-1))/Para.Mass-0.5*Para.dt*(Wing(1).Thrust(tk-2)+Wing(2).Thrust(tk-1))/Para.Mass;

    wing(ib).dotcy = Wing(ib).dotCy(1);
    wing(ib).domega = Wing(ib).dOmega(1);
        end
        
%         cx=Sys.CenterX(tk-1)+1.5*Para.dt*Sys.dCenterX(tk-1)-0.5*Para.dt*Sys.dCenterX(tk-2);
%         dcx=Sys.dCenterX(tk-1)+1.5*Para.dt*(Sys.Thrust(tk-1)+Para.g)/Para.Mass-0.5*Para.dt*(Sys.Thrust(tk-2)+Para.g)/Para.Mass;
%         
%         cy=0;
%         dcy=0;
%         
%         omega=0;
%         domega=0;

%        cy=Sys.CenterY(tk-1)+1.5*Para.dt*Sys.dCenterY(tk-1)-0.5*Para.dt*Sys.dCenterY(tk-2);
%        dcy=Sys.dCenterY(tk-1)+1.5*Para.dt*Sys.Lift(tk-1)/Para.Mass-0.5*Para.dt*Sys.Lift(tk-2)/Para.Mass;
%        
%        omega=Sys.Omega(tk-1)+1.5*Para.dt*Sys.dOmega(tk-1)-0.5*Para.dt*Sys.dOmega(tk-2);
%        domega=Sys.dOmega(tk-1)+1.5*Para.dt*Sys.Torque(tk-1)/Para.I-0.5*Para.dt*Sys.Torque(tk-2)/Para.I;
    end
    
    % Set initial guess for broyden's method
    X=zeros(2,M+2);
    if tk==2
        for ib=1:Nw
            X(ib,1:M+1)=0;
            X(ib,M+2)=0;
        end
        X=reshape(X(1:Nw,1:M+2),1,(M+2)*Nw);
        X((M+2)*Nw+1)=wing(1).cx;
        X((M+2)*Nw+2)=wing(1).dotcx;
        X((M+2)*Nw+3)=wing(2).cx;
        X((M+2)*Nw+4)=wing(2).dotcx;
%         X((M+2)*Nw+5)=wing(1).omega;
%         X((M+2)*Nw+6)=wing(1).domega;
    else if tk==3
            for ib=1:Nw
                X(ib,1:M+1)=2*Wing(ib).Nu(tk-1,:)-Wing(ib).Nu(tk-2,:);
                X(ib,M+2)=2*Wing(ib).GammaMax(tk-1)-Wing(ib).GammaMax(tk-2);
            end
            X=reshape(X(1:Nw,1:M+2),1,(M+2)*Nw);
%             X((M+2)*Nw+1)=cx;
%             X((M+2)*Nw+2)=dcx;
%             X((M+2)*Nw+3)=cy;
%             X((M+2)*Nw+4)=dcy;
%             X((M+2)*Nw+5)=omega;
%             X((M+2)*Nw+6)=domega;
        X((M+2)*Nw+1)=wing(1).cx;
        X((M+2)*Nw+2)=wing(1).dotcx;
        X((M+2)*Nw+3)=wing(2).cx;
        X((M+2)*Nw+4)=wing(2).dotcx;
        else
            for ib=1:Nw
                X(ib,1:M+1)=3*Wing(ib).Nu(tk-1,:)-3*Wing(ib).Nu(tk-2,:)+Wing(ib).Nu(tk-3,:);
                X(ib,M+2)=3*Wing(ib).GammaMax(tk-1)-3*Wing(ib).GammaMax(tk-2)+Wing(ib).GammaMax(tk-3);
            end
            X=reshape(X(1:Nw,1:M+2),1,(M+2)*Nw);
%             X((M+2)*Nw+1)=cx;
%             X((M+2)*Nw+2)=dcx;
%             X((M+2)*Nw+3)=cy;
%             X((M+2)*Nw+4)=dcy;
%             X((M+2)*Nw+5)=omega;
%             X((M+2)*Nw+6)=domega;
        X((M+2)*Nw+1)=wing(1).cx;
        X((M+2)*Nw+2)=wing(1).dotcx;
        X((M+2)*Nw+3)=wing(2).cx;
        X((M+2)*Nw+4)=wing(2).dotcx;
        end
    end
    
    a = X((M+2)*Nw+2)
% implicit solver
    broyden(X,@f,(M+2)*Nw+4,1e-10);
    toc
    
    tic
    
    %%
    %     fig=figure(1);
    %     for ib=Para.Nw
    % %         plot(1:81,wing(1).nu,'r.-',40,wing(1).gammaMax,'r*',1:81,wing(2).nu,'b.-',40,wing(2).gammaMax,'b*');hold on;
    %         plot(1:81,wing(1).nu+wing(2).nu,'r.-',40,wing(1).gammaMax+wing(2).gammaMax,'r*');hold on;
    %     end
    %     title(num2str(tk));
    %     hold off;
    
    %% --record the solutions
    for ib=1:Nw
        Wing(ib).Nu(tk,:)=wing(ib).nu;
        Wing(ib).ZetaBody(tk,:)=wing(ib).zetabody;
        %         for j=1:wing(ib).K
        Wing(ib).ZetaF{tk}(1:wing(ib).K)=wing(ib).zetaf(1:wing(ib).K);
        Wing(ib).GammaF{tk}(1:wing(ib).K)=wing(ib).gammaf(1:wing(ib).K);
        %         end
        Wing(ib).GammaMax(tk)=wing(ib).gammaMax;
        Wing(ib).Lenf(tk)=length(wing(ib).zetaf);
        Wing(ib).PointVX{tk}=wing(ib).PointVX;
        Wing(ib).PointVCirl{tk}=wing(ib).PointVCirl;
        Wing(ib).Cx(tk) = wing(ib).cx;
        Wing(ib).dotCx(tk) = wing(ib).dotcx;
        Wing(ib).Cy(tk) = wing(ib).cy;
        Wing(ib).dotCy(tk) = wing(ib).dotcy;
        Wing(ib).Omega(tk) = wing(ib).omega;
        Wing(ib).dOmega(tk) = wing(ib).domega;
    end
    
    gammaMaxLeft(tk)=-Wing(2).GammaMax(tk);
    gammaMaxRight(tk)=-Wing(1).GammaMax(tk);
    
%     Sys.CenterX(tk)=cx;
%     Sys.CenterY(tk)=cy;
%     Sys.Omega(tk)=omega;
%     
%     Sys.dCenterX(tk)=dcx;
%     Sys.dCenterY(tk)=dcy;
%     Sys.dOmega(tk)=domega;
    
    %record the current terminate time
    tem_at=tk;
    % write data into files
    w = whos('Wing');
    fprintf('mem: %6.2f (MB)\n', w.bytes/1e6)
    
%     if SaveAllData==1
        if mod(tk,20) == 0
            fprintf('save data!!!');
        if w.bytes/1e6 >2000
            fprintf('mem: %6.2f (MB)\n', w.bytes/1e6)
            save( [str,'/data.mat'],'tem_at','Sys','Wing','Vp','PointVVp','StartNew','-v7.3');
        else
            save( [str,'/data.mat'],'tem_at','Sys','Wing','Vp','PointVVp','StartNew');
        end
        end
%     end
    %%
    if ifSaveDividedData==1
        save([strd,'/tem.mat'],'tem_at','Sys');
        
        strs=[strd,'T',num2str(tk)];
        for ib=1:Para.Nw
            %wing(ib).nu=Wing(ib).Nu(tk,:);
            %wing(ib).zetabody=Wing(ib).ZetaBody(tk,:);
            wing(ib).Lenf=Wing(ib).Lenf(tk);
            %wing(ib).zetaf=Wing(ib).ZetaF{tk}(1:wing(ib).Lenf);
            %wing(ib).gammaf=Wing(ib).GammaF{tk}(1:wing(ib).Lenf);
            %wing(ib).gammaMax=Wing(ib).GammaMax(tk);
        end
        
%         sys.CenterX=Sys.CenterX(tk);
%         sys.CenterY=Sys.CenterY(tk);
%         sys.Omega=Sys.Omega(tk);
%         
%         sys.dCenterX=Sys.dCenterX(tk);
%         sys.dCenterY=Sys.dCenterY(tk);
%         sys.dOmega=Sys.dOmega(tk);
        save([strs,'.mat'],'tem_at','wing');
    end
    
    %%
    save([str,'/GMAX.mat'],'gammaMaxLeft','gammaMaxRight');
    save( [str,'/parameter.mat'], 'Para','ifSaveDividedData');
    
%     fK1=wing(1).K
%     fK2=wing(2).K
    Tstep=toc(step)
end

TimeTotal=toc(Ttotal)
end
