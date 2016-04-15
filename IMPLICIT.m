function IMPLICIT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   IMPLICIT.m
%
%
%   Written by Fang Fang - 03/19/2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global tk wing_p1 wing_p2 wing_p3
global wing Para
global Nw s M
%% --implicit nonlinear solver for bounded sheet
for ib=1:Para.Nw
    wing(ib).K=wing(ib).K+1; %add one node on free sheet
end

%% --use explicit as initial guess
if tk<4
    for ib = 1:Nw
        wing(ib).cx = wing_p1(ib).cx;
        wing(ib).cy = wing_p1(ib).cy;
        wing(ib).omega = wing_p1(ib).omega;
        wing(ib).dotcx = wing_p1(ib).dotcx;
        wing(ib).dotcy = wing_p1(ib).dotcy;
        wing(ib).domega = wing_p1(ib).domega;
    end
else
    for ib = 1:Nw
        wing(ib).cx = wing_p1(ib).cx + 1.5*Para.dt*wing_p1(ib).dotcx - 0.5*Para.dt*wing_p2(ib).dotcx;
        wing(ib).cy = wing_p1(ib).cy(1);
        wing(ib).omega = wing_p1(ib).omega(1);
        
        %%% wing1 and wing2 are separate bodies
        wing(ib).dotcx = wing_p1(ib).dotcx + 1.5*Para.dt*wing_p1(ib).thrust/Para.Mass-0.5*Para.dt*wing_p2(ib).thrust/Para.Mass;
        
        %%% wing1 and wing2 is a whole body
        %     wing(ib).dotcx = wing_p1(ib).dotcx+1.5*Para.dt*(wing_p1(1).thrust+wing_p1(2).thrust)/Para.Mass-0.5*Para.dt*(wing_p2(1).thrust+wing_p2(2).thrust)/Para.Mass;
        
        wing(ib).dotcy = wing_p1(ib).dotcy(1);
        wing(ib).domega = wing_p1(ib).domega(1);
    end
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
            X(ib,1:M+1)=2*wing_p1(ib).nu-wing_p2(ib).nu;
            X(ib,M+2)=2*wing_p1(ib).gammaMax-wing_p2(ib).gammaMax(tk-2);
        end
        X=reshape(X(1:Nw,1:M+2),1,(M+2)*Nw);
        X((M+2)*Nw+1)=wing(1).cx;
        X((M+2)*Nw+2)=wing(1).dotcx;
        X((M+2)*Nw+3)=wing(2).cx;
        X((M+2)*Nw+4)=wing(2).dotcx;
    else
        X(ib,1:M+1) = 3*wing_p1(ib).nu - 3*wing_p2(ib).nu + wing_p3(ib).nu;
        X(ib,M+2)   = 3*wing_p1(ib).gammaMax - 3*wing_p2(ib).gammaMax + wing_p3(ib).gammaMax;
        X=reshape(X(1:Nw,1:M+2),1,(M+2)*Nw);
        X((M+2)*Nw+1)=wing(1).cx;
        X((M+2)*Nw+2)=wing(1).dotcx;
        X((M+2)*Nw+3)=wing(2).cx;
        X((M+2)*Nw+4)=wing(2).dotcx;
    end
end

a = X((M+2)*Nw+2)
% implicit solver
broyden(X,@f,(M+2)*Nw+4,1e-10);
end
