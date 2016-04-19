function [result]=f(x)
% Implicit solver, evaluate f(x), x is a vector of length M+2

% # called in broyden.m
% # Function CHEBYSHEVF.m is called
global wing Para tk Sys cx dcx cy dcy omega domega gamma1 rest bdelta fang Wing wing_p1 wing_p2 wing_p3
global gammaf zetaf K nu zetabody
M=Para.M; s=Para.s;   Nw=Para.Nw;

%% read in x
y(1:Nw,1:M+2)=reshape(x(1:(M+2)*Nw),Nw,M+2);

wing(1).cx = x((M+2)*Nw+1);
wing(1).dotcx = x((M+2)*Nw+2);
wing(2).cx=x((M+2)*Nw+3);
wing(2).dotcx=x((M+2)*Nw+4);
wing(2).cy=x((M+2)*Nw+5);
wing(2).dotcy=x((M+2)*Nw+6);
% omega=x((M+2)*Nw+5);
% domega=x((M+2)*Nw+6);

%% get the updated wing location info
wing = FLATPOSITION(wing, Para, tk);


for ib = 1:Nw
    % add the trailing edge in free sheet
    wing(ib).zetaf(wing(ib).K) = wing(ib).zetabody(1);
end

for ib = 1:Nw
    wing(ib).nu = y(ib,1:M+1);
    wing(ib).gammaMax = y(ib,M+2);
end

% add the trailing edge in free sheet
for ib=1:Nw
    wing(ib).gammaf(wing(ib).K)=wing(ib).gammaMax;
end

% Calculate the cos expansion of bounded part of f
for ibb=1:Nw
    K(ibb)=wing(ibb).K;
    zetaf(ibb,1:K(ibb))=wing(ibb).zetaf;
    gammaf(ibb,1:K(ibb))=wing(ibb).gammaf;
    nu(ibb,:)=wing(ibb).nu;
    zetabody(ibb,:)=wing(ibb).zetabody;
end

rest=zeros(2,M+1); bdelta=zeros(2,M+1);
for ib=1:Nw
    ChebyshevF(ib,:)=CHEBYSHEVF(ib);
end

%%
Sys.Torque(tk)=0; Sys.Lift(tk)=0; Sys.Thrust(tk)=0;
for ib=1:Nw
    [torque,lift,thrust]=TorqueForce(ib);
    %     Wing(ib).Torque(tk) = torque;
    %     Wing(ib).Lift(tk) = lift;
    %     Wing(ib).Thrust(tk) = thrust;
    %     Sys.Torque(tk)=Sys.Torque(tk)+torque;
    %     Sys.Lift(tk)=Sys.Lift(tk)+lift;
    %     Sys.Thrust(tk)=Sys.Thrust(tk)+thrust;
end

% a = Wing(1).Thrust(2)
%% evaluate the functions
result=zeros(Nw,M+1);
for ib=1:Nw
    f_SinSeriesHat = (2*M)/(2*1i)*[0, ChebyshevF(ib,2:M), ChebyshevF(ib,M+1), -ChebyshevF(ib,M:-1:2)];
    f_SinSeries = (ifft(f_SinSeriesHat));
    
    result(ib,1:M+1)=wing(ib).nu-2*real(f_SinSeries(1:M+1).*sqrt(1-Para.s.^2))+...
        ChebyshevF(ib,2)+2*ChebyshevF(ib,1)*Para.s-wing(ib).gammaMax/pi+...
        gamma1(ib)*(1-sqrt(1-Para.s.^2).*(pi-acos(Para.s))+ Para.s*log(2))/pi;
    
    result(ib,M+2)=-ChebyshevF(ib,2)-2*ChebyshevF(ib,1)+wing(ib).gammaMax/pi-gamma1(ib)*(1+log(2))/pi;
    
    fS(ib,:)=f_SinSeriesHat;
end

result=reshape(result(1:Nw,1:M+2),1,(M+2)*Nw);


%% wing1 and wing2 are seperate bodies
result((M+2)*Nw+1)=wing(1).cx-wing_p1(1).cx-Para.dt/2*(wing(1).dotcx+wing_p1(1).dotcx);

result((M+2)*Nw+2)=wing(1).dotcx-wing_p1(1).dotcx-Para.dt/(2*Para.Mass)*(wing(1).thrust+wing_p1(1).thrust);

result((M+2)*Nw+3)=wing(2).cx-wing_p1(2).cx-Para.dt/2*(wing(2).dotcx+wing_p1(2).dotcx);

result((M+2)*Nw+4)=wing(2).dotcx-wing_p1(2).dotcx-Para.dt/(2*Para.Mass)*(wing(2).thrust+wing_p1(2).thrust);

result((M+2)*Nw+5)=wing(2).cy-wing_p1(2).cy-Para.dt/2*(wing(2).dotcy+wing_p1(2).dotcy);

result((M+2)*Nw+6)=wing(2).dotcy-wing_p1(2).dotcy-Para.dt/(2*Para.Mass)*(wing(2).lift + wing_p1(2).lift);

%%%
%result((M+2)*Nw+1)=wing(1).cx-Wing(1).Cx(tk-1)-Para.dt/2*(wing(1).dotcx+Wing(1).dotCx(tk-1));
%
%result((M+2)*Nw+2)=wing(1).dotcx-Wing(1).dotCx(tk-1)-Para.dt/(2*Para.Mass)*(Wing(1).Thrust(tk)+Wing(1).Thrust(tk-1));
%
%result((M+2)*Nw+3)=wing(2).cx-Wing(2).Cx(tk-1)-Para.dt/2*(wing(2).dotcx+Wing(2).dotCx(tk-1));
%
%result((M+2)*Nw+4)=wing(2).dotcx-Wing(2).dotCx(tk-1)-Para.dt/(2*Para.Mass)*(Wing(2).Thrust(tk)+Wing(2).Thrust(tk-1));

%%% wing1 and wing2 is a whole body
%result((M+2)*Nw+1)=wing(1).cx-Wing(1).Cx(tk-1)-Para.dt/2*(wing(1).dotcx+Wing(1).dotCx(tk-1));
%
%result((M+2)*Nw+2)=wing(1).dotcx-Wing(1).dotCx(tk-1)-Para.dt/(2*Para.Mass)*(Wing(1).Thrust(tk)+Wing(1).Thrust(tk-1)+Wing(2).Thrust(tk)+Wing(2).Thrust(tk-1));
%
%result((M+2)*Nw+3)=wing(2).cx-Wing(2).Cx(tk-1)-Para.dt/2*(wing(2).dotcx+Wing(2).dotCx(tk-1));
%
%result((M+2)*Nw+4)=wing(2).dotcx-Wing(2).dotCx(tk-1)-Para.dt/(2*Para.Mass)*(Wing(1).Thrust(tk)+Wing(1).Thrust(tk-1)+Wing(2).Thrust(tk)+Wing(2).Thrust(tk-1));

end
