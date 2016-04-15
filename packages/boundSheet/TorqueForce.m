function [torque,lift,thrust]=TorqueForce(ib)
global Para wing tk Rt bt wing_p1 wing_p2 wing_p3

M=Para.M; s=Para.s; FreeVelocity=Para.FreeVelocity(tk-1); dt=Para.dt;

%% compute tau(j), mu(j), j=1,...M+1 and d\nu/dt, d\Gamma/dt
mu=zeros(1,M+1);    dotnu=zeros(1,M+1);
tau=real(wing(ib).tangential*conj(wing(ib).dzetabody));

for j=1:M+1
    mu(j)=real(wing(ib).tangential*(bt(ib,j)+Rt(ib,j)+conj(FreeVelocity)));
    
    % 2nd order BDF (1st order when tk==2)
    if (tk==2) %| tk-input_tk<2
        dotnu(j)=DERIVATIVE(wing(ib).nu(j),wing_p1(ib).nu(j),0,0,dt,111);
    else
        dotnu(j)=DERIVATIVE(wing(ib).nu(j),wing_p1(ib).nu(j),wing_p2(ib).nu(j),0,dt,10);
    end
end

if (tk==2) %| tk-input_tk<2
    dotgammaMax=DERIVATIVE(wing(ib).gammaMax,wing_p1(ib).gammaMax,0,0,dt,111);
else
    dotgammaMax=DERIVATIVE(wing(ib).gammaMax,wing_p1(ib).gammaMax,wing_p2(ib).gammaMax,0,dt,10);
end

%% sum to get the integral of pressure jump over chord
integratep=sum((mu-tau).*wing(ib).nu-(1+s).*dotnu)-...
    .5*((mu(1)-tau(1))*wing(ib).nu(1)-(1+s(1))*dotnu(1))-...
    .5*((mu(M+1)-tau(M+1))*wing(ib).nu(M+1)-(1+s(M+1))*dotnu(M+1));
integratepS=sum((mu-tau).*wing(ib).nu.*s+(1-s.^2).*dotnu/2)-...
    .5*((mu(1)-tau(1))*wing(ib).nu(1)*s(1)+(1-s(1)^2)*dotnu(1)/2)-...
    .5*((mu(M+1)-tau(M+1))*wing(ib).nu(M+1)*s(M+1)+(1-s(M+1)^2)*dotnu(M+1)/2);
integratep=-(pi*integratep/M+2*dotgammaMax);
integratepS=-pi*integratepS/M;


%% calculate net force
wing(ib).integratepS=integratepS;
wing(ib).integrateP=integratep;
wing(ib).les=-pi*wing(ib).nu(M+1)^2/8;
wing(ib).skin = .5*Para.skin*skindrag(wing(ib),mu,tau); % .5*Para.skin is the real coefficient...
%wing(ib).skin = Para.skin*power(abs(wing(ib).dotcx),1.5);

thrust=wing(ib).integrateP*real(wing(ib).normal)+wing(ib).les*real(wing(ib).tangential);
skin_thrust = wing(ib).skin*real(wing(ib).tangential);

lift=wing(ib).integrateP*imag(wing(ib).normal)+wing(ib).les*imag(wing(ib).tangential);
skin_lift = wing(ib).skin*imag(wing(ib).tangential);

%% calculate torque
c=wing(ib).cx+wing(ib).cy*1i;
PART1=(wing(ib).L+wing(ib).tangential-c)*wing(ib).integrateP;
PART2=wing(ib).tangential*integratepS;
torque=real(conj(wing(ib).tangential)*(PART1+PART2))+wing(ib).les*real((wing(ib).zetabody(M+1)-c)*conj(wing(ib).normal));

skin_torque = wing(ib).skin*real((wing(ib).zetabody(M+1)-c)*conj(wing(ib).normal));

if ib == 2
    Para.exforce(tk) = 0; % add external force
    thrust = thrust + Para.exforce(tk);
end
wing(ib).thrust=thrust + skin_thrust;
wing(ib).lift=lift + skin_lift;
wing(ib).torque=torque + skin_torque;

end
