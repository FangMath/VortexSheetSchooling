function wing = FLATPOSITION(wing, Para, tk)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Find the motion of fish
%   Written by Fang Fang - 03/19/2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% define the instantaneous angle and period
t=(tk-1)*Para.dt;
period=2*pi;

for ib = 1:Para.Nw;
    % angle of fish, position of center mass -> the leading edge
    wing(ib).theta=Para.PAmp*cos(period*(t+Para.phase(ib)))+wing(ib).omega;
    wing(ib).L=Para.HAmp(ib)*cos(period*t)*1i+wing(ib).cx+wing(ib).cy*1i;
    
    % find the position of plate
    wing(ib).zetabody=wing(ib).L+(Para.s+1)*exp(wing(ib).theta*1i);
    
    % velocities of body
    wing(ib).dtheta=-Para.PAmp*period*sin(period*t)+wing(ib).domega;
    wing(ib).dotL = -Para.HAmp(ib)*period*sin(period*t)*1i+wing(ib).dotcx+wing(ib).dotcy*1i;
    wing(ib).dzetabody=wing(ib).dotL+(Para.s+1)*1i*wing(ib).dtheta*exp(wing(ib).theta*1i);
    
    % find the tangential and normal vector of plate
    wing(ib).tangential=exp(1i*wing(ib).theta);     % tangential vector
    wing(ib).normal=1i*exp(1i*wing(ib).theta);      % normal vector
end

