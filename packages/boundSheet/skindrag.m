function skin = skindrag(twing,mu, tau)
% evaluate the skin friction drag based on mean fluid velocity on sides of wing.
% note that skin friction coefficient has not been multiplied

M = length(mu)-1; % M = Para.M

%%%% Integrate tangential average velocity mu %%%%%%
[xs,whts]=chebpts(M+1); % build cheby nodes and weights
AvgVel = sum(whts.*(mu-tau)); % mu-tau is the two-side average of fluid-body relative tangential velocity difference
AvgVel = 0.5*AvgVel;

%%%% Integrate tangential difference velocity nu/sqrt(1-s^2) %%%%%%
DiffVel = sum(twing.nu)*pi/M;
DiffVel = 0.25*DiffVel;

%%%% Calculate mean velocity on up side and bottom side %%%%
VelUp = AvgVel + DiffVel;
VelBtm = AvgVel - DiffVel;

%%%% Calculate skin drag on up side and bottom side %%%%
SkinUp = power(abs(VelUp),1.5);
SkinBtm = power(abs(VelBtm),1.5);

%%%% add up up/bottom drags %%%%
skin = SkinUp + SkinBtm;
end
