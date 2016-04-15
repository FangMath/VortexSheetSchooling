function pot=BFMM2(targetb, zetaf, gammaf, ib)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% An alternative of old BFMM.m file. 
% Evaluate the potential induced by

global wing s gamma1

tangential=wing(ib).tangential;
zetaf_ = zetaf(1:end-1);
gammaf_ = gammaf(1:end-1);

nsource = length(zetaf_);
source(1,:)=real(zetaf_);
source(2,:)=imag(zetaf_);

ntarget = length(targetb);
target = zeros(2,ntarget);

target(1,:)=real(targetb);
target(2,:)=imag(targetb);

N=nsource;

if N==1
    X=0;
else
    X(1)=gammaf_(2)-gammaf_(1);
    if N>2
        X(2:N-1)=gammaf_(3:N)-gammaf_(1:N-2);
    end
    X(N)=gammaf_(N)-gammaf_(N-1);
end
X=-X/(4*1i*pi);
dipstr=X;


    [J,K] = ndgrid(1:ntarget,1:nsource);
    dx = bsxfun(@minus,target(1,1:ntarget)',source(1,1:nsource));
    dy = bsxfun(@minus,target(2,1:ntarget)',source(2,1:nsource));
    F = 1./(dx+dy*1i);
    pot = F*dipstr.';
    pot = pot.';

% %% -add the integral over the last interval, which contains a log singularity.
distance=abs(zetaf(end-1)-zetaf(end));
gamma1(ib)=(gammaf(end-1)-gammaf(end))/distance; % approx \gamma(1)

% if j==1 (zetaj is the trailing edge), use adjacent node s(2) instead
pot(2:end)=pot(2:end)+gamma1(ib)*log((s(2:end)-1)./(s(2:end)-1-distance))./(2*pi*1i*tangential);
pot(1)=pot(1)+gamma1(ib)*log((s(2)-1)/(s(2)-1-distance))/(2*pi*1i*tangential);
end

