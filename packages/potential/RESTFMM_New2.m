function pot=RESTFMM_New2(targetp, zetaf, gammaf, zetabodyw, nuw, PointVX, PointVCirl,outputX)

nsource = length(zetaf)+length(zetabodyw)+length(PointVX);
ntarget = length(targetp);
if nsource == 0
    pot = zeros(1, ntarget);
    return;
end
source(1,:)=real([zetaf,PointVX,zetabodyw]);
source(2,:)=imag([zetaf,PointVX,zetabodyw]);

target = zeros(2,ntarget);

target(1,:)=real(targetp);
target(2,:)=imag(targetp);

M=nsource;
MN=length(zetaf)+length(PointVX);
N=length(zetaf);


if N==1
    X=0;
else if N == 0
    X = [];
else
    X(1)=gammaf(2)-gammaf(1);
    if N>2
        X(2:N-1)=gammaf(3:N)-gammaf(1:N-2);
    end
    X(N)=gammaf(N)-gammaf(N-1);
end
X=-X/(4*1i*pi);
X(N+1:MN)=PointVCirl/(2*1i*pi);
if M ~= MN
X(MN+1:M)=nuw/(-2*(M-MN-1)*1i);
X([MN+1,M])=0.5*X([MN+1,M]);
end
dipstr=X;


    [J,K] = ndgrid(1:ntarget,1:nsource);
    dx = bsxfun(@minus,target(1,1:ntarget)',source(1,1:nsource));
    dy = bsxfun(@minus,target(2,1:ntarget)',source(2,1:nsource));
    F = 1./(dx+dy*1i);
    pot = F*dipstr.';
    pot = pot.';
end
