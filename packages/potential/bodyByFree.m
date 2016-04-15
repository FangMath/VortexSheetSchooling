function pot=bodyByFree(targetb, zetaf, gammaf)

nsource = length(zetaf);
source(1,:)=real(zetaf);
source(2,:)=imag(zetaf);

ntarget = length(targetb);
target = zeros(2,ntarget);

target(1,:)=real(targetb);
target(2,:)=imag(targetb);

N=nsource;

distance = abs(zetaf(1:N-1)-zetaf(2:N));
if N==1
    X=0;
else
    X(1:N-1)=(gammaf(2:N)-gammaf(1:N-1))./(zetaf(2:N)-zetaf(1:N-1));
end
X=-X/(2*1i*pi);
dipstr=X;

if N > 1

    [J,K] = ndgrid(1:ntarget,1:nsource);
    dx = bsxfun(@minus,target(1,1:ntarget).',source(1,1:nsource));
    dy = bsxfun(@minus,target(2,1:ntarget).',source(2,1:nsource));
    dz = dx+dy*1i;
    F = log(dz(:,1:N-1)./dz(:,2:N));
%    F = 1./(dx+dy*1i);
    pot = F*dipstr.';
    pot = pot.';
else
    pot = zeros(1,nsource);
end

end
