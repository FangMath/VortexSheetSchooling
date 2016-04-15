function [pot] = VField2(ib, tk, target)
global tem_at wing Sys Vp Para left right down up meshref 

        iib=Para.Nw-ib+1;

for Ib=1:Para.Nw
    K(Ib)=wing(Ib).K;
    zetaf(Ib,1:K(Ib))=wing(Ib).zetaf;
    gammaf(Ib,1:K(Ib))=wing(Ib).gammaf;
    nu(Ib,:)=wing(Ib).nu;
    zetabody(Ib,:)=wing(Ib).zetabody;
end

%%%%% calculate smoothing parameter delta_b and delta_f

    % vel conj of f-sht by other wings and point vortices, by 1/z summation
        byPoints = RESTFMM_New2(target, [], [], [], [], [wing(iib).PointVX,wing(ib).PointVX],[wing(iib).PointVCirl,wing(ib).PointVCirl]);

    % vel conj of f-sht by other f-shts, by log(z) summation
        bySheets1 =  bodyByFree(target, zetaf(ib,1:K(ib)), gammaf(ib,1:K(ib)));
        bySheets2 =  bodyByFree(target, zetaf(iib,1:K(iib)), gammaf(iib,1:K(iib)));

    % vel conj of f-sht by other free wings, by hlbt-trans(de-sing)
        byWings1 = freeByBodyHil(target, zetabody(ib,:), nu(ib,:), Para);
        byWings2 = freeByBodyHil(target, zetabody(iib,:), nu(iib,:), Para);

    VelocityF{ib}=conj(byPoints + bySheets1 + bySheets2) + byWings1 + byWings2;
    pot = VelocityF{ib};
end
