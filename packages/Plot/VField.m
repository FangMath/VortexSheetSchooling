function [pot, target] = VField(ib, tk)
global tem_at wing Sys Vp Para left right down up meshref 

        iib=Para.Nw-ib+1;

for Ib=1:Para.Nw
    K(Ib)=wing(Ib).K;
    zetaf(Ib,1:K(Ib))=wing(Ib).zetaf;
    gammaf(Ib,1:K(Ib))=wing(Ib).gammaf;
    nu(Ib,:)=wing(Ib).nu;
    zetabody(Ib,:)=wing(Ib).zetabody;
end

[x,y] = meshgrid(-(0:.2:20), (-1:.2:1));
%z = bsxfun(@plus,1i*(-5:.01:5)',(0:1:100));
z = x + 1i*y;
z = reshape(z,1, []);
size(z)

target(1,:) = real(z);
target(2,:) = imag(z);
%target(1,:) = real(zetaf(ib, 1:K(ib)));
%target(2,:) = imag(zetaf(ib, 1:K(ib)));

%%%%% calculate smoothing parameter delta_b and delta_f
        byOwnWing = freeByBodyHil(z, zetabody(ib,:), nu(ib,:), Para);
        byOwnSheet =  bodyByFree(z, zetaf(ib,1:K(ib)), gammaf(ib,1:K(ib)));
        byPoints = RESTFMM_New2(z, [], [], [], [], [wing(ib).PointVX],[wing(ib).PointVCirl]);

        %byPoints = RESTFMM_New2(z, [], [], [], [], [wing(iib).PointVX,wing(ib).PointVX],[wing(iib).PointVCirl,wing(ib).PointVCirl]);

        %byOtherSheet =  bodyByFree(z, zetaf(iib,1:K(iib)), gammaf(iib,1:K(iib)));

        %byOtherWing = freeByBodyHil(z, zetabody(iib,:), nu(iib,:), Para);

    %VelocityF{ib}=conj(byPoints + byOtherSheet + byOwnSheet)+ byOtherWing + byOwnWing;

    %VelocityF{ib}=conj(byPoints + byOwnSheet)+ byOwnWing;
    VelocityF{ib}=conj(byPoints + byOwnSheet)+ byOwnWing - Para.FishV;

    %VelocityF{ib}=conj(Bf.'+ Rt1 + byOtherSheets)+Para.FreeVelocity(tk) + byOtherWings;
    pot = VelocityF{ib};
end
