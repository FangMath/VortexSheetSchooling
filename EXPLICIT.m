function EXPLICIT(Para)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   EXPLICIT.m
%
%   This is the main function of the project FishSchool.
%   Free sheet advection by explicit method
%   Advection velocity induced by both bounded sheet and free sheet
%
%   # Function BFreeDelta.m is called
%   - output: VelocityF(tk,j);
%           zetaf_new (new free sheet, except for trailing edge)
%
%   Written by Fang Fang - 03/19/2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global wing tk zetaf gammaf K nu zetabody
global delta_f Zf Gf takeindex fang Zf_new Gf_new VelocityF VelocityF_new

FreeVelocity = Para.FreeVelocity(tk); dt = Para.dt;

% compute the velocity (Velocity(tk,j)) on free sheet at time tk
for ib = 1:Para.Nw
    K(ib) = wing(ib).K;
    zetaf(ib,1:K(ib)) = wing(ib).zetaf;
    gammaf(ib,1:K(ib)) = wing(ib).gammaf;
    nu(ib,:) = wing(ib).nu;
    zetabody(ib,:) = wing(ib).zetabody;
end

temp_wing_zetaf = [];

for ib = 1:Para.Nw
    %-------------------------------------------------------------------------
    % calculate smoothing parameter delta_b and delta_f
    %-------------------------------------------------------------------------
    % compute smoothing parameter \delta_b on bounded sheet
    delta_b = Para.Delta1*exp(-(abs(Para.s-1)/Para.Delta1).^2);
    
    % arclength of free sheet from the trailing edge (length K vector)
    distance = zeros(1,K(ib));
    distance(K(ib)) = 0;
    for k = K(ib)-1:-1:1
        distance(k) = abs(zetaf(ib,k)-zetaf(ib,k+1)) + distance(k+1);
    end
    delta_f{ib} = Para.Delta1+(Para.Delta0-Para.Delta1)*((distance/Para.eta).^Para.p./(1+(distance./Para.eta).^Para.p));
    
    %-------------------------------------------------------------------------
    % compute the velocity (Velocity(tk,j)) on free sheet at time tk
    %-------------------------------------------------------------------------
    % vel conj of f-sht by f-sht itself and it's wing, by Ken's FMM
    Bf = [];
    Bf = FMMken(zetaf(ib,1:K(ib)), gammaf(ib,1:K(ib)), delta_f{ib}, zetabody(ib,:), nu(ib,:), delta_b, 32, 64, 1e-12, 0);
    
    tocBf = toc
    if Para.Nw==2
        iib = Para.Nw-ib+1;
        % Calculate the velocity conj on sheet induced by wings and the other sheet
        tic
        % vel conj of f-sht by other wings and point vortices, by 1/z summation
        Rt1 = [];
        Rt1 = RESTFMM_New2(zetaf(ib,1:K(ib)), [], [], [], [], [wing(iib).PointVX,wing(ib).PointVX],[wing(iib).PointVCirl,wing(ib).PointVCirl]);
        toc
        
        tic
        % vel conj of f-sht by other f-shts, by log(z) summation
        byOtherSheets = [];
        byOtherSheets =  bodyByFree(zetaf(ib,1:K(ib)), zetaf(iib,1:K(iib)), gammaf(iib,1:K(iib)));
        toc
        tic
        % vel conj of f-sht by other free wings, by hlbt-trans(de-sing)
        byOtherWings = [];
        byOtherWings = freeByBodyHil(zetaf(ib,1:K(ib)), zetabody(iib,:), nu(iib,:), Para);
        toc
    end
    
    VelocityF{ib} = conj(Bf.'+ Rt1 + byOtherSheets)+FreeVelocity + byOtherWings;
    
    %-------------------------------------------------------------------------
    %  calculate velocity of  Point Vortices
    %-------------------------------------------------------------------------
    if ~isempty(wing(ib).PointVX)
        iib = Para.Nw-ib+1;
        
        % vel conj of self f-sht and self point vortices except selfpoint, by 1/z summation
        Rp1 = [];
        for pk = 1:length(wing(ib).PointVX)
            Rp1(pk) = RESTFMM_New2(wing(ib).PointVX(pk), zetaf(ib,1:K(ib)), gammaf(ib,1:K(ib)),...
                zetabody(ib,:), nu(ib,:), wing(ib).PointVX(1:end~=pk), wing(ib).PointVCirl(1:end~=pk));
        end
        
        % vel conj of other f-sht and other point vortices except self, by 1/z summation
        Rp2 = [];
        Rp2 = RESTFMM_New2(wing(ib).PointVX, zetaf(3-ib,1:K(3-ib)), gammaf(3-ib,1:K(3-ib)),...
            zetabody(3-ib,:), nu(3-ib,:), wing(3-ib).PointVX, wing(3-ib).PointVCirl);
        
        VelocityPF{ib} = conj(Rp1+Rp2)+FreeVelocity;
        
        %px = wing(ib).PointVX;
        index_piece(wing,ib,1);hold on;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Update point vortices
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for ib = 1:Para.Nw
    if ~isempty(wing(ib).PointVX)
        temp_wing_PointVX{ib} = zeros(1,length(wing(ib).PointVX));
        if wing(ib).StartNew
            temp_wing_PointVX{ib}(end-wing(ib).StartNew+1:end) = ...
                wing(ib).PointVX(end-wing(ib).StartNew+1:end)+dt*VelocityPF{ib}(end-wing(ib).StartNew+1:end);
        end
        
        if length(wing(ib).PointVX)-wing(ib).StartNew>0
            temp_wing_PointVX{ib}(1:end-wing(ib).StartNew) = wing(ib).PointVX(1:end-wing(ib).StartNew) + dt*(1.5*VelocityPF{ib}(1:end-wing(ib).StartNew)-.5*wing(ib).PointVVp(1:end));
        end
        wing(ib).PointVX = temp_wing_PointVX{ib};
        wing(ib).PointVVp = VelocityPF{ib};
        
    end
    wing(ib).StartNew = 0;
    
    
    %-------------------------------------------------------------------------
    % Update the position of free sheets
    % newly born node advected by forward Euler (j = K)
    % Other nodes advected by Adam-Bashforth (j = 1,2,...,K-1)
    %-------------------------------------------------------------------------
    if K(ib)>1
        temp_wing_zetaf(ib,1:K(ib)-1) = zetaf(ib,1:K(ib)-1) + dt*(1.5*VelocityF{ib}(1:K(ib)-1)-.5*wing(ib).Vp(1:K(ib)-1));
    end
    temp_wing_zetaf(ib,K(ib)) = zetaf(ib,K(ib))+dt*VelocityF{ib}(K(ib));
    
    %-------------------------------------------------------------------------
    % Dynamically add/delete points on free sheets
    %-------------------------------------------------------------------------
    Zf = temp_wing_zetaf(ib,1:K(ib));
    takeindex = logical(ones(1,length(Zf)));
    fang = 0;
    DelVortex(ib, Para.perc);
    Zf = Zf(takeindex);
    Gf = gammaf(ib,takeindex);
    VelocityF{ib} = VelocityF{ib}(takeindex);
    
    Zf_new = Zf; Gf_new = Gf; VelocityF_new = VelocityF{ib};
    AddVortex(ib, dt, Para.Delta0);
    wing(ib).zetaf = Zf_new;
    wing(ib).gammaf = Gf_new;
    wing(ib).K = length(wing(ib).zetaf);
    wing(ib).Vp  =  VelocityF_new;
end

%-------------------------------------------------------------------------
% Turn far sheet pieces into point vortices
%-------------------------------------------------------------------------
if mod(tk,50)==0 && tk >= 100
    %     Para.farR = 5;
    Para.farR(2) = -wing(2).dotcx+2; Para.farR(1) = Para.farR(2)+wing(2).cx-wing(1).cx;
    farR = Para.farR
    for ib = 1:Para.Nw
        TurnToPointVortices(ib,Para.farR(ib));
    end
end

end
% global
% vector initialize, change size
% clean up sub routines
