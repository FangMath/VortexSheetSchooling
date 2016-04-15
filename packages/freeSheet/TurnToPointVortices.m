function TurnToPointVortices(ib,criticalL)
global Cirl_vortex piece_n piece_p wing NumPiece
pause(1);fprintf('turn to point vortices');
ib

ploton=1;
index_piece(wing,ib,ploton);
PointXN{ib}=[]; PointCirlN{ib}=[];
PointXP{ib}=[]; PointCirlP{ib}=[];

ContinueTurn=1; kPiece=1; % # of dipoles to be turned

while ContinueTurn==1&&kPiece<=min(NumPiece(ib))
    kPiece
    
    SheetPieceN=wing(ib).zetaf(piece_n{ib}(kPiece,1):piece_n{ib}(kPiece,2));
    SheetCirlN=Cirl_vortex{ib}(piece_n{ib}(kPiece,1):piece_n{ib}(kPiece,2));
    
    PointXN{ib}(kPiece)=sum(SheetCirlN.*SheetPieceN)/sum(SheetCirlN);
    PointCirlN{ib}(kPiece)=sum(SheetCirlN);
    
    SheetPieceP=wing(ib).zetaf(piece_p{ib}(kPiece,1):piece_p{ib}(kPiece,2));
    SheetCirlP=Cirl_vortex{ib}(piece_p{ib}(kPiece,1):piece_p{ib}(kPiece,2));
    
    PointXP{ib}(kPiece)=sum(SheetCirlP.*SheetPieceP)/sum(SheetCirlP);
    PointCirlP{ib}(kPiece)=sum(SheetCirlP);

a = norm(min([PointXN{ib}(kPiece),PointXP{ib}(kPiece)]-wing(ib).cx))
criticalL
kPiece
min(NumPiece(ib))
    
    if norm(min([PointXN{ib}(kPiece),PointXP{ib}(kPiece)]-wing(ib).cx))>criticalL &&kPiece<=min(NumPiece(ib))-2
        kPiece=kPiece+1;
    else
        ContinueTurn=0;
        PointXN{ib}=PointXN{ib}(1:kPiece-1);
        PointXP{ib}=PointXP{ib}(1:kPiece-1);
        PointCirlN{ib}=PointCirlN{ib}(1:kPiece-1);
        PointCirlP{ib}=PointCirlP{ib}(1:kPiece-1);
    end
end

PointNum=kPiece-1 % # of dipoles turned

wing(ib).StartNew=PointNum*2; % # of new vortices (2*dipoles)
wingStartNew = wing(ib).StartNew
OldLength=length(wing(ib).PointVX)
if ib==1
    wing(ib).PointVX(OldLength+1:2:OldLength+wing(ib).StartNew)=PointXN{ib};
    wing(ib).PointVX(OldLength+2:2:OldLength+wing(ib).StartNew)=PointXP{ib};
    wing(ib).PointVCirl(OldLength+1:2:OldLength+wing(ib).StartNew)=PointCirlN{ib};
    wing(ib).PointVCirl(OldLength+2:2:OldLength+wing(ib).StartNew)=PointCirlP{ib};
else
    wing(ib).PointVX(OldLength+1:2:OldLength+wing(ib).StartNew)=PointXP{ib};
    wing(ib).PointVX(OldLength+2:2:OldLength+wing(ib).StartNew)=PointXN{ib};
    wing(ib).PointVCirl(OldLength+1:2:OldLength+wing(ib).StartNew)=PointCirlP{ib};
    wing(ib).PointVCirl(OldLength+2:2:OldLength+wing(ib).StartNew)=PointCirlN{ib};
end

wing(ib).zetaf=wing(ib).zetaf(min(piece_p{ib}(kPiece,1),piece_n{ib}(kPiece,1)):end);
wing(ib).gammaf=wing(ib).gammaf(min(piece_p{ib}(kPiece,1),piece_n{ib}(kPiece,1)):end);
wing(ib).K=length(wing(ib).zetaf);

wing(ib).Vp=wing(ib).Vp(min(piece_p{ib}(kPiece,1),piece_n{ib}(kPiece,1)):end);

index_piece(wing,ib,ploton);
end
