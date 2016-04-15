function DivideData(str_input)
% global tem_at Wing Sys
if nargin~=0
    str=str_input;
else
    str='./output/';
end

load( [str,'data.mat'], 'tem_at','Wing','Sys');
str=[str,'Ddata/'];
if exist(str,'dir')
    rmdir(str,'s');
end
mkdir(str);
save([str,'tem.mat'],'tem_at');
for tk=1:tem_at
    tk
    strs=[str,'T',num2str(tk)];
for ib=1:2
    wing(ib).nu=Wing(ib).Nu(tk,:);
    wing(ib).zetabody=Wing(ib).ZetaBody(tk,:);
    wing(ib).Lenf=Wing(ib).Lenf(tk);
    wing(ib).zetaf=Wing(ib).ZetaF{tk}(1:wing(ib).Lenf);
    wing(ib).gammaf=Wing(ib).GammaF{tk}(1:wing(ib).Lenf);
    
    wing(ib).gammaMax=Wing(ib).GammaMax(tk);

end

sys.CenterX=Sys.CenterX(tk);
sys.CenterY=Sys.CenterY(tk);
sys.Omega=Sys.Omega(tk);

sys.dCenterX=Sys.dCenterX(tk);
sys.dCenterY=Sys.dCenterY(tk);
sys.dOmega=Sys.dOmega(tk);
save([strs,'.mat'],'tem_at','sys','wing');

gammaMaxLeft(tk)=-Wing(2).GammaMax(tk);
gammaMaxRight(tk)=-Wing(1).GammaMax(tk);
end
save([str_input,'GMAX.mat'],'gammaMaxLeft','gammaMaxRight');

end