function [wing, Para, tk] = INITIALIZE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   INITIALIZE.m
%   start from initial time
%
%   Written by Fang Fang - 03/19/2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%-------------------------------------------------------------------------
Para = parameter;
fileDir = Para.dir

if exist(fileDir,'dir')
    rmdir(fileDir,'s');
end
mkdir(fileDir);


fileDird=[fileDir,'/Ddata/'];
if exist(fileDird,'dir')
    rmdir(fileDird,'s');
end
mkdir(fileDird);

save( [fileDir,'/parameter.mat'], 'Para');

%-------------------------------------------------------------------------
%% --Initialize data at tk=1, i.e. initial time t=0:
tk=1;
wing(1).cx = 0;
wing(1).cy = 0;
wing(1).omega = 0;

wing(2).cx = wing(1).cx + Para.dx;
wing(2).cy = wing(1).cy + Para.dy;
wing(2).omega = 0;

wing(1).dotcx = Para.FishV(1);
wing(1).dotcy = 0;
wing(1).domega = 0;
wing(2).dotcx = Para.FishV(2);
wing(2).dotcy = 0;
wing(2).domega = 0;

wing = FLATPOSITION(wing, Para, tk);
for ib=1:Para.Nw
    wing(ib).nu=zeros(1,Para.M+1);
    wing(ib).K=1;
    wing(ib).zetaf=wing(ib).zetabody(1);
    wing(ib).gammaf(1)=0;
    wing(ib).gammaMax(1)=0;
    wing(ib).PointVX=[];
    wing(ib).PointVCirl=[];
    
    wing(ib).thrust  = 0;
    wing(ib).torque = 0;
    wing(ib).lift = 0;
    
    wing(ib).Vp = 0;
    wing(ib).PointVVp = [];
    
    wing(ib).StartNew=0;
end

%-------------------------------------------------------------------------
% save data at tk=1

tem_at=tk;
save([fileDird,'/tem.mat'],'tem_at');

fileDirs=[fileDir,'/Ddata/T',num2str(tk)];
for ib=1:Para.Nw
    wing(ib).Lenf=wing(ib).K;
end

save([fileDirs,'.mat'],'tk','wing');
end
