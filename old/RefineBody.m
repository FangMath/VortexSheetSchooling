function  RefineBody( str_input,M,tem_at_input )
global Para
load([str_input,'data.mat'],'Wing','Sys','tem_at','Vp');
load([str_input,'parameter.mat'],'Para');
OldM=Para.M
Olds=Para.s;
Para.M=M;
Para.s=cos((0:Para.M)*pi/Para.M);
s=Para.s;
% load([str_input,'data.mat'],'Wing');
if tem_at_input<tem_at
    tem_at=tem_at_input;
end
for ib=1:Para.Nw
    for k=1:tem_at-3
        NewWing(ib).ZetaBody(k,1:OldM+1)=Wing(ib).ZetaBody(k,1:OldM+1);
        NewWing(ib).ZetaBody(k,OldM+2:M+1)=Wing(ib).ZetaBody(k,1);
        
        NewWing(ib).Nu(k,1:OldM+1)=Wing(ib).Nu(k,1:OldM+1);
        NewWing(ib).Nu(k,OldM+2:M+1)=0;
    end
    size(Olds)
    size(Wing(ib).ZetaBody(k,1:OldM+1))
    for k=tem_at-2:tem_at
        NewWing(ib).ZetaBody(k,1:M+1)=interp1(Olds,Wing(ib).ZetaBody(k,1:OldM+1),s,'pchip');
        NewWing(ib).Nu(k,1:M+1)=interp1(Olds,Wing(ib).Nu(k,1:OldM+1),s,'pchip');
    end
    Wing(ib).ZetaBody=NewWing(ib).ZetaBody;
    Wing(ib).Nu=NewWing(ib).Nu;
end
% Para.M=80;Para.s=cos((0:Para.M)*pi/Para.M);
str=[str_input(1:end-1),'_refine/'];
if exist(str,'dir')
    rmdir(str,'s');
end
mkdir(str);
save( [str,'parameter.mat'], 'Para');
if tem_at>700
    save( [str,'data.mat'],'tem_at','Sys','Wing','Vp','-v7.3');
else
    save( [str,'data.mat'],'tem_at','Sys','Wing','Vp');
end
end
