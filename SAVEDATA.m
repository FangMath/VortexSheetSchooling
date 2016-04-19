function SAVEDATA(fileDir)
global tk wing_p1 wing_p2 wing_p3
global Wing wing Para
%% --record the solutions
for ib=1:Para.Nw
    Wing(ib).Cx(tk) = wing(ib).cx;
    Wing(ib).Cy(tk) = wing(ib).cy;
    Wing(ib).dotCx(tk) = wing(ib).dotcx;
    Wing(ib).dotCy(tk) = wing(ib).dotcy;
    Wing(ib).L(tk) = wing(ib).L;
    Wing(ib).dotL(tk) = wing(ib).dotL;

    Wing(ib).IntegrateP(tk)=wing(ib).integrateP;
    Wing(ib).Skin(tk)=wing(ib).skin;
    Wing(ib).Les(tk)=wing(ib).les;

    Wing(ib).Spring(tk)=wing(ib).spring;
    Wing(ib).Hydrolift(tk)=wing(ib).hydrolift;

    Wing(ib).Thrust(tk)=wing(ib).thrust;
    Wing(ib).Lift(tk)=wing(ib).lift;
    Wing(ib).Torque(tk)=wing(ib).torque;
end

gammaMaxLeft(tk)=-wing(2).gammaMax;
gammaMaxRight(tk)=-wing(1).gammaMax;

%record the current terminate time
tem_at=tk;
% write data into files
w = whos('Wing');
fprintf('mem: %6.2f (MB)\n', w.bytes/1e6)

if w.bytes/1e6 >2000
    fprintf('mem: %6.2f (MB)\n', w.bytes/1e6)
    save( [fileDir,'/data.mat'],'tem_at', 'Wing','-v7.3');
else
    save( [fileDir,'/data.mat'],'tem_at', 'Wing');
end

%%
save([fileDir,'/Ddata/tem.mat'],'tem_at');

fileDirs=[fileDir,'/Ddata/T',num2str(tk)];
for ib=1:Para.Nw
    wing(ib).Lenf = length(wing(ib).zetaf);
end
save([fileDirs,'.mat'], 'tk', 'wing');

%%
save([fileDir,'/GMAX.mat'],'gammaMaxLeft','gammaMaxRight');
save( [fileDir,'/parameter.mat'], 'Para');
end
