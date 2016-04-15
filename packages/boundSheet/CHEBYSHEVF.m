function [chebyshevF]=CHEBYSHEVF(ib)
% Calculate the cos expansion of bounded part of f, i.e. the RHS of integral equation f(j)

% # called in f.m
% # Function B.m is called
global Para wing tk gamma1 Rt bt fang
global Nw M zetaf gammaf K nu zetabody
M=Para.M; s=Para.s; FreeVelocity=Para.FreeVelocity(tk-1);
%% calculate the RHS of integral equation f(j)
% PART1=zeros(M+1,1); PART2=zeros(M+1,1); f=zeros(1,M+1);

if Nw==2
    iib=Nw-ib+1;
    Rt(ib,1:M+1)=RESTFMM_New2(wing(ib).zetabody, [], [], zetabody(iib,:), nu(iib,:), [wing(iib).PointVX,wing(ib).PointVX],[wing(iib).PointVCirl,wing(ib).PointVCirl]) +...
        bodyByFree(wing(ib).zetabody, zetaf(iib,1:K(iib)), gammaf(iib,1:K(iib))); % induced by other wing, other sheet, and all point vortices
end
% tic
bt(ib,1:M+1)=BFMM2(wing(ib).zetabody, zetaf(ib,1:K(ib)), gammaf(ib, 1:K(ib)), ib); % induced by self sheet
% cbt=toc
PART11=real(wing(ib).normal*conj(wing(ib).dzetabody-FreeVelocity));

f=PART11-real(wing(ib).normal*(bt(ib,:)+Rt(ib,:)));

%% calculate \tilde{f}, the bounded part of f
f_tilde(1)=f(1)+gamma1(ib)*log(1-s(2))/(2*pi);
f_tilde(2:M+1)=f(2:M+1)+gamma1(ib)*log(1-s(2:M+1))/(2*pi);
fang(ib,:)=wing(ib).dzetabody(1);
real(wing(ib).normal*conj(wing(ib).dzetabody(1)-FreeVelocity));
%% use fft to get the cos expansion
f_tilde(M+2:2*M)=f_tilde(M:-1:2);
chebyshevF=fft(f_tilde);
chebyshevF=real(chebyshevF(1:M+1));

chebyshevF=chebyshevF/M;
chebyshevF([1,M+1])=.5*chebyshevF([1,M+1]);

end
