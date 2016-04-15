clear all; close all;
ns = 80; nt = 80;
x = []; y = [];

%%%%%%% set target points %%%%%%%%
st = 1;
for k = 0+st:8+st
    yy(k)=1/10^k;
r = yy(k);
y = [y, ones(1,nt)*r];
x = [x, -1:2/(nt-1):1];

y = [y, r*sin([pi/2 : -pi/nt: -pi/2])];
x = [x, r*cos([pi/2 : -pi/nt: -pi/2])+1];
w = x + 1i*y;

%y = [y, r*sin([pi/2, -pi/2])];
%x = [x, r*cos([pi/2, -pi/2])+1];

end
yy1 = yy(1)
yyn = yy(end)
nt = length(w)

%% set test polynomial coefficients
a = [ 2, 3, 0, 1];

%% set chebyshev nodes and function values at chebyshev nodes
[xs,whts]=chebpts(ns+1);
val = exp(xs');
%val = exp(xs'.^2);
%val = cos(xs');
%val = sin(xs');
val = polyval(a, xs');

%%%%% by Michael's scheme %%%%%%
%[gxs,gwhts]=legewhts(ns,1);
[gxs2,gwhts2]=lgwt(ns,-1,1); gxs = fliplr(gxs2');
tic;
fxs = lgrang(xs', val, gxs);
[pot_r,pot_i]=hilbert(ns,fxs,zeros(1,ns),nt,x,y);
pr_mk = pot_r; pi_mk = pot_i;
toc

%quiver(x,y,pr_mk,pi_mk,'r'); hold on;

%%%%% by exact Qn %%%%%%
for k = 1:nt
    pot_exact(k) = 2*hilbert_poly3(a,x(k)+1i*y(k));
end
pr_exact = real(pot_exact); pi_exact = imag(pot_exact);

%%%%% by Fang's scheme %%%%%%
tic
pot_gq = hilbert_sig(val, w, 1);
pr_gq = real(pot_gq); pi_gq = imag(pot_gq);
toc

%% by clenshaw_curtis quadrature
tic
pot_ccq = hilbert_sig(val, w, 2);
pr_ccq = real(pot_ccq); pi_ccq = imag(pot_ccq);
toc


%%% error by michael's scheme
fprintf('Error of gauss and mk is:');
er = pr_mk - pr_gq;
ei = pi_mk - pi_gq;
err = norm(er)/sqrt(nt) + 1i* norm(ei)/sqrt(nt)

%% error by cc quadrature
fprintf('Error of gauss and cc is:');
er = pr_ccq - pr_gq;
ei = pi_ccq - pi_gq;
err = norm(er)/sqrt(nt) + 1i* norm(ei)/sqrt(nt)

%% error by michael's scheme
fprintf('Error of mk and exact is:');
er = pr_mk - pr_exact;
ei = pi_mk - pi_exact;
err = norm(er)/sqrt(nt) + 1i* norm(ei)/sqrt(nt)

%% error by gauss quadrature
fprintf('Error of gauss and exact is:');
er = pr_gq - pr_exact;
ei = pi_gq - pi_exact;
err = norm(er)/sqrt(nt) + 1i* norm(ei)/sqrt(nt)

%% error by cc quadrature
fprintf('Error of and exact is:');
er = pr_ccq - pr_exact;
ei = pi_ccq - pi_exact;
err = norm(er)/sqrt(nt) + 1i* norm(ei)/sqrt(nt)
