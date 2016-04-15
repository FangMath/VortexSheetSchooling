function pot = hilbert_sig(val, w, qflag)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   hilbert_sig.m
%   This function is to evaluate the hilbert transform of a function 
%   f(x) where x is a real number. 
%   
%       I(f,z) = \int_{-1}^1 \frac{f(x)}/{z-x} dx
%   
%   The function f(x) is intepolated by N+1 chebyshev nodes cos(k\pi/N), 
%   singluarity is pulled out artificially. The regular polynomial 
%   g(x)=\frac{f(x)-f(z)}/{x-z} is evaluated by gauss quadrature (qflag=1) 
%   or clenshaw_curtis quadrature (qflag=2) 
%   
%   val: function values at chebyshev nodes
%   x + 1i*y: target points
%   pot: I(f,z), potential at target points z
%   
%   Written by Fang Fang - 03/11/2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ns = length(val) - 1; % # of source points
nt = length(w); % # of target points

[xs,whts]=chebpts(ns+1); % build cheby nodes and weights
funw = lgrang(xs',val, w); % f_n(w) polynomial approx at target points

%% by gauss quadrature
if qflag == 1

[gxs,gwhts]=lgwt(ns,-1,1); % gauss points and weights
gval = lgrang(xs', val, gxs'); % interpolation from cheby points to gauss points
for k = 1:nt
    pot_sig = funw(k)*log((1+w(k))./(w(k)-1));
    % use gaussian quadrature of g(z)=(f(z)-f(w))/(z-w) 
    g_xs = (gval-funw(k))./(gxs'-w(k));
    pot_gau = sum(g_xs*gwhts);
    % add singular part and quadrature part
    pot(k) = pot_sig - pot_gau;
end

%% by clenshaw_curtis quadrature
else if qflag == 2

for k = 1:nt
    pot_sig = funw(k)*log((1+w(k))./(w(k)-1));
    % use clenshaw_curtis quadrature of g(z)=(f(z)-f(w))/(z-w) 
    g_xs = (val-funw(k))./(xs'-w(k));
    pot_gau = sum(whts.*g_xs);
    % add singular part and quadrature part
    pot(k) = pot_sig - pot_gau;
end

end
end
