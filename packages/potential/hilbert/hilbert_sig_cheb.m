function pot = hilbert_sig_cheb(val, w)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   hilbert_sig_cheb.m
%   This function is to evaluate the hilbert transform of a function with 
%   inverse square root singularity 
%   Here the function f(x) is a real function.
%   
%       I(f,w) = \int_{-1}^1 \frac{f(x)}/{\sqrt{1-x^2}(w-x)} dx
%   
%   The function f(x) is intepolated by N+1 chebyshev nodes cos(k\pi/N), 
%   singluarity is pulled out artificially. The regular polynomial 
%   g(x)=\frac{f(x)-f(w)}/{x-w} is evaluated by chebyshev quadrature 
%   
%   val: function values at chebyshev nodes
%   w : target points
%   pot: I(f,w), potential at target points w
%   
%   Written by Fang Fang - 03/17/2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ns = length(val) - 1; % # of source points
nt = length(w); % # of target points

[xs,whts]=chebpts(ns+1); % build cheby nodes and weights
funw = lgrang(xs',val, w); % f_n(w) polynomial approx at target points

for k = 1:nt
    pot_sig = funw(k)*pi./(sqrt(w(k)+1)*sqrt(w(k)-1));
    %pot_sig = funw(k)*pi*log((1+w(k))./(w(k)-1));

    % use clenshaw_curtis quadrature of g(z)=(f(z)-f(w))/(z-w) 
    g_xs = (val-funw(k))./(xs'-w(k));
    g_xs(1) = g_xs(1)/2;
    g_xs(ns+1) = g_xs(ns+1)/2;

    pot_gau = sum(g_xs)*pi/ns;
    % add singular part and quadrature part
    pot(k) = pot_sig - pot_gau;
end

end
