function p = lgrang( x, f, xt )
n = length(x); m = length(xt);
% Evaluates Lagrange interpolating poly
% for the data vectors x and f,
% at the points xt(1), ..., xt(m),
% returning results in p(1), ..., p(m)
for i = 1 : n % determine weights
w(i) = prod( x(i) - x([1:i-1 i+1:n]) );
end
for l = 1 : m % determine poly values at xt’s
wt = w.*( xt(l) - x );
p(l) = sum(f./wt)/sum(1./wt);
end
