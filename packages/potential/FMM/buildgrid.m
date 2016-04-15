function AX = buildgrid(zeta, X, delta, N, radius)
% global delta
% size(X)

x=[real(zeta);imag(zeta)];
[d,n] = size(x);
xmin = min(x,[],2)';
xmax = max(x,[],2)';
l = xmax - xmin;
ngrid=ceil(l/radius)+3;

% tic
xindex_2d = floor(bsxfun(@minus, x, xmin')/radius)+1;
% toc
xindex = [1,ngrid(1)]*xindex_2d + 1; %xindex(i) is the index of box in which xi locates
m = histc(xindex,1:ngrid(1)*ngrid(2));
% toc

% tic
for in=1:length(m)
    a{in}=find(bsxfun(@eq, xindex,in));
end
% toc

AK=zeros(1,N);
% TT=zeros(1,N);
for xk=1:N
    xin=xindex(xk);
    T=[a{xin-1},a{xin},a{xin+1},...
        a{xin+ngrid(1)-1},a{xin+ngrid(1)},a{xin+ngrid(1)+1},...
        a{xin-ngrid(1)-1},a{xin-ngrid(1)},a{xin-ngrid(1)+1}];
    T(logical(T==xk))=[]; 
%     TT(xk,:)=T;
    G=conj(zeta(xk)-zeta(T))./(abs(zeta(xk)-zeta(T)).^2+delta(T).^2)-1./(zeta(xk)-zeta(T));
%     G=-1./(zeta(xk)-zeta(T));
%     xk
%     T
    AX(xk)=G*X(T);
end
% ccc2=conj(zeta(1)-zeta(TT(1,:)))./(abs(zeta(1)-zeta(TT(1,:))).^2+delta(TT(1,:)).^2);
% ddd2=X(TT(1,:));
% ddd21=ddd2(1:3)
% ccc21=ccc2(1:3)
% TTT=X(2:4)
