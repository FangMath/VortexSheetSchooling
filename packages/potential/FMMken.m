function pot = FMMken(zetaf, gammaf, delta_f, zetabodyw, nuw, delta_b, p,occ, tol, store)
% velocity conj of free sheet induced by free sheet itself and the wing attached to this sheet

x = [real(zetaf), real(zetabodyw); imag(zetaf), imag(zetabodyw)];
M = size(x,2);
theta = (1:p)*2*pi/p;
proxy = [cos(theta); sin(theta)];
proxy_ = proxy;
for i = 1:3
    proxy_ = 2*proxy_;
    proxy = [proxy proxy_];
end
clear theta proxy_

delta=[delta_f,delta_b];
mean_delta = mean(delta);
F = ifmm(@Afun,x,occ,@Kfun,proxy,tol,store);
% w = whos('F');
% fprintf([repmat('-',1,80) '\n'])
% fprintf('mem: %6.2f (MB)\n', w.bytes/1e6)

X=zeros(M,1);
N=length(zetaf);
if N==1
    X(1)=0;
else
    X(1)=gammaf(2)-gammaf(1);
    if N>2
        X(2:N-1)=gammaf(3:N)-gammaf(1:N-2);
    end
    X(N)=gammaf(N)-gammaf(N-1);
end
X=X*(M/(2*1i));
X(N+1:M)=-nuw*(M*pi)/((M-N-1)*1i);
X([N+1,M])=0.5*X([N+1,M]);

% tic
% A=Afun(1:M,1:M);
% toc

% tic
pot = ifmm_mv(F,X,@Afun);
% t = toc
% Z = A*X;
% e = norm(pot - Z)/norm(Z)
pot = pot(1:N);
% pot=conj(pot');

    function A = Afun(i,j)
        [I,J] = ndgrid(i,j);
        dx = bsxfun(@minus,x(1,i)',x(1,j));
        dy = bsxfun(@minus,x(2,i)',x(2,j));
        %          ddelta = repmat(delta_f(j),length(i),1);
        %          A = -1/(N*2*pi)*(dx-dy*1i)./dz;
        %          A = -1/(2*pi*N)./(dx + dy*1i);
        %          A = -1/(N*2*pi)*(dx-dy*1i)./(dx.^2 + dy.^2+0.2^2);
        %          A = -1/(N*2*pi)*(dx-dy*1i)./(dx.^2 + dy.^2+ddelta.^2);
        
        %          A = -1/(2*pi)*log(sqrt(dx.^2 + dy.^2))/N; %% Ken's
        
        dz = dx.^2 + dy.^2;
        ddelta = repmat(delta(j),length(i),1);
        if j
            dz = bsxfun(@plus,dz,ddelta.^2);
        end
        A = -1/(M*2*pi)*(dx-dy*1i)./dz;
        A(I == J) = 0;
    end

    function K = Kfun(pxy,j)
        %         dx = bsxfun(@minus,x(1,i)',x(1,j));
        %         dy = bsxfun(@minus,x(2,i)',x(2,j));
        %         K = (dx-dy*1i)./(dx.^2 + dy.^2+delta.^2);
        %         K = (dx-dy*1i)./(dx.^2 + dy.^2+0.2^2);
        %         K = 1./(dx + dy*1i);
        
        %         K = log(sqrt(dx.^2 + dy.^2)); %% Ken's
        
        dx = bsxfun(@minus,pxy(1,:)',x(1,j));
        dy = bsxfun(@minus,pxy(2,:)',x(2,j));
        dz = dx.^2 + dy.^2;
        ddelta = repmat(delta(j),size(dz,1),1); 
        %if isempty(j) % bugs here?
        %    ddelta = 0*dz;
        %end

        K1 = (dx - dy*1i)./bsxfun(@plus,dz,ddelta.^2);
        K2 = (dx - dy*1i)./(dz + mean_delta);
        K = [K1; K2];
    end
end
