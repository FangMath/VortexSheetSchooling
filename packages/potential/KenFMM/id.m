% ID    Interpolative decomposition.
%
%    [SK,RD,T] = ID(A) produces skeleton and redundant indices SK and RD,
%    respectively, and an interpolation matrix T so that A(:,RD) = A(:,SK)*T.
%
%    [SK,RD,T] = ID(A,TOL) produces an interpolative decomposition to relative
%    precision TOL. The approximation error is estimated using one step of the
%    randomized power method.
%
%    See also QR, SVD.

function [sk,rd,T] = id(A,tol,iter)

  % set default parameters
  if nargin < 2 || isempty(tol)
    tol = 1e-15;
  end

  % get matrix size
  [m,n] = size(A);
  l = min(m,n);

  % quick return if matrix is empty
  if isempty(A)
    sk = [];
    rd = 1:n;
    T = zeros(0,n);
    return
  end

  % compute ID
  [~,R,E] = qr(A,0);
  S = abs(diag(R));
  M = abs(R(1))*tol;
  k = sum(S > M);
  sk = E(1:k);
  rd = E(k+1:end);
  T = R(1:k,1:k)\R(1:k,k+1:end);
%    while k < l && rerr() > rnrm()*tol
%      M = 0.5*M;
%      k = sum(S > M);
%      sk = E(1:k);
%      rd = E(k+1:end);
%      T = R(1:k,1:k)\R(1:k,k+1:end);
%    end
%
%    % estimate redundant matrix norm
%    function s = rnrm()
%      x = rand(length(rd),1);
%      x = x/norm(x);
%      y = A(:,rd)'*(A(:,rd)*x);
%      s = sqrt(norm(y));
%    end
%
%    % estimate approximation error
%    function s = rerr()
%      x = rand(length(rd),1);
%      x = x/norm(x);
%      y = A(:,sk)*(T*x) - A(:,rd)*x;
%      y = T'*(A(:,sk)'*y) - A(:,rd)'*y;
%      s = sqrt(norm(y));
%    end
end