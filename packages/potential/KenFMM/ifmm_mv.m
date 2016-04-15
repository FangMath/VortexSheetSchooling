% IFMM_MV       Multiply using matrix factors from the interpolative fast
%               multipole method.
%
%    Y = IFMM_MV(F,X) produces the matrix Y resulting from applying the factored
%    matrix F to the matrix X. This form requires that F have all near-field
%    interactions stored.
%
%    Y = IFMM_MV(F,X,A) applies the factorization F with all near-field
%    interactions generated from A.
%
%    See also IFMM, RSKEL, GRSKEL2, GRSKEL3.

function Y = ifmm_mv(F,X,A)

  % set default parameters
  if nargin < 3
    A = [];
  end

  % initialize
  N = F.N;
  nlvl = F.nlvl;
  store = F.store;
  rem  = logical(ones(N,1));
  rem_ = logical(zeros(N,1));
  Z = cell(nlvl,1);
  Y = cell(nlvl,1);

  % upward sweep
  Z{1} = X;
  for lvl = 1:nlvl-1
    prem1 = cumsum(rem);
    for i = F.lvpu(lvl)+1:F.lvpu(lvl+1)
      rem(F.U(i).rd) = 0;
    end
    prem2 = cumsum(rem);
    Z{lvl+1} = Z{lvl}(prem1(rem),:);
    for i = F.lvpu(lvl)+1:F.lvpu(lvl+1)
      rd  = prem1(F.U(i).rd);
      sk1 = prem1(F.U(i).sk);
      sk2 = prem2(F.U(i).sk);
      Z{lvl+1}(sk2,:) = Z{lvl+1}(sk2,:) + F.U(i).T*Z{lvl}(rd,:);
    end
  end

  % downward sweep
  Y{nlvl} = zeros(size(Z{nlvl}));
  for lvl = nlvl-1:-1:1
    Y{lvl} = zeros(size(Z{lvl}));
    prem2 = cumsum(rem);
    for i = F.lvpu(lvl)+1:F.lvpu(lvl+1)
      rem(F.U(i).rd) = 1;
    end
    prem1 = cumsum(rem);
    Y{lvl}(prem1(rem_),:) = Y{lvl+1};
    for i = F.lvpu(lvl)+1:F.lvpu(lvl+1)
      rem_(F.U(i).rd) = 1;
    end
    for i = F.lvpu(lvl)+1:F.lvpu(lvl+1)
      rd  = prem1(F.U(i).rd);
      sk1 = prem1(F.U(i).sk);
      sk2 = prem2(F.U(i).sk);
      Y{lvl}(rd,:) = F.U(i).T'*Y{lvl+1}(sk2,:);
      Y{lvl}(sk1,:) = Y{lvl+1}(sk2,:);
    end
    for i = F.lvpd(lvl)+1:F.lvpd(lvl+1)
      if store
        j = prem1(F.D(i).i);
        k = prem1(F.D(i).j);
        Y{lvl}(j,:) = Y{lvl}(j,:) + F.D(i).D*Z{lvl}(k,:);
      else
        j = F.D(i).i;
        k = F.D(i).j;
        D = A(j,k);
        j = prem1(j);
        k = prem1(k);
        Y{lvl}(j,:) = Y{lvl}(j,:) + D*Z{lvl}(k,:);
      end
    end
  end

  % extract output
  Y = Y{1};
end