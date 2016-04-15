% IFMM  Interpolative fast multipole method.
%
%    F = IFMM(A,X,OCC,KFUN,PROXY) produces a factorization F of the matrix A,
%    corresponding to interactions between the points X, by using the
%    interpolative fast multipole method with tree occupancy parameter OCC and
%    proxy function KFUN to capture the far field. Far-field interactions are
%    compressed by evaluating KFUN against artificial proxy points PROXY, which
%    are given for a reference box of unit size centered at the origin and
%    translated and scaled as appropriate. For example, if A is the 2D Laplace
%    Green's function, i.e., A(I,J) = LOG(NORM(X(:,I) - X(:,J))), then we might
%    use KFUN(X,Y) = LOG(NORM(X - Y)) with PROXY being the circle centered about
%    the origin with radius 1.5.
%
%    F = IFMM(A,X,OCC,KFUN,PROXY,TOL) produces a factorization using local
%    relative precision TOL.
%
%    F = IFMM(A,X,OCC,KFUN,PROXY,TOL,STORE) produces a factorization with all
%    near-field interactions stored if STORE = 1; they are not stored if
%    STORE = 0.
%
%    See also IFMM_MV, RSKEL, GRSKEL2, GRSKEL3, HYPOCT.

function F = ifmm(A,x,occ,Kfun,proxy,tol,store)
  start = tic;

  % set default parameters
  if nargin < 6
    tol = 1e-15;
  end
  if nargin < 7
    store = 0;
  end

  % build tree
  N = size(x,2);
  tic
  t = hypoct(x,occ);
  fprintf([repmat('-',1,80) '\n'])
  fprintf('%3s | %63.2e (s)\n', '-', toc)

  % count nonempty boxes at each level
  pblk = zeros(t.nlvl+1,1);
  for lvl = 1:t.nlvl
    pblk(lvl+1) = pblk(lvl);
    for i = t.lvp(lvl)+1:t.lvp(lvl+1)
      if ~isempty(t.nodes(i).xi)
        pblk(lvl+1) = pblk(lvl+1) + 1;
      end
    end
  end

  % initialize
  nbox = t.lvp(end);
  md = 128;
  mn = 128;
  e = cell(md,1);
  D = struct('i',e,'j',e,'D',e);
  e = cell(md,1);
  U = struct('rd',e,'sk',e,'T',e);
  F = struct('N',N,'nlvl',t.nlvl+1,'lvpd',zeros(1,t.nlvl+2), ...
             'lvpu',zeros(1,t.nlvl+2),'D',D,'U',U,'store',store);
  nd = 0;
  nu = 0;
  rem = ones(N,1);
  nrem1 = sum(rem);
  mnz = 128;
  I = zeros(mnz,1);
  J = zeros(mnz,1);
  S = zeros(mnz,1);

  % process self interactions
  nz = 0;
  for i = 1:nbox
    slf = t.nodes(i).xi;
    if ~isempty(slf)
      while md < nd + 1
        e = cell(md,1);
        s = struct('i',e,'j',e,'D',e);
        F.D = [F.D; s];
        md = 2*md;
      end
      nd = nd + 1;
      F.D(nd).i = slf;
      F.D(nd).j = slf;
      if store
        F.D(nd).D = A(slf,slf);
      end

      % find neighbors
      if t.nlvl > 1
        nbr = t.nodes(i).nbor;
        nnbr = length(nbr);
        m = 2*nnbr;
        while mnz < nz + m
          e = zeros(mnz,1);
          I = [I; e];
          J = [J; e];
          S = [S; e];
          mnz = 2*mnz;
        end
        I_ = zeros(m,1);
        J_ = zeros(m,1);
        idx = 1:nnbr;
        I_(idx) = i;
        J_(idx) = nbr;
        idx = idx + nnbr;
        I_(idx) = nbr;
        J_(idx) = i;
        S_ = ones(size(I_));
        I(nz+1:nz+m) = I_(:);
        J(nz+1:nz+m) = J_(:);
        S(nz+1:nz+m) = S_(:);
        nz = nz + m;
      end
    end
  end
  F.lvpd(2) = nd;

  % compress near field
  M = sparse(I(1:nz),J(1:nz),S(1:nz),nbox,nbox);
%   for lvl = 1:t.nlvl
%     for i = t.lvp(lvl)+1:t.lvp(lvl+1)
%       slf = t.nodes(i).xi;
%       if ~isempty(slf)
%         while md < nd + 1
%           e = cell(md,1);
%           s = struct('rd',e,'sk',e,'T',e);
%           F.U = [F.U; s];
%           md = 2*md;
%         end
% 
%         % add neighbor interactions
%         nbr = [t.nodes(find(M(:,i))).xi];
%         K1 = full(A(nbr,slf));
%         K2 = full(A(slf,nbr));
%         K = [K1; K2'];
% 
%         % add proxy interactions
%         if ~isempty(Kfun)
%           l = t.lrt / 2^(lvl - 1);
%           pxy = bsxfun(@plus,proxy*l,t.nodes(i).ctr');
%           K3 = Kfun(pxy,slf);
%           K3 = K3*mean(sqrt(sum(K.^2,2)))/mean(sqrt(sum(K3.^2,2)));
%           K = [K; K3];
%         end
% 
%         % skeletonize
%         [sk,rd,T] = id(K,tol);
% 
%         % store factors
%         nu = nu + 1;
%         F.U(nu).rd = slf(rd);
%         F.U(nu).sk = slf(sk);
%         F.U(nu).T = T;
% 
%         % restrict to skeletons
%         rem(slf(rd)) = 0;
%         t.nodes(i).xi = slf(sk);
%       end
%     end
%   end
%   F.lvpu(2) = nu;

  % process neighbor interactions
  if t.nlvl > 1
    [~,Jm] = find(M);
    Jm = unique(Jm)';
    m = length(Jm);
    while md < nd + m
      e = cell(md,1);
      s = struct('i',e,'j',e,'D',e);
      F.D = [F.D; s];
      md = 2*md;
    end
    for j = Jm
      i = find(M(:,j))';
      ixi = [t.nodes(i).xi];
      jxi = t.nodes(j).xi;
      nd = nd + 1;
      F.D(nd).i = ixi;
      F.D(nd).j = jxi;
      if store
        F.D(nd).D = A(ixi,jxi);
      end
    end
  end
  F.lvpd(3) = nd;

  % print summary
  nrem2 = sum(rem);
  lvl = t.nlvl;
  nblk = pblk(lvl) + t.lvp(lvl+1) - t.lvp(lvl);
  fprintf('%3d | %6d | %8d | %8d | %8.2f | %8.2f | %10.2e (s)\n', ...
          lvl,nblk,nrem1,nrem2,nrem1/nblk,nrem2/nblk,toc)

  % loop over tree levels
  for lvl = t.nlvl:-1:1
    tic
    clvl = t.nlvl - lvl + 1;
    nrem1 = sum(rem);

    % pull up skeletons from children
    if lvl < t.nlvl
      for i = t.lvp(lvl)+1:t.lvp(lvl+1)
        t.nodes(i).xi = [t.nodes(i).xi [t.nodes(t.nodes(i).chld).xi]];
      end
    end

    % loop over nodes
    for i = t.lvp(lvl)+1:t.lvp(lvl+1)
      slf = t.nodes(i).xi;

      % skeletonize
      if lvl > 2 && ~isempty(Kfun)
        l = t.lrt / 2^(lvl - 1);
        pxy = bsxfun(@plus,proxy*l,t.nodes(i).ctr');
        K = Kfun(pxy,slf);
      else
        K = zeros(0,length(slf));
      end
      [sk,rd,T] = id(K,tol);

      % store factors
      nu = nu + 1;
      F.U(nu).rd = slf(rd);
      F.U(nu).sk = slf(sk);
      F.U(nu).T = T;

      % restrict to skeletons
      rem(slf(rd)) = 0;
      t.nodes(i).xi = slf(sk);
    end
    F.lvpu(clvl+2) = nu;

    % process interaction lists
    if lvl > 1
      m = t.lvp(lvl+1) - t.lvp(lvl);
      while md < nd + m
        e = cell(md,1);
        s = struct('i',e,'j',e,'D',e);
        F.D = [F.D; s];
        md = 2*md;
      end
      for i = t.lvp(lvl)+1:t.lvp(lvl+1)
        ilst = t.nodes(t.nodes(i).prnt).nbor;
        ilst = ilst(ilst > t.lvp(lvl-1));
        ilst = [t.nodes(ilst).chld];
        ilst = setdiff(ilst,t.nodes(i).nbor);
        ixi = [t.nodes(ilst).xi];
        jxi = t.nodes(i).xi;
        nd = nd + 1;
        F.D(nd).i = ixi;
        F.D(nd).j = jxi;
        if store
          F.D(nd).D = A(ixi,jxi);
        end
      end
      F.lvpd(clvl+3) = nd;
    end
    nrem2 = sum(rem);

    % print summary
    nblk = pblk(lvl) + t.lvp(lvl+1) - t.lvp(lvl);
    fprintf('%3d | %6d | %8d | %8d | %8.2f | %8.2f | %10.2e (s)\n', ...
            lvl,nblk,nrem1,nrem2,nrem1/nblk,nrem2/nblk,toc)
  end
  F.D = F.D(1:nd);
  F.U = F.U(1:nu);
  fprintf([repmat('-',1,80) '\n'])
  toc(start)
end