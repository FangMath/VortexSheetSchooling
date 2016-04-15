% HYPOCT        Build hyperoctree.
%
%    T = HYPOCT(X) builds a hyperoctree T (a generalization of a binary tree in
%    1D, a quadtree in 2D, and an octree in 3D) over a set of points X such that
%    each nonempty hypercube node in T is recursively subdivided until it
%    contains a unique point. The tree T is structured as follows:
%
%        NLVL  - tree depth
%        LVP   - tree level pointer array
%        LRT   - size of tree root
%
%    It also contains the tree node data array NODES, with structure:
%
%        CTR  - node center
%        XI   - node point indices
%        PRNT - node parent
%        CHLD - node children
%        NBOR - node neighbors
%
%    Some examples of how to access the tree data are given below:
%
%      - The nodes on level I are T.NODES(T.LVP(I)+1:T.LVP(I+1)).
%      - The size of each node on level I is T.LRT / 2^(I-1).
%      - The points in node index I are X(:,T.NODES(I).XI).
%      - The parent of node index I is T.NODES(T.NODES(I).PRNT).
%      - The children of node index I are [T.NODES(T.NODES(I).CHLD)].
%      - The neighbors of node index I are [T.NODES(T.NODES(I).NBOR)].
%
%    T = HYPOCT(X,OCC) builds a hyperoctree such that each nonempty node
%    contains at most OCC points.

function T = hypoct(x,occ)

% set default parameters
if nargin < 2 || isempty(occ)
    occ = 1;
end

% initialize
[d,n] = size(x);
xmin = min(x,[],2)';
xmax = max(x,[],2)';
l = max(xmax - xmin);
s = struct('ctr',0.5*(xmin+xmax),'xi',1:n,'prnt',[],'chld',[],'nbor',[]);
T = struct('nlvl',1,'lvp',[0 1],'lrt',l,'nodes',s);
nlvl = 1;
nbox = 1;
mlvl = 1;
mbox = 1;
flag = n > occ;

% loop over all boxes in the tree
while flag
    l = 0.5*l;
    flag = 0;
    for prnt = T.lvp(nlvl)+1:T.lvp(nlvl+1)
        xi = T.nodes(prnt).xi;
        xn = length(xi);
        
        % subdivide box if it contains too many points
        if xn > occ
            flag = 1;
            ctr = T.nodes(prnt).ctr;
            idx = bsxfun(@gt, x(:,xi), ctr');
            idx = 2.^((1:d) - 1)*idx + 1;
            m = histc(idx,1:2^d);
            for i = 1:2^d
                if m(i) > 0
                    nbox = nbox + 1;
                    while mbox < nbox
                        e = cell(mbox,1);
                        s = struct('ctr',e,'xi',e,'prnt',e,'chld',e,'nbor',e);
                        T.nodes = [T.nodes; s];
                        mbox = 2*mbox;
                    end
                    s = struct( 'ctr', ctr + l*(bitget(i-1,1:d) - 0.5), ...
                        'xi', xi(idx == i),                    ...
                        'prnt', prnt,                            ...
                        'chld', [],                              ...
                        'nbor', []);
                    T.nodes(nbox) = s;
                    T.nodes(prnt).chld = [T.nodes(prnt).chld nbox];
                end
            end
            T.nodes(prnt).xi = [];
        end
    end
    
    % update if any new boxes added
    if flag
        nlvl = nlvl + 1;
        T.nlvl = nlvl;
        while mlvl < nlvl
            T.lvp = [T.lvp zeros(1,mlvl)];
            mlvl = 2*mlvl;
        end
        T.lvp(nlvl+1) = nbox;
    end
end
T.lvp = T.lvp(1:nlvl+1);
T.nodes = T.nodes(1:nbox);

% find neighbors of each box
l = T.lrt;
for lvl = 2:nlvl
    l = 0.5*l;
    for i = T.lvp(lvl)+1:T.lvp(lvl+1)
        prnt = T.nodes(i).prnt;
        j = T.nodes(prnt).chld;
        T.nodes(i).nbor = j(j ~= i);
        for j = T.nodes(prnt).nbor
            if ~isempty([T.nodes(j).xi])
                T.nodes(i).nbor = [T.nodes(i).nbor j];
            end
        end
        idx = [T.nodes(T.nodes(prnt).nbor).chld];
        c = reshape([T.nodes(idx).ctr],d,[])';
        dist = round(abs(bsxfun(@minus,T.nodes(i).ctr,c)) / l);
        j = idx(max(dist,[],2) <= 1);
        T.nodes(i).nbor = [T.nodes(i).nbor j];
    end
end
end