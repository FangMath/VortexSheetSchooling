function PreMOVIE(radius,tk)
% global xindex_2d xindex a Wing Para zeta
global W Wing Nw
clearvars W m n

Cirl_vortex=[];
for ib=1:Nw
    
    n(ib)=Wing(ib).Lenf(tk)
    gammaf(ib,1:n(ib))=Wing(ib).GammaF{tk};
    Cirl_vortex(ib,1)=0;
    Cirl_vortex(ib,2:n(ib))=gammaf(ib,1:n(ib)-1)-gammaf(ib,2:n(ib));
    
    zeta=Wing(ib).ZetaF{tk};
    x=[real(zeta);imag(zeta)];
    xmin = min(x,[],2)';
    xmax = max(x,[],2)';
    l = xmax - xmin;
    ngrid=ceil(l/radius)+3;
    
    tic
    xindex_2d = floor(bsxfun(@minus, x, xmin')/radius)+1;
    toc
    xindex = [1,ngrid(1)]*xindex_2d + 1; %xindex(i) is the index of box in which xi locates
    m = histc(xindex,1:ngrid(1)*ngrid(2));
    toc
    length(m)
    tic
    for in=1:length(m)
        W(ib).Box{in}=find(bsxfun(@eq, xindex,in)); % which particles are contained in box{in} 
        W(ib).Strength(in)=sum(Cirl_vortex(ib,W(ib).Box{in}));
        W(ib).Position(in)=sum(Wing(ib).ZetaF{tk}(W(ib).Box{in}))/length(W(ib).Box{in});
    end
    toc
    tic
    W(ib).plusPos=W(ib).Position(bsxfun(@gt,W(ib).Strength,0));
    W(ib).minusPos=W(ib).Position(bsxfun(@lt,W(ib).Strength,0));
    toc
msize=5;
    plot(imag(W(ib).plusPos),-real(W(ib).plusPos),'r.','MarkerSize',msize);hold on;
    plot(imag(W(ib).minusPos),-real(W(ib).minusPos),'b.','MarkerSize',msize);hold on;
%         axis([left right down up]);
%     grid on;
end
% hold off;
end
