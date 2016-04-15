function index_piece(wing,ib,ploton)
global Cirl_vortex piece_p piece_n NumPiece
gammaf=[];

n=length(wing(ib).gammaf);
% gammaf{ib}(1:n)=wing(ib).gammaf;
gammaf{ib}=wing(ib).gammaf;
Cirl_vortex{ib}(1)=0;
%note that Gamma we defined is the minus of true circulation
Cirl_vortex{ib}(2:n)=gammaf{ib}(1:n-1)-gammaf{ib}(2:n);
pos=0;neg=0; %index of pos or neg piece
piece_p{ib}=[];piece_n{ib}=[];
for k=1:n
    if k==1
        if Cirl_vortex{ib}(k)==0
            if Cirl_vortex{ib}(k+1)>=0 %determine the first piece
                pos=1; piece_p{ib}(pos,1)=k; %start of pos
            else
                neg=1; piece_n{ib}(neg,1)=k; %start of neg
            end
            
        else
            if Cirl_vortex{ib}(k)>0 %determine the first piece
                pos=1; piece_p{ib}(pos,1)=k; %start of pos
            else
                neg=1; piece_n{ib}(neg,1)=k; %start of neg
            end
        end
    end
    if k<n&&Cirl_vortex{ib}(k)>0&&Cirl_vortex{ib}(k+1)<0
        neg=neg+1; %pos piece ends, neg piece begins
        piece_p{ib}(pos,2)=k;
        piece_n{ib}(neg,1)=k+1;  %end of pos, start of neg
    else if k<n&&Cirl_vortex{ib}(k)<0&&Cirl_vortex{ib}(k+1)>0
            pos=pos+1; %neg piece ends, pos piece begins
            piece_n{ib}(neg,2)=k; piece_p{ib}(pos,1)=k+1;  %end of pos, start of neg
        else if k==n&&Cirl_vortex{ib}(k)>0 %final end is pos
                piece_p{ib}(pos,2)=k;
            else if k==n&&Cirl_vortex{ib}(k)<0 %final end is neg
                    piece_n{ib}(neg,2)=k;
                end
            end
        end
    end
end

NumPiece(ib)=min(pos,neg) % # of dipoles

msize=2;
for k=1:pos
    if ploton==1
        if ib == 1
        plot(real(wing(ib).zetaf(piece_p{ib}(k,1):piece_p{ib}(k,2))),imag(wing(ib).zetaf(piece_p{ib}(k,1):piece_p{ib}(k,2))),...
            'r-','LineWidth',msize);hold on;
        end
        if ib == 2
        plot(real(wing(ib).zetaf(piece_p{ib}(k,1):piece_p{ib}(k,2))),imag(wing(ib).zetaf(piece_p{ib}(k,1):piece_p{ib}(k,2))),...
            'r-','LineWidth',msize);hold on;
        end
    end
end

for k=1:neg
    if ploton==1
        if ib == 1
        plot(real(wing(ib).zetaf(piece_n{ib}(k,1):piece_n{ib}(k,2))),imag(wing(ib).zetaf(piece_n{ib}(k,1):piece_n{ib}(k,2))),...
            'b-','LineWidth',msize);hold on;
        end
        if ib == 2
        plot(real(wing(ib).zetaf(piece_n{ib}(k,1):piece_n{ib}(k,2))),imag(wing(ib).zetaf(piece_n{ib}(k,1):piece_n{ib}(k,2))),...
            'b-','LineWidth',msize);hold on;
        end
    end
end
if ploton==1
if ib == 1
    if ~isempty(wing(ib).PointVX)
        if wing(ib).PointVCirl(1)>0
            plot(real(wing(ib).PointVX(1:2:end)),imag(wing(ib).PointVX(1:2:end)),'r*');hold on;
            plot(real(wing(ib).PointVX(2:2:end)),imag(wing(ib).PointVX(2:2:end)),'b*');hold on;
        else
            plot(real(wing(ib).PointVX(1:2:end)),imag(wing(ib).PointVX(1:2:end)),'b*');hold on;
            plot(real(wing(ib).PointVX(2:2:end)),imag(wing(ib).PointVX(2:2:end)),'r*');hold on;
        end
    end
end

if ib == 2
    if ~isempty(wing(ib).PointVX)
        if wing(ib).PointVCirl(1)>0
            plot(real(wing(ib).PointVX(1:2:end)),imag(wing(ib).PointVX(1:2:end)),'ro');hold on;
            plot(real(wing(ib).PointVX(2:2:end)),imag(wing(ib).PointVX(2:2:end)),'bo');hold on;
        else
            plot(real(wing(ib).PointVX(1:2:end)),imag(wing(ib).PointVX(1:2:end)),'bo');hold on;
            plot(real(wing(ib).PointVX(2:2:end)),imag(wing(ib).PointVX(2:2:end)),'ro');hold on;
        end
    end
end

    axis equal;
end
