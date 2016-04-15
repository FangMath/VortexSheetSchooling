% delete nodes on free sheet according to delta_f
% if abs(zetaf(k)-zetaf(k+2))<perc*delta_f(k), then delete zetaf(k+1)
function AddVortex(ib, dt, Delta0)
global Zf Gf tk Zf_new Gf_new VelocityF VelocityF_new

t=(1:length(Zf))*dt;
space=abs(diff(Zf));
count=0;
R=0;Zf_interp={}; Gf_interp={}; tnew={}; index=[]; 
% Zf_new=Zf; Gf_new=Gf; VelocityF_new=VelocityF;
for k=1:length(space)
    if space(k)>Delta0
        count=count+1;
        drag=1;
        R(count)=floor(space(k)/Delta0); % numer of nodes inserted
        tspace=dt/(R(count)+1);
        tnew{count}=tspace*(1:R(count))+t(k);
        
        Zf_interp{count}=interp1(t(max(1,k-drag):k+1),Zf(max(1,k-drag):k+1),tnew{count},'pchip');
        Gf_interp{count}=interp1(t(max(1,k-drag):k+1),Gf(max(1,k-drag):k+1),tnew{count},'pchip');
        VelocityF_interp{count}=interp1(t(max(1,k-drag):k+1),VelocityF{ib}(max(1,k-drag):k+1),tnew{count},'pchip');
        
        index(count)=k;
    end
end

head_new=1;
t_new=t;
for ct=1:count
    if ct==1
        tail_new=index(1)+head_new-1;
        head=1; tail=index(1);
    else
        tail_new=head_new+index(ct)-index(ct-1)-1;
        head=index(ct-1)+1; tail=index(ct);
    end
    t_new(head_new:tail_new)=t(head:tail);
    Zf_new(head_new:tail_new)=Zf(head:tail);
    Gf_new(head_new:tail_new)=Gf(head:tail);
    VelocityF_new(head_new:tail_new)=VelocityF{ib}(head:tail);

    t_new(tail_new+1:tail_new+R(ct))=tnew{ct};
    Zf_new(tail_new+1:tail_new+R(ct))=Zf_interp{ct};
    Gf_new(tail_new+1:tail_new+R(ct))=Gf_interp{ct};
    VelocityF_new(tail_new+1:tail_new+R(ct))=VelocityF_interp{ct};
    
    head_new=tail_new+R(ct)+1;
end

if count>0
head=index(count)+1;
tail=length(t);
tail_new=head_new+tail-head;

t_new(head_new:tail_new)=t(head:tail);
Zf_new(head_new:tail_new)=Zf(head:tail);
Gf_new(head_new:tail_new)=Gf(head:tail);
VelocityF_new(head_new:tail_new)=VelocityF{ib}(end);
end

end
