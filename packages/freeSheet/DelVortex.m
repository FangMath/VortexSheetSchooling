% delete nodes on free sheet according to delta_f
% if abs(zetaf(k)-zetaf(k+2))<perc*delta_f(k), then delete zetaf(k+1)

% #called in main.m
function DelVortex(ib, perc)
global Para Zf fang delta_f takeindex
fang=fang+1;

a=1:length(Zf);
takef=Zf(takeindex);
index=a(takeindex); %index of taken points in orinial array
space=abs(takef(3:2:2*ceil(end/2)-1)-takef(1:2:2*ceil(end/2)-3));
do=0;
for k=1:length(space)
    kindex=index(2*k);
    if space(k)<perc*delta_f{ib}(kindex)
        do=1;
        takeindex(kindex)=0;
    end
end
if do==1
    DelVortex(ib, perc);
end

end

