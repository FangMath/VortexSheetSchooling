function MOVIEshow(str_inputF)
close all;
global vidObj

load([str_inputF, '.mat'],'F');
vidObj = VideoWriter([str_inputF, '.avi']);

vidObj.FrameRate=20;
open(vidObj);
tem_at=length(F)
for i=1:tem_at
    i
    currFrame = F(i);
    writeVideo(vidObj,currFrame);
end
close(vidObj);
end
