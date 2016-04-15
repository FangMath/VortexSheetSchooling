function MOVIEfromSnap(str_input,apstr,tem_at)
close all;
global Cirl_vortex Wing Sys vidObj
if nargin~=0
    str=str_input;
else
    str='./output/';
end

load( [str,'parameter.mat'],'Para');
dt=Para.dt;M=Para.M;Nw=Para.Nw;
% load( [str,'data.mat'], 'tem_at','Wing','Sys');

% Initialize data at tk=tem_at(1):
start=1;

% set axis frame
left=-10;right=10; down=-5; up=5;

%set x, y for the central line
x=[0,0]; y=[down,up];

%  plot body and free sheet at each time step
if exist([str,'movie_01',apstr,'.avi'],'file')
    delete([str,'movie_01',apstr,'.avi']);
end

%calculate circulation for each vortex

tem_at

vidObj = VideoWriter([str,'movie_01',apstr,'.avi']);
vidObj.FrameRate=50;
open(vidObj);

for i=2:tem_at-1
    % h=figure(1);
    i
    h=openfig([str,'Snapshot/T',num2str(i),'.fig'],'new','invisible');
    currFrame = getframe(h)
%     currFrame.cdata=currFrame.cdata(1:420,1:520,:);
    writeVideo(vidObj,currFrame);
    close(h);
end
close(vidObj);
end
