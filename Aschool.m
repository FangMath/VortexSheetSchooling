function Aschool(input_tk)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Aschool.m
%
%   This is the main function of the project FishSchool.
%
%   It can be divided in parts of
%   Initialize;
%   Explicit free sheet evolution;
%   Implicit bound sheet solver;
%   Data saving;
%
%   fileDir: directory of data files
%
%   inputTk: starting timestep.
%             If start from initial point inputTk=1; if start from the end
%             of last calculation, inputTk=0; start from inputTk else.
%
%   Written by Fang Fang - 03/19/2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all;
clearvars -except input_tk;

Ttotal = tic; % total run time

global tk cx cy omega dcx dcy domega wing_p1 wing_p2 wing_p3
global wing Para
global Nw s M

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% initialize
if input_tk == 1 % start from initialization
    
    [wing, Para, tk] = INITIALIZE;
    Nw = Para.Nw; s = Para.s; M = Para.M;
    
else % start from input_tk, or last end

    Para = parameter;
    fileDir = Para.dir;
    
    % if input_tk == 0 then start from last end
    load( [fileDir,'/data.mat'],'tem_at', 'Wing');
    load( [fileDir,'/Ddata/tem.mat'],'tem_at');
    if input_tk ~= 0
        tk = input_tk
    else
        tk = tem_at;
    end

    if input_tk == -1 % start from tk=1 but with given initial conditions
        tk = 1;
    end
    
    % modify parameters here
    load( [fileDir,'/parameter.mat'], 'Para');
    Para.dir = fileDir;
    %Para.Mass = mass;
    Nw = Para.Nw; s = Para.s; M = Para.M;
    
    % load last 3 time steps
    load([fileDir,'/Ddata/T',num2str(tk-2),'.mat'], 'wing'); wing_p2 = wing;
    load([fileDir,'/Ddata/T',num2str(tk-1),'.mat'], 'wing'); wing_p1 = wing;
    load([fileDir,'/Ddata/T',num2str(tk),'.mat'], 'wing');
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Para.dir
pause(3)
if tk == 1
    for ib = 1:Para.Nw
    wing(ib).integrateP=0;
    wing(ib).skin=0;
    wing(ib).les=0;

    wing(ib).spring=0;
    wing(ib).hydrolift=0;

    wing(ib).thrust=0;
    wing(ib).lift=0;
    wing(ib).torque=0;
end
    SAVEDATA(Para.dir)
end

%% start the main solver
MaxT = 5005;
while tk <= MaxT
    
    % update previous time steps
    wing_p3 = wing_p2;
    wing_p2 = wing_p1;
    wing_p1 = wing;
    
    step = tic;
    
    % set background flow initially
    %    Para.FreeVelocity(tk) = max(Para.iniV-Para.iniV/Para.numV*tk*Para.dt,0);
    Para.FreeVelocity(tk) = Para.iniV;
    
    % (result would be at (tk+1)^th time step)
    tk
    
    %% explicit solver for free sheet
    tex = tic;
    EXPLICIT(Para);
    Texplicit = toc(tex)
    
    tk = tk+1; %update to next time step
    
    %% implicit solver for bound sheet
    tim = tic
    IMPLICIT
    Timplicit = toc(tim)
    
    %% save result data
    SAVEDATA(Para.dir)
    Tstep = toc(step)
end

TimeTotal = toc(Ttotal)
end
