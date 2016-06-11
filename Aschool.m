function Aschool(input_tk)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Aschool.m
%
%   This is the main function of the project FishSchool.
%
%   It basically includes the following parts:
%   * Initialization;
%   * Explicit free sheets update;
%   * Implicit bound sheets solver;
%   * Data saving;
%
%   fileDir: directory of data files
%
%   inputTk: starting timestep.
%            - If start from the beginning, use inputTk=1;
%            - If resume from where we left, use inputTk=0;
%            - Otherwise resume from timestep inputTk.
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
%%                       Initialization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if input_tk == 1 % start from the beginning
    
    [wing, Para, tk] = INITIALIZE;
    Nw = Para.Nw; s = Para.s; M = Para.M;
    
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
    
    
else % resume from timestep input_tk, or where we stopped (input_tk=0)
    
    % set the directory of the simulation to resume
    Para = parameter; fileDir = Para.dir;
    
    load( [fileDir,'/data.mat'],'tem_at', 'Wing');
    load( [fileDir,'/Ddata/tem.mat'],'tem_at');
    
    if input_tk ~= 0
        tk = input_tk
    else % if input_tk == 0 then resume from where we stopped
        tk = tem_at;
    end
    
    % start from tk=1 but with given initial conditions (rarely used)
    if input_tk == -1
        tk = 1;
    end
    
    % you could modify some parameters in the file parameter.mat before loaded
    load( [fileDir,'/parameter.mat'], 'Para');
    Para.dir = fileDir;
    %Para.Mass = mass;
    Nw = Para.Nw; s = Para.s; M = Para.M;
    
    % load date of the last 3 timesteps
    load([fileDir,'/Ddata/T',num2str(tk-2),'.mat'], 'wing'); wing_p2 = wing;
    load([fileDir,'/Ddata/T',num2str(tk-1),'.mat'], 'wing'); wing_p1 = wing;
    load([fileDir,'/Ddata/T',num2str(tk),'.mat'], 'wing');
    
end

%%--------------------------------------------------------------------%%
Para.dir
pause(3)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                  Start the main solver
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
MaxT = 5005; % maximum of timesteps

while tk <= MaxT
    
    step = tic; % run time of each step
    
    % update last three time steps
    wing_p3 = wing_p2;
    wing_p2 = wing_p1;
    wing_p1 = wing;
    
    % set background flow velocity
    %    Para.FreeVelocity(tk) = max(Para.iniV-Para.iniV/Para.numV*tk*Para.dt,0);
    Para.FreeVelocity(tk) = Para.iniV;
    
    % (result would be at (tk+1)^th time step)
    tk % current timestep
    
    %% explicit solver for free vortex sheets
    tex = tic;
    EXPLICIT(Para);
    Texplicit = toc(tex)
    
    tk = tk+1; %update to next timestep
    
    %% implicit solver for bound vortex sheets
    tim = tic;
    IMPLICIT
    Timplicit = toc(tim)
    
    %% save simulation data
    SAVEDATA(Para.dir)
    Tstep = toc(step)
end

TimeTotal = toc(Ttotal)
end
