function Para = parameter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   parameter.m
%   set parameters
%
%   Written by Fang Fang - 03/19/2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%-------------------------------------------------------------------------
%% physical parameters
Para.Nw=2;           % number of fish
Para.Mass=2;       % fish mass
Para.I=10;           % fish moment of inertia
Para.g=0;            % gravity

%****** vary skin friction *%
Para.skin = 2*0.02;      % skin friction coefficient
%2*[0.01, 0.02, 0.06, 0.1]
%*********************** 
%***** heaving amplitude
Para.HAmp(1) = .2;       
Para.HAmp(2) = 0;  % same fish
%*********************** 
%***** fish iniV *******
Para.FishV(1) = -13.25; 
Para.FishV(2) = Para.FishV(1); 
% amp = [.1, .15, .2, .25, .3];
% C_f=0.02: 5.7; 12.7; 23.04; 35.86; 50.01
% C_f=0.03: 4.15; 8.33; 15; 24.02; 34.32
% C_f=0.04: 3.44; 6.35; 11.16; 17.83; 25.92
% Para.skin = 2*0.02; Para.HAmp(1) = 0.2; (good parameters)

%% single wing velocity
% C_f=0.02:    7.22;   14.77;   25.4;    38.28;    52.6
% C_f=0.03:    5.22;   10.2;    17.3;    26.4;     36.9
% C_f=0.04:    4.19;   7.96;  **13.25;   20.13;    28.4
%*********************** 
%***** initial schooling number
Para.iniSch = 1;
% background flow
Para.iniV=3;
Para.numV=3; % # of period to decay

%------------------ add spring constant ---------------------
Para.spring = 1;

%-------------------------------------------------------------------------
%% fish movement
Para.PAmp = 0;     % pitching amplitude

Para.dx = -Para.FishV(1)*Para.iniSch - 2;        % initial position of fish2
Para.dy = 0;
Para.phase = 0*ones(1,Para.Nw);    % phase difference of two fish
%Para.phase = [0, .5];    % phase difference of two fish

%-------------------------------------------------------------------------
%% numerical parameters
TimeSteps=100;
Para.dt=1/TimeSteps;
Para.M=80;
Para.s=cos((0:Para.M)*pi/Para.M);
Para.tem_at=1;
%% desingularizing parameters
Para.Delta0=0.2;
% Para.Delta1=1.6*Para.FreeVelocity*Para.dt;
Para.Delta1=0.1;
Para.eta=2*Para.Delta0;
Para.p=2;
Para.perc=.2;

%-------------------------------------------------------------------------
% trancation distance -> point vortices approximation
Para.farR(2) = -2*Para.FishV(2)+2;
%Para.farR(2) = 18;
Para.farR(1) = Para.farR(2)+Para.dx;

%%%%%%%%%%% directory of data files %%%%%%%%%%%%%%%%
Para.dir = ['../data_twowing_schooling/amp', num2str(Para.HAmp(1)), 'skin', num2str(Para.skin), 'sch', num2str(Para.iniSch)];

end
