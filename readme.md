# Vortex sheet simulation program for flapping-wing schooling (MATLAB)
The program uses a vortex sheet algorithm to simulate the hydrodynamic interactions between self-propelled flapping wings. For details of the algorithm, please contact Fang Fang (ff559@nyu.edu) for a paper manuscript.

## Program first touch
- Download the code. Open MATLAB and go to the program directory. Run `startup.m` (add necessary directories to the program path).
- Run `Aschool(1)`. You will see a log of things in MATLAB. Stop the program anytime by `CTRL+c`, and resume by running `Aschool(0)`.
- Find your simulation data files at `../data_twowing_schooling/amp0.2skin0.04sch1/`.
- Plot some simulation results: `plotcxwh('../data_twowing_schooling/amp0.2skin0.04sch1')`.
- View the movies: `MOVIEframeD('../data_twowing_schooling/amp0.2skin0.04sch1/',1)`.

## Basic structure of the program
- `Aschool.m`: the main function. Run
  * `Aschool(1)`: Start a new simulation (from timestep tk=1);
  * `Aschool(0)`: Resume the simulation from where it stopped;
  * `Aschool(tk)`: Resume the simulation from timestep tk;

  Note:  To resume a simulation, you need to specify the directory of that simulation in the file `parameter.m`. The program would then use the parameters used in that simulation. 

- `INITIALIZE.m`: initialize the simulation with specified initial condition (initialize from static by default).
- `parameter.m`: specify physical and numerical paramters.
- `EXPLICIT.m`: update the free vortex sheets using an explicit method.
- `IMPLICIT.m`: update the bound vortex sheets (flapping wings) through an implicit solver.
- `SAVEDATA.m`: save simulation data.
 
  

