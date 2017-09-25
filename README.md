# IEC_Optimization
An MATLAB optimizer that runs CUDA simulation code through the MEX interface.  The simulation is a nuclear fusor.
CUDA functions are located in the "headers" directory.

## CUDA Overview
This is a 2D3V simulation of ions in the IEC and contains the following functions.  (1), (2), and (5) are the most computationally intensive, and all but (7) are CUDA kernels

1. Find Density
Interpolates the charge of each particle to the grid, using `atomicAdd` to avoid race conditions
2. Calculate Potential
Iterative linear solver for Poisson's equation
3. Calculate Electric Field
Gradient of the Potential
4. Interpolate Acceleration to Particles
Reverse of (1) except race conditions not possible, so it runs quite quickly
5. Perform Collisions
Randomly pairs off particles within each cell, performs the collision method outlined in my paper [Coulomb collision model for use in nonthermal plasma simulation] (https://journals.aps.org/pre/abstract/10.1103/PhysRevE.95.063209)
6. Move Particles
Updates particles using the Boris method for updating velocity in an electric field and the Birdsall method for moving particles in a 2D3V curvilinear coordinate system
7. Kill/Create Particles
Kills particles when they leave the domain and replaces them with new particles when the new period starts.

## How to use
1. Run "IEC_CUDA_2D3V_Optimizer.m", setting "Optimize = true" in line 13.
2. A new folder will automatically be generated called "Opt_datetime" where datetime will be the current date and time.  In this folder you will find (1) a file called "RUN_ME_Opt.m" (2) a directory called "data" storing the data output of each evaluation of the cost function, (3) a directory named "optimizedVoltages" and (4) text file whos name starts with "notes".
3. For a visualization of the optimization process, run "RUN_ME_Opt.m", turning on video recording if desired.
4. To test the optimization results, open "IEC_CUDA_2D3V_Optimizer.m" with "Optimize = false" in line 13, and in line 30 enter the path of the optimized voltages (from 2.(3) above) and run.
5. A new folder will automatically be generated (as in 2 above) containing "RUN_ME_Run.m" which can be run to visualize the performance of the optimized voltages.
