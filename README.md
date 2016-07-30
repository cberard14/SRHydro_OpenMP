# Special Relativistic 1.5 Dimensional Lagrangian Hydrocode

## Introduction: what is it, why was it made?
Welcome to the page for my special relativistic Hydrocode! It is a Lagrangian Hydrocode that uses a conservative scheme for updating cell values. It uses a special relativistic Riemann solver which works with sideflow -- hence, it is "1.5" dimensional. The conservative scheme was written following the PPM (piecewise parabolic method) paper from Collela & Woodward 1984, while the exact Riemann solver was written following the approach of Pons, Marti and Muller 2000.

This project began as a part of my 4th year thesis project (advised by Chris Matzner at University of Toronto). At the end of the term, we had a working python code and a working Fortran version (with the added feature of PPM) that could do some pretty cool things -- but this did not accomplish the original goal! The original goal was to develop a hydrocode that was relativistic, and that would be used to model a particular kind of shock breakout in a core-collapse supernova, which required the added exotic property of handling movement perpendicular to the 1D grid. Over the course of the summer (2015), a special relativistic version in Python was written, and ultimately put away after it was found that it was far too slow to solve the problem we had set out to study. Unfortunately, a slower code meant the capability of dealing with fewer particle cells, which meant insufficient accuracy.

This C++ version was born out of a labour of love, and a necessity to see what kind of speed (and therefore, accuracy) improvement we might see writing in a compiled language with better coding practices and the further improvement of OpenMP. 
Plots of this speed/accuracy improvement are to come!

## Dependencies
I tested this code using Scinet, on the GPC. It has been tested and works there. The modules I've loaded to get the code to work are in the "setup" file. Simply by typing "source setup" the appropriate modules will be loaded. 

Module list:

gcc/4.8.1

openmpi/gcc/1.6.4

boost/1.54.0-gcc4.8.1

gsl/1.16-gcc

## Documentation: what does everything do?
For now, please see the comments in the header files for each function; these give a summary of what the functions are supposed to do.

## How to Run
Type "make" to compile the code. This gives 3 executables; "Hydro", the main hydrocode, "RiemannTester" which allows one to test problematic points for the Riemann solver, and "ConservativeToPrimitiveTester" which does the same but tests the routine which maps conservative variables to their primitive counterparts. Conservative variables are the variables that are updated directly (energy (Em), x- and y-momentum (Smx, Smy, where x is along the grid, y is perpendicular), volume (Vm)), while primitive variables are pressure (P), density (rho), x- and y-velocity (vx, vy) and are what the Riemann solver takes as input and returns as output in order to update the conserved variables. 

To run:
choose a number of threads, n, to run the code on. This depends on your hardware.

export OMP_NUM_THREADS=n

To run the full simulation, give the Hydro executable and the parameters file. Right now, only Sod shocktube conditions are supported, so use SodParamts.txt, like so:

./Hydro SodParams.txt 

The code will then run, outputting the timestep and percentage completed, if you're impatient like me and like to watch the progress! You will notice a slight slow-down as the code progresses. This is because the code is designed not to update cells which are not being acted upon (because it's pointless to waste time on the Riemann solver and variable mapping if nothing has changed in the conditions!) and so as the solution evolves, more cells are affected, and more cells therefore need to be updated by the expensive root-finding and differential-equation solving routines.

Similarly, to test problematic points, you can type:

./RiemannTester RiemannValues.txt

where RiemannValues is a set of problematic Riemann solver points. You would do the same for the ConservedToPrimitive variable mapper. As of yet, the code is not written to determine which points these are -- the hope is that I have identified places where the Riemann solver breaks (returns NaNs, infs by division by 0, sqrt of negative number, etc) and that there should be little need to use these executables. My testing has so far revealed this to be true. 

PLEASE NOTE, the header "inifile.h" was not written by me -- it was used in a software course I took one semester and I can't find the attribution fot it! If you know who this file belongs to, I would like to credit the writer!

## To do list:
To do:
- write up official documentation!
- add error catching for NaNs and Infs coming out of the Riemann solver and conserved to primitive variable mapper
- add Python code for the analytic solution
- add initial conditions for Sedov Shock tube and fully custom initial conditions (ex. read initial pressure, density, velocities from file)
- add reflective, free boundary conditions
- add different symmetries (when sideflow disabled -- sideflow only works in planar symmetry)
- add pictures/visualizations of the code at work


