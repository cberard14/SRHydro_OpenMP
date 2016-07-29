# SRHydro_OpenMP

Welcome to the page for my special relativistic Hydrocode! It is a Lagrangian Hydrocode that uses a conservative scheme for updating cell values. It uses a special relativistic Riemann solver which works with sideflow -- hence, it is "1.5" dimensional. 

This project began as a part of my 4th year thesis project (advised by Chris Matzner at University of Toronto), and at the end of that thesis project we had a working python code and a working Fortran version (with the added feature of PPM) that could do some pretty cool things -- but this did not accomplish the original goal! The original goal was to develop a hydrocode that was relativistic, and that could be used to model a particular kind of shock breakout in a core-collapse supernova, which required the added exotic property of being capable of modelling movement perpendicular to the 1D grid. Over the course of the summer (2015) a special relativistic version in Python was written, and ultimately put away after it was found that it was far too slow to study the problem we had set out to study. Unfortunately, a slower code meant the capability of dealing with fewer particle cells, which meant insufficient accuracy.

This C++ version was born out of a labour of love, and a necessity to see what kind of speed (and therefore, accuracy) improvement we might see writing in a compiled language with better coding practices and the further improvement of OpenMP. 
Plots of this speed/accuracy improvement are to come!

To do:
- write up official documentation!
- further modularize functions in "Hydro.cpp" -- should get their own files and headers to allow easier editing...
- add error catching for NaNs and Infs coming out of the Riemann solver and conserved to primitive variable mapper
- add Python code for the analytic solution
- add easy initial conditions for Sedov Shock tube and fully custom initial conditions
- add reflective, free boundary conditions
- add different symmetries (when sideflow disabled -- sideflow only works in planar symmetry)
- add pictures/visualizations of the code at work

