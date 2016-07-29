# SRHydro_OpenMP

Welcome to the page for my special relativistic Hydrocode! It is a Lagrangian Hydrocode that uses a conservative scheme for updating cell values.

It uses a special relativistic Riemann solver which works with sideflow -- hence, it is "1.5" dimensional.

More documentation to follow...

To do:
- write up official documentation!
- further modularize functions in "Hydro.cpp" -- should get their own files and headers to allow easier editing...
- add error catching for NaNs and Infs coming out of the Riemann solver and conserved to primitive variable mapper
- add Python code for the analytic solution
- add easy initial conditions for Sedov Shock tube and fully custom initial conditions
- add reflective, free boundary conditions
- add different symmetries (when sideflow disabled -- sideflow only works in planar symmetry)
- add pictures/visualizations of the code at work

