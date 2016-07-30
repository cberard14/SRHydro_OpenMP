#ifndef BOUNDARYCONDITIONS_H
#define BOUNDARYCONDITIONS_H

// function calls the boundary conditions at each timestep
// effectively just updates the ghost cells (cells indexed at -1 and N)
// with the appropriate values depending on whether BCs are fixed, reflective, free

void ApplyBCs(int BC, std::vector<double>& P, std::vector<double>& rho, \
			  std::vector<double>& vx, std::vector<double>& vy, \
			  std::vector<double>& Plong, std::vector<double>& rholong, \
			  std::vector<double>& vxlong, std::vector<double>& vylong);

#endif