#ifndef SODINITIALCONDITIONS_H
#define SODINITIALCONDITIONS_H

// Holds the functions meant to generate Sod shocktube initial conditions 
// Sod Initial condition parameters given by SodParams.txt

void GenerateSodICs( double boundary, double grid_length, \
  				double PL, double PR, double rhoL, double rhoR, \
  				double vxL, double vxR, double vyL, double vyR, \
  				std::vector<double>& P, std::vector<double>& rho, \
  				std::vector<double>& vx, std::vector<double>& vy, \
  				std::vector<double>& x_l, std::vector<double>& x_c, \
  				std::vector<double>& x_r, std::vector<double>& vgrid);

#endif