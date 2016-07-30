#ifndef CALLCONSERVEDTOPRIMITIVE_H
#define CALLCONSERVEDTOPRIMITIVE_H

#include "./ConservedToPrimitiveDriver.h"

// function calls ConservedToPrimitive mapper which maps conserved variables
// to their pressures, densities, velocities (because the Riemann solver needs these!)
// parallelized with OpenMP
// threads do not need to communicate with one another, only to know that they
// are not calculating anything for the same cell!

void CallConservedToPrimitive(double gamma, std::vector<double>& E_m, \
				std::vector<double>& S_mx,std::vector<double>& S_my, \
				std::vector<double>& V_m, \
				std::vector<double>& P, std::vector<double>& rho, \
			      std::vector<double>& vx, std::vector<double>& vy,\
			      std::vector<double>& Smxold,std::vector<double>& Vmold);

#endif