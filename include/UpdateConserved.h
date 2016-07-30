#ifndef UPDATECONSERVED_H
#define UPDATECONSERVED_H
#include <vector>
// this function updates the conserved variables (energy, volume, momenta)
// with a Godunov advective scheme
// see Collela & Woodward 1984 for the details

void UpdateConserved(double dt,std::vector<double>& x_l,std::vector<double>& x_c, \
				std::vector<double>& x_r,std::vector<double>& dtdm, \
				std::vector<double>& E_m,std::vector<double>& Sm_x, \
				std::vector<double>& Sm_y,std::vector<double>& V_m, \
		     std::vector<double>& Pstar,std::vector<double>& vxstar, \
		     std::vector<double>& Smxold,std::vector<double>& Vmold);

#endif
