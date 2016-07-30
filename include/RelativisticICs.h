#ifndef RELATIVISTICICS_H
#define RELATIVISTICICS_H

// this function converts the initial conditions (P, rho, vx vy)
// to their corresponding energies, momenta, volume
// can be called regardless of initial condition type

void GetRelativisticICs(double gamma, std::vector<double>& P, \
				std::vector<double>& rho,std::vector<double>& vx, \
				std::vector<double>& vy,std::vector<double>& vgrid, \
				std::vector<double>& vol,std::vector<double>& E_m,
				std::vector<double>& V_m,std::vector<double>& dm, \
				std::vector<double>& S_mx,std::vector<double>& S_my, \
				std::vector<double>& e_internal,std::vector<double>& h, \
			std::vector<double>& lorentz,std::vector<double>& Smxold,std::vector<double>& Vmold);

#endif