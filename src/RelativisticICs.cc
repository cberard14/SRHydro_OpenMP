#include <iostream>
#include <cmath>
#include <vector>

void GetRelativisticICs(double gamma, std::vector<double>& P, \
				std::vector<double>& rho,std::vector<double>& vx, \
				std::vector<double>& vy,std::vector<double>& vgrid, \
				std::vector<double>& vol,std::vector<double>& E_m,
				std::vector<double>& V_m,std::vector<double>& dm, \
				std::vector<double>& S_mx,std::vector<double>& S_my, \
				std::vector<double>& e_internal,std::vector<double>& h, \
			std::vector<double>& lorentz,std::vector<double>& Smxold,std::vector<double>& Vmold)
{
  	// take ICs and convert to relativistic quantities
	int n = P.size();

	for (int i=0;i<n;i++)
	{
		lorentz[i] = 1.0/sqrt(1.0-pow(vx[i],2)-pow(vy[i],2));
		e_internal[i] = P[i]/((gamma-1.0)*rho[i]);
		h[i] = 1.0+e_internal[i]+P[i]/rho[i];
		S_mx[i] = vx[i]*h[i]*lorentz[i];
		S_my[i] = vy[i]*h[i]*lorentz[i];
		E_m[i] = h[i]*lorentz[i]-(1.0/lorentz[i])*(P[i]/rho[i]);
		vol[i] = vgrid[i]*lorentz[i];
		dm[i] = rho[i]*vol[i];
		V_m[i] = vgrid[i]/dm[i];
		Smxold[i]=S_mx[i];
		Vmold[i]=V_m[i];
	}
}
