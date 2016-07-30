#include <iostream>
#include "../include/ConservedToPrimitiveDriver.h"
#include <cmath>
#include <vector>

void CallConservedToPrimitive(double gamma, std::vector<double>& E_m, \
				std::vector<double>& S_mx,std::vector<double>& S_my, \
				std::vector<double>& V_m, \
				std::vector<double>& P, std::vector<double>& rho, \
			      std::vector<double>& vx, std::vector<double>& vy,\
			      std::vector<double>& Smxold,std::vector<double>& Vmold)
{
	// Find corresponding pressure, density, velocities to the conserved variables
	// fill an array 
	int n = E_m.size();
	double Pval,rhoval,vxval,vyval,Em,Smx,Smy,Vm;

#pragma omp parallel for shared(P,rho,vx,vy,Smxold,S_mx,Vmold,V_m)
	for (int i=0;i<n;i++)
	{
	  Em=E_m[i];
	  Smx=S_mx[i];
	  Smy=S_my[i];
	  Vm=V_m[i];
   
      	if ( fabs(Smxold[i]-S_mx[i])<1.0e-4 || abs(Vmold[i]-V_m[i])<1.0e-4 )
	    {
	      ConservedToPrimitiveDriver(gamma,Em,Smx,Smy,Vm,P[i],rho[i],vx[i],vy[i]);
	    }
	}
}
