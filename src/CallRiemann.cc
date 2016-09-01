#include <iostream>
#include "../include/RiemannDriver.h"

void CallRiemann(double gamma, std::vector<double>& Plong, \
				std::vector<double>& rholong,std::vector<double>& vxlong, \
				std::vector<double>& vylong,std::vector<double>& Pstar, \
				std::vector<double>& vxstar)
{
	// fill Pstar, vxstar (boundary values for P, vx)
  	int n = Plong.size();

#pragma omp parallel for
	for (int i=0;i<n-1;i++)
	{
		RiemannDriver(gamma,Plong[i],Plong[i+1],rholong[i],rholong[i+1], \
			      vxlong[i],vxlong[i+1],vylong[i],vylong[i+1],Pstar[i],vxstar[i]);
	}
}
