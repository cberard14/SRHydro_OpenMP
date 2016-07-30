#include <iostream>
#include <vector>

void ApplyBCs(int BC, std::vector<double>& P, std::vector<double>& rho, \
			  std::vector<double>& vx, std::vector<double>& vy, \
			  std::vector<double>& Plong, std::vector<double>& rholong, \
			  std::vector<double>& vxlong, std::vector<double>& vylong)
{
	int n=Plong.size();
	// Update ghost cell values
	for (int i=1;i<n-1;i++)
	{
		Plong[i] = P[i-1];
		rholong[i] = rho[i-1];
		vxlong[i] = vx[i-1];
		vylong[i] = vy[i-1];
	}

	// 0 --> fixed BCs
	if (BC == 0)
	{
		Plong[n-1] = P[n-3];
		Plong[0] = P[0];
		rholong[n-1] = rho[n-3];
		rholong[0] = rho[0];
		vxlong[n-1] = vx[n-3];
		vxlong[0] = vx[0];
		vylong[n-1] = vy[n-3];
		vylong[0] = vy[0];
	}
}
