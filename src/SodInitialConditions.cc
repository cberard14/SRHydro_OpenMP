#include <iostream>
#include <vector>

void GenerateSodICs( double boundary, double grid_length,		\
  				double PL, double PR, double rhoL, double rhoR, \
  				double vxL, double vxR, double vyL, double vyR, \
  				std::vector<double>& P, std::vector<double>& rho, \
  				std::vector<double>& vx, std::vector<double>& vy, \
  				std::vector<double>& x_l, std::vector<double>& x_c, \
  				std::vector<double>& x_r, std::vector<double>& vgrid)
{
	// Read right and left P, rho, vx, vy
	// (pressure, density, x-velocity, y-velocity)
	// where x-velocity is along the grid, y is perpendicular
	// left and right conditions alotted before and after the boundary, resp.
	// Output the filled vectors
	// allot cell positions first

	int n = P.size();
	double cell_length = grid_length/n;
	for (int i=0;i<n;i++)
  	{
		x_l[i] = i*cell_length;
		x_c[i] = (i+0.5)*cell_length;
		x_r[i] = (i+1.0)*cell_length;
		vgrid[i] = x_r[i]-x_l[i];
  	}

  	for (int i=0;i<n;i++)
  	{
		if (x_c[i] < boundary)
		{
	  		P[i] = PL;
	  		rho[i] = rhoL;
	  		vx[i] = vxL;
	  		vy[i] = vyL;
		}
		else
		{
			P[i] = PR;
  			rho[i] = rhoR;
  			vx[i] = vxR;
  			vy[i] = vyR;
		}
  	}
}
