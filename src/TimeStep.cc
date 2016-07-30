#include "../include/UpdateConserved.h"
#include <iostream>
#include <cmath>

double GetTimestep(double CFL,double& tstamp,double t_end,std::vector<double>& x_l, \
		   std::vector<double>& x_r,std::vector<double>& vx, \
		   std::vector<double>& dtdm, std::vector<double>& dm)
{
	// calculate the timestep to be used
	int n = vx.size();
	double dt,min_dx,max_dvxL,max_dvxR;
	std::vector<double> dvxL(n),dvxR(n),dx(n);

    for (int i=1;i<n-1;i++)
    {
        dvxL[i] = abs(vx[i]-vx[i-1]);
        dvxR[i] = abs(vx[i]-vx[i+1]);
    }
 	
    min_dx = x_r[1]-x_l[1];
    max_dvxL=0.1;
    max_dvxR=0.1;
    for (int i=0;i<n;i++)
    {
      dx[i] = x_r[i]-x_l[i]; 
      if (dx[i] < min_dx)
	{
	  min_dx=dx[i]; // find minimum cell width
	}
      if (i>1 && i < n-1)
	{
	  if (dvxL[i]>max_dvxL){max_dvxL=dvxL[i];}
	  if (dvxR[i]>max_dvxR){max_dvxR=dvxR[i];}
	}
    }

    dt = CFL*min_dx/(1.0/sqrt(3.0)+max_dvxL+max_dvxR);
    
    if ( (dt+tstamp) > t_end)
    {
      dt = t_end-tstamp;
    }
    tstamp+=dt; // update tstamp directly
    for (int i=0;i<n;i++){dtdm[i]=dt/dm[i];}
    return dt;
}

void TimeStep(double& tstamp, double t_end, double CFL, \
			std::vector<double>& x_l,std::vector<double>& x_r, \
			std::vector<double>& vx,std::vector<double>& x_c, \
	                std::vector<double>& dtdm,std::vector<double>& dm,\
			std::vector<double>& E_m,std::vector<double>& Sm_x, \
			std::vector<double>& Sm_y,std::vector<double>& V_m, \
	      std::vector<double>& Pstar,std::vector<double>& vxstar, \
	      std::vector<double>& Smxold,std::vector<double>& Vmold)
{
	int n = Sm_x.size();
	double dt;

	dt=GetTimestep(CFL,tstamp,t_end,x_l,x_r,vx,dtdm,dm);
	UpdateConserved(dt,x_l,x_c,x_r,dtdm,E_m,Sm_x,Sm_y,V_m,Pstar,vxstar,Smxold,Vmold);
}
