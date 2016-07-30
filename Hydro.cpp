// Readtest.cpp

#include <iostream>
#include <fstream>
#include <cmath>
#include "./include/inifile.h"
#include "./include/ConservedToPrimitiveDriver.h"
#include "./include/RiemannDriver.h"
#include <iomanip>

/*
*****************************************************************************
*                1.5D SPECIAL RELATIVISTIC LAGRANGIAN HYDROCODE             *
*****************************************************************************
> Exact relativistic Riemann solver is an implementation of Pons, Marti and Muller (2000)
> Advection scheme based on Colella & Woodward (1984) 
*/

void GenerateSodICs( double boundary, double grid_length, \
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


void CallRiemann(double gamma, std::vector<double>& Plong, \
				std::vector<double>& rholong,std::vector<double>& vxlong, \
				std::vector<double>& vylong,std::vector<double>& Pstar, \
				std::vector<double>& vxstar)
{
	// fill Pstar, vxstar (boundary values for P, vx)
  	int n = Plong.size();

#pragma omp parallel for shared(Pstar,vxstar)
	for (int i=0;i<n-1;i++)
	{
		RiemannDriver(gamma,Plong[i],Plong[i+1],rholong[i],rholong[i+1], \
			      vxlong[i],vxlong[i+1],vylong[i],vylong[i+1],Pstar[i],vxstar[i]);
	}
}


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

void UpdateConserved(double dt,std::vector<double>& x_l,std::vector<double>& x_c, \
				std::vector<double>& x_r,std::vector<double>& dtdm, \
				std::vector<double>& E_m,std::vector<double>& Sm_x, \
				std::vector<double>& Sm_y,std::vector<double>& V_m, \
		     std::vector<double>& Pstar,std::vector<double>& vxstar, \
		     std::vector<double>& Smxold,std::vector<double>& Vmold)
{
	// update conserved quantities via Godunov advection scheme
	int n=E_m.size();
	std::vector<double> x_lnew(n), x_rnew(n), x_cnew(n);
	std::vector<double> dtdmnew(n), Sm_ynew(n), Sm_xnew(n);
	std::vector<double> E_mnew(n), V_mnew(n);

    
    for (int i=0;i<n+1;i++)
    {	
        if (i <= n-1)
	{
            x_lnew[i] = x_l[i]+vxstar[i]*dt; // fills all x_rnew
        }
        else
	{
            x_rnew[i-1] = x_r[i-1]+vxstar[i]*dt; // fills last x_rnew element
        }
    }

    for (int i=0;i<n-1;i++)
    {
        x_rnew[i] = x_l[i+1]; // fills first n-2 elements of x_rnew
    }

    for (int i=0;i<n;i++)
    {
      Smxold[i]=Sm_x[i];
      Vmold[i]=V_m[i];
      //std::cout<<"SMx old="<<Smxold[i]<<", new="<<Sm_x[i]<<std::endl;
      //std::cout<<"Vm old="<<Vmold[i]<<", new="<<V_m[i]<<std::endl;
      
        dtdmnew[i] = dtdm[i];
        Sm_ynew[i] = Sm_y[i];
        V_mnew[i] = V_m[i]+dtdm[i]*(vxstar[i+1]-vxstar[i]);
        Sm_xnew[i] = Sm_x[i]-dtdm[i]*(Pstar[i+1]-Pstar[i]);
        E_mnew[i] = E_m[i]-dtdm[i]*(Pstar[i+1]*vxstar[i+1]-Pstar[i]*vxstar[i]);
        x_cnew[i] = 0.5*(x_lnew[i]+x_rnew[i]);

        // now swap old = new!
        dtdm[i] = dtdmnew[i];
        Sm_y[i] = Sm_ynew[i];
        Sm_x[i] = Sm_xnew[i];
        E_m[i] = E_mnew[i];
        V_m[i] = V_mnew[i];
        x_c[i] = x_cnew[i];
        x_l[i] = x_lnew[i];
        x_r[i] = x_rnew[i];
    }
}

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

int main(int argc, char* argv[])
{
	Inifile parameter(argv[1]);

	// Simulation parameters 
	int ncells = parameter.get<int>("ncells", 100); 
	double gamma = parameter.get<double>("gamma", 5.0/3.0); 
	double  CFL = parameter.get<double>("CFL", 0.7);  
	int simulation_type = parameter.get<int>("simulation_type", 0); 
	int BC_type = parameter.get<int>("BC_type", 0); 
	int write_interval = parameter.get<int>("write_interval", 20);           
	bool do_PPM = parameter.get<bool>("do_PPM", false); 
	double  t_end = parameter.get<double>("t_end", 0.1); 

	// Initial Conditions 
	// Sod Shock Tube 
	// If simulation_type == 0
	double grid_length = parameter.get<double>("grid_length", 1.0); 
	double boundary_center = parameter.get<double>("boundary_center", 0.5); 
	double PL = parameter.get<double>("PL", 1000.0); 
	double PR = parameter.get<double>("PR", 0.01); 
	double rhoL = parameter.get<double>("rhoL", 1.0); 
	double rhoR = parameter.get<double>("rhoR", 1.0); 
	double vxL = parameter.get<double>("vxL", 0.0); 
	double vxR = parameter.get<double>("vxR", 0.0); 
	double vyL = parameter.get<double>("vyL", 0.0); 
	double vyR = parameter.get<double>("vyR", 0.0); 

	// Initiate arrays
	// Primitive arrays
	std::vector<double> P(ncells);
	std::vector<double> rho(ncells);
	std::vector<double> vx(ncells);
	std::vector<double> vy(ncells);
	std::vector<double> Plong(ncells+2);
	std::vector<double> rholong(ncells+2);
	std::vector<double> vxlong(ncells+2);
	std::vector<double> vylong(ncells+2);
	std::vector<double> x_l(ncells);
	std::vector<double> x_c(ncells);
	std::vector<double> x_r(ncells);
	std::vector<double> vgrid(ncells);
	std::vector<double> vol(ncells);
	std::vector<double> dm(ncells);
	std::vector<double> dtdm(ncells);
	std::vector<double> h(ncells);
	std::vector<double> cs(ncells);

	// Conserved arrays
	std::vector<double> E_m(ncells);
	std::vector<double> D_m(ncells);
	std::vector<double> Sm_x(ncells);
	std::vector<double> Sm_y(ncells);
	std::vector<double> lorentz(ncells);
	std::vector<double> V_m(ncells);
	std::vector<double> e_internal(ncells);
	std::vector<double> Smxold(ncells);
	std::vector<double> Vmold(ncells);

	// Riemann arrays (P, vx at cell boundaries)
	std::vector<double> Pstar(ncells+1);
	std::vector<double> vxstar(ncells+1);

	double tstamp;


	//******** GET INITIAL CONDITIONS ********\\
	if (simulation_type == 0)
	{
		GenerateSodICs( boundary_center, grid_length, PL, PR, rhoL, rhoR, \
  					vxL, vxR, vyL, vyR, P, rho, vx, vy, x_l, x_c, \
  					x_r, vgrid);
	}

	GetRelativisticICs(gamma,P,rho,vx,vy,vgrid,vol,E_m,V_m,dm,Sm_x,Sm_y,e_internal,h,lorentz,Smxold,Vmold);

	//******** BEGIN TIMESTEPPER ********\\
	double tstamp = 0.0;
	while (tstamp<t_end)
	{
       	        std::cout<<"tstamp="<<tstamp<<" % done="<<tstamp*100.0/t_end<<std::endl;
		ApplyBCs(BC_type, P, rho, vx, vy, Plong, rholong, vxlong, vylong);
		CallRiemann(gamma, Plong, rholong, vxlong, vylong, Pstar, vxstar);
		TimeStep(tstamp,t_end,CFL,x_l,x_r,vx,x_c,dtdm,dm,E_m,Sm_x,Sm_y,V_m,Pstar,vxstar,Smxold,Vmold);
		CallConservedToPrimitive(gamma,E_m,Sm_x,Sm_y,V_m,P,rho,vx,vy,Smxold,Vmold);
	}
	std::cout<<"All done! t="<<tstamp<<std::endl;
	std::cout<<"i        xc       P      rho      vx        vy"<<std::endl;
	for (int i=0;i<ncells;i++)
	  {
	    std::cout<<i<<" "<<P[i]<<" "<<x_c[i]<<" "<<rho[i]<<" "<<vx[i]<<" "<<vy[i]<<std::endl;
	  }
}

