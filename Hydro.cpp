#include <iostream>
#include <cmath>
#include "./include/inifile.h"
#include "./include/SodInitialConditions.h"
#include "./include/RelativisticICs.h"
#include "./include/BoundaryConditions.h"
#include "./include/RiemannDriver.h"
#include "./include/CallRiemann.h"
#include "./include/ConservedToPrimitiveDriver.h"
#include "./include/CallConservedToPrimitive.h"
#include "./include/UpdateConserved.h"
#include "./include/TimeStep.h"
#include <iomanip>


/*
*****************************************************************************
*                1.5D SPECIAL RELATIVISTIC LAGRANGIAN HYDROCODE             *
*****************************************************************************
> Exact relativistic Riemann solver is an implementation of Pons, Marti and Muller (2000)
> Advection scheme based on Colella & Woodward (1984) 
*/


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

