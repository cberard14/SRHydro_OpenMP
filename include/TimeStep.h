#ifndef TIMESTEP_H
#define TIMESTEP_H

// contains functions needed to step the code forwards in time
// GetTimestep calculates timestep using CFL condition (ensures
// that fluid flow doesn't pass cell boundary during a timestep)
// function Timestep just passes the timestep along to UpdateConservedVariables

double GetTimestep(double CFL,double& tstamp,double t_end,std::vector<double>& x_l, \
		   std::vector<double>& x_r,std::vector<double>& vx, \
		   std::vector<double>& dtdm, std::vector<double>& dm);

void TimeStep(double& tstamp, double t_end, double CFL, \
			std::vector<double>& x_l,std::vector<double>& x_r, \
			std::vector<double>& vx,std::vector<double>& x_c, \
	                std::vector<double>& dtdm,std::vector<double>& dm,\
			std::vector<double>& E_m,std::vector<double>& Sm_x, \
			std::vector<double>& Sm_y,std::vector<double>& V_m, \
	      std::vector<double>& Pstar,std::vector<double>& vxstar, \
	      std::vector<double>& Smxold,std::vector<double>& Vmold);

#endif