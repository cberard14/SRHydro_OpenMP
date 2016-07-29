// contains Riemann class
#ifndef RIEMANNDRIVER_H
#define RIEMANNDRIVER_H
#include <boost/array.hpp>

// This class handles the Riemann solver
// The Riemann solver follows the algorithm of Pons, Marti, and Muller 2000
// solves for the pressure and velocity at a cell's boundary
// given left and right pressures, densities, x-velocities, y-velocities, adiabatic index

class Riemann
{
public:
	Riemann(double leftP, double rightP, double leftrho, double rightrho, \
		double leftvx, double rightvx, double leftvy, double rightvy, double gam);

	double PL,PR,rhoL,rhoR,vxL,vxR,vyL,vyR,gamma;

	double shockvelocity(double Pg, int sign);   		// return velocity if shock wave
	double rarefactionvelocity(double Pg, int sign);	// return velocity if rarefaction wave

};

// Helper functions: 
// RfODE is the ordindary differential equation that obtains the velocity for a rarefaction wave
// FRiemann is the function whos root we want to find (returns current root estimate)
// Riemann driver is called by the main Hydrocode and returns the values for pressure and velocity
void RfODE(const boost::array<double,6> &vrhovt, \
	boost::array<double,6> &dvrhovtdp, double t);
double FRiemann(double Pg, void* Riemannvalues);
void RiemannDriver(double gamma, double &PL, double &PR, \
				double &rhoL, double &rhoR, double &vxL, \
				double &vxR, double &vyL, double &vyR, \
		   double &P_star, double &vx_star);

#endif
