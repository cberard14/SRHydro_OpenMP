#ifndef CONSERVEDTOPRIMITIVEDRIVER_H
#define CONSERVEDTOPRIMITIVEDRIVER_H

// This class handles the conversion between conserved and primitive variables
// conserved: Em, Smx, Smy, Vm (energy, x momentum, y momentum, volume)
// primitive: P, rho, vx vy (pressure, density, x velocity, y velocity)

class ConservedtoPrimitive {
public:
	ConservedtoPrimitive(double xmomentum, double ymomentum, double energy, \
			   double volume, double adiabatindex);

	// member variables
	double Sx, Sy, E, V, gamma;

	// member functions 
	double lors(double Pg); 	// return lorentz factor
	double rhos(double Pg);  	// return density
	double es(double Pg);   	// return internal energy
	double Pvalue(double Pg);  	// return pressure

};

// goal is to solve F(P)=0, function FP returns current estimate for the root
// ConservedToPrimitiveDriver is called by the Hydrocode and returns the primitive variables
double FP(double Pg, void* PrimitiveValues);   
void ConservedToPrimitiveDriver(double gamma, double Em, \
						double Smx, double Smy, double Vm, \
				double& Pval,double& rhoval, double& vxval, double& vyval);

#endif
