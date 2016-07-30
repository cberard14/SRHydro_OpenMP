#include <cmath>
#include <iostream>
#include <gsl/gsl_roots.h>

class ConservedtoPrimitive {
public:
	ConservedtoPrimitive(double xmomentum, double ymomentum, double energy, \
			   double volume, double adiabatindex);

	double Sx, Sy, E, V, gamma;

	double lors(double Pg); 
	double rhos(double Pg);  
	double es(double Pg);  
	double Pvalue(double Pg);  

};

ConservedtoPrimitive::ConservedtoPrimitive(double xmomentum, double ymomentum, double energy, \
					   double volume, double adiabatindex)
{
	Sx=xmomentum;
	Sy=ymomentum;
	E=energy;
	V=volume;
	gamma=adiabatindex;
}

double ConservedtoPrimitive::lors(double Pg)
{
	// get lorentz factor
	double vx,vy,lorentz;
	vx = Sx/(E+Pg*V);
	vy = Sy/(E+Pg*V);
	lorentz = 1.0/sqrt(1.0-pow(vx,2)-pow(vy,2));
	return lorentz;
}

double ConservedtoPrimitive::rhos(double Pg)
{
	// get primitive density
	double rho,lorentz;
	lorentz = lors(Pg);
	rho = 1.0/(V*lorentz);
	return rho;
}

double ConservedtoPrimitive::es(double Pg)
{
	// get internal energy
	double e_int,lorentz;
	lorentz = lors(Pg);
	e_int = (E-lorentz+V*Pg*(1.0-pow(lorentz,2) ))/lorentz;
	return e_int;
}


double ConservedtoPrimitive::Pvalue(double Pg)
{
	return (gamma-1.0)*(rhos(Pg))*(es(Pg));
}

double FP(double Pg, void* PrimitiveValues)
{
	ConservedtoPrimitive* Pvals = (ConservedtoPrimitive*)PrimitiveValues;
	double Pressure;
	Pressure=Pvals->Pvalue(Pg);
	return Pressure-Pg;
}

void ConservedToPrimitiveDriver(double gamma, double Em, \
						double Smx, double Smy, double Vm, \
				double& Pval,double& rhoval, double& vxval, double& vyval )
{
	double A = pow(Vm,2);
	double B = 2.0*Em*Vm;
	double CC = (pow(Em,2)-pow(Smx,2)-pow(Smy,2));
	double Plow = 0.00001; 
	double Phigh = 1500.0;
	double Proot;
	double precision=1.0e-6;

	ConservedtoPrimitive PrimitiveValues(Smx, Smy, Em,Vm, gamma);

	gsl_root_fsolver* solver;
	gsl_function fwrapper;
	solver = gsl_root_fsolver_alloc(gsl_root_fsolver_brent);
	fwrapper.function = FP;
	fwrapper.params = &PrimitiveValues;
	gsl_root_fsolver_set(solver,&fwrapper,Plow,Phigh);

	int status = 1;
	for (int iter=0; status and iter < 100; ++iter) 
	{
		gsl_root_fsolver_iterate(solver);
		Proot = gsl_root_fsolver_root(solver);
		Plow = gsl_root_fsolver_x_lower(solver);
		Phigh = gsl_root_fsolver_x_upper(solver);
		status = gsl_root_test_interval(Plow,Phigh,precision,precision);
	} 
	Pval = Proot;
	rhoval = PrimitiveValues.rhos(Proot);
	vxval = Smx/(Em+Proot*Vm);
	vyval = Smy/(Em+Proot*Vm);
	gsl_root_fsolver_free(solver); 

}
