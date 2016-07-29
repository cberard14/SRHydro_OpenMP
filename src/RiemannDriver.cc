#include <iostream>
#include <cmath>
#include <complex>
#include <boost/array.hpp>
#include <boost/numeric/odeint.hpp>
#include <gsl/gsl_roots.h>

class Riemann
{
public:
	Riemann(double leftP, double rightP, double leftrho, double rightrho, \
		double leftvx, double rightvx, double leftvy, double rightvy, double gam);

	double PL,PR,rhoL,rhoR,vxL,vxR,vyL,vyR,gamma;

	double shockvelocity(double Pg, int sign);
	double rarefactionvelocity(double Pg, int sign);

};

Riemann::Riemann(double leftP, double rightP, double leftrho, double rightrho, \
		double leftvx, double rightvx, double leftvy, double rightvy, double gam)
{
	PL=leftP;
	PR=rightP;
	rhoL=leftrho;
	rhoR=rightrho;
	vxL=leftvx;
	vxR=rightvx;
	vyL=leftvy;
	vyR=rightvy;
	gamma=gam;
}

double Riemann::shockvelocity(double Pg,int sign)
{
	double vshock,e,h,hb,W,A,B,C;
	double j,VSS,WSS;
	double P,rho,vx,vy;

	// Assign P, rho, vx, vy from
	// appropriate subset of public Riemann variables
	// variable 'sign' is determined in function FP
	// a non-Riemann class member function which
	// determines the form of the equation
	// whose root we needs to solve.

	if (sign==-1)
	{
		P=PL;
		rho=rhoL;
		vx=vxL;
		vy=vyL;
	}
	else if (sign==1)
	{
		P=PR;
		rho=rhoR;
		vx=vxR;
		vy=vyR;
	}
	
	// get velocity behind shock as function of Pguess 
	e = P/((gamma-1.0)*rho);
	h = 1.0+e+P/rho;
	W = 1.0/(sqrt(1.0-(pow(vx,2)+pow(vy,2))));
	
	A = (1.0+(gamma-1.0)*(P-Pg)/(gamma*Pg));
	B = -(gamma-1.0)*(P-Pg)/(gamma*Pg);
	C = (P-Pg)*h/(rho)-pow(h,2);

	// Marti & Muller 2003, eq 31 - positive root is unique, so this approach works
	if ((-B+(sqrt(pow(B,2)-4.0*A*C)))/(2.0*A) > 0.0) 
	{
	  hb = (-B+(sqrt(pow(B,2)-4.0*A*C)))/(2.0*A);
	}
	else if ((-B-(sqrt(pow(B,2)-4.0*A*C)))/(2.0*A) > 0.0 ) 
	{
	  hb = (-B-(sqrt(pow(B,2)-4.0*A*C)))/(2.0*A);
	}
	else
	  hb = 0.0;

	if ( ((2.0*h/rho*(P-Pg))-(pow(h,2)-pow(hb,2))) >= 0.0){j = sign*(sqrt(-pow((P-Pg),2)/-(1.0e-8)));}
	else{
	  j = sign*(sqrt(-pow((P-Pg),2)/  ((2.0*h/rho*(P-Pg))-(pow(h,2)-pow(hb,2)))  ) );
	}
	
	VSS = ( pow((rho*W),2)*vx+sign*fabs(j)*(sqrt(pow(j,2)+
							 pow(rho,2)*pow(W,2)*(1.0-pow(vx,2)))))/(pow(rho,2)*pow(W,2)+pow(j,2));
	
	WSS = 1.0/(sqrt(1.0-pow(VSS,2)));
	vshock = (h*W*vx+(WSS*(Pg-P))/j)*pow((h*W+(Pg-P)*(WSS*vx/j+1.0/(rho*W))),-1);
	return vshock;
}


void RfODE(const boost::array<double,6> &vrhovt, \
	boost::array<double,6> &dvrhovtdp, double t){
	// this function outputs dvrhovtdp (an array of 3)
	// dvrhovtdp[0], vrhovt = dvdp, v
	// dvrhovtdp[1], vrhovt = drhodp, rho
	// dvrhovtdp[2], vrhovt = dvtdp, vt
	// we want v (vrhovt[0]) at Pguess (pressure for which FP=0)
	double Prf, vxrf, rhorf, vyrf, gamrf;
	double dprfdp, dvxdp, drhodp, dvydp, dgamrfdp;
	int sign, dsigndp;
	double maxv=0.99999;
	// intermediate variables
	double h,e,vtot,Wtot,Wx,ksi,cs,g;

	// fill variables
	vxrf=vrhovt[0];
	rhorf=vrhovt[1];
	vyrf=vrhovt[2];

	// dummy variables for purpose of importing [vx,rho,vy,sign,P,gamma]
	// only v (vrhovt[0]) is actually wanted from the ODE
	sign=vrhovt[3];    // constant
	Prf=vrhovt[4];     // vary linearly by ODE solver (independent variable)
	gamrf=vrhovt[5];   // constant

	// ensure rho > 0.0 else terrible things happen
	if (rhorf<0.0)
	{
		rhorf*=-1.0;
	}

	// Fill intermediate variables
	e = Prf/((gamrf-1.0)*rhorf);
	h = 1.0+e+Prf/rhorf;

	if (rhorf <= 1.0e-15 || Prf <= 1.0e-15)
	{
		cs = sqrt(gamrf-1.0);
	}
	else
	{
		cs = sqrt(gamrf*Prf/(h*rhorf));
	}

	vtot= sqrt(pow(vxrf,2)+pow(vyrf,2));
	if (vtot>=1.0){vtot=maxv;}
	if (vxrf>=1.0 && vyrf==0.0){vxrf=maxv;}

	Wtot = 1.0/(sqrt(1.0-pow(vtot,2)));
	Wx = 1.0/(sqrt(1.0-pow(vxrf,2)));
	if ( (pow(vtot*cs,2)+pow(vxrf,2)*(1.0-pow(cs,2))) > 1.0) {ksi = (vxrf*(1.0-pow(cs,2))+sign*cs*(  sqrt(1.0e-10)  ))/(1.0-pow(vtot*cs,2)); }
	else {ksi = (vxrf*(1.0-pow(cs,2))+sign*cs*(  sqrt((1.0-pow(vtot,2))*(1.0-pow(vtot*cs,2)-pow(vxrf,2)*(1.0-pow(cs,2))))  ))/(1.0-pow(vtot*cs,2));} // x/t PMM3.6 

	g = pow(vyrf,2)*(pow(ksi,2)-1.0)/pow((1.0-ksi*vxrf),2); //PMM 3.11
	if (g<-1.0){g=-0.999999;}

	dsigndp = 0.0;
	dprfdp = 1.0;
	dgamrfdp = 0.0;
	
	dvxdp = sign/(rhorf*h*pow(Wtot,2)*cs*(sqrt(1.0+g))); //3.10 Pons Marti Muller 2000 (PMM 3.10)
	drhodp = 1.0/(h*pow(cs,2));  //PMM 3.5
	dvydp = -(vyrf/h*(1.0+1.0/(gamrf-1.0))*(1.0/rhorf-Prf*drhodp/pow(rhorf,2))+(pow(Wtot,2)*vxrf*dvxdp*vyrf))/(1.0+pow(vyrf,2)*pow(Wtot,2));

	// pass derivatives to derivative array --> OUT VARIABLES
	dvrhovtdp[0]=dvxdp;
	dvrhovtdp[1]=drhodp;
	dvrhovtdp[2]=dvydp;

	// dummy variables to match dimension of [vx,rho,vy,sign,P,gamma] array
	dvrhovtdp[3]=dsigndp;  // 0
	dvrhovtdp[4]=dprfdp;   // 1
	dvrhovtdp[5]=dgamrfdp; // 0
}

double Riemann::rarefactionvelocity(double Pg, int sign){
  
        double dpsize;
	double P,rho,vx,vy,vtest;
	// get public variables from Riemann class instance
	if (sign==-1)
	{
		P=PL;
		rho=rhoL;
		vx=vxL;
		vy=vyL;
		if (Pg>PL){dpsize=0.0001;}
		else{dpsize=-0.0001;}
	}
	else if (sign==1)
	{
		P=PR;
		rho=rhoR;
		vx=vxR;
		vy=vyR;
		if (Pg>PR){dpsize=0.0001;}
		else{dpsize=-0.0001;}
	}
	boost::array<double,6> vrhovt = {vx,rho,vy,double(sign),P,gamma}; 
	boost::numeric::odeint::integrate(RfODE,vrhovt,P,Pg,dpsize);
	return vrhovt[0];
}

double FRiemann(double Pg, void* Riemannvalues){
	Riemann* Rvals = (Riemann*)Riemannvalues;
	// get propagation signs for Shock, Rf waves
	// shock(left data)-shock(right data)+Rf(left data)-Rf(right data)=0
	double vLwave,vRwave,Pleft,Pright;
	Pleft = Rvals->PL;
	Pright = Rvals->PR;
	// get left side wave form
	if ( Pg < Pleft)
	{
	  vLwave = Rvals->rarefactionvelocity(Pg,-1);
	}
	else
	{
	  vLwave = Rvals->shockvelocity(Pg,-1);
	}

	// get right side wave form
	if (Pg < Pright)
	{
	  vRwave = Rvals->rarefactionvelocity(Pg,1);
	}
	else
	{
		vRwave = Rvals->shockvelocity(Pg,1);
	}
	return vLwave-vRwave;
}


void RiemannDriver(double gamma, double &PL, double &PR, \
				double &rhoL, double &rhoR, double &vxL, \
				double &vxR, double &vyL, double &vyR, \
		   double &P_star, double &vx_star)
{
	double Phigh = 1200.0;
	double Plow = 0.0001;
	double Proot;
	double precision = 1.0e-5;

	Riemann RiemannValues(PL,PR,rhoL,rhoR,vxL,vxR,vyL,vyR,gamma);

	if (fabs(PL-PR) < 1.0e-5 && fabs(vxL-vxR) < 1.0e-5)
	  {
	    P_star=PL;
	    vx_star=vxL;
	  }
	else
	  {
	    gsl_root_fsolver* solver;
	    gsl_function fwrapper;
	    solver = gsl_root_fsolver_alloc(gsl_root_fsolver_bisection);
	    fwrapper.function = FRiemann;
	    fwrapper.params = &RiemannValues;
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
	    
	    P_star = Proot;
	    if (Proot <= PL)
	      {
		vx_star = RiemannValues.rarefactionvelocity(Proot,-1);
	      }
	    else
	      {
		vx_star = RiemannValues.shockvelocity(Proot,-1);
	      }
	    
	    gsl_root_fsolver_free(solver);
	    
	  }
}
