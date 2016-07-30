#include <iostream>
#include <vector>

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
