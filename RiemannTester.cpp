#include "./include/RiemannDriver.h"
#include "./include/inifile.h"

int main(int argc, char* argv[])
{
  Inifile parameter(argv[1]);
  
  double gamma = parameter.get<double>("gamma", 5.0/3.0);
  double PL = parameter.get<double>("PL", 1000.0);
  double PR = parameter.get<double>("PR", 0.01);
  double rhoL = parameter.get<double>("rhoL", 1.0);
  double rhoR = parameter.get<double>("rhoR", 1.0);
  double vxL = parameter.get<double>("vxL", 0.0);
  double vxR = parameter.get<double>("vxR", 0.0);
  double vyL = parameter.get<double>("vyL", 0.0);
  double vyR = parameter.get<double>("vyR", 0.0);

  double Pstar,vxstar;
  
  RiemannDriver(gamma,PL,PR,rhoL,rhoR,vxL,vxR,vyL,vyR,Pstar,vxstar);

  std::cout<<"P*="<<Pstar<<" , v*="<<vxstar<<std::endl;
}
