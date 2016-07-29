#include "./include/ConservedToPrimitiveDriver.h"
#include "./include/inifile.h"

int main(int argc, char* argv[])
{
  Inifile parameter(argv[1]);
  
  double gamma = parameter.get<double>("gamma", 5.0/3.0);
  double Em = parameter.get<double>("Em", 1400.0);
  double Smx = parameter.get<double>("Smx", 0.0);
  double Smy = parameter.get<double>("Smy", 0.0);
  double Vm = parameter.get<double>("Vm", 1.0);

  std::vector<double> Prhovxvy(4);

  ConservedToPrimitiveDriver(gamma, Em, Smx, Smy, Vm, Prhovxvy[0], Prhovxvy[1], Prhovxvy[2], Prhovxvy[3]);

  std::cout<<"P="<<Prhovxvy[0]<<std::endl;
  std::cout<<"rho="<<Prhovxvy[1]<<std::endl;
  std::cout<<"vx="<<Prhovxvy[2]<<std::endl;
  std::cout<<"vy="<<Prhovxvy[3]<<std::endl;
}