#ifndef CALLRIEMANN_H
#define CALLRIEMANN_H

#include "../include/RiemannDriver.h"
// function calls the Riemann solver function "Riemann Driver"
// this function is parallelized with OpenMP
// threads allotted over the grid, have no need to communicate with 
// one another just in that they should not compute for the same cell!

void CallRiemann(double gamma, std::vector<double>& Plong, \
				std::vector<double>& rholong,std::vector<double>& vxlong, \
				std::vector<double>& vylong,std::vector<double>& Pstar, \
				std::vector<double>& vxstar);

#endif