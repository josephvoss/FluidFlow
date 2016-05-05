#include "Simulation.h"
#include <iostream>

//Pressure function can't be solved independently - how solve?
//Solve for velocity iteratively - then use points? Issue is if
//no stable solution exists

int main(int argc, char** argv)
{
	MPI::Init();
	
	std::cout<<"Hello?"<<std::endl;
	Simulation* workBench = new Simulation;
	workBench->init();
	workBench->iterate();

	std::cout<<workBench->solvedVelData[99][5].u<<std::endl;

	
	return 0;
}
