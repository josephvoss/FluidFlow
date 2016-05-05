#include "Simulation.h"
#include <iostream>
//#include "H5Cpp.h"

//Pressure function can't be solved independently - how solve?
//Solve for velocity iteratively - then use points? Issue is if
//no stable solution exists

int main(int argc, char** argv)
{
	MPI::Init();
	
//	std::cout<<"Hello?"<<std::endl;

//	H5::H5File file("./output.txt", HDF_ACC_RDWR);
//	H5::Dataset dataSet = file.openDataSet("dset");

	Simulation* workBench = new Simulation;
	workBench->iterate();
//	std::cout<<workBench->solvedVelData[50][50].u<<std::endl;
//	dataSet.write(workBench->solvedVelData,H5::PredType::NATIVE_DOUBLE);

	
	return 0;
}
