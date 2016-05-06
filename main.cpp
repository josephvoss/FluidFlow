#include "Simulation.h"
#include <iostream>
//#include "H5Cpp.h"

//Pressure function can't be solved independently - how solve?
//Solve for velocity iteratively - then use points? Issue is if
//no stable solution exists

int main(int argc, char** argv)
{
	MPI_Init(&argc, &argv);
	
//	std::cout<<"Hello?"<<std::endl;

	Simulation* workBench = new Simulation;
	workBench->iterate();

	if (workBench->getRank() == 0)
	{
		H5::H5File file("./output.txt", HDF_ACC_RDWR);
		hsize_t dimsf[2];
		dimsf[0] = workBench->getNt();
		dimsf[1] = workBench->getProblemSize();
		H5::Dataspace dataSpace(2, dimsf);
	 
		H5::Dataset dataSet = file.createDataSet("Data", H5::PredType::NATIVE_DOUBLE,dataSpace);

		dataSet.write(workBench->solvedVelData, H5::PredType::NATIVE_DOUBLE);
	}
/*
	int t,i,j;
	for(t=0;t<workBench->getNt();t++)
		for(i=0;i<ny;i++)
			for(j=0;j<nx;j++)
				

	if (workBench->myRank == 0)
	{
		dataSet.write(workBench->solvedVelData,H5::PredType::NATIVE_DOUBLE);
	}
*/
	
	return 0;
}
