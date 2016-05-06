#include "Simulation.h"
#include <iostream>
#include "H5Cpp.h"

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
		H5::H5File* file = new H5::H5File("output.hdf5",H5F_ACC_RDWR);
		hsize_t dimsf[2];
		dimsf[0] = workBench->getNt();
		dimsf[1] = workBench->getProblemSize();
		H5::DataSpace dataSpace(2, dimsf);
	 
		H5::DataSet* dataSet = new H5::DataSet(file->createDataSet("Data", H5::PredType::NATIVE_DOUBLE,dataSpace));

		dataSet->write(workBench->solvedVelData, H5::PredType::NATIVE_DOUBLE);
		std::cout<<"Hellllo"<<std::endl;
	}
	
	return 0;
}
