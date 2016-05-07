#include "Simulation.h"
#include <iostream>
#include "hdf5.h"

//Pressure function can't be solved independently - how solve?
//Solve for velocity iteratively - then use points? Issue is if
//no stable solution exists

int main(int argc, char** argv)
{
	MPI_Init(&argc, &argv);
	
//	std::cout<<"Hello?"<<std::endl;

	Simulation* workBench = new Simulation;
	workBench->iterate();

	std::cout<<"Running"<<std::endl;
	if (workBench->getRank() == 0)
	{
		hid_t file, dataset;
		hid_t datatype, dataspace;
		file = H5Fcreate("output.hdf5",H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
		std::cout<<"Made a file"<<std::endl;
		hsize_t dimsf[2];
		dimsf[0] = workBench->getNt();
		dimsf[1] = workBench->getProblemSize();
		dataspace = H5Screate_simple(2, dimsf, NULL);
		datatype = H5Tcopy(H5T_NATIVE_DOUBLE);
		dataset = H5Dcreate2(file, "dset", datatype, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
		H5Dwrite(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, workBench->solvedVelData);
		H5Sclose(dataspace);
		H5Tclose(datatype);
		H5Dclose(dataset);
		H5Fclose(file); 

		std::cout<<"Hellllo"<<std::endl;
	}

	MPI_Finalize();
	
	return 0;
}
