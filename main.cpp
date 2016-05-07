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
		hid_t file, datasetP, datasetU, datasetV;
		hid_t datatype, dataspace;
		file = H5Fcreate("output.hdf5",H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
		std::cout<<"Made a file"<<std::endl;

		hid_t offsetP[3] = {0,0,0};
		hid_t offsetU[3] = {1,1,1};
		hid_t offsetV[3] = {2,2,2};
		hid_t stride[3] = {3,3,3};

		hsize_t dimsf[3];
		dimsf[0] = workBench->getNt();
		dimsf[1] = workBench->getNy();
		dimsf[2] = workBench->getNx();
		dataspace = H5Screate_simple(3, dimsf, NULL);
		datatype = H5Tcopy(H5T_NATIVE_DOUBLE);

//		H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, offsetP, stride_3d, count_3d, NULL)


		datasetP = H5Dcreate2(file, "P", datatype, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
		H5Dwrite(datasetP, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, workBench->solvedPMat);
		datasetU = H5Dcreate2(file, "U", datatype, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
		H5Dwrite(datasetU, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, workBench->solvedUMat);
		datasetV = H5Dcreate2(file, "V", datatype, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
		H5Dwrite(datasetV, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, workBench->solvedVMat);

		H5Sclose(dataspace);
		H5Tclose(datatype);
		H5Dclose(datasetP);
		H5Dclose(datasetU);
		H5Dclose(datasetV);
		H5Fclose(file); 

		std::cout<<"Hellllo"<<std::endl;
	}

	MPI_Finalize();
	
	return 0;
}
