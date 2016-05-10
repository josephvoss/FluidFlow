#include "Simulation.h"
#include <iostream>
#include "hdf5.h"

//Pressure function can't be solved independently - how solve?
//Solve for velocity iteratively - then use points? Issue is if
//no stable solution exists

int main(int argc, char** argv)
{
	MPI_Init(&argc, &argv);
	
	Simulation* workBench = new Simulation;
	workBench->iterate();

	std::cout<<"Running"<<std::endl;
	if (workBench->getRank() == 0)
	{
		hid_t file, datasetP, datasetU, datasetV;
		hid_t datatype, dataspace, memspaceP, memspaceU, memspaceV;
		file = H5Fcreate("./logFiles/output.hdf5",H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
		std::cout<<"Made a file"<<std::endl;

		hsize_t offsetP[3] = {0,0,0};
		hsize_t offsetU[3] = {0,0,1};
		hsize_t offsetV[3] = {0,0,2};
		hsize_t stride1[3] = {1,1,1};
		hsize_t stride2[3] = {1,1,3};

		hsize_t dimsf[3], count[3];
		dimsf[0] = workBench->getNt();
		dimsf[1] = workBench->getNy();
		dimsf[2] = workBench->getNx()*3;
		count[0] = dimsf[0];
		count[1] = dimsf[1];
		count[2] = dimsf[2]/3;
		dataspace = H5Screate_simple(3, count, NULL);
		datatype = H5Tcopy(H5T_NATIVE_DOUBLE);

		//P
		H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, offsetP, stride1, count, NULL); //start stride count;
		memspaceP = H5Screate_simple(3, dimsf, NULL);
		H5Sselect_hyperslab(memspaceP, H5S_SELECT_SET, offsetP, stride2, count, NULL); //start stride count;
		datasetP = H5Dcreate2(file, "P", datatype, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
		H5Dwrite(datasetP, H5T_NATIVE_DOUBLE, memspaceP, dataspace, H5P_DEFAULT, &(workBench->solvedVelData[0][0].p));
		//U
		memspaceU = H5Screate_simple(3, dimsf, NULL);
		H5Sselect_hyperslab(memspaceU, H5S_SELECT_SET, offsetP, stride2, count, NULL); //start stride count;
		datasetU = H5Dcreate2(file, "U", datatype, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
		H5Dwrite(datasetU, H5T_NATIVE_DOUBLE, memspaceU, dataspace, H5P_DEFAULT, &(workBench->solvedVelData[0][0].u));

		//V
		memspaceV = H5Screate_simple(3, dimsf, NULL);
		H5Sselect_hyperslab(memspaceV, H5S_SELECT_SET, offsetP, stride2, count, NULL); //start stride count;
		datasetV = H5Dcreate2(file, "V", datatype, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
		H5Dwrite(datasetV, H5T_NATIVE_DOUBLE, memspaceV, dataspace, H5P_DEFAULT, &(workBench->solvedVelData[0][0].v));

		H5Sclose(dataspace);
		H5Sclose(memspaceP);
		H5Sclose(memspaceU);
		H5Tclose(datatype);
		H5Dclose(datasetP);
		H5Fclose(file); 

	}

	MPI_Finalize();
	
	return 0;
}
