#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <mpi.h>

//Pressure function can't be solved independently - how solve?

typedef struct
{
	double p;
	double u;
	double v;
} datumPoint;

class Simulation
{
	private:
		//Size values
		const static int nx = 41;
		const static int ny = 41;
		const static int nt = 100;
		const static int nit = 50;
		double dx;
		double dy;
		double dt;

		//Physical Values
		int rho=1;
		int nu=.1; 
		int F=10; //needed?

		int numCells;
		const static int problemSize = nx*ny;
		int leftOvers;
		int startingLocation;
		int xLocation, yLocation;

		//Solved data storage
		datumPoint* localVelData;
		double* localPrePresData;
		double* solvedPrePresData[nit];
		double* localB;
	
		double uijn,vijn,uim1jn,uip1jn,uijm1n,uijp1n,vim1jn,vip1jn,vijm1n,vijp1n,pip1jn,pim1jn,pijp1n,pijm1n;

		int* recCounts;

		//MPI data
		int counter = 1;
		int subCounter = 1;
		int myRank;
		int size;

	public:
		Simulation();
		double buildUpB(int xLocation, int yLocation);
		double pressureSolve(int xLocation, int yLocation, int i);
		double xMomentumSolve(int xLocation, int yLocation);
		double yMomentumSolve(int xLocation, int yLocation);
		void iterate(void);

		int getNt(void)
			{ 	return nt; 	}
		int getNx(void)
			{ 	return nx; 	}
		int getNy(void)
			{ 	return ny; 	}
		int getProblemSize(void)
			{ 	return problemSize; 	}
		int getRank(void)
			{ 	return myRank; 	}

		double solvedPMat[nt][ny][nx];
		double solvedUMat[nt][ny][nx];
		double solvedVMat[nt][ny][nx];
		datumPoint solvedVelData[nt][problemSize];
};
